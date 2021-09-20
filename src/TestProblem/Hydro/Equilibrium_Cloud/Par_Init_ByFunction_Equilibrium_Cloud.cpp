/***Total deep copy***/
#include "GAMER.h"

#include"Particle_IC_Constructor.h"
#define DEBUG
#ifdef PARTICLE

static RandomNumber_t *RNG = NULL;

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Equilibrium_Cloud
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------

void Par_Init_ByFunction_Equilibrium_Cloud( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] )
{
   // Define particles' attributes array
   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

   // Define the Particle IC Constructor
   Particle_IC_Constructor Filename_Loader;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   // only the master rank will construct the initial condition
   if ( MPI_Rank == 0 ){

      // Allocate memory for particles' attributes array
      Mass_AllRank = new real [NPar_AllRank];
      for (int d=0; d<3; d++)
      {
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      }

      // Input filenames as parameters into Filename_Loader
      Filename_Loader.Read_Filenames("Input__TestProb");
      int Par_Idx_Last = 0;
      double Par_Ratio = 0.0;

      for(int k=0;k<Filename_Loader.filenames.Cloud_Num;k++){
         // initialize Particle_IC_Constructor for each cloud
         Particle_IC_Constructor Cloud_Constructor;
         Cloud_Constructor.Load_Physical_Params(Filename_Loader.filenames,k);
         Cloud_Constructor.Init();

         // calculate particle number for each cloud
         if((Par_Ratio + Cloud_Constructor.params.Cloud_Par_Num_Ratio) > 1.0){
            Aux_Error( ERROR_INFO, "The sum of particle number ratios of all clouds exceeds 1!! Please check!");
         }
         if((k==Filename_Loader.filenames.Cloud_Num-1)&&((Par_Ratio + Cloud_Constructor.params.Cloud_Par_Num_Ratio) != 1.0)){
            Aux_Error( ERROR_INFO, "The sum of particle number ratios of all clouds is not equal to 1.0!! Please check!");
         }
         int Par_Idx_Start = Par_Idx_Last;
         Par_Idx_Last += Cloud_Constructor.params.Cloud_Par_Num_Ratio*NPar_AllRank;
         Par_Ratio += Cloud_Constructor.params.Cloud_Par_Num_Ratio;
         if(k==Filename_Loader.filenames.Cloud_Num-1){
            Par_Idx_Last = NPar_AllRank;
         }

         // set equilibrium initial conditions for each cloud
         // Input Mass_AllRank, Pos_AllRank, Vel_AllRank into this function for particles' masses, positions, and velocities.
         // NPar_AllRank is the total number (including all clouds) of particles.
         // Par_Idx_Start is the starting index for the particle in this cloud, and Par_Idx_Last is the last index.
         Cloud_Constructor.Par_SetEquilibriumIC(Mass_AllRank, Pos_AllRank, Vel_AllRank,NPar_AllRank,Par_Idx_Start,Par_Idx_Last);
         
      }//for(int k=0;k<Filename_Loader.filenames.Cloud_Num;k++)

   }//if ( MPI_Rank == 0 )
   
   
   // synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)   ParTime[p] = Time[0];
   
   // get the number of particles in each rank and set the corresponding offsets
   if ( NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
               NPar_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];
   int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                             // --> (ii) is necessary for OpenMPI version < 1.7

   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }

   // send particle attributes from the master rank to all ranks
   real *Mass   =   ParMass;
   real *Pos[3] = { ParPosX, ParPosY, ParPosZ };
   real *Vel[3] = { ParVelX, ParVelY, ParVelZ };

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   } 
#  endif
if ( MPI_Rank == 0 )
   {
      delete RNG;
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
} // FUNCTION : Par_Init_ByFunction_Equilibrium_Cloud

#endif // #ifdef PARTICLE
