/* Treat black hole as a partilce and put it in the soliton, with given initial position and velocity */
/* If BH_AddParForRestart is set, particle position will be reset; If BH_AddParForRestart is not set, treat the simulation as normal restart. */
#include "GAMER.h"
#include "TestProb.h"

// problem-specific global variables
// =======================================================================================
static double   Soliton_CM_MaxR;                         // maximum radius for determining Soliton CM
static double   Soliton_CM_TolErrR;                      // maximum allowed errors for determining Soliton CM
static double   SolitonCenter[3];                        // user defined center for placing the compressed soliton
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the compressed soliton density profile
static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static double  *Soliton_DensProf   = NULL;               // soliton density profile [radius/density]

       double   NSDHalfMassRadius;                       // NSD half mass radius; used for compute potential profile; will be called by extern
       double   NSDPotCenter[3];                         // user defined center for external NSD potential center, will be called by extern
static bool     first_run_flag;                          // flag suggesting first run (for determining whether write header in log file or not )
static bool     Fluid_Periodic_BC_Flag;                  // flag for checking the fluid boundary condtion is setup to periodic (0: user defined; 1: periodic)
#ifdef PARTICLE
static int      NewParAttTracerIdx = Idx_Undefined;      // particle attribute index for labelling particles
static int      WriteDataInBinaryFlag;                   // flag for determining output data type (0:text 1:binary Other: Do not write)

static bool     ParRefineFlag;                           // flag for refinement based on particles
static bool     BH_AddParForRestart;                     // flag for adding new particle after restart
static long     BH_AddParForRestart_NPar;                // number for particle will be added after restart
static char     Particle_Data_Filename[MAX_STRING];      // filename of the particles mass, initial position, initial velocity data
static char     Particle_Log_Filename[MAX_STRING];       // filename for recording particle data
static double  *Particle_Data_Table = NULL;              // particle data table [mass/position/velocity]
#endif
// =======================================================================================
//
// external potential routines
extern void Init_ExtPot_ELBDM_NSDPot();
//



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

   // errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif

# ifdef PARTICLE
   if ( ( OPT__INIT == INIT_BY_FUNCTION ) && ( amr->Par->Init != PAR_INIT_BY_FUNCTION ) )
      Aux_Error( ERROR_INFO, "must set PAR_INIT == PAR_INIT_BY_FUNCTION for OPT__INIT == INIT_BY_FUNCTION !!\n" );
# endif

// only accept OPT__INIT == INIT_BY_FUNCTION or OPT__INIT == INIT_BY_RESTART
   if ( ( OPT__INIT != INIT_BY_FUNCTION ) && ( OPT__INIT != INIT_BY_RESTART ) )
      Aux_Error( ERROR_INFO, "enforced to accept only OPT__INIT == INIT_BY_FUNCTION or OPT__INIT == INIT_BY_RESTART !!\n" );

// only accept OPT__RESTART_RESET == 1 or OPT__INIT == INIT_BY_FUNCTION for OPT__EXT_POT == EXT_POT_FUNC
//   if ( ( OPT__EXT_POT == EXT_POT_FUNC ) && ( OPT__RESTART_RESET != 1 ) && ( OPT__INIT != INIT_BY_FUNCTION ) )
//      Aux_Error( ERROR_INFO, "must set OPT__RESTART_RESET == 1 or OPT__INIT == INIT_BY_FUNCTION for OPT__EXT_POT == EXT_POT_FUNC !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == ELBDM && defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//                3. Must NOT call any EoS routine here since it hasn't been initialized at this point
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb_ELBDM_Black_Hole_in_Soliton";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",          &VARIABLE,               DEFAULT,          MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Soliton_CM_MaxR",          &Soliton_CM_MaxR,        NoMax_double,     Eps_double,       NoMax_double      );
   ReadPara->Add( "Soliton_CM_TolErrR",       &Soliton_CM_TolErrR,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "SolitonCenter_x",          &SolitonCenter[0],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "SolitonCenter_y",          &SolitonCenter[1],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "SolitonCenter_z",          &SolitonCenter[2],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  NoDef_str,     Useless_str,      Useless_str     );
   ReadPara->Add( "Fluid_Periodic_BC_Flag",   &Fluid_Periodic_BC_Flag,  false,           Useless_bool,     Useless_bool      );
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   {
      ReadPara->Add( "NSDHalfMassRadius",        &NSDHalfMassRadius,      Eps_double,       Eps_double,       NoMax_double      );
      ReadPara->Add( "NSDPotCenter_x",           &NSDPotCenter[0],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "NSDPotCenter_y",           &NSDPotCenter[1],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "NSDPotCenter_z",           &NSDPotCenter[2],           0.0,           NoMin_double,     NoMax_double      );
   }

#ifdef PARTICLE
   ReadPara->Add( "ParRefineFlag",            &ParRefineFlag,            false,         Useless_bool,      Useless_bool      );
   ReadPara->Add( "WriteDataInBinaryFlag",    &WriteDataInBinaryFlag,         -1,          NoMin_int,      NoMax_int         );
   if ( ( amr->Par->Init == PAR_INIT_BY_FUNCTION ) && ( OPT__INIT == INIT_BY_FUNCTION ) )
      ReadPara->Add( "Particle_Data_Filename",   Particle_Data_Filename,  Useless_str,     Useless_str,       Useless_str       );
   ReadPara->Read( FileName );
   delete ReadPara;

   ReadPara  = new ReadPara_t;
   if ( ( WriteDataInBinaryFlag == 0 ) || ( WriteDataInBinaryFlag == 1 ) )
      ReadPara->Add( "Particle_Log_Filename",    Particle_Log_Filename,   Useless_str,     Useless_str,       Useless_str       );
   if ( ( OPT__RESTART_RESET == 1 ) || ( OPT__INIT == INIT_BY_RESTART ) )
   {
      ReadPara->Add( "BH_AddParForRestart",      &BH_AddParForRestart,      false,         Useless_bool,      Useless_bool      );
      ReadPara->Read( FileName );
      delete ReadPara;

      if ( BH_AddParForRestart == 1 )
      {
         ReadPara  = new ReadPara_t;
         ReadPara->Add( "BH_AddParForRestart_NPar", &BH_AddParForRestart_NPar,  -1L,          NoMin_long,        NoMax_long        );
         ReadPara->Add( "Particle_Data_Filename",   Particle_Data_Filename,  Useless_str,     Useless_str,       Useless_str       );
         ReadPara->Read( FileName );
         delete ReadPara;
      }
   }
   else
      ReadPara->Read( FileName );
      delete ReadPara;
#else
   ReadPara->Read( FileName );
   delete ReadPara;
#endif

// (1-2) set the default values
   if ( Soliton_CM_TolErrR < 0.0 )           Soliton_CM_TolErrR = 1.0*amr->dh[MAX_LEVEL];

// (1-3) check the runtime parameters
#ifdef PARTICLE
   if ( ( BH_AddParForRestart == 1 ) &&  ( OPT__RESTART_RESET != 1 ) && ( OPT__INIT != INIT_BY_RESTART ) )  
      Aux_Error( ERROR_INFO, "must set OPT__RESTART_RESET == 1 or OPT__INIT == INIT_BY_RESTART if BH_AddParForRestart is enabled !!\n" );
#endif
   if ( Fluid_Periodic_BC_Flag )  // use periodic boundary condition
   {
      for ( int direction = 0; direction < 6; direction++ )
      {
         if ( OPT__BC_FLU[direction] != BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "must set periodic BC for fluid --> reset OPT__BC_FLU[%d] to 1 !!\n", direction );
      }
   }
   else  // use user define boundary condition
   {
      for ( int direction = 0; direction < 6; direction++ )
      {
         if ( OPT__BC_FLU[direction] != BC_FLU_USER )
            Aux_Error( ERROR_INFO, "must adopt user defined BC for fluid --> reset OPT__BC_FLU[%d] to 4 !!\n", direction );
      }
   }


// (2) load the reference soliton density profile
   if ( OPT__INIT == INIT_BY_FUNCTION )
   {
//    load the reference profile
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for Soliton_DensProf
      const int  NCol         = 2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

      Soliton_DensProf_NBin = Aux_LoadTable( Soliton_DensProf, Soliton_DensProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );
   }
   
// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=================================================================================\n"  );
      Aux_Message( stdout, "  test problem ID                              = %d\n",     TESTPROB_ID               );
      Aux_Message( stdout, "  soliton CM max radius                        = %13.6e\n", Soliton_CM_MaxR           );
      Aux_Message( stdout, "  soliton CM tolerated error                   = %13.6e\n", Soliton_CM_TolErrR        );
      Aux_Message( stdout, "  soliton center_x                             = %13.6e\n", SolitonCenter[0]          );
      Aux_Message( stdout, "  soliton center_y                             = %13.6e\n", SolitonCenter[1]          );
      Aux_Message( stdout, "  soliton center_z                             = %13.6e\n", SolitonCenter[2]          );
      Aux_Message( stdout, "  soliton density profile filename             = %s\n",     Soliton_DensProf_Filename );
      Aux_Message( stdout, "  number of bins of soliton density profile    = %d\n",     Soliton_DensProf_NBin     );
      Aux_Message( stdout, "  fluid periodic boundary condition flag       = %d\n",     Fluid_Periodic_BC_Flag    );
      if ( OPT__EXT_POT == EXT_POT_FUNC )
      {
         Aux_Message( stdout, "  NSD half mass radius                         = %13.6e\n", NSDHalfMassRadius         );
         Aux_Message( stdout, "  NSD potential center_x                       = %13.6e\n", NSDPotCenter[0]           );
         Aux_Message( stdout, "  NSD potential center_y                       = %13.6e\n", NSDPotCenter[1]           );
         Aux_Message( stdout, "  NSD potential center_z                       = %13.6e\n", NSDPotCenter[2]           );
      }

#ifdef PARTICLE
      Aux_Message( stdout, "  refine grid based on particles               = %d\n",     ParRefineFlag              );
      Aux_Message( stdout, "  write particle data in binary format         = %d\n",     WriteDataInBinaryFlag      );
      if ( (WriteDataInBinaryFlag == 0) || (WriteDataInBinaryFlag == 1) )
         Aux_Message( stdout, "  particle log filename                        = %s\n",     Particle_Log_Filename     );
      if ( ( amr->Par->Init == PAR_INIT_BY_FUNCTION ) && ( OPT__INIT == INIT_BY_FUNCTION ) )
         Aux_Message( stdout, "  particle data filename                       = %s\n",     Particle_Data_Filename    );
      if ( ( OPT__RESTART_RESET == 1 ) && ( OPT__INIT == INIT_BY_RESTART ) )
      {
         Aux_Message( stdout, "  add particles after restart                  = %d\n",     BH_AddParForRestart        );
         if ( BH_AddParForRestart == 1 )
         {
            Aux_Message( stdout, "  number of particles to be added              = %ld\n",    BH_AddParForRestart_NPar  );
            Aux_Message( stdout, "  particle data filename                       = %s\n",     Particle_Data_Filename    );
         }
      }
#endif
      Aux_Message( stdout, "=================================================================================\n"   );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



#ifdef PARTICLE

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByRestart_Black_Hole_in_Soliton()
// Description :  User-specified initialization
//
// Note        :  1. Add particles after restart
//                2. Set the central coordinates of tidal field
//
// Parameter   :  None
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime
//-------------------------------------------------------------------------------------------------------
static void Par_Init_ByRestart_Black_Hole_in_Soliton() 
{
   const bool RowMajor_No_particle_data             = false;                 // load data into the column-major order
   const bool AllocMem_Yes_particle_data            = true;                  // allocate memory for Soliton_DensProf
   const int  NCol_particle_data                    = 7;                     // total number of columns to load for particle data
   const int  Col_particle_data[NCol_particle_data] = {0, 1, 2, 3, 4, 5, 6}; // target columns: (mass, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z)
   long       BH_AddParForRestart_Check;         

   if ( MPI_Rank == 0 )
   {
       Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

       if ( !Aux_CheckFileExist(Particle_Data_Filename) )
          Aux_Error( ERROR_INFO, "Error!! Initial particle data file %s does not exists!!\n", Particle_Data_Filename);

       BH_AddParForRestart_Check = Aux_LoadTable( Particle_Data_Table, Particle_Data_Filename, NCol_particle_data, Col_particle_data, RowMajor_No_particle_data, AllocMem_Yes_particle_data );

       if ( BH_AddParForRestart_Check != BH_AddParForRestart_NPar )
          Aux_Error( ERROR_INFO, "BH_AddParForRestart_Check(%ld) != BH_AddParForRestart_Npar(%ld) !!\n", BH_AddParForRestart_Check, BH_AddParForRestart_NPar );
   }
   MPI_Bcast(&BH_AddParForRestart_NPar, 1, MPI_LONG, 0, MPI_COMM_WORLD);

   const long   NNewPar        = ( MPI_Rank == 0 ) ? BH_AddParForRestart_NPar : 0;
   const long   NPar_AllRank   = NNewPar;

   real *NewParAtt[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)   NewParAtt[v] = new real [NNewPar];

// set particle attributes
// ============================================================================================================
   real *Time_AllRank      = NewParAtt[PAR_TIME];
   real *Mass_AllRank      = NewParAtt[PAR_MASS];
   real *Type_AllRank      = NewParAtt[PAR_TYPE];
   real *Pos_AllRank[3]    = { NewParAtt[PAR_POSX], NewParAtt[PAR_POSY], NewParAtt[PAR_POSZ] };
   real *Vel_AllRank[3]    = { NewParAtt[PAR_VELX], NewParAtt[PAR_VELY], NewParAtt[PAR_VELZ] };
   real *TracerIdx_AllRank = NewParAtt[NewParAttTracerIdx];


// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      double  TotM;
      const double *Mass_table        = Particle_Data_Table;
      const double *Position_table[3];
      const double *Velocity_table[3];

      for (int d=0; d<3; d++)
      {
         Position_table[d] = Particle_Data_Table+(1+d)*NPar_AllRank;
         Velocity_table[d] = Particle_Data_Table+(4+d)*NPar_AllRank;
      }


//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       time
         Time_AllRank[p] = Time[0];


//       mass
         Mass_AllRank[p] = Mass_table[p];
         TotM += Mass_AllRank[p];


//       position
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = Position_table[d][p];


//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }

//       velocity
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = Velocity_table[d][p];

//       particle type
         Type_AllRank[p] = PTYPE_GENERIC_MASSIVE;   // use root rank to declare type and MPI_Scatter to other ranks, for generality such that particle type might be different for different particles

//       particle tracer index
         TracerIdx_AllRank[p] = (real)p;

      } // for (long p=0; p<NPar_AllRank; p++)

//      Aux_Message( stdout, "=====================================================================================================\n" );
//      Aux_Message( stdout, "Total mass = %13.7e\n",  TotM );
//      for (long p=0; p<NPar_AllRank; p++)
//           Aux_Message( stdout, "Pisition for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Pos_AllRank[0][p], Pos_AllRank[1][p], Pos_AllRank[2][p] );
//      for (long p=0; p<NPar_AllRank; p++)
//           Aux_Message( stdout, "Velocity for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Vel_AllRank[0][p], Vel_AllRank[1][p], Vel_AllRank[2][p] );
//      Aux_Message( stdout, "=====================================================================================================\n" );

   } // if ( MPI_Rank == 0 )

// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAtt );


// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] NewParAtt[v];

// refine the grids
   if ( ParRefineFlag )
   {
#  ifdef LOAD_BALANCE
      const bool   Redistribute_Yes = true;
      const bool   SendGridData_Yes = true;
      const bool   ResetLB_Yes      = true;
      const double Par_Weight       = amr->LB->Par_Weight;
      const UseLBFunc_t UseLB       = USELB_YES;
      const bool   SortRealPatch_No = false;
#  else
      const UseLBFunc_t UseLB       = USELB_NO;
#  endif

      for (int lv=0; lv<MAX_LEVEL; lv++)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ... ", lv );

         Flag_Real( lv, UseLB );

         Refine( lv, UseLB );

#     ifdef LOAD_BALANCE
//       no need to exchange potential since we haven't calculated it yet
         Buf_GetBufferData( lv,   amr->FluSg[lv  ], NULL_INT, NULL_INT, DATA_AFTER_REFINE, _TOTAL, _NONE, Flu_ParaBuf, USELB_YES );

         Buf_GetBufferData( lv+1, amr->FluSg[lv+1], NULL_INT, NULL_INT, DATA_AFTER_REFINE, _TOTAL, _NONE, Flu_ParaBuf, USELB_YES );

         LB_Init_LoadBalance( Redistribute_Yes, SendGridData_Yes, Par_Weight, ResetLB_Yes, SortRealPatch_No, lv+1 );
#     endif

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=0; lv<MAX_LEVEL; lv++)
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByRestart_Black_Hole_in_Soliton



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Black_Hole_in_Soliton
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
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Black_Hole_in_Soliton( const long NPar_ThisRank, const long NPar_AllRank,
                                                real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                                real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                                real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{
   const bool RowMajor_No_particle_data             = false;                 // load data into the column-major order
   const bool AllocMem_Yes_particle_data            = true;                  // allocate memory for Soliton_DensProf
   const int  NCol_particle_data                    = 7;                     // total number of columns to load for particle data
   const int  Col_particle_data[NCol_particle_data] = {0, 1, 2, 3, 4, 5, 6}; // target columns: (mass, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z)
   long       BH_AddParForRestart_Check;         

   if ( MPI_Rank == 0 )
   {
       Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

       if ( !Aux_CheckFileExist(Particle_Data_Filename) )
          Aux_Error( ERROR_INFO, "Error!! Initial particle data file %s does not exists!!\n", Particle_Data_Filename);

       BH_AddParForRestart_Check = Aux_LoadTable( Particle_Data_Table, Particle_Data_Filename, NCol_particle_data, Col_particle_data, RowMajor_No_particle_data, AllocMem_Yes_particle_data );

       if ( BH_AddParForRestart_Check != NPar_AllRank )
          Aux_Error( ERROR_INFO, "BH_AddParForRestart_Check(%ld) != NPar_AllRank(%ld) !!\n", BH_AddParForRestart_Check, NPar_AllRank );
   }


   real *NewParAtt[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)   NewParAtt[v] = new real [NPar_AllRank];

// set particle attributes
// ============================================================================================================
   real *Time_AllRank      = NewParAtt[PAR_TIME];
   real *Mass_AllRank      = NewParAtt[PAR_MASS];
   real *Type_AllRank      = NewParAtt[PAR_TYPE];
   real *Pos_AllRank[3]    = { NewParAtt[PAR_POSX], NewParAtt[PAR_POSY], NewParAtt[PAR_POSZ] };
   real *Vel_AllRank[3]    = { NewParAtt[PAR_VELX], NewParAtt[PAR_VELY], NewParAtt[PAR_VELZ] };
   real *TracerIdx_AllRank = NewParAtt[NewParAttTracerIdx];


// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      double  TotM;
      const double *Mass_table        = Particle_Data_Table;
      const double *Position_table[3];
      const double *Velocity_table[3];

      for (int d=0; d<3; d++)
      {
         Position_table[d] = Particle_Data_Table+(1+d)*NPar_AllRank;
         Velocity_table[d] = Particle_Data_Table+(4+d)*NPar_AllRank;
      }


//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       time
         Time_AllRank[p] = Time[0];


//       mass
         Mass_AllRank[p] = Mass_table[p];
         TotM += Mass_AllRank[p];


//       position
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = Position_table[d][p];


//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }

//       velocity
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = Velocity_table[d][p];

//       particle type
         Type_AllRank[p] = PTYPE_GENERIC_MASSIVE;   // use root rank to declare type and MPI_Scatter to other ranks, for generality such that particle type might be different for different particles

//       particle tracer index
         TracerIdx_AllRank[p] = (real)p;

      } // for (long p=0; p<NPar_AllRank; p++)

//      Aux_Message( stdout, "=====================================================================================================\n" );
//      Aux_Message( stdout, "Total mass = %13.7e\n",  TotM );
//      for (long p=0; p<NPar_AllRank; p++)
//           Aux_Message( stdout, "Pisition for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Pos_AllRank[0][p], Pos_AllRank[1][p], Pos_AllRank[2][p] );
//      for (long p=0; p<NPar_AllRank; p++)
//           Aux_Message( stdout, "Velocity for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Vel_AllRank[0][p], Vel_AllRank[1][p], Vel_AllRank[2][p] );
//      Aux_Message( stdout, "=====================================================================================================\n" );

   } // if ( MPI_Rank == 0 )

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
   real *Mass      =   ParMass;
   real *Type      =   ParType;
   real *Pos[3]    = { ParPosX, ParPosY, ParPosZ };
   real *Vel[3]    = { ParVelX, ParVelY, ParVelZ };
   real *TracerIdx = AllAttribute[NewParAttTracerIdx];

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank,      NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Scatterv( Type_AllRank,      NSend, SendDisp, MPI_DOUBLE, Type, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Scatterv( TracerIdx_AllRank, NSend, SendDisp, MPI_DOUBLE, TracerIdx, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank,      NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   MPI_Scatterv( Type_AllRank,      NSend, SendDisp, MPI_FLOAT,  Type, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   MPI_Scatterv( TracerIdx_AllRank, NSend, SendDisp, MPI_FLOAT,  TracerIdx, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   }
#  endif

// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] NewParAtt[v];

} // FUNCTION : Par_Init_ByFunction_Black_Hole_in_Soliton



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Particle_Data_Text
// Description :  Output the particle position and velocity in text file
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Record_Particle_Data_Text( char *FileName )
{
// check
//   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
//      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File;
   static bool first_enter_flag_particle = true;
   
// header
   if ( MPI_Rank == 0 )
   {
      if ( first_enter_flag_particle )
      {
         if ( !Aux_CheckFileExist(FileName) )
         {
            File = fopen( FileName, "w" );
            fprintf( File, "# Time                    Step                    Active_Particles   ");
            for (int v=0; v<PAR_NATT_TOTAL; v++)
            for (int v=0; v<PAR_NATT_TOTAL; v++)
                fprintf( File, "  %*s", (v==0)?20:21, ParAttLabel[v] );
            fprintf( File, "\n" );
         }
         else if ( first_run_flag )
         {
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
            File = fopen( FileName, "a" );
            fprintf( File, "# Time                    Step                    Active_Particles   ");

            for (int v=0; v<PAR_NATT_TOTAL; v++)
               fprintf( File, "  %*s", (v==0)?20:21, ParAttLabel[v] );
            fprintf( File, "\n" );
            first_run_flag = false;
         }
      }
      fclose( File );
      first_enter_flag_particle = false;
   }

// data
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         File = fopen( FileName, "a" );

         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
//          skip inactive particles
            if ( amr->Par->Mass[p] < 0.0 )   continue;

            fprintf( File, "%20.14e    %13ld    %13ld          ", Time[0], Step, amr->Par->NPar_Active_AllRank );
            for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %21.14e", amr->Par->Attribute[v][p] );

            fprintf( File, "\n" );
         }

         fclose( File );
      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
} // FUNCTION : Record_Particle_Data_Text



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Particle_Data_Binary
// Description :  Output the particle position and velocity in binary file
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Record_Particle_Data_Binary( char *FileName )
{
// check
//   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
//      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File;
   int par_natt_total = PAR_NATT_TOTAL;
   
// open the file by root rank
   if ( MPI_Rank == 0 )
   {
       File = fopen( FileName, "w" );                  // overwrite the file no matter how
      if ( first_run_flag )
      {
//          File = fopen( FileName, "w" );                  // overwrite the file no matter how
//          int par_natt_total = PAR_NATT_TOTAL;
//          fwrite(&par_natt_total, sizeof(int), 1, File);  // write number of attribute only at simulation start and first_run_flag == true
          first_run_flag = false;
      }
      fclose( File );
   }
   MPI_Barrier( MPI_COMM_WORLD );

// data
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         File = fopen( FileName, "a" );
         if ( MPI_Rank == 0 )
         {
             fwrite(&(Time[0])                      , sizeof(double), 1, File);  // write level 0 time by rank == 0 for each time step
             fwrite(&Step                           , sizeof(long)  , 1, File);  // write step index by rank == 0 for each time step  
             fwrite(&par_natt_total                 , sizeof(int)   , 1, File);  // write number of attribute by rank == 0 for each time step 
             fwrite(&(amr->Par->NPar_Active_AllRank), sizeof(long)  , 1, File);  // write active particle number  by rank ==0 for each time step
         }

         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
//          skip inactive particles
            if ( amr->Par->Mass[p] < 0.0 )   continue;
            for (int v=0; v<PAR_NATT_TOTAL; v++)   fwrite( &(amr->Par->Attribute[v][p]),  sizeof(real), 1, File );  // write all attribute for selected particle for each time step
         }

         fclose( File );
      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
} // FUNCTION : Record_Particle_Data_Binary



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_Black_Hole_in_Soliton
// Description :  Add the problem-specific particle attributes: used for labelling particle
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void AddNewParticleAttribute_Black_Hole_in_Soliton(void)
{
   if ( NewParAttTracerIdx == Idx_Undefined )
      NewParAttTracerIdx = AddParticleAttribute( "ParticleTracerIdx" );
}

#endif // end of ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  GetPhase
// Description :  Calculate the wave function phase based on real/imaginary part
//
// Note        :  1.  Will be called whenever phase is needed
//
// Parameter   :  real dens_sqrt: square root of wave function 
//                real real_part: real part of wave function
//                real imag_part: imaginary part of wave function 
//
// Return      :  phase
//-------------------------------------------------------------------------------------------------------
static double GetPhase(real dens_sqrt, real real_part, real imag_part)
{
   double phase;
   if ( fabs(real_part) < fabs(imag_part) )  // use acos when abs(real_part) < abs(imag_part) since it will be more precise in that regeion
   {
      phase = acos (real_part/dens_sqrt);
      if ( imag_part < 0. )
         phase = 2.*M_PI-phase;
   }
   else // use asin when abs(real_part) >= abs(imag_part) since it will be more precise in that region
   {
      phase = asin (imag_part/dens_sqrt);
      if ( real_part < 0. )
         phase = 1.*M_PI-phase;
   }

   // makes the phase always between [0,2\pi)
   if ( phase < 0.)
      phase += 2.*M_PI;
   else if ( phase>=2.*M_PI )
      phase -= 2.*M_PI;
   //

   if ( ( phase < 0. ) || ( phase >= 2.*M_PI ) )
      Aux_Error( ERROR_INFO, "Phase %.8e is not in range [0,2*pi) !!\n", phase );
   if ( phase!=phase )
      Aux_Error( ERROR_INFO, "Phase is NaN !! Square root of density is %.8e ; real part = %.8e ; imaginary part = %.8e \n", dens_sqrt, real_part, imag_part );
   return phase;
} // FUNCTION : GetPhase




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_ELBDM_Black_Hole_in_Soliton
// Description :  Set the particle IC if BH_AddParForRestart is enabled; erase the soliton initial velocity by phase modulation scheme if EraseSolVelFlag is enabled; treated as normal restart if neither of them is enabled 
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Init_User_ELBDM_Black_Hole_in_Soliton(void)
{
   if ( ( OPT__INIT == INIT_BY_FUNCTION ) || ( OPT__RESTART_RESET == 1 ) )
      first_run_flag = true;
   else
      first_run_flag = false;
#ifdef PARTICLE
   if ( BH_AddParForRestart == 1 )
      Par_Init_ByRestart_Black_Hole_in_Soliton();
#endif

} // FUNCTION : Init_User_ELBDM_Black_Hole_in_Soliton



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
   const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

   double r_tar    = sqrt( SQR(x-SolitonCenter[0]) + SQR(y-SolitonCenter[1]) + SQR(z-SolitonCenter[2]) );
   double dens_tar = Mis_InterpolateFromTable( Soliton_DensProf_NBin, Table_Radius, Table_Density, r_tar );

// linear interpolation
   if ( dens_tar == NULL_REAL )
   {
      if      ( r_tar <  Table_Radius[0] )
         dens_tar = Table_Density[0];

      else if ( r_tar >= Table_Radius[Soliton_DensProf_NBin-1] )
         dens_tar = Table_Density[Soliton_DensProf_NBin-1];

      else
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                    r_tar, Table_Radius[0], Table_Radius[Soliton_DensProf_NBin-1] );
   }
   fluid[DENS] = dens_tar;
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Soliton
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//                2. Set the BC as isolated
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
////-------------------------------------------------------------------------------------------------------
static void BC_Soliton( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_HALO



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass
// Description :  Record the center of mass (CM)
//
// Note        :  1. Invoked by  recursively
//                2. Only include cells within CM_MaxR from CM_Old[] when updating CM
//
// Parameter   :  CM_Old[] : Previous CM
//                CM_New[] : New CM to be returned
//-------------------------------------------------------------------------------------------------------
static void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR, const long DensMode )
{

   const double CM_MaxR2          = SQR( CM_MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;

   int   *PID0List = NULL;
   double M_ThisRank, MR_ThisRank[3], M_AllRank, MR_AllRank[3];
   real (*TotalDens)[PS1][PS1][PS1];

   M_ThisRank = 0.0;
   for (int d=0; d<3; d++)    MR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {

//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensMode, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, 0.0, DE_Consistency_No );

      delete [] PID0List;

//    calculate the center of mass
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         double x, y, z, dx, dy, dz;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;  dz = z - CM_Old[2];
                                       if ( Periodic ) {
                                          if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                          else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                       }
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;  dy = y - CM_Old[1];
                                       if ( Periodic ) {
                                          if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                          else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                       }
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;  dx = x - CM_Old[0];
                                       if ( Periodic ) {
                                          if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                          else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                       }

//          only include cells within CM_MaxR
            const double R2 = SQR(dx) + SQR(dy) + SQR(dz);
            if ( R2 < CM_MaxR2 )
            {
               const double dm = TotalDens[PID][k][j][i]*dv;

               M_ThisRank     += dm;
               MR_ThisRank[0] += dm*x;
               MR_ThisRank[1] += dm*y;
               MR_ThisRank[2] += dm*z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks to calculate the CM
// --> note that all ranks will get CM_New[]
   MPI_Allreduce( &M_ThisRank, &M_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( MR_ThisRank, MR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)    CM_New[d] = MR_AllRank[d] / M_AllRank;

// map the new CM back to the simulation domain
   if ( Periodic )
   for (int d=0; d<3; d++)
   {
      if      ( CM_New[d] >= amr->BoxSize[d] )  CM_New[d] -= amr->BoxSize[d];
      else if ( CM_New[d] < 0.0              )  CM_New[d] += amr->BoxSize[d];

   }

   for (int d=0; d<3; d++)
      if ( CM_New[d] >= amr->BoxSize[d]  ||  CM_New[d] < 0.0 )
         Aux_Error( ERROR_INFO, "CM_New[%d] = %14.7e lies outside the domain !!\n", d, CM_New[d] );

} // FUNCTION : GetCenterOfMass



//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//                4. When simulation starts, this function will be called to calculate center of whole halo for calculating initial density 
//                   profile, which will be used to calculate correlation function, if ComputeCorrelation is true.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Record_CenterOfMass( void )
{
   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];
   long   DensMode;
   DensMode  = _DENS;  // find peak densit for FDM component only

   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;

// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    get the total density on grids
      real (*TotalDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensMode, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, 0.0, DE_Consistency_No );

      delete [] PID0List;

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
         for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
         for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

            dens = TotalDens[PID][k][j][i];
            pote = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];

            if ( dens > max_dens_loc )
            {
               max_dens_loc        = dens;
               real_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i];
               imag_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i];
               max_dens_pos_loc[0] = x;
               max_dens_pos_loc[1] = y;
               max_dens_pos_loc[2] = z;
            }

            if ( pote < min_pote_loc )
            {
               min_pote_loc        = pote;
               min_pote_pos_loc[0] = x;
               min_pote_pos_loc[1] = y;
               min_pote_pos_loc[2] = z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = real_loc;
   send[2] = imag_loc;
   send[3] = max_dens_pos_loc[0];
   send[4] = max_dens_pos_loc[1];
   send[5] = max_dens_pos_loc[2];
   send[6] = min_pote_loc;
   send[7] = min_pote_pos_loc[0];
   send[8] = min_pote_pos_loc[1];
   send[9] = min_pote_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );

// record the maximum density and center coordinates
   double max_dens      = -__DBL_MAX__;
   double min_pote      = +__DBL_MAX__;
   int    max_dens_rank = -1;
   int    min_pote_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }

         if ( recv[r][6] < min_pote )
         {
            min_pote      = recv[r][6];
            min_pote_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );

      if ( min_pote_rank < 0  ||  min_pote_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect min_pote_rank (%d) !!\n", min_pote_rank );

      static bool first_enter_flag_center = true;
      FILE       *file_center;
      if ( first_enter_flag_center )
      {
         if ( !Aux_CheckFileExist(filename_center) )
         {
            file_center = fopen( filename_center, "w" );
            fprintf( file_center, "# %s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Real", "Imag", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter_s", "CM_x_s", "CM_y_s", "CM_z_s");
            fclose( file_center );
         }
         else if ( first_run_flag )
         {
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
            file_center = fopen( filename_center, "a" );
            fprintf( file_center, "# %s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Real", "Imag", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter_s", "CM_x_s", "CM_y_s", "CM_z_s");
            fclose( file_center );
#ifndef PARTICLE
            first_run_flag = false;   // if #define PARTICLE, first_run_flag will be turned to false after recording the data for first time step, so in that case no need to turn it to false here
#endif
         }
         first_enter_flag_center = false;
      }

      file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2], recv[max_dens_rank][3],
                              recv[max_dens_rank][4], recv[max_dens_rank][5], recv[min_pote_rank][6], recv[min_pote_rank][7],
                              recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   double TolErrR2 = SQR( Soliton_CM_TolErrR );
   const int    NIterMax = 20;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];
   
   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   
   while ( true )
   {
      GetCenterOfMass( CM_Old, CM_New, Soliton_CM_MaxR, _DENS );      // for soliton center of mass, use FDM density
   
      dR2 = SQR( CM_Old[0] - CM_New[0] )
          + SQR( CM_Old[1] - CM_New[1] )
          + SQR( CM_Old[2] - CM_New[2] );
      NIter ++;
   
      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( CM_Old, CM_New, sizeof(double)*3 );
   }
   
   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > Soliton_CM_TolErrR (%13.7e) !!\n", sqrt(dR2), Soliton_CM_TolErrR );
   
      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }
   delete [] recv;

} // FUNCTION : Record_CenterOfMass 



//-------------------------------------------------------------------------------------------------------
// Function    :  Do_COM_and_CF
// Description :  Do record center of mass and calculate correlation function 
//
// Note        :  1. It will call center of mass routine 
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Do_COM_and_CF( void )
{
   Record_CenterOfMass();
#ifdef PARTICLE
   char Particle_Log_Filename_Full[MAX_STRING];
   if ( WriteDataInBinaryFlag == 0 )
   {
       sprintf(Particle_Log_Filename_Full, "%s.txt", Particle_Log_Filename);
       Record_Particle_Data_Text(Particle_Log_Filename_Full);
   }
   else if ( WriteDataInBinaryFlag == 1 )
   {
       sprintf(Particle_Log_Filename_Full, "%s_StepIdx=%06ld.bin", Particle_Log_Filename, Step);
       Record_Particle_Data_Binary(Particle_Log_Filename_Full);
   }
   else
       first_run_flag = false;
#endif
}



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Black_Hole_in_Soliton
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
static void End_Black_Hole_in_Soliton()
{
#ifdef PARTICLE
   delete [] Particle_Data_Table;
#endif
} // FUNCTION : End_Black_Hole_in_Soliton
#endif // end of if ( MODEL == ELBDM && defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Black_Hole_in_Soliton()
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Black_Hole_in_Soliton()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

//#  if ( NCOMP_PASSIVE_USER > 0 )
//   if ( MPI_Rank == 0)     Aux_Message( stdout, "NCOMP_PASSIVE_USER = %d \n", NCOMP_PASSIVE);
//#  endif

// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();

   Init_Function_User_Ptr      = SetGridIC;
   BC_User_Ptr                 = BC_Soliton;
   Aux_Record_User_Ptr         = Do_COM_and_CF;
   Init_User_Ptr               = Init_User_ELBDM_Black_Hole_in_Soliton;
   Init_ExtPot_Ptr             = Init_ExtPot_ELBDM_NSDPot;
   End_User_Ptr                = End_Black_Hole_in_Soliton;
#  ifdef PARTICLE
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_Black_Hole_in_Soliton;
   Par_Init_ByFunction_Ptr     = Par_Init_ByFunction_Black_Hole_in_Soliton;
#  endif // #ifdef PARTICLE
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )

// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
   Src_Init_User_Ptr           = NULL; // option: SRC_USER;                example: SourceTerms/User_Template/CPU_Src_User_Template.cpp

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Black_Hole_in_Soliton
