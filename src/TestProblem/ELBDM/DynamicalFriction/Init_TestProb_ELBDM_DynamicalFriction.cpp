#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
double DynFri_CM_MaxR;           // maximum radius for determining CM
double DynFri_CM_TolErrR;        // maximum allowed errors for determining CM
double DynFri_CM_Init[3];        // initial guess of CM

double DynFri_ProfM_MaxR;        // maximum radius for computing the mass profile
double DynFri_ProfM_BinSize;     // bin size for computing the mass profile

double DynFri_CM_Old[3];         // previous CM
double DynFri_CM_New[3];         // new CM
// =======================================================================================




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

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
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
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,                 DEFAULT,             MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "DynFri_CM_InitX",      &DynFri_CM_Init[0],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_InitY",      &DynFri_CM_Init[1],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_InitZ",      &DynFri_CM_Init[2],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_MaxR",       &DynFri_CM_MaxR,           -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_TolErrR",    &DynFri_CM_TolErrR,        -1.0,                NoMin_double,     NoMax_double      );

   ReadPara->Add( "DynFri_ProfM_MaxR",    &DynFri_ProfM_MaxR,        -1.0,                Eps_double,       NoMax_double      );
   ReadPara->Add( "DynFri_ProfM_BinSize", &DynFri_ProfM_BinSize,     -1.0,                Eps_double,       NoMax_double      );
// ********************************************************************************************************************************

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   for (int d=0; d<3; d++)
   {
      if ( DynFri_CM_Init[d] < 0.0 )
      {
         if ( OPT__INIT == INIT_BY_RESTART )
            Aux_Error( ERROR_INFO, "Must specify DynFri_CM_Init* during restart !!\n" );

         else
            DynFri_CM_Init[d] = amr->BoxCenter[d];
      }
   }

   if ( DynFri_CM_MaxR    < 0 )  DynFri_CM_MaxR    = 20.0*amr->dh[MAX_LEVEL];
   if ( DynFri_CM_TolErrR < 0 )  DynFri_CM_TolErrR =  1.0*amr->dh[MAX_LEVEL];

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )   Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported !!\n" );


// (2) set the problem-specific derived parameters
   for (int d=0; d<3; d++)    DynFri_CM_New[d] = DynFri_CM_Init[d];


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID      = %d\n",     TESTPROB_ID          );
      Aux_Message( stdout, "  CM_MaxR              = %13.7e\n", DynFri_CM_MaxR       );
      Aux_Message( stdout, "  CM_TolErrR           = %13.7e\n", DynFri_CM_TolErrR    );
      Aux_Message( stdout, "  CM_InitX             = %13.7e\n", DynFri_CM_Init[0]    );
      Aux_Message( stdout, "  CM_InitY             = %13.7e\n", DynFri_CM_Init[1]    );
      Aux_Message( stdout, "  CM_InitZ             = %13.7e\n", DynFri_CM_Init[2]    );
      Aux_Message( stdout, "  DynFri_ProfM_MaxR    = %13.7e\n", DynFri_ProfM_MaxR    );
      Aux_Message( stdout, "  DynFri_ProfM_BinSize = %13.7e\n", DynFri_ProfM_BinSize );
      Aux_Message( stdout, "======================================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
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

   Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported !!\n" );

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass
// Description :  Record the center of mass (CM)
//
// Note        :  1. Invoked by Aux_Record_DynamicalFriction() recursively
//                2. Only include cells within CM_MaxR from CM_Old[] when updating CM
//
// Parameter   :  CM_Old[] : Previous CM
//                CM_New[] : New CM to be returned
//                CM_MaxR  : Maximum radius to compute CM
//
// Return      :  CM_New[]
//-------------------------------------------------------------------------------------------------------
void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR )
{

   const double CM_MaxR2          = SQR( CM_MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

   int   *PID0List = NULL;
   double M_ThisRank, MR_ThisRank[3], M_AllRank, MR_AllRank[3];
   real (*TotalDens)[PS1][PS1][PS1];

   M_ThisRank = 0.0;
   for (int d=0; d<3; d++)    MR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                    TimingSendPar_No );
#     endif

//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], 0, amr->NPatchComma[lv][1]/8, PID0List, _TOTAL_DENS,
                         OPT__RHO_INT_SCHEME, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif


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
// --> note that all ranks need to set CM_New[] since it is required for finding the next CM
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
// Function    :  Aux_Record_DynamicalFriction
// Description :  Record the center of mass
//
// Note        :  1. Linked to the function pointer "Aux_Record_User_Ptr"
//                2. Invoke GetCenterOfMass() recursively until the center of mass converges
//                3. Results are stored in the file "Record__CM"
//
// Parameter   :  None
//
// Return      :  DynFri_CM_Old[], DynFri_CM_New[]
//-------------------------------------------------------------------------------------------------------
void Aux_Record_DynamicalFriction()
{

// update the center of mass until convergence
   const double TolErrR2 = SQR( DynFri_CM_TolErrR );
   const int    NIterMax = 10;

   double dR2;
   int NIter = 0;

   memcpy( DynFri_CM_Old, DynFri_CM_New, sizeof(double)*3 );

   while ( true )
   {
      GetCenterOfMass( DynFri_CM_Old, DynFri_CM_New, DynFri_CM_MaxR );

      dR2 = SQR( DynFri_CM_Old[0] - DynFri_CM_New[0] )
          + SQR( DynFri_CM_Old[1] - DynFri_CM_New[1] )
          + SQR( DynFri_CM_Old[2] - DynFri_CM_New[2] );
      NIter ++;

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( DynFri_CM_Old, DynFri_CM_New, sizeof(double)*3 );
   }


// compute the mass profile
   Profile_t Prof;

   const bool LogBin_No          = false;
   const bool RemoveEmptyBin_Yes = true;

   Aux_ComputeProfile( &Prof, DynFri_CM_New, DynFri_ProfM_MaxR, DynFri_ProfM_BinSize,
                       LogBin_No, NULL_REAL, RemoveEmptyBin_Yes );

   if ( MPI_Rank == 0 )
   {
      Prof.Data[0] *= Prof.Weight[0];
      for (int b=1; b<Prof.NBin; b++)  Prof.Data[b] = Prof.Data[b-1] + Prof.Data[b]*Prof.Weight[b];

                     FILE *File = fopen( "Profile.txt", "w" );
                     fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
                     for (int b=0; b<Prof.NBin; b++)
                        fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
                                 Prof.Radius[b], Prof.Data[b], Prof.Weight[b], Prof.NCell[b] );
                     fclose( File );
   }


// record the center of mass and the mass profile
   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > TolErrR (%13.7e) !!\n", sqrt(dR2), DynFri_CM_TolErrR );

      const char FileName[] = "Record__CM";
      static bool FirstTime = true;

      FILE *File = NULL;

//    header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         File = fopen( FileName, "a" );
         fprintf( File, "#%12s   %9s   %5s   %13s   %13s   %13s\n", "Time", "Step", "NIter", "CM-x", "CM-y", "CM-z" );
         fclose( File );

         FirstTime = false;
      }

      File = fopen( FileName, "a" );
      fprintf( File, "%13.7e   %9ld   %5d   %13.7e   %13.7e   %13.7e\n",
               Time[0], Step, NIter, DynFri_CM_New[0], DynFri_CM_New[1], DynFri_CM_New[2] );
      fclose( File );
   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_Record_DynamicalFriction



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_DynamicalFriction
// Description :  Flag cells around the center of mass
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr"
//                2. Use the runtime table Input__Flag_User to specify the target refinement radius
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the
//                              file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_DynamicalFriction( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

   bool Flag = false;

// flag cells within a distance of "Threshold" from the center of mass
   const double dh   = amr->dh[lv];
   const double r[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                         amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                         amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh };

   double dr[3], dr2;

   for (int d=0; d<3; d++)    dr[d] = fabs( r[d] - DynFri_CM_New[d] );

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC )
   {
      for (int d=0; d<3; d++)
         if ( dr[d] > 0.5*amr->BoxSize[d] )  dr[d] -= amr->BoxSize[d];
   }

   dr2 = SQR( dr[0] ) + SQR( dr[1] ) + SQR( dr[2] );

   Flag = dr2 < SQR( Threshold );

   return Flag;

} // FUNCTION : Flag_DynamicalFriction
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_DynamicalFriction
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_DynamicalFriction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_DynamicalFriction;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = Aux_Record_DynamicalFriction;
   End_User_Ptr             = NULL;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_DynamicalFriction
