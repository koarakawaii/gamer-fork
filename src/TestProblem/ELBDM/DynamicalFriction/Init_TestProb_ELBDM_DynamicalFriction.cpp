#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
double DynFri_CM_MaxR;           // maximum radius for determining CM
double DynFri_CM_TolErrR;        // maximum allowed errors for determining CM
double DynFri_CM_Init[3];        // initial guess of CM

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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,                 DEFAULT,             MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "DynFri_CM_InitX",   &DynFri_CM_Init[0],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_InitY",   &DynFri_CM_Init[1],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_InitZ",   &DynFri_CM_Init[2],        -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_MaxR",    &DynFri_CM_MaxR,           -1.0,                NoMin_double,     NoMax_double      );
   ReadPara->Add( "DynFri_CM_TolErrR", &DynFri_CM_TolErrR,        -1.0,                NoMin_double,     NoMax_double      );
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

   if ( DynFri_CM_MaxR    < 0 )  DynFri_CM_MaxR    = 10.0*amr->dh[MAX_LEVEL];
   if ( DynFri_CM_TolErrR < 0 )  DynFri_CM_TolErrR =  1.0*amr->dh[MAX_LEVEL]; 

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )   Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported !!\n" );


// (2) set the problem-specific derived parameters
   for (int d=0; d<3; d++)    DynFri_CM_Old[d] = DynFri_CM_Init[d];


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID = %d\n",     TESTPROB_ID       );
      Aux_Message( stdout, "  CM_MaxR         = %13.7e\n", DynFri_CM_MaxR    );
      Aux_Message( stdout, "  CM_TolErrR      = %13.7e\n", DynFri_CM_TolErrR );
      Aux_Message( stdout, "  CM_InitX        = %13.7e\n", DynFri_CM_Init[0] );
      Aux_Message( stdout, "  CM_InitY        = %13.7e\n", DynFri_CM_Init[1] );
      Aux_Message( stdout, "  CM_InitZ        = %13.7e\n", DynFri_CM_Init[2] );
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


// set the periodic CM_Old[]
   double CM_Old_Periodic[3];

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC )
      for (int d=0; d<3; d++)
         CM_Old_Periodic[d] = ( CM_Old[d] > amr->BoxCenter[d] ) ? CM_Old[d]-2.0*amr->BoxSize[d] :
                                                                  CM_Old[d]+2.0*amr->BoxSize[d];
   else
      for (int d=0; d<3; d++)
         CM_Old_Periodic[d] = CM_Old[d];


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

         for (int k=0; k<PS1; k++)  {  const double zz = z0 + k*dh;
                                       const double z1 = zz - CM_Old         [2];
                                       const double z2 = zz - CM_Old_Periodic[2];
                                       const double z  = ( fabs(z1) <= fabs(z2) ) ? z1 : z2;
         for (int j=0; j<PS1; j++)  {  const double yy = y0 + j*dh;
                                       const double y1 = yy - CM_Old         [1];
                                       const double y2 = yy - CM_Old_Periodic[1];
                                       const double y  = ( fabs(y1) <= fabs(y2) ) ? y1 : y2;
         for (int i=0; i<PS1; i++)  {  const double xx = x0 + i*dh;
                                       const double x1 = xx - CM_Old         [0];
                                       const double x2 = xx - CM_Old_Periodic[0];
                                       const double x  = ( fabs(x1) <= fabs(x2) ) ? x1 : x2;

//          only include cells within CM_MaxR
            const double R2 = x*x + y*y + z*z;
            if ( R2 < CM_MaxR2 )
            {
               const double dm = TotalDens[PID][k][j][i]*dv;

               M_ThisRank     += dm;
               MR_ThisRank[0] += dm*xx;
               MR_ThisRank[1] += dm*yy;
               MR_ThisRank[2] += dm*zz;
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

   while ( true )
   {
      GetCenterOfMass( DynFri_CM_Old, DynFri_CM_New, DynFri_CM_MaxR );

      dR2 = SQR( DynFri_CM_Old[0] - DynFri_CM_New[0] )
          + SQR( DynFri_CM_Old[1] - DynFri_CM_New[1] )
          + SQR( DynFri_CM_Old[2] - DynFri_CM_New[2] );
      NIter ++;

      memcpy( DynFri_CM_Old, DynFri_CM_New, sizeof(double)*3 );

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )   break;
   }


// record the center of mass
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
   Flag_User_Ptr            = NULL;
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
