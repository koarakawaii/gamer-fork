#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
// parameters of the external density array
static char   Stream_ExtDens_Filename[MAX_STRING];    // filename of the external density array
static int    Stream_ExtDens_N;                       // number of cells along each direction of the external density array
static double Stream_ExtDens_L;                       // physical box size of the external density array
static int    Stream_ExtDens_IntOrder;                // interpolation order on the external density array

static real  *Stream_ExtDens = NULL;                  // external density array
static double Stream_ExtDens_dh;                      // cell size of the external density array

// parameters of the Plummer model
       int    Stream_Plummer_RSeed;                   // random seed for setting particle position and velocity
       double Stream_Plummer_M;                       // total mass
       double Stream_Plummer_R0;                      // scale radius
       double Stream_Plummer_MaxR;                    // maximum radius of particles
       double Stream_Plummer_Center[3];               // central coordinates
       double Stream_Plummer_BulkVel[3];              // bulk velocity
       int    Stream_Plummer_MassProfNBin;            // number of radial bins in the mass profile table

       double Stream_Plummer_Rho0;                    // peak density
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_Stream( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                 real *AllAttribute[PAR_NATT_TOTAL] );
#endif




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

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic --> reset NX0_TOT_* !!\n" );


// warnings
   /*
   if ( MPI_Rank == 0 )
   {
   } // if ( MPI_Rank == 0 )
   */


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
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",             &VARIABLE,                      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Stream_ExtDens_Filename",      Stream_ExtDens_Filename,       Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Stream_ExtDens_N",            &Stream_ExtDens_N,             -1,             1,                NoMax_int         );
   ReadPara->Add( "Stream_ExtDens_L",            &Stream_ExtDens_L,             -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Stream_ExtDens_IntOrder",     &Stream_ExtDens_IntOrder,       1,             0,                1                 );
   ReadPara->Add( "Stream_Plummer_RSeed",        &Stream_Plummer_RSeed,          123,           0,                NoMax_int         );
   ReadPara->Add( "Stream_Plummer_M",            &Stream_Plummer_M,             -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Stream_Plummer_R0",           &Stream_Plummer_R0,            -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Stream_Plummer_MaxR",         &Stream_Plummer_MaxR,          -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Stream_Plummer_CenterX",      &Stream_Plummer_Center[0],     -1.0,           0.0,              amr->BoxSize[0]   );
   ReadPara->Add( "Stream_Plummer_CenterY",      &Stream_Plummer_Center[1],     -1.0,           0.0,              amr->BoxSize[1]   );
   ReadPara->Add( "Stream_Plummer_CenterZ",      &Stream_Plummer_Center[2],     -1.0,           0.0,              amr->BoxSize[2]   );
   ReadPara->Add( "Stream_Plummer_BulkVelX",     &Stream_Plummer_BulkVel[0],     0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Stream_Plummer_BulkVelY",     &Stream_Plummer_BulkVel[1],     0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Stream_Plummer_BulkVelZ",     &Stream_Plummer_BulkVel[2],     0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Stream_Plummer_MassProfNBin", &Stream_Plummer_MassProfNBin,   1000,          2,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( Stream_ExtDens_L <= 0.0 )   Stream_ExtDens_L = amr->BoxSize[0];

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   Stream_ExtDens_dh   = Stream_ExtDens_L / Stream_ExtDens_N;
   Stream_Plummer_Rho0 = Stream_Plummer_M / ( 4.0/3.0*M_PI*CUBE(Stream_Plummer_R0) );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID = %d\n", TESTPROB_ID );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Parameters of the external density array:\n" );
      Aux_Message( stdout, "     Stream_ExtDens_Filename     = %s\n",     Stream_ExtDens_Filename );
      Aux_Message( stdout, "     Stream_ExtDens_N            = %d\n",     Stream_ExtDens_N        );
      Aux_Message( stdout, "     Stream_ExtDens_L            = %13.7e\n", Stream_ExtDens_L        );
      Aux_Message( stdout, "     Stream_ExtDens_IntOrder     = %d\n",     Stream_ExtDens_IntOrder );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Parameters of the Plummer model:\n" );
      Aux_Message( stdout, "     random seed                 = %d\n",     Stream_Plummer_RSeed         );
      Aux_Message( stdout, "     total mass                  = %14.7e\n", Stream_Plummer_M             );
      Aux_Message( stdout, "     peak density                = %14.7e\n", Stream_Plummer_Rho0          );
      Aux_Message( stdout, "     scale radius                = %14.7e\n", Stream_Plummer_R0            );
      Aux_Message( stdout, "     maximum radius of particles = %14.7e\n", Stream_Plummer_MaxR          );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "     central coordinate [%d]     = %14.7e\n", d, Stream_Plummer_Center[d]  );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "     bulk velocity [%d]          = %14.7e\n", d, Stream_Plummer_BulkVel[d] );
      Aux_Message( stdout, "     mass profile bins           = %d\n",     Stream_Plummer_MassProfNBin  );
      Aux_Message( stdout, "=============================================================================\n" );
   }


// (5) load the external density array
   if ( OPT__GRAVITY_EXTRA_MASS )
   {
//    check file existence
      if ( !Aux_CheckFileExist(Stream_ExtDens_Filename) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Stream_ExtDens_Filename );

//    check file size
      FILE *File = fopen( Stream_ExtDens_Filename, "rb" );
      fseek( File, 0, SEEK_END );
      const long N3         = CUBE( long(Stream_ExtDens_N) );
      const long ExpectSize = N3*sizeof(real);
      const long FileSize   = ftell( File );
      if ( FileSize != ExpectSize )
         Aux_Error( ERROR_INFO, "size of the file <%s> (%ld) != expect (%ld) !!\n", Stream_ExtDens_Filename, FileSize, ExpectSize );
      MPI_Barrier( MPI_COMM_WORLD );

//    allocate array
      Stream_ExtDens = new real [N3];

//    load data
      rewind( File );
      fread( Stream_ExtDens, sizeof(real), N3, File );

      fclose( File );
   } // if ( OPT__GRAVITY_EXTRA_MASS )


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

   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Stream
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Stream()
{

   delete [] Stream_ExtDens;

} // FUNCTION : End_Stream



//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_AddExtraMassForGravity_Stream
// Description :  Add extra mass source when computing gravity
//
// Note        :  1. Invoked by several functions using the function pointer "Poi_AddExtraMassForGravity_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Mass introduced here will only be used for computing gravity (i.e., when solving the Poisson eq.)
//                   --> It will NOT be used when solving other eqs. (e.g., hydro/MHD/Schroedinger)
//                   --> It will NOT be stored in the output data
//                   --> It WILL be included when computing the average density in Poi_GetAverageDensity(),
//                       which will then be used as DC in the periodic Poisson solver
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  Mass density at (x, y, z, t)
//-------------------------------------------------------------------------------------------------------
real Poi_AddExtraMassForGravity_Stream( const double x, const double y, const double z, const double Time,
                                        const int lv, double AuxArray[] )
{

// return zero if (x,y,z) lies outside the external density array
   const double EdgeL = 0.5*Stream_ExtDens_dh;
   const double EdgeR = Stream_ExtDens_L - 0.5*Stream_ExtDens_dh;

   if ( x < EdgeL  ||  y < EdgeL  ||  z < EdgeL  ||  x > EdgeR  ||  y > EdgeR  ||  z > EdgeR )  return (real)0.0;


   const double _dh = 1.0/Stream_ExtDens_dh;
   const long   N   = Stream_ExtDens_N;
   real ExtDens = 0.0;

   switch ( Stream_ExtDens_IntOrder )
   {
      case 0:
      {
         const long i   = x*_dh;
         const long j   = y*_dh;
         const long k   = z*_dh;
         const long idx = IDX321( i, j, k, N, N );

#        ifdef GAMER_DEBUG
         if ( i < 0  ||  i >= N )   Aux_Error( ERROR_INFO, "incorrect i = %ld (x %14.7e, N %ld) !!\n", i, x, N );
         if ( j < 0  ||  j >= N )   Aux_Error( ERROR_INFO, "incorrect j = %ld (y %14.7e, N %ld) !!\n", j, y, N );
         if ( k < 0  ||  k >= N )   Aux_Error( ERROR_INFO, "incorrect k = %ld (z %14.7e, N %ld) !!\n", k, z, N );
#        endif

         ExtDens = Stream_ExtDens[idx];
      } // case 0
      break;


      case 1:
      {
         int    idxLR[2][3];  // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
         double dr      [3];  // distance to the center of the left cell
         double Frac [2][3];  // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells
         double xyz     [3] = { x, y, z };

         for (int d=0; d<3; d++)
         {
//          calculate the array index of the left and right cells
            dr      [d] = xyz[d]*_dh - 0.5;
            idxLR[0][d] = int( dr[d] );
            idxLR[1][d] = idxLR[0][d] + 1;

            if ( idxLR[0][d] < 0 )
            {
               idxLR[0][d] = 0;
               idxLR[1][d] = 1;
            }

            else if ( idxLR[1][d] >= N )
            {
               idxLR[0][d] = N-2;
               idxLR[1][d] = N-1;
            }

//          get the weighting of the nearby 8 cells
            dr     [d] -= (double)idxLR[0][d];
            Frac[0][d]  = 1.0 - dr[d];
            Frac[1][d]  =       dr[d];
         }

//       calculate acceleration
         for (int k=0; k<2; k++)
         for (int j=0; j<2; j++)
         for (int i=0; i<2; i++)
         {
            const long idx = IDX321( idxLR[i][0], idxLR[j][1], idxLR[k][2], N, N );

#           ifdef GAMER_DEBUG
            if ( idxLR[i][0] < 0  ||  idxLR[i][0] >= N )   Aux_Error( ERROR_INFO, "incorrect i = %ld (LR %d, x %14.7e, N %ld) !!\n", idxLR[i][0], i, x, N );
            if ( idxLR[j][1] < 0  ||  idxLR[j][1] >= N )   Aux_Error( ERROR_INFO, "incorrect j = %ld (LR %d, y %14.7e, N %ld) !!\n", idxLR[j][1], j, y, N );
            if ( idxLR[k][2] < 0  ||  idxLR[k][2] >= N )   Aux_Error( ERROR_INFO, "incorrect k = %ld (LR %d, z %14.7e, N %ld) !!\n", idxLR[k][2], k, z, N );
#           endif

            ExtDens += Stream_ExtDens[idx]*Frac[i][0]*Frac[j][1]*Frac[k][2];
         }
      } // case 1
      break;


      default:
         Aux_Error( ERROR_INFO, "unsupported order of interpolation (%d) !!\n", Stream_ExtDens_IntOrder );
   } // switch ( Stream_ExtDens_IntOrder )


   return ExtDens;

} // FUNCTION : Poi_AddExtraMassForGravity_Stream



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Stream
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_Stream( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_Stream
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Stream
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Stream()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr         = SetGridIC;
   Init_Field_User_Ptr            = NULL;
   Flag_User_Ptr                  = NULL;
   Mis_GetTimeStep_User_Ptr       = NULL;
   BC_User_Ptr                    = BC_Stream;
   Flu_ResetByUser_Func_Ptr       = NULL;
   Output_User_Ptr                = NULL;
   Aux_Record_User_Ptr            = NULL;
   Init_User_Ptr                  = NULL;
   End_User_Ptr                   = End_Stream;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr           = NULL;
   Init_ExternalPot_Ptr           = NULL;
   Poi_AddExtraMassForGravity_Ptr = Poi_AddExtraMassForGravity_Stream;
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr        = Par_Init_ByFunction_Stream;
   Par_Init_Attribute_User_Ptr    = NULL;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Stream
