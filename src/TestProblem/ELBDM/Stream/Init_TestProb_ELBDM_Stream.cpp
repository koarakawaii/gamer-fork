#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static char  Stream_ExtDens_Filename[MAX_STRING];  // filename of the external density array
static int   Stream_ExtDens_N;                     // Grid size of the external density array

static real *Stream_ExtDens = NULL;                // external density array
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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif


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
// ReadPara->Add( "KEY_IN_THE_FILE",            &VARIABLE,                 DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Stream_ExtDens_Filename",    Stream_ExtDens_Filename,   Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Stream_ExtDens_N",          &Stream_ExtDens_N,         -1,             1,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


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
      Aux_Message( stdout, "  test problem ID         = %d\n", TESTPROB_ID             );
      Aux_Message( stdout, "  Stream_ExtDens_Filename = %s\n", Stream_ExtDens_Filename );
      Aux_Message( stdout, "  Stream_ExtDens_N        = %d\n", Stream_ExtDens_N        );
      Aux_Message( stdout, "=============================================================================\n" );
   }


// (5) load the external density array
// check file existence
   if ( !Aux_CheckFileExist(Stream_ExtDens_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Stream_ExtDens_Filename );

// check file size
   FILE *File = fopen( Stream_ExtDens_Filename, "rb" );
   fseek( File, 0, SEEK_END );
   const long N3         = CUBE( long(Stream_ExtDens_N) );
   const long ExpectSize = N3*sizeof(real);
   const long FileSize   = ftell( File );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> (%ld) != expect (%ld) !!\n", Stream_ExtDens_Filename, FileSize, ExpectSize );
   MPI_Barrier( MPI_COMM_WORLD );

// allocate array
   Stream_ExtDens = new real [N3];

// load data
   rewind( File );
   fread( Stream_ExtDens, sizeof(real), N3, File );

   fclose( File );


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

   const double dh  = amr->dh[0];    // replace it by an input parameter later
   const long   N   = Stream_ExtDens_N;
   const long   i   = x/dh;
   const long   j   = y/dh;
   const long   k   = z/dh;
   const long   idx = IDX321( i, j, k, N, N );

#  ifdef GAMER_DEBUG
   if ( i < 0  ||  i >= N )   Aux_Error( ERROR_INFO, "incorrect i = %ld (x %14.7e, N %ld) !!\n", i, x, N );
   if ( j < 0  ||  j >= N )   Aux_Error( ERROR_INFO, "incorrect j = %ld (y %14.7e, N %ld) !!\n", j, y, N );
   if ( k < 0  ||  k >= N )   Aux_Error( ERROR_INFO, "incorrect k = %ld (z %14.7e, N %ld) !!\n", k, z, N );
#  endif

   return Stream_ExtDens[idx];

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
   Par_Init_ByFunction_Ptr        = NULL;
   Par_Init_Attribute_User_Ptr    = NULL;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Stream
