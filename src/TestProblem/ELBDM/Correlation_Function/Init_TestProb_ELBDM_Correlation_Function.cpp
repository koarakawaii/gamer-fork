#include "GAMER.h"
#include "TestProb.h"

// extern functions

// This function computes desnity profil, with standare deviation 
void Aux_ComputeProfile_with_Sigma( Profile_with_Sigma_t *Prof[], const double Center[], const double r_max_input, const double dr_min,
                                    const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const long TVarBitIdx[],
                                    const int NProf, const int MinLv, const int MaxLv, const PatchType_t PatchType,
                                    const double PrepTime );
// This function computes correlation function
//void Aux_ComputeCorrelation( Profile_t *Correlation[], FieldIdx_t *Passive_idx[], const Profile_t *prof_init[], const double Center[],
void Aux_ComputeCorrelation( Profile_t *Correlation[], const Profile_with_Sigma_t *prof_init[], const double Center[],
                             const double r_max_input, const double dr_min, const bool LogBin, const double LogBinRatio,
                             const bool RemoveEmpty, const long TVarBitIdx[], const int NProf, const int MinLv, const int MaxLv,
                             const PatchType_t PatchType, const double PrepTime, const double dr_min_prof );
//
// intern functions
static void Record_CenterOfMass( bool record_flag );

// problem-specific global variables
// =======================================================================================
static FieldIdx_t Idx_Dens0 = Idx_Undefined;  // field index for storing the **initial** density
static double   System_CM_MaxR;               // maximum radius for determining System CM
static double   System_CM_TolErrR;            // maximum allowed errors for determining System CM
static double   Soliton_CM_MaxR;              // maximum radius for determining Soliton CM
static double   Soliton_CM_TolErrR;           // maximum allowed errors for determining Soliton CM
static double   Center[3];                    // use CoM coordinate of the whole halo as center
static double   dr_min_prof;                  // bin size of correlation function statistics (minimum size if logarithic bin) (profile)
static double   LogBinRatio_prof;             // ratio of bin size growing rate for logarithmic bin (profile)
static double   RadiusMax_prof;               // maximum radius for correlation function statistics (profile)
static double   dr_min_corr;                  // bin size of correlation function statistics (minimum size if logarithic bin) (correlation)
static double   LogBinRatio_corr;             // ratio of bin size growing rate for logarithmic bin (correlation)
static double   RadiusMax_corr;               // maximum radius for correlation function statistics (correlation)
//static double   PrepTime;                     // time for doing statistics
static bool     ComputeCorrelation;           // flag for compute correlation
static bool     ReComputeCorrelation;         // flag for recompute correlation for restart; use the simulation time of RESTART as initial time for computing time correlation; only available for RESTART
static bool     LogBin_prof;                  // logarithmic bin or not (profile)
static bool     RemoveEmpty_prof;             // remove 0 sample bins; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (profile)
static bool     LogBin_corr;                  // logarithmic bin or not (correlation)
static bool     RemoveEmpty_corr;             // remove 0 sample bins; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (correlation)
static int      MinLv;                        // do statistics from MinLv to MaxLv
static int      MaxLv;                        // do statistics from MinLv to MaxLv
static int      OutputCorrelationMode;        // output correlation function mode=> 0: constant interval 1: by table
static int      StepInitial;                  // inital step for recording correlation function (OutputCorrelationMode = 0) 
static int      StepInterval;                 // interval for recording correlation function (OutputCorrelationMode = 0)
static int      *StepTable;                   // step index table for output correlation function (OutputCorrelationMode = 1)
static bool     Fluid_Periodic_BC_Flag;       // flag for checking the fluid boundary condtion is setup to periodic (0: user defined; 1: periodic)
static char     FilePath_corr[MAX_STRING];    // output path for correlation function text files

static int step_counter;                             // counter for caching consumed step indices
static Profile_with_Sigma_t Prof_Dens_initial;                      // pointer to save initial density profile
static Profile_with_Sigma_t *Prof[] = { &Prof_Dens_initial };
static Profile_t            Correlation_Dens;                       // pointer to save density correlation function
static Profile_t            *Correlation[] = { &Correlation_Dens };       
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

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
   #  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
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
   const char FileName[] = "Input__TestProb_ELBDM_Correlation_Function";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",          &VARIABLE,               DEFAULT,          MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara  = new ReadPara_t;
   ReadPara->Add( "dr_min_corr",              &dr_min_corr,            Eps_double,       Eps_double,       NoMax_double      );
   ReadPara->Add( "LogBinRatio_corr",         &LogBinRatio_corr,          1.0,           Eps_double,      NoMax_double       );
   ReadPara->Add( "RadiusMax_corr",           &RadiusMax_corr,         Eps_double,       Eps_double,      NoMax_double       );
   ReadPara->Add( "LogBin_corr",              &LogBin_corr,             false,          Useless_bool,     Useless_bool       );
   ReadPara->Add( "RemoveEmpty_corr",         &RemoveEmpty_corr,        false,          Useless_bool,     Useless_bool       );
   ReadPara->Add( "dr_min_prof",              &dr_min_prof,            Eps_double,       Eps_double,       NoMax_double      );
   ReadPara->Add( "MinLv",                    &MinLv,                       0,                     0,        MAX_LEVEL       );
   ReadPara->Add( "MaxLv",                    &MaxLv,               MAX_LEVEL,                     0,        MAX_LEVEL       );
   ReadPara->Add( "OutputCorrelationMode",    &OutputCorrelationMode,       0,                     0,             1          );
   ReadPara->Add( "StepInitial",              &StepInitial,                 0,                     0,       NoMax_int        );
   ReadPara->Add( "StepInterval",             &StepInterval,                1,                     1,        NoMax_int       );
   ReadPara->Add( "FilePath_corr",            FilePath_corr,       Useless_str,           Useless_str,      Useless_str      );
   if ( OPT__INIT == INIT_BY_RESTART )
      ReadPara->Add( "ReComputeCorrelation",     &ReComputeCorrelation,     false,         Useless_bool,      Useless_bool      );
   ReadPara->Read( FileName );
   delete ReadPara;


// (1-2) set the default values

   if ( dr_min_corr <=Eps_double )          dr_min_corr = 1e-3*0.5*amr->BoxSize[0];
   if ( RadiusMax_corr<=Eps_double )        RadiusMax_corr = 0.5*amr->BoxSize[0];
   if ( LogBinRatio_corr<=1.0 )             LogBinRatio_corr = 2.0;

   if ( dr_min_prof <=Eps_double )          dr_min_prof = dr_min_corr;
   RadiusMax_prof                           = RadiusMax_corr * 1.05;   // assigned by Test Problem
   LogBinRatio_prof                         = 1.0;                     // hard-coded by Test Problem (no effect)
   LogBin_prof                              = false;                   // hard-coded by Test Problem
   RemoveEmpty_prof                         = false;                   // hard-coded by Test Problem

   if ( MinLv < 0 ) MinLv = 0;
   if ( MaxLv <= MinLv ) MaxLv = MAX_LEVEL;
   if ( FilePath_corr == "\0" )  sprintf( FilePath_corr, "./" );
   else
   { 
      FILE *file_checker = fopen(FilePath_corr, "r");
      if (!file_checker)
         Aux_Error( ERROR_INFO, "File path %s for saving correlation function text files does not exist!! Please create!!\n", FilePath_corr );
      else
         fclose(file_checker);
   }
   

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );

// (2) set the problem-specific derived parameters


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
      Aux_Message( stdout, "=================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                             = %d\n",     TESTPROB_ID               );
      Aux_Message( stdout, "  histogram bin size  (correlation)           = %13.6e\n", dr_min_corr            );
      Aux_Message( stdout, "  log bin ratio       (correlation)           = %13.6e\n", LogBinRatio_corr       );
      Aux_Message( stdout, "  radius maximum      (correlation)           = %13.6e\n", RadiusMax_corr         );
      Aux_Message( stdout, "  use logarithmic bin (correlation)           = %d\n"    , LogBin_corr            );
      Aux_Message( stdout, "  remove empty bin    (correlation)           = %d\n"    , RemoveEmpty_corr       );
      Aux_Message( stdout, "  histogram bin size  (profile)               = %13.6e\n", dr_min_prof            );
      Aux_Message( stdout, "  log bin ratio       (profile, no effect)    = %13.6e\n", LogBinRatio_prof       );
      Aux_Message( stdout, "  radius maximum      (profile, assigned)     = %13.6e\n", RadiusMax_prof         );
      Aux_Message( stdout, "  use logarithmic bin (profile, assigned)     = %d\n"    , LogBin_prof            );
      Aux_Message( stdout, "  remove empty bin    (profile, assigned)     = %d\n"    , RemoveEmpty_prof       );
//      Aux_Message( stdout, "  prepare time                                = %13.6e\n", PrepTime               );
      Aux_Message( stdout, "  minimum level                               = %d\n"    , MinLv                  );
      Aux_Message( stdout, "  maximum level                               = %d\n"    , MaxLv                  );
      Aux_Message( stdout, "  output correlation function mode            = %d\n"    , OutputCorrelationMode  );
      Aux_Message( stdout, "  file path for correlation text file         = %s\n"    , FilePath_corr          );
      if ( OPT__INIT == INIT_BY_RESTART )
         Aux_Message( stdout, "  re-compute correlation using restart time as initial time = %d\n", ReComputeCorrelation );
      Aux_Message( stdout, "=================================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_StepTable
// Description :  Load the dump table from the file "Input__StepTable"
//-------------------------------------------------------------------------------------------------------
static void Init_Load_StepTable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ...\n" );


   const char FileName[] = "Input__StepTable";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   const int MaxLine = 10000;
   char *input_line = NULL;
   size_t len = 0;
   int Trash, line, n;


// allocate the step table
   StepTable = new int [MaxLine];


// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (line=0; line<MaxLine; line++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
         Aux_Error( ERROR_INFO, "incorrect reading at line %d of the file <%s> !!\n", line+2, FileName );

      sscanf( input_line, "%d%d", &Trash, &StepTable[line] );

//    stop the reading
      if ( input_line[0] == 42 )                   // '*' == 42
      {

//       ensure that at least one step index is loaded
         if ( line == 0 )
            Aux_Error( ERROR_INFO, "please provide at least one step index in the step table !!\n" );

         int StepTable_NDump   = line;             // record the number of step indices
         StepTable[line]   = __INT_MAX__;          // set the next step as an extremely large number

         if ( StepTable[line-1] < END_STEP )
         {
            END_STEP          = StepTable[line-1];    // reset the ending time as the time of the last step

            if ( MPI_Rank == 0 )
               Aux_Message( stdout, "NOTE : the END_STEP is reset to the time of the last step index = %de\n",
                            END_STEP );
         }


//       verify the loaded dump table
         for (int t=1; t<=line; t++)
         {
            if ( StepTable[t] < StepTable[t-1] )
               Aux_Error( ERROR_INFO, "values recorded in \"%s\" must be monotonically increasing !!\n",
                          FileName );
         }

         break;

      } // if ( input_line[0] == 42 )
   } // for (line=0; line<MaxLine; line++)


   if ( line == MaxLine )
      Aux_Error( ERROR_INFO, "please prepare a symbol * in the end of the file <%s> !!\n", FileName );


   fclose( File );

   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ... done\n" );

} // FUNCTION : Init_Load_StepTable


//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_ELBDM_Correlation_Function
// Description :  Add the problem-specific fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void AddNewField_ELBDM_Correlation_Function(void)
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   Idx_Dens0 = AddField( "Dens0", NORMALIZE_NO, INTERP_FRAC_NO );
//   if ( MPI_Rank == 0 )   printf("Idx_Dens0 = %d \n", Idx_Dens0);
#  endif

} // FUNCTION : AddNewField_ELBDM_Correlation_Function


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_ELBDM_Correlation_Function
// Description :  Store the initial density
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Init_User_ELBDM_Correlation_Function(void)
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
//    store the initial density in both Sg so that we don't have to worry about which Sg to be used
//    a. for restart and ReComputeCorrelation disabled, the initial density has already been loaded and we just need to copy the data to another Sg
      if ( ( OPT__INIT == INIT_BY_RESTART ) && ( !ReComputeCorrelation ) ) {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i];

         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }

//    b. for starting a new simulation, we must copy the initial density to both Sg
      else {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

         amr->patch[   amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }
   }

   step_counter = 0;  // initialize global variable step counter
   const double InitialTime = Time[0];
   if ( OutputCorrelationMode == 0)
   {
      if ( MPI_Rank==0 ) Aux_Message( stdout, "StepInitial = %d ; StepInterval = %d \n", StepInitial, StepInterval);
   }
   else if ( OutputCorrelationMode == 1)
      Init_Load_StepTable();
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "InitialTime = %13.6e \n", InitialTime );

   // compute the enter position for passive field
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate halo center for passive field:\n");

   Center = ...;
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "Center of passive field is ( %14.11e,%14.11e,%14.11e )\n", Center[0], Center[1], Center[2] );
   // commpute density profile for passive field;
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate density profile for passive field:\n");

   const long TVar[] = {BIDX(Idx_Dens0)};
   Aux_ComputeProfile_with_Sigma( Prof, Center, RadiusMax_prof, dr_min_prof, LogBin_prof, LogBinRatio_prof, RemoveEmpty_prof, TVar, 1, MinLv, MaxLv, PATCH_LEAF, InitialTime );

   if ( MPI_Rank == 0 )
   {
      char Filename[MAX_STRING];
      sprintf( Filename, "%s/initial_profile_with_Sigma.txt", FilePath_corr );
      FILE *output_initial_prof = fopen(Filename, "w");
      fprintf( output_initial_prof, "#%19s  %21s  %21s  %21s  %11s\n", "Radius", "Dens", "Dens_Sigma" , "Weighting", "Cell_Number");
      for (int b=0; b<Prof[0]->NBin; b++)
         fprintf( output_initial_prof, "%20.14e  %21.14e  %21.14e  %21.14e  %11ld\n",
                   Prof[0]->Radius[b], Prof[0]->Data[b], Prof[0]->Data_Sigma[b], Prof[0]->Weight[b], Prof[0]->NCell[b] );
      fclose(output_initial_prof);
   }
#  endif

} // FUNCTION : Init_User_ELBDM_Correlation_Function

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
//void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
//                const int lv, double AuxArray[] )
//{
//
//// HYDRO example
//   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
//
//   Dens = 1.0;
//   MomX = 0.0;
//   MomY = 0.0;
//   MomZ = 0.0;
//   Pres = 2.0;
//   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
//                                    EoS_AuxArray_Int, h_EoS_Table, NULL ); // assuming EoS requires no passive scalars
//   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );         // do NOT include magnetic energy here
//
//// set the output array
//   fluid[DENS] = Dens;
//   fluid[MOMX] = MomX;
//   fluid[MOMY] = MomY;
//   fluid[MOMZ] = MomZ;
//   fluid[ENGY] = Etot;
//
//} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_HALO
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
static void BC_HALO( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_HALO



//-------------------------------------------------------------------------------------------------------
// Function    :  Do_CF
// Description :  Do correlation function calculation
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
static void Do_CF( void )
{
   if ( ((OutputCorrelationMode==1) && (Step==StepTable[step_counter])) || ((OutputCorrelationMode==0) && (Step>=StepInitial) && (((Step-StepInitial)%StepInterval)==0)) )
   {
      const long TVar[] = {_DENS};
      Aux_ComputeCorrelation( Correlation, (const Profile_with_Sigma_t**)Prof, Center, RadiusMax_corr, dr_min_corr, LogBin_corr, LogBinRatio_corr,
                              RemoveEmpty_corr, TVar, 1, MinLv, MaxLv, PATCH_LEAF, Time[0], dr_min_prof);

      if ( MPI_Rank == 0 )
      {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s/correlation_function_t=%.4e.txt", FilePath_corr, Time[0] );
         FILE *output_correlation = fopen(Filename, "w");
         fprintf( output_correlation, "#%19s  %21s  %21s  %11s\n", "Radius", "Correlation_Function", "Weighting", "Cell_Number");
         for (int b=0; b<Correlation[0]->NBin; b++)
             fprintf( output_correlation, "%20.14e  %21.14e  %21.14e  %11ld\n",
                      Correlation[0]->Radius[b], Correlation[0]->Data[b], Correlation[0]->Weight[b], Correlation[0]->NCell[b] );
         fclose(output_correlation);
      }
      // accumulate the step counter
      step_counter ++;
   }
}



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Correlation_Function
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
static void End_Correlation_Function()
{
   
} // FUNCTION : End_Correlation_Function
#endif // end of if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Correlation_Function
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Correlation_Function()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


//   Init_Function_User_Ptr = SetGridIC;
   Init_Field_User_Ptr    = AddNewField_ELBDM_Correlation_Function;
   Init_User_Ptr          = Init_User_ELBDM_Correlation_Function;
   Aux_Record_User_Ptr    = Do_COM_and_CF;
   BC_User_Ptr            = BC_HALO;
   End_User_Ptr           = End_Correlation_Function;
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )

// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
   Src_Init_User_Ptr      = NULL; // option: SRC_USER;                example: SourceTerms/User_Template/CPU_Src_User_Template.cpp

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Correlation_Function
