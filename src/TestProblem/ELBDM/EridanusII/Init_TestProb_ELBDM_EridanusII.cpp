#include "GAMER.h"
#include "TestProb.h"



// soliton-specific global variables
// =======================================================================================
static double   Soliton_CoreRadius;                      // soliton core radius
static int      Soliton_InputMode;                       // soliton input mode: 1/2 -> table/approximate analytical form
static double   Soliton_OuterSlope;                      // soliton outer slope (only used by Soliton_InputMode=2)
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the reference soliton density profile

static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static double  *Soliton_DensProf   = NULL;               // soliton density profile [radius/density]
static double   Soliton_ScaleL     = NULL;               // L/D: length/density scale factors of each soliton
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference soliton profiles)
static double   Soliton_ScaleD     = NULL;
// =======================================================================================


// particle-specific global variables
// =======================================================================================
       int    Star_RSeed;           // random seed for setting particle position and velocity
       double Star_Rho0;            // peak density
       double Star_R0;              // scale radius
       double Star_MaxR;            // maximum radius for particles
       int    Star_MassProfNBin;    // number of radial bins in the mass profile table

static double Star_FreeT;           // free-fall time at Star_R0
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_EridanusII( const long NPar_ThisRank, const long NPar_AllRank,
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

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
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
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Soliton_CoreRadius",        &Soliton_CoreRadius,        -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Soliton_InputMode",         &Soliton_InputMode,          1,             1,                2                 );
   ReadPara->Add( "Soliton_OuterSlope",        &Soliton_OuterSlope,        -8.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Star_RSeed",                &Star_RSeed,                 123,           0,                NoMax_int         );
   ReadPara->Add( "Star_Rho0",                 &Star_Rho0,                 -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_R0",                   &Star_R0,                   -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_MaxR",                 &Star_MaxR,                 -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_MassProfNBin",         &Star_MassProfNBin,          1000,          2,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
#  ifdef GRAVITY
   Star_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Star_Rho0) );
#  endif


// (3) load the reference soliton density profile and evaluate the scale factors
   if ( OPT__INIT != INIT_BY_RESTART  &&  Soliton_InputMode == 1 )
   {
//    load the reference profile
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for Soliton_DensProf
      const int  NCol         = 2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

      Soliton_DensProf_NBin = Aux_LoadTable( Soliton_DensProf, Soliton_DensProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );


//    get the core radius of the reference profile
      const double *RadiusRef = Soliton_DensProf + 0*Soliton_DensProf_NBin;
      const double *DensRef   = Soliton_DensProf + 1*Soliton_DensProf_NBin;
      const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

      double CoreRadiusRef = NULL_REAL;

      for (int b=1; b<Soliton_DensProf_NBin-1; b++)
      {
         if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
         {
            CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
            break;
         }
      }

      if ( CoreRadiusRef == NULL_REAL )
         Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );


//    evaluate the scale factors of each soliton
      Soliton_ScaleL = Soliton_CoreRadius / CoreRadiusRef;
      Soliton_ScaleD = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(Soliton_ScaleL) );
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 12.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (5) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID                );
      Aux_Message( stdout, "  soliton core radius                       = %13.6e\n", Soliton_CoreRadius         );
      Aux_Message( stdout, "  soliton input mode                        = %d\n",     Soliton_InputMode          );
      if      ( Soliton_InputMode == 2 )
      Aux_Message( stdout, "  soliton outer slope                       = %13.6e\n", Soliton_OuterSlope         );
      else if ( Soliton_InputMode == 1 ) {
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Soliton_DensProf_Filename  );
      Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Soliton_DensProf_NBin      ); }
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  star cluster properties:\n" );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     Star_RSeed );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", Star_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", Star_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", Star_MaxR );
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     Star_MassProfNBin );
      Aux_Message( stdout, "  free-fall time at the scale radius        = %13.7e\n", Star_FreeT );
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

   const double Soliton_Center[3] = { amr->BoxCenter[0],
                                      amr->BoxCenter[1],
                                      amr->BoxCenter[2] };
   const double r_tar             = sqrt( SQR(x-Soliton_Center[0]) +
                                          SQR(y-Soliton_Center[1]) +
                                          SQR(z-Soliton_Center[2]) );

   if ( Soliton_InputMode == 1 )
   {
      const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
      const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

      double r_ref, dens_ref;

//    rescale radius (target radius --> reference radius)
      r_ref = r_tar / Soliton_ScaleL;

//    linear interpolation
      dens_ref = Mis_InterpolateFromTable( Soliton_DensProf_NBin, Table_Radius, Table_Density, r_ref );

      if ( dens_ref == NULL_REAL )
      {
         if      ( r_ref <  Table_Radius[0] )
            dens_ref = Table_Density[0];

         else if ( r_ref >= Table_Radius[Soliton_DensProf_NBin-1] )
            dens_ref = Table_Density[Soliton_DensProf_NBin-1];

         else
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                       r_ref, Table_Radius[0], Table_Radius[Soliton_DensProf_NBin-1] );
      }

//    rescale density (reference density --> target density) and add to the fluid array
      fluid[DENS] = dens_ref*Soliton_ScaleD;
   }

   else if ( Soliton_InputMode == 2 )
   {
      const double m22    = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
      const double rc_kpc = Soliton_CoreRadius*UNIT_L/Const_kpc;

      fluid[DENS] = 1.945e7/SQR( m22*rc_kpc*rc_kpc )*pow( 1.0+9.06e-2*SQR(r_tar/rc_kpc), Soliton_OuterSlope );
   }

   else
      Aux_Error( ERROR_INFO, "Unsupported Soliton_InputMode (%d) !!\n", Soliton_InputMode );


// set the real and imaginary parts
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_EridanusII
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_EridanusII()
{

   delete [] Soliton_DensProf;

} // FUNCTION : End_EridanusII



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
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
void BC_EridanusII( real fluid[], const double x, const double y, const double z, const double Time,
                    const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_EridanusII
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_EridanusII()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = BC_EridanusII;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = End_EridanusII;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr  = Par_Init_ByFunction_EridanusII;
#  endif
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_EridanusII
