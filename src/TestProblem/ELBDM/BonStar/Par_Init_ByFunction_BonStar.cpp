#include "GAMER.h"
#ifdef SUPPORT_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#else
#error : ERROR : please turn on SUPPORT_GSL for the cluster merger test problem !!
#endif

#ifdef PARTICLE

extern int    BonStar_RSeed;
extern double BonStar_MaxR;
extern double BonStar_Rzero;
extern int    BonStar_ProfNBin;

static RandomNumber_t *RNG = NULL;


static double DensProf( const double r );
static double MassProf( const double r );
static void GetSigmaProf( const int NBin, const double *r, double *Sigma );
static void   RanVec_FixRadius( const double r, double RanVec[] );
static double GSL_IntFunc_SigmaProf( double r, void *IntPara );


// density profile parameters
#define ANAL_FORM 2

// analytical form 1
#if ( ANAL_FORM == 1 )
const double r_core = 5.0;
const double a0     = +1.99761e-04;
const double a1     = +1.71915e-07;
const double a2     = -1.81062e-08;
const double a3     = +9.13137e-11;
const double a4     = -1.33860e-13;

// analytical form 2
#elif ( ANAL_FORM == 2 )
const double rho0   = 2.0e-4;
const double r_core = 2.5e2;
const double pow1   = 2.0;
const double pow2   = 4.0;

#else
#  error: ERROR : unsupported ANAL_FORM !!
#endif


#define POW5( a )       ( (a)*(a)*(a)*(a)*(a) )
#define POW6( a )       ( (a)*(a)*(a)*(a)*(a)*(a) )
#define POW7( a )       ( (a)*(a)*(a)*(a)*(a)*(a)*(a) )
#define POW8( a )       ( (a)*(a)*(a)*(a)*(a)*(a)*(a)*(a) )

// GSL parameters
#ifdef SUPPORT_GSL
const int    GSL_IntRule  = GSL_INTEG_GAUSS41;     // Gauss-Kronrod integration rules (options: 15,21,31,41,51,61)
const size_t GSL_WorkSize = 1000000;               // work size used by GSL

gsl_integration_workspace *GSL_WorkSpace = NULL;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_BonStar
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
void Par_Init_ByFunction_BonStar( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      double *Prof_R = NULL;
      double *Prof_M = NULL;
      double *Prof_S = NULL;
      double *Table_SigmaProf_M = NULL;
      double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3], Sigma;

      Mass_AllRank = new real [NPar_AllRank];
      for (int d=0; d<3; d++)
      {
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      }


//    initialize the random number generator
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, BonStar_RSeed );


//    initialize GSL
      gsl_rng *GSL_RNG = gsl_rng_alloc( gsl_rng_mt19937 );
      gsl_rng_set( GSL_RNG, BonStar_RSeed );
      GSL_WorkSpace = gsl_integration_workspace_alloc( GSL_WorkSize );


//    determine the total enclosed mass within the maximum radius
      TotM = MassProf( BonStar_MaxR );
      ParM = TotM / NPar_AllRank;


//    construct the profile tables
      Prof_R = new double [BonStar_ProfNBin];
      Prof_M = new double [BonStar_ProfNBin];
      Prof_S = new double [BonStar_ProfNBin];

      dr = BonStar_MaxR / (BonStar_ProfNBin-1);

      for (int b=0; b<BonStar_ProfNBin; b++)
      {
         Prof_R[b] = dr*b;
         Prof_M[b] = MassProf( Prof_R[b] );
      }

      GetSigmaProf( BonStar_ProfNBin, Prof_R, Prof_S );


//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       (1) mass
         Mass_AllRank[p] = ParM;


//       (2) position
//       --> sample from the cumulative mass profile with linear interpolation
         RanM = RNG->GetValue( 0, 0.0, 1.0 )*TotM;
         RanR = Mis_InterpolateFromTable( BonStar_ProfNBin, Prof_M, Prof_R, RanM );

//       record the maximum error
         EstM     = MassProf( RanR );
         ErrM     = fabs( (EstM-RanM)/RanM );
         ErrM_Max = fmax( ErrM, ErrM_Max );

//       randomly set the position vector with a given radius
         RanVec_FixRadius( RanR, RanVec );
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + amr->BoxCenter[d];

//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }


//       (3) velocity
//       interpolate velocity dispersion at the particle position
         Sigma = Mis_InterpolateFromTable( BonStar_ProfNBin, Prof_R, Prof_S, RanR );

//       RanR must already lie within the interpolation table since it is returned from the table of mass profile
         if ( Sigma == NULL_REAL )
            Aux_Error( ERROR_INFO, "cannot determine velocity dispersion (target radius = %20.14e) !!\n", RanR );

//       assume velocity distribution is isotropic and Gaussian with a standard deviation of Sigma
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = gsl_ran_gaussian( GSL_RNG, Sigma );

      } // for (long p=0; p<NPar_AllRank; p++)

      Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
      Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max );


//    free memory
      delete RNG;
      delete [] Prof_R;
      delete [] Prof_M;
      delete [] Prof_S;
      gsl_rng_free( GSL_RNG );
      gsl_integration_workspace_free( GSL_WorkSpace );
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
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_BonStar



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r      : Input radius
//                RanVec : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = RNG->GetValue( 0, -1.0, +1.0 );
         RanR2     += SQR( RanVec[d] );
      }
   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf
// Description :  Mass profile
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//-------------------------------------------------------------------------------------------------------
double MassProf( const double r )
{

   double M;

#  if   ( ANAL_FORM == 1 )
   if ( r < r_core)  M = 4.0/3.0*M_PI*CUBE(r)*DensProf(r_core);
   else              M = 4.0/3.0*M_PI*CUBE(r_core)*DensProf(r_core) +
                         4.0*M_PI*( 1.0/7.0*a4*( POW7(r) - POW7(r_core) ) +
                                    1.0/6.0*a3*( POW6(r) - POW6(r_core) ) +
                                    1.0/5.0*a2*( POW5(r) - POW5(r_core) ) +
                                    1.0/4.0*a1*( POW4(r) - POW4(r_core) ) +
                                    1.0/3.0*a0*( CUBE(r) - CUBE(r_core) ) );
#  elif ( ANAL_FORM == 2 )
   M = rho0*4.0*M_PI*(  CUBE(r_core)*atan(r/r_core)/16.0 +
                     ( 3.0*POW4(r_core)*POW5(r) + 8.0*POW6(r_core)*CUBE(r) - 3.0*POW8(r_core)*r ) /
                     ( 48.0*POW6(r) + 144.0*SQR(r_core)*POW4(r) + 144.0*POW4(r_core)*SQR(r) + 48.0*POW6(r_core) )  );
#  endif

   return M;

} // FUNCTION : MassProf



//-------------------------------------------------------------------------------------------------------
// Function    :  DensProf
// Description :  Density profile
//
// Note        :  Calculate the mass density for a given radius
//
// Parameter   :  r  : Input radius
//-------------------------------------------------------------------------------------------------------
double DensProf( const double r )
{

   double rho;

#  if   ( ANAL_FORM == 1 )
   double rr;
   if ( r <= r_core )   rr = r_core;
   else                 rr = r;

   rho = a4*POW4(rr) + a3*CUBE(rr) + a2*SQR(rr) + a1*rr + a0;

#  elif ( ANAL_FORM == 2 )
   rho = rho0/pow( 1.0 + pow(r/r_core,pow1), pow2 );
#  endif

   return rho;

} // FUNCTION : DensProf



//-------------------------------------------------------------------------------------------------------
// Function    :  GetSigmaProf
// Description :  Calculate the interpolation table of velocity dispersion profile
//
// Note        :  1. Calculate velocity dispersion from Jeans equation
//                   --> grad( rho(r)*sigma(r)^2 ) = -rho(r)*grad( phi(r) ) = -G*rho(r)*M_tot(r)/r^2
//                   --> rho(r1)*sigma(r1)^2 = rho(r2)*sigma(r2)^2 + G*Integrate( rho(r)*M_tot(r)/r^2, [r, r1, r2] )
//                   --> It's essentially the same as calculating gas pressure
//                       --> Just replace pressure_gas(r) by rho_dm(r)*sigma(r)^2
//                2. Normalized to rho(Rzero)*sigma(Rzero)^2 = 0
//                3. r[] must be set in advance
//
// Parameter   :  NBin  : Number of radial bins in the table
//                r     : Radius at each bin
//                Sigma : Dark matter velocity dispersion at each bin
//
// Return      :  r, Sigma
//-------------------------------------------------------------------------------------------------------
void GetSigmaProf( const int NBin, const double *r, double *Sigma )
{

// set up GSL
   const double GSL_MaxAbsErr = 0.0;      // maximum allowed absolute error (0 --> disable)
   const double GSL_MaxRelErr = 1.0e-6;   // maximum allowed relative error (0 --> disable)
   gsl_function GSL_Func;

   GSL_Func.function = &GSL_IntFunc_SigmaProf;
   GSL_Func.params   = NULL;


// calculate rho*sigma^2 at Rmax by having rho(Rzero)*sigma(Rzero)^2 = 0
   double GSL_Result, GSL_AbsErr;
   gsl_integration_qag( &GSL_Func, r[NBin-1], BonStar_Rzero, GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                        GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

// gravitational constant is not included in the GSL integrand for better performance
   Sigma[NBin-1] = NEWTON_G*GSL_Result;


// calculate rho*sigma^2 at r < Rmax
   for (int b=NBin-2; b>=0; b--)
   {
      gsl_integration_qag( &GSL_Func, r[b], r[b+1], GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                           GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

      Sigma[b] = Sigma[b+1] + NEWTON_G*GSL_Result;
   }


// convert rho*sigma^2 to sigma
   for (int b=0; b<NBin; b++)    Sigma[b] = sqrt(  Sigma[b] / DensProf( r[b] )  );

} // FUNCTION : GetSigmaProf



//-------------------------------------------------------------------------------------------------------
// Function    :  GSL_IntFunc_SigmaProf
// Description :  GSL integrand for calculating the velocity dispersion profile
//
// Note        :  1. Return rho(r)*M(r)/r^2
//                2. Gavitational constant is not included in this integrand and must be multiplied later
//
// Parameter   :  r        : Radius
//                IntPara  : Integration parameters
//
// Return      :  See note above
//-------------------------------------------------------------------------------------------------------
double GSL_IntFunc_SigmaProf( double r, void *IntPara )
{

   return DensProf( r ) * MassProf( r ) / SQR(r);

} // FUNCTION : GSL_IntFunc_SigmaProf



#endif // #ifdef PARTICLE
