#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
extern double SolitonSubCenter[3];
extern double SolitonCoreRadius;
extern double SolitonPotScale;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_ELBDM_SolitonPot
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_ELBDM_SolitonPot()
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_SolitonPot()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_ELBDM_SolitonPot( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = SolitonSubCenter[0];
   AuxArray_Flt[1] = SolitonSubCenter[1];
   AuxArray_Flt[2] = SolitonSubCenter[2];
   AuxArray_Flt[3] = SolitonCoreRadius;
   AuxArray_Flt[4] = (real)SolitonPotScale;

} // FUNCTION : SetExtPotAuxArray_ELBDM_SolitonPot
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================


//-------------------------------------------------------------------------------------------------------
// Function    :  SolitonMass
// Description :  Compute the soliton mass at given radius
//
// Note        :  1. Invoked by SetParameter in Init_TestProb_ELBDM_Black_Hole_in_Halo.cpp 
//                2. "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  sooliton_mass_scale : proportional factor for scaling soliton mass
//                r_normalized        : radius normalized by soliton core radius
//
// Return      :  soliton mass at a given normalized radius
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
double SolitonMass(double soliton_mass_scale, double r_normalized)
{
    double factor       = pow( (pow(2.,1./8.) -1.) , 0.5 );
    double soliton_mass = soliton_mass_scale/pow(1.+pow(factor*r_normalized,2),7)*(3465.*pow(factor*r_normalized,13.)+23100.*pow(factor*r_normalized,11.)+65373.*pow(factor*r_normalized,9.)+101376.*pow(factor*r_normalized,7.)+92323.*pow(factor*r_normalized,5.)+48580.*pow(factor*r_normalized,3.)-3465.*(factor*r_normalized)+3465.*pow(pow(factor*r_normalized,2)+1.,7.)*atan(factor*r_normalized));
   return soliton_mass;
}
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  SolitonPot
// Description :  Compute the soliton potential at given radius
//
// Note        :  1. Invoked by ExtPot_ELBDM_SolitonPot
//
// Parameter   :  soliton_pot_scale : proportional factor for scaling soliton potential
//                r_normalized      : radius normalized by soliton core radius
//
// Return      :  soliton potential at a given normalized radius
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real SolitonPot(real soliton_pot_scale, real r_normalized)
{
   real factor            = POW( (POW(2.,1./8.) - 1.) , 0.5 );  
   real X                 = factor*r_normalized;

   real soliton_potential = soliton_pot_scale*(factor*(-1732.5/(X*X+1.)-6641.25/POW(X*X+1.,2.)+3927./POW(X*X+1.,3.)-5915.25/POW(X*X+1.,4.)+324.5/POW(X*X+1.,5.)-1568.75/POW(X*X+1.,6.)+288.75*POW(X,12.)/POW(X*X+1.,6.)+3465.*LOG(factor))-3465.*ATAN(X)/r_normalized); // in code unit

   return soliton_potential;
}



//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_ELBDM_SolitonPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_ELBDM_SolitonPot(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = core radius
//                      UserArray_Flt[3] = soliton potential proportinoal factor
//                3. GenePtr has the size of EXT_POT_NGENE_MAX defined in Macro.h (default = 6)
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//                GenePtr           : Array of pointers for general potential tables
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_ELBDM_SolitonPot( const double x, const double y, const double z, const double Time,
                                     const double UserArray_Flt[], const int UserArray_Int[],
                                     const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

   const double Cen[3]                  = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const double r_core                  = UserArray_Flt[3];
   const double soliton_potential_scale = UserArray_Flt[4];
   const real   dx                      = (real)(x - Cen[0]);
   const real   dy                      = (real)(y - Cen[1]);
   const real   dz                      = (real)(z - Cen[2]);
   const real   r_normalized            = SQRT( dx*dx + dy*dy + dz*dz )/r_core;

   real soliton_pot = SolitonPot(soliton_potential_scale, r_normalized);
    
   return soliton_pot;
} // FUNCTION : ExtPot_ELBDM_SolitonPot



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_ELBDM_SolitonPot;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_ELBDM_SolitonPot
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_SolitonPot()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_ELBDM_SolitonPot( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_ELBDM_SolitonPot( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_ELBDM_SolitonPot( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_ELBDM_SolitonPot( double [], int [] );
void SetCPUExtPot_ELBDM_SolitonPot( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_ELBDM_SolitonPot( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_ELBDM_SolitonPot
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_ELBDM_SolitonPot()
{

   SetExtPotAuxArray_ELBDM_SolitonPot( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_ELBDM_SolitonPot( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_ELBDM_SolitonPot( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_ELBDM_SolitonPot

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
