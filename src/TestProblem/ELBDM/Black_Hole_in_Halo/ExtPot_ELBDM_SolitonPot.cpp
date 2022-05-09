#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY


//-------------------------------------------------------------------------------------------------------
// Function    :  SolitonPot
// Description :  Compute the soliton potential
//
// Note        :  1. Invoked by ExtPot_ELBDM_SolitonPot
//
// Parameter   :  r_normalized: radius normalized by soliton core radius
//                PotFactor   : proportional factor for scaling soliton potential
//
// Return      :  soliton potential at a given normalized radius
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real SolitonPot(real r_normalized, real PotFactor)
{
   real soliton_potential;
   return soliton_potential;
}



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
extern double ELBDM_SolitonPot_M;
extern double SolitonSubCenter[3];

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
   AuxArray_Flt[3] = ELBDM_SolitonPot_Factor;

} // FUNCTION : SetExtPotAuxArray_ELBDM_SolitonPot
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================
extern double CoreRadius

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_ELBDM_SolitonPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_ELBDM_SolitonPot(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = soliton potential factor
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

   const double Cen[3]        = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const real   PotFactor     = (real)UserArray_Flt[3];
   const real   dx            = (real)(x - Cen[0]);
   const real   dy            = (real)(y - Cen[1]);
   const real   dz            = (real)(z - Cen[2]);
   const real   r_normaliized = SQRT( dx*dx + dy*dy + dz*dz )/CoreRadius;

   return SolitonPot(r_normalized, PotFactor);

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
