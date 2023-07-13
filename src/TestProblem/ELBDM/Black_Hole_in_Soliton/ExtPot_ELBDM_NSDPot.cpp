#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__
extern double NSDPotCenter[3];
extern double NEWTON_G;
extern double NSDHalfMassRadius;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_ELBDM_NSDPot
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_ELBDM_NSDPot()
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_NSDPot()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_ELBDM_NSDPot( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = NSDPotCenter[0];
   AuxArray_Flt[1] = NSDPotCenter[1];
   AuxArray_Flt[2] = NSDPotCenter[2];
   AuxArray_Flt[3] = NSDHalfMassRadius;
   AuxArray_Flt[4] = NEWTON_G;

} // FUNCTION : SetExtPotAuxArray_ELBDM_NSDPot
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================



//-------------------------------------------------------------------------------------------------------
// Function    :  NSDPot
// Description :  Compute the NSD potential at given radius
//
// Note        :  1. Invoked by ExtPot_ELBDM_NSDPot
//
// Parameter   :  r_normalized      : radius normalized by half mass radius
//                newton_g          : NEWTON_G in code unit
//
// Return      :  soliton potential at a given normalized radius
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real NSDPot(real r_normalized, real newton_g)
{
   const real r_normalized_criteria1 = 2.6719e-04;        // criteria 1 for normalized radius
   const real r_normalized_criteria2 = 1.0254e+01;        // criteria 2 for normalized radius
   const real nsd_psuedo_mass        = 7.7273e+01;
   real nsd_pot;

   if ( r_normalized < r_normalized_criteria1 )
      nsd_pot = -4.3868e+00;
   else if ( r_normalized < r_normalized_criteria2 )
      nsd_pot =  ( 2.95066303e-06*pow(log10(r_normalized),20) +\
                   7.52118430e-05*pow(log10(r_normalized),19) +\
                   8.23946187e-04*pow(log10(r_normalized),18) +\
                   4.95235688e-03*pow(log10(r_normalized),17) +\
                   1.67569639e-02*pow(log10(r_normalized),16) +\
                   2.47868162e-02*pow(log10(r_normalized),15) +\
                  -2.93701826e-02*pow(log10(r_normalized),14) +\
                  -1.86084959e-01*pow(log10(r_normalized),13) +\
                  -2.25567594e-01*pow(log10(r_normalized),12) +\
                   2.56302237e-01*pow(log10(r_normalized),11) +\
                   8.89869540e-01*pow(log10(r_normalized),10) +\
                   2.88839423e-01*pow(log10(r_normalized), 9) +\
                  -1.43347858e+00*pow(log10(r_normalized), 8) +\
                  -1.37561221e+00*pow(log10(r_normalized), 7) +\
                   1.30999880e+00*pow(log10(r_normalized), 6) +\
                   2.28694285e+00*pow(log10(r_normalized), 5) +\
                  -7.56088670e-01*pow(log10(r_normalized), 4) +\
                  -2.68798838e+00*pow(log10(r_normalized), 3) +\
                   3.35077389e-01*pow(log10(r_normalized), 2) +\
                   3.37026028e+00*pow(log10(r_normalized), 1) +\
                   2.00631692e+00) - 4.3868e+00;
  else
     nsd_pot = -newton_g*nsd_psuedo_mass/r_normalized;

   return nsd_pot;
}



//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_ELBDM_NSDPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_ELBDM_NSDPot(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = NSD half mass radius
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
static real ExtPot_ELBDM_NSDPot( const double x, const double y, const double z, const double Time,
                                     const double UserArray_Flt[], const int UserArray_Int[],
                                     const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{

   const double Cen[3]                  = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
   const double r_half_mass             = UserArray_Flt[3];
   const real   dx                      = (real)(x - Cen[0]);
   const real   dy                      = (real)(y - Cen[1]);
   const real   dz                      = (real)(z - Cen[2]);
   const real   r_normalized            = SQRT( dx*dx + dy*dy + dz*dz )/r_half_mass;
   const real   newton_g                = (real)UserArray_Flt[4];

   real nsd_pot = NSDPot(r_normalized, newton_g);
    
   return nsd_pot;
} // FUNCTION : ExtPot_ELBDM_NSDPot



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_ELBDM_NSDPot;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_ELBDM_NSDPot
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_ELBDM_NSDPot()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_ELBDM_NSDPot( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_ELBDM_NSDPot( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_ELBDM_NSDPot( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_ELBDM_NSDPot( double [], int [] );
void SetCPUExtPot_ELBDM_NSDPot( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_ELBDM_NSDPot( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_ELBDM_NSDPot
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
void Init_ExtPot_ELBDM_NSDPot()
{

   SetExtPotAuxArray_ELBDM_NSDPot( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_ELBDM_NSDPot( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_ELBDM_NSDPot( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_ELBDM_NSDPot

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
