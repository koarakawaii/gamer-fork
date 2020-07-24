#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_EridanusII
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by Init_ExtPotAuxArray_EridanusII(), where
//                      UserArray[0] = x coordinate of the external potential center
//                      UserArray[1] = y ...
//                      UserArray[2] = z ..
//                      UserArray[3] = gravitational_constant*point_source_mass
//                3. Currently it does not support the soften length
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_EridanusII( const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   const bool RotatingFrame = ( UserArray[9] > 0.0 ) ? true : false;

   if ( RotatingFrame )
   {
      const double CM[3]       = { UserArray[0], UserArray[1], UserArray[2] };
      const real   GM          = UserArray[3];
      const real   R           = UserArray[4];
      const real   Vrot        = UserArray[5];
      const bool   FixedPos    = ( UserArray[6] > 0.0 ) ? true : false;
      const bool   Centrifugal = ( UserArray[7] > 0.0 ) ? true : false;
      const real   Angle0      = UserArray[8];

      real dx, dy, dz, dr2, _R, theta, Rx, Ry, Rz, tmp, phi;

//    calculate the relative coordinates between the target cell and the center of mass of the satellite
      dx  = (real)( x - CM[0] );
      dy  = (real)( y - CM[1] );
      dz  = (real)( z - CM[2] );
      dr2 = SQR(dx) + SQR(dy) + SQR(dz);

//    calculate the relative coordinates between the MW center and the center of mass of the satellite
//    --> assuming a circular orbit on the xy plane with an initial azimuthal angle of theta=0, where theta=acos(Rx/R)
      _R    = (real)1.0 / R;
      theta = ( FixedPos ) ? Angle0 : Vrot*Time*_R + Angle0;
      Rx    = R*COS( theta );
      Ry    = R*SIN( theta );
      Rz    = (real)0.0;
      tmp   = ( dx*Rx + dy*Ry + dz*Rz )*_R;
      phi   = (real)0.5*GM*CUBE(_R)*( dr2 - (real)3.0*SQR(tmp) );

      if ( Centrifugal )
      phi  -= (real)0.5*GM*CUBE(_R)*( SQR(dx) + SQR(dy) );

      return phi;
   } // if ( RotatingFrame )

   else
   {
      const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
      const real   GM     = (real)UserArray[3];
      const real   dx     = (real)(x - Cen[0]);
      const real   dy     = (real)(y - Cen[1]);
      const real   dz     = (real)(z - Cen[2]);
      const real   _r     = 1.0/SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

      return -GM*_r;
   } // if ( RotatingFrame ) ... else ...

} // FUNCTION : ExtPot_EridanusII



// =================================
// get the CPU/GPU function pointers
// =================================

#ifdef __CUDACC__
__device__
#endif
static ExtPot_t ExtPot_Ptr = ExtPot_EridanusII;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_EridanusII
// Description :  Return the function pointers to the CPU/GPU external potential routines
//
// Note        :  1. To enable this routine, link to the function pointers "SetCPU/GPUExtPot_Ptr"
//                   in a test problem initializer as follows:
//
//                      void SetCPUExtPot_EridanusII( ExtPot_t &CPUExtPot_Ptr );
//                      # ifdef GPU
//                      void SetGPUExtPot_EridanusII( ExtPot_t &GPUExtPot_Ptr );
//                      # endif
//
//                      ...
//
//                      SetCPUExtPot_Ptr = SetCPUExtPot_EridanusII;
//                      # ifdef GPU
//                      SetGPUExtPot_Ptr = SetGPUExtPot_EridanusII;
//                      # endif
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                2. Must obtain the CPU and GPU function pointers by separate routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_EridanusII( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_EridanusII( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_EridanusII( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...




#endif // #ifdef GRAVITY
