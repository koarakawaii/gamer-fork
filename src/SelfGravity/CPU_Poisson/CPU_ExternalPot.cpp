#ifndef __EXTERNALPOT__
#define __EXTERNALPOT__



#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExternalPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by "Init_ExternalPot_Ptr", which
//                   points to Init_ExternalPot() by default but may be overwritten by various
//                   test problem initializers
//                3. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                   --> But one can easily modify this file to change the default behavior
//                4. Currently it does not support the soften length
//
// Return      :  External potential
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   const double CM[3]       = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM          = UserArray[3];
   const real   R           = UserArray[4];
   const real   Vrot        = UserArray[5];
   const bool   FixedPos    = ( UserArray[6] > 0.0 ) ? true : false;
   const bool   Centrifugal = ( UserArray[7] > 0.0 ) ? true : false;

   real dx, dy, dz, dr2, _R, theta, Rx, Ry, Rz, tmp, phi;

// calculate the relative coordinates between the target cell and the center of mass of the satellite
   dx  = (real)( x - CM[0] );
   dy  = (real)( y - CM[1] );
   dz  = (real)( z - CM[2] );
   dr2 = SQR(dx) + SQR(dy) + SQR(dz);

// calculate the relative coordinates between the MW center and the center of mass of the satellite
// --> assuming a circular orbit on the xy plane with an initial azimuthal angle of theta=0, where theta=acos(Rx/R)
   _R    = (real)1.0 / R;
   theta = ( FixedPos ) ? (real)0.0 : Vrot*Time*_R;
   Rx    = R*COS( theta );
   Ry    = R*SIN( theta );
   Rz    = (real)0.0;
   tmp   = ( dx*Rx + dy*Ry + dz*Rz )*_R;
   phi   = (real)0.5*GM*CUBE(_R)*( dr2 - (real)3.0*SQR(tmp) );

   if ( Centrifugal )
   phi  -= (real)0.5*GM*CUBE(_R)*( SQR(dx) + SQR(dy) );

   return phi;


} // FUNCTION : ExternalPot



#endif // #ifdef GRAVITY



#endif // #ifndef __EXTERNALPOT__
