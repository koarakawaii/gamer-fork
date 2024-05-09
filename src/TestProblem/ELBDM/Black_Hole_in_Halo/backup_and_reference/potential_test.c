#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
   double r_normalized;
   printf("Enter the value for r_normalized (in unit of r_c).\n");
   scanf("%lf", &r_normalized);
   printf("r_normalized is %.4e .\n", r_normalized);

   double mass_of_sun       = 1.98892e33;             // in unit of g
   double newton_g_cgs      = 6.6743e-8;              // in cgs
   double ELBDM_MASS        = 1e-23;                  // in unit of eV/c^2
   double Const_eV          = 1.602176634e-12; 	      // in erg
   double Const_c           = 299792458.0*100;        // cm/2
   double UNIT_M            = 2.57725326366036e+44;   // in unit of g
   double UNIT_L            = 4.58351745817846e+24;   // in unit of cm
   double UNIT_V            = 1e7;                    // in unit of cm/s
   double ScaleFactor       = 1.;                     // scaling factor
   double CoreRadius        = 1.07713872e-03;         // in unit of Mpc/h
   double h_0               = 0.6732117;

   ELBDM_MASS               = ELBDM_MASS*Const_eV/pow(Const_c,2.)/UNIT_M;
   double m_a_22            = ELBDM_MASS*UNIT_M/(Const_eV/pow(Const_c,2.))/1e-22;        // eV/c^2 -> 10^-22 ev/c^2
   double factor            = pow( (pow(2.,1./8.) - 1.) , 0.5 );
   double X                 = factor*r_normalized;

//   double soliton_potential = newton_g_cgs*4.077703890131877e6/ScaleFactor/pow(m_a_22*10.,2.)/(CoreRadius*1000./h_0)*(factor*(-1732.5/(X*X+1.)-6641.25/pow(X*X+1.,2.)+3927./pow(X*X+1.,3.)-5915.25/pow(X*X+1.,4.)+324.5/pow(X*X+1.,5.)-1568.75/pow(X*X+1.,6.)+288.75*pow(X,12.)/pow(X*X+1.,6.)+3465.*log(factor))-3465.*atan(X)/r_normalized); // in unit of GM_sun/r_c, where r_c in unit of Mpc/h

//   printf("Potential at %.4e r_c is %.8e\n",r_normalized, soliton_potential*mass_of_sun/(CoreRadius*UNIT_L)/pow(UNIT_V,2.));

   double potential_scale   = newton_g_cgs*4.077703890131877e6/ScaleFactor/pow(m_a_22*10.,2.)/(CoreRadius*1000./h_0); // to here is in unit of GM_sun/r_c, where  r_c is in unit of Mpc/h
          potential_scale  *= mass_of_sun/(CoreRadius*UNIT_L)/pow(UNIT_V,2.);  // convert to code unit
//   printf("%.8e\n", potential_scale);
   double soliton_potential = potential_scale*(factor*(-1732.5/(X*X+1.)-6641.25/pow(X*X+1.,2.)+3927./pow(X*X+1.,3.)-5915.25/pow(X*X+1.,4.)+324.5/pow(X*X+1.,5.)-1568.75/pow(X*X+1.,6.)+288.75*pow(X,12.)/pow(X*X+1.,6.)+3465.*log(factor))-3465.*atan(X)/r_normalized); // in unit of GM_sun/r_c, where r_c in unit of Mpc/h
   printf("Potential at %.4e r_c is %.8e\n",r_normalized, soliton_potential);

   return 0;
}
