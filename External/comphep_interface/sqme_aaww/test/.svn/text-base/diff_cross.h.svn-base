#ifndef __DIFF_CROSS_H
#define __DIFF_CROSS_H
#include <stdio.h>
#include <math.h>

namespace anom_aaww {
extern double sqme(double, double, double);
extern double sqme_deg(double, double, double);
}; // namespace anom_aaww

#define debug 0

const double  GEV2toB = 0.3894*1E-3;   // 1 Gev^2 = GeV2b b
const double GEV2toNB = 0.3894*1E-3*1E9;    // 1 Gev^2 = GeV2b b


double ds_dom_form(double s, double theta_deg, double mw){

   if(debug > 99) printf("Entering herwig_factor ...\n");

   const double pifac = acos(-1.);
   const double alpha = 1.0 / 137.0359895;

   //  change to rad      
   const double theta_rad = pifac/180.*theta_deg;

   // if negative mass, use omega settings
   if(mw < 0) mw =  80.41;
             

   //     u,t gamma gamma -> w w  
   const double mwsq = mw*mw;
   const double t = -s/2.*(1-sqrt(1-4*mwsq/s)*cos(theta_rad))+mwsq;
   const double u = -s-t+2*mwsq;

   // check s,t,u
   if(debug > 99){
      printf("\n s = %f\n", s);
      printf("t = %f\n", t);
      printf("u = %f\n\n", u);
   }   
   
   // matrix factor evaluation
   const double beta = sqrt(1-4*mwsq/s);
   const double ss = s*s;
   const double mwsq_t = mwsq-t;
   const double mwsq_u = mwsq-u;
   const double mwsq_tt = mwsq_t*mwsq_t;
   const double mwsq_uu = mwsq_u*mwsq_u;

   double factor = 3*alpha*alpha/(2*s)*beta*
      ( 1 - 2*s*(2*s+3*mwsq)/(3*mwsq_t*mwsq_u) + 2*ss*(ss+3*mwsq*mwsq)/(3*mwsq_tt*mwsq_uu) );
   
   // convert to nb
   factor *= GEV2toNB;

   return(factor);
}

// triangular functions
double triang_funct(double x, double y, double z){
   const double result = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;

   return(result);
}

// phase space factor to multiply |ME|^2 to get dsigma/domega
double phase_space_factor(double s, double mw) {
   const double pifac = acos(-1.);
   const double mwsq = mw*mw;

   double result = 1./(64.*pifac*pifac) /s
         *sqrt( triang_funct(s, mwsq, mwsq) /triang_funct(s, 0, 0) );
   
   // convert to nb
   result *= GEV2toNB;

   return result;
}

double ds_dom_comp(double s, double angle_deg, double mw) {
   const double pi = acos(-1);
//   return sqme*phase_space_factor(s, mw);
  
   double beta = sqrt(1-4*mw*mw/s);
   return anom_aaww::sqme_deg(s,angle_deg, mw)*.3894*1E+6*beta*1/s/(64.*pi*pi); //[nb]
}

#endif
