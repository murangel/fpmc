#include <stdio.h>
#include <math.h>

#include "diff_cross.h"

namespace anom_aaww {
extern int asgn_(int,double);
extern double sqme_deg(double, double, double);
extern double va[11];
}; // namespace anom_aaww

#include "diff_cross.h"

int main() {
   printf("Call to comphem anomalous coupling ME^2\n");

   double sqrts = 167000.82;//4000;
//   double angle_deg = 0;

   double s = sqrts*sqrts;
   const double alpha = 1.0 / 137.0359895;

   double mw = 80.4;

   double dkappa = 0;
   double lambda = 0;
   double a0w = 0.0;
   double a0c = 0.0;
//   double mz = va[3];
   double sw = anom_aaww::va[2];
   double cw = sqrt(1-sw*sw);
   double mz = mw/cw;

   anom_aaww::asgn_(1, alpha);
   anom_aaww::asgn_(2, sw);
   anom_aaww::asgn_(3, mz);
   anom_aaww::asgn_(4, dkappa);
   anom_aaww::asgn_(5, lambda);
   anom_aaww::asgn_(6, a0w);
   anom_aaww::asgn_(7, a0c);

   {
      int i;
      sqrts = 2*mw+.0001; s = sqrts*sqrts;
      printf(" test for energy s= %f\n", sqrts);
      for(i=0; i <= 10; i++) {
         double ds_dom_c = ds_dom_comp(s, 18.*i, mw);
         double ds_dom_f = ds_dom_form(s, 18.*i, mw);
         double sqme_c = anom_aaww::sqme_deg(s, 18.*i, mw);

         printf("deg= %f  ds/do(comp) = %e sqme = %e rel_diff = %f\n", 
               18.*i, ds_dom_c, sqme_c, ds_dom_c/ds_dom_f );
      }   

      sqrts = 14000;s = sqrts*sqrts;
      printf(" test for energy s= %f\n", sqrts);
      for(i=0; i <= 10; i++) {
         double ds_dom_c = ds_dom_comp(s, 18.*i, mw);
         double ds_dom_f = ds_dom_form(s, 18.*i, mw);
         double sqme_c = anom_aaww::sqme_deg(s, 18.*i, mw);

         printf("deg= %f  ds/do(comp) = %e ds/do(form) = %f sqme = %e rel_diff = %f\n", 
               18.*i, ds_dom_c, ds_dom_f, sqme_c, ds_dom_c/ds_dom_f );
      }
   }   

   /*
   printf("ds/domega(comphep) = %f\n", ds_dom_c);
   printf("ds/domega(form)    = %f\n", ds_dom_f);
   printf("comphem/form       = %f\n", ds_dom_c/ds_dom_f);
   */
   return 0;
}
