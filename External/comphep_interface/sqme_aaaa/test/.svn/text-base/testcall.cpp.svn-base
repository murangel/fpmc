#include <stdio.h>
#include <math.h>

#include "diff_cross.h"

namespace anom_aaaa {
extern int asgn_(int,double);
extern double sqme_deg(double, double, double);
extern double va[11];
}; // namespace anom_aaaa

#include "diff_cross.h"

int main() {

   printf("Call to comphep anomalous coupling ME^2\n");

   double sqrts;
   double s;

   double mw = 987654321; // aaaa !!!!!!!!! should NOT matter

   printf("mw SHOULD NOT MATTER (gg->gg). It is just a remanent from aaWW\n");
   printf("mw = %f\n",mw);

   //assigning aA1/aA2
   double aA1 = 1;  //TeV-4
   double aA2 = 0;



   anom_aaaa::asgn_(1, aA1);
   anom_aaaa::asgn_(2, aA2);

   printf("aA1 = %f, aA2 = %f\n",aA1,aA2);

   {
      int i;
      sqrts = 7000; s = sqrts*sqrts;
      printf(" test for energy sqrt(s)= %f\n", sqrts);
      for(i=0; i <= 10; i++) {
         double ds_dom_c = ds_dom_comp(s, 18.*i, mw);
         double sqme_c = anom_aaaa::sqme_deg(s, 18.*i, mw);

         printf("deg= %f  ds/do(comp) = %enb sqme = %e\n", 
               18.*i, ds_dom_c, sqme_c );
      }   

      sqrts = 14000;s = sqrts*sqrts;
      printf(" test for energy sqrt(s)= %f\n", sqrts);
      for(i=0; i <= 10; i++) {
         double ds_dom_c = ds_dom_comp(s, 18.*i, mw);
         double sqme_c = anom_aaaa::sqme_deg(s, 18.*i, mw);

         printf("deg= %f  ds/do(comp) = %enb sqme = %e\n", 
               18.*i, ds_dom_c,sqme_c);
      }
   }   

   return 0;
}
