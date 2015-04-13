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
   printf("Call to comphem anomalous coupling ME^2\n");

   double mw = 987654321; // aaaa !!!!!!!!! should NOT matter

   printf("mw SHOULD NOT MATTER (gg->gg). It is just a remanent from aaWW\n");
   printf("mw = %f\n",mw);

   double sqrts;
   double tests;
   double coup;

   sqrts = 14000;
   tests = sqrts*sqrts;
   coup = 0.1;
   printf("Testing at energy %f and mw = %f\n", sqrts, mw);
   for(int j=0; j<9; j++) {
      printf("coupling aA1 = %f\n", coup);
      anom_aaaa::asgn_(1, coup);
      for(int i=0; i <= 10; i++) {
         double sqme_c = anom_aaaa::sqme_deg(tests, 18.*i, mw);
         printf("angle %f\t sqme = %f\n", 18.*i, sqme_c); 
      }
      coup/=10;
   }

   sqrts = 1000; 
   tests = sqrts*sqrts;
   coup = 0.1;
   printf("----------------------\n");
   printf("Testing at energy %f and mw = %f\n", sqrts, mw);
   for(int j=0; j<9; j++) {
      printf("coupling aA1 = %f\n", coup);
      anom_aaaa::asgn_(1, coup);
      for(int i=0; i <= 10; i++) {
         double sqme_c = anom_aaaa::sqme_deg(tests, 18.*i, mw);
         printf("angle %f\t sqme = %f\n", 18.*i, sqme_c); 
      }
      coup/=10;
   }


   return 0;
}
