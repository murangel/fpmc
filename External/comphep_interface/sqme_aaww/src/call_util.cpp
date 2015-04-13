#include <stdio.h>
#include <math.h>

namespace anom_aaww {

extern double sqme_(int nsub, double *momenta, int *err);


//////////////////////////////////////////////////////////////////// 
// sqme wrapper for comphep
//////////////////////////////////////////////////////////////////// 
double sqme(double s, double t, double mw) {

   double sqrts = sqrt(s);   
   const double phi = 0;

   double u = -t - s + 2*mw*mw;
   double beta=sqrt(1-4*mw*mw/s);
   //double costheta=(1+2*(t-mw*mw)/s)/beta;
   // this holds in general case too!
   double costheta = (t-u)/beta/s;

   double momenta[4*4];

   momenta[0] = sqrts/2;
   momenta[1] = 0;
   momenta[2] = 0;
   momenta[3] = sqrts/2;

   momenta[4+0] = sqrts/2;
   momenta[4+1] = 0;
   momenta[4+2] = 0;
   momenta[4+3] = -sqrts/2;

   double e = sqrts/2;
   double p = sqrt(e*e - mw*mw);

   double pz = p*costheta;
   double pt = p*sqrt(1-costheta*costheta);
   double px = pt*cos(phi);
   double py = pt*sin(phi);

   momenta[8+0] = e;
   momenta[8+1] = px;
   momenta[8+2] = py;
   momenta[8+3] = pz;

   momenta[12+0] = e;
   momenta[12+1] = -px;
   momenta[12+2] = -py;
   momenta[12+3] = -pz;

   int err;
   return sqme_(1, momenta, &err);
}

double sqme_deg(double s, double angle_deg, double mw) {
   const double pi = acos(-1);
   double cs = cos(angle_deg/180.*pi);

   double beta=sqrt(1-4*mw*mw/s);
   double t = -s/2.*(1.-beta*cs)+mw*mw;

   return sqme(s, t, mw);
}

} //namespace anom_aaww


