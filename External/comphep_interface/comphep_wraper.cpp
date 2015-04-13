#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEBUG 0



 #ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
 extern "C" {
 #endif
// // routines called by fpmc - note the difference masses of mw and mz 
 void sqme_aaww_c__(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a0w, double *_aCw, double *_cutoff);

 void sqme_aaww_c_(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a0w, double *_aCw, double *_cutoff);
 
 void sqme_aazz_c__(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mz, double *_sw, double *_a0z, double *_aCz, double *_cutoff);

 void sqme_aazz_c_(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mz, double *_sw, double *_a0z, double *_aCz, double *_cutoff);

 void sqme_aaaa_c__(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a1a, double *_a2a, double *_cutoff);

 void sqme_aaaa_c_(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a1a, double *_a2a, double *_cutoff);
 #ifdef __cplusplus
 }
 #endif
  
// forward declaration
   namespace anom_aaww {
      extern int asgn_(int,double);
      extern double sqme_deg(double, double, double);
      extern double sqme(double, double, double);
   }; // namespace anom_aaww
   
   namespace anom_aazz {
      extern int asgn_(int,double);
      extern double sqme_deg(double, double, double);
      extern double sqme(double, double, double);
   }; // namespace anom_aazz

   namespace anom_aaaa {
      extern int asgn_(int,double);
      extern double sqme_deg(double, double, double);
      extern double sqme(double, double, double);
   }; // namespace anom_aaaa



// wrapper for g77
void sqme_aaww_c__(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a0w, double *_aCw, double *_cutoff) {

    sqme_aaww_c_(_amp2, _s, _t, _alpha,
      _mw, _sw, _dkappa, _lambda,
     _a0w, _aCw, _cutoff);
 }


void sqme_aaww_c_(double *_amp2, double *_s, double *_t, double *_alpha,
     double* _mw, double *_sw, double *_dkappa, double *_lambda,
     double *_a0w, double *_aCw, double *_cutoff) {
//////////////////////////////////////////////////////////////////// 
// anomalous aaww coupling:
//    _cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 

   // apply cutoff  for anomalous couplings
   double fact = (*_cutoff > 0 ) ? 1/pow(1+*_s/ *_cutoff / *_cutoff, 2) : 1;
   
   double dkappa = *_dkappa*fact;
   double lambda = *_lambda*fact;
   double a0w = *_a0w*fact;
   double aCw = *_aCw*fact;

// mw in comphep can be set through sw and mz only. see service.c for 
// details 

//   double sw = anom_aaww::va[2];
   double sw = *_sw;
   double cw = sqrt(1-sw*sw);
   double mz = *_mw/cw;

   anom_aaww::asgn_(1, *_alpha);
   anom_aaww::asgn_(2, sw);
   anom_aaww::asgn_(3, mz);
   anom_aaww::asgn_(4, dkappa);
   anom_aaww::asgn_(5, lambda);
   anom_aaww::asgn_(6, a0w);
   anom_aaww::asgn_(7, aCw);

   double amp2 = anom_aaww::sqme(*_s, *_t, *_mw);
   *_amp2 = amp2;
}   

//////////////////////////////////////////////////////////////////// 
//wrapper for g77
void sqme_aazz_c__(double *_amp2, double *_s, double *_t, double *_alpha,
     double* _mz, double *_sw, double *_a0z, double *_aCz, double *_cutoff) {

sqme_aazz_c_(_amp2, _s, _t, _alpha,
     _mz, _sw, _a0z, _aCz, _cutoff);
}

void sqme_aazz_c_(double *_amp2, double *_s, double *_t, double *_alpha,
     double* _mz, double *_sw, double *_a0z, double *_aCz, double *_cutoff) {
//////////////////////////////////////////////////////////////////// 
// anomalous aazz coupling:
//////////////////////////////////////////////////////////////////// 

   // apply cutoff  for dkappa and lambda
   double fact = (*_cutoff > 0 ) ? 1/pow(1+*_s/ *_cutoff / *_cutoff, 2) : 1;
   double a0z = *_a0z*fact;
   double aCz = *_aCz*fact;

//   double sw = anom_aazz::va[2];
   double sw = *_sw;
   double mz = *_mz;

   anom_aazz::asgn_(1, *_alpha);
   anom_aazz::asgn_(2, sw);
   anom_aazz::asgn_(3, mz);
   anom_aazz::asgn_(4, a0z);
   anom_aazz::asgn_(5, aCz);
//   anom_aazz::asgn_(6, ); //ee
//   anom_aazz::asgn_(7, ); //cw

   double amp2 = anom_aazz::sqme(*_s, *_t, mz);
   *_amp2 = amp2;
}   


//////////////////////////////////////////////////////////////////// 
//wrapper for g77
void sqme_aaaa_c__(double *_amp2, double *_s, double *_t, double *_alpha,
      double* _mw, double *_sw, double *_dkappa, double *_lambda,
      double *_a1a, double *_a2a, double *_cutoff) {

    sqme_aaaa_c_(_amp2, _s, _t, _alpha,
      _mw, _sw, _dkappa, _lambda,
     _a1a, _a2a, _cutoff);
 }

void sqme_aaaa_c_(double *_amp2, double *_s, double *_t, double *_alpha,
     double* _mw, double *_sw, double *_dkappa, double *_lambda,
     double *_a1a, double *_a2a, double *_cutoff) {
//////////////////////////////////////////////////////////////////// 
// anomalous aaaa coupling:
//    _cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 

   // apply cutoff  for anomalous couplings 
   //R.S. Gupta arXiv:1111.3354 [hep-ph] 
   double fact = (*_cutoff > 0 ) ? 1/(1+pow(*_s/ *_cutoff / *_cutoff, 2)) : 1;

   double dkappa = *_dkappa*fact;
   double lambda = *_lambda*fact;
   double a1a = *_a1a*fact;
   double a2a = *_a2a*fact;

// mw in comphep can be set through sw and mz only. see service.c for 
// details 

//   double sw = anom_aaaa::va[2];
   double sw = *_sw;
   double cw = sqrt(1-sw*sw);
   double mz = *_mw/cw;

   anom_aaaa::asgn_(1, a1a);
   anom_aaaa::asgn_(2, a2a);

   double amp2 = anom_aaaa::sqme(*_s, *_t, *_mw);
   *_amp2 = amp2;
} 
