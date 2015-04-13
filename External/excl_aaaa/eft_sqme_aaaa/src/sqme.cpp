#include<iostream>
#include<fstream>
#include<math.h>
#include"helicity_amplitudes.h"
using namespace std;

namespace eft_aaaa {
    
const double alpha_em = 1./137.036; // EM coupling at zero momentum (on shell scheme)



void me_SM(void (*me)(double,double ,double *, double *, int),
	  double s,double t, double *re, double*im, int exclude_loops){
  // This routine computes the complex SM amplitude
  // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm


  // SM fermion content: (e,mu,tau,u,c,t,d,s,b)
  // SM_weight equals (number of colors) * (el. charge)^4  
  // SM masses in GeV
  
  const double SM_weight [9] = {1, 1, 1, 16./27., 16./27., 16./27., 1./27.,  1./27., 1./27.};
  const double SM_masses [9] = {0.5e-3,0.105,1.77,0.0023,1.28,173.07,0.0048,0.095,4.18};

  double d_re;
  double d_im;
  

  *re=0;
  *im=0;

  for (int i=0;i<=8;i++){
    me(s/(4*SM_masses[i]*SM_masses[i]),t/(4*SM_masses[i]*SM_masses[i]), &d_re, &d_im, exclude_loops);
    *re += d_re * SM_weight[i];
    *im += d_im * SM_weight[i];
  };

  //  cout<<*re<<"  "<<*im<<endl;


  // Add also the W contribution

  const double mW=80.385;  // W mass in GeV
 
 
  if (me==Mpppp_fermion){
    Mpppp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mppmm_fermion){
    Mppmm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpmpm_fermion){
    Mpmpm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpmmp_fermion){
    Mpmmp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpppm_fermion){
    Mpppm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);  
  }
   
  *re += d_re; 
  *im += d_im;
    
  *re *= 8 * alpha_em*alpha_em;
  *im *= 8 * alpha_em*alpha_em;
   
  // the factor of 8 is needed because of the conventions in
  // Costantini, DeTollis, Pistoni

  
  return;
    
};


// Computes the  squared matrix element and the SM interference from free zeta_1, zeta_2
double sqme (double s,double t, int exclude_loops_SM, double z1, double z2){

  double zeta1 = z1; // are in GeV^-4 units
  double zeta2 = z2;

  if (s<0 || t >0 || t<-s ){
    cout<<"Invalid domain. Valid range is s>=0 and -s<=t<=0"<<endl;
    return 0;
  }
 
// Now parameters 
/* 
  // read in input data

  ifstream data;
  
  data.open ("./zetas");
  
  data >> zeta1;
  data >> zeta2;


  data.close();
*/
  
  double re_ex;
  double im_ex;
  double re_SM;
  double im_SM;

  double value=0;
    
  // Mpppp:

  // the exotic matrix element:
  Mpppp_eft(zeta1,zeta2,s,t,&re_ex,&im_ex);
  re_ex *= 8;
  im_ex *= 8;
  // the factor of 8 is needed because of the conventions in
  // Costantini, DeTollis, Pistoni

  // the SM matrix element:
  me_SM(Mpppp_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);
    
  value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

  
  // repeat for the other helicities
  
  // Mppmm:

  Mppmm_eft(zeta1,zeta2,s,t,&re_ex,&im_ex);
  re_ex *= 8;
  im_ex *= 8;

  me_SM(Mppmm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);
    
  value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

  // Mpmmp:

  Mpmmp_eft(zeta1,zeta2,s,t,&re_ex,&im_ex);
  re_ex *= 8;
  im_ex *= 8;

  me_SM(Mpmmp_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);
    
  value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

  // Mpmpm:

  Mpmpm_eft(zeta1,zeta2,s,t,&re_ex,&im_ex);
  re_ex *= 8;
  im_ex *= 8;

  me_SM(Mpmpm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);
    
  value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

  // Mpppm

  Mpppm_eft(zeta1,zeta2,s,t,&re_ex,&im_ex);
  re_ex *= 8;
  im_ex *= 8;

  me_SM(Mpppm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);
    
  value += 4* (  re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) );
 
  return 0.5*value;

};


}  //namespace eft_aaaa
