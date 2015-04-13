// Computes different helicity amplitudes as defined in 
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

#include<math.h>
#include"./functions.h"
#include<iostream>
#include<fstream>


const double PI = 4*atan(1);
const int low=1;
const int high=2;
const int forward=3;
const int backward=4;




int limits(double sred,double tred, double ured){

  double const shigh=pow(10.,9.);

  if ( sred <= 0.001 ) return low; // EFT limit
  else if ( ( sred <= 10. && -tred < 0.0001*sred ) ||
	    ( sred >  10. && sred <= shigh && -tred < 0.001)  ||
	    ( sred > shigh && -tred < 1.) ) return forward; // forward limit
  else if ( ( sred <= 10.  && -ured < 0.0001*sred ) ||
	    ( sred >  10.  && -ured < 0.001 ) ||
	    ( sred > shigh && -ured < 1.)) return backward;  // backward limit
  else if ( sred > shigh ) return high; // high energy limit
  else return 0; // no limit

  // explanation: 
  // for sred>shigh, optimal value to switch from HE limit to forward limit is |tred|=1
  // (at these value both limits are somewhat bad but quickly converge at either side)
  // for sred<shigh, only switch from exact result to forward limit at |t| < 0.001 for better accuracy  
};


void Mxxxx_fermion(double x, double y, double * re, double * im){
  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=1;
  *im=0;

  double z = - x - y;
  
  double temp;


 
  temp = 2 * ( y*y + z*z ) / ( x * x ) - 2/x;
  *re += temp * ( ReT(y) + ReT(z) );
  *im += temp * ( ImT(y) + ImT(z) );

  temp =  1/(2 * x * y) - 1/y ; 
  *re += temp * ReI(x,y); 
  *im += temp * ImI(x,y);

  temp =  1/(2 * x * z) - 1/z ; 
  *re += temp * ReI(x,z); 
  *im += temp * ImI(x,z);

  temp =  4/x +1/y +1/z + 1/(2*z*y)  -  2 * ( y*y + z*z ) / ( x * x );
  *re += temp * ReI(y,z);
  *im += temp * ImI(y,z);

  temp = 2* (y-z)/x;
  *re += temp * ( ReB(y) - ReB(z) );
  *im += temp * ( ImB(y) - ImB(z) );
     


  return;

};

void Mpppp_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
  // M++++ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region=limits(sred,tred,ured);
    if( region==low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *sred*sred ;  
      *im=0;  }
    else if( region == forward || region == backward )
      {            // Forward and backward limit 
	*re= 1./(2.* sred*sred)*( 2.* sred*sred+(-2.*sred+4.*sred*sred)*ReB(sred)
				  +(2.*sred-8.*sred*sred)*ReB(-sred) +(-1.+2.*sred)*ReT(sred)+(-1.-2.*sred+4.*sred*sred)*ReT(-sred) )  ;
	*im= 1./(2.* sred*sred)*(               (-2.*sred+4.*sred*sred)*ImB(sred)
						+(2.*sred-8.*sred*sred)*ImB(-sred) +(-1.+2.*sred)*ImT(sred)+(-1.-2.*sred+4.*sred*sred)*ImT(-sred) )  ;
      }
    else if( region== high ) 
      {            // high energy limit
	*re = 1.+(tred-ured)/sred * log(tred/ured)+(tred*tred+ured*ured) / (2.*sred*sred)*(pow( log(tred/ured) , 2 )+PI*PI);
	*im = 0;
      }
    else{   
      Mxxxx_fermion(sred,tred,re,im);
    };
    
  } 
  
  return;

};

void Mpmmp_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
  // M+--+ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;
  
  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else
    {
      int region = limits(sred,tred,ured); 
      if( region == low ){ // EFT limit
	*re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *tred*tred ;
	*im=0;  
      }
      else if( region==forward )
	{                // Forward limit 
	  *re=0.; *im=0.;
	}
      else if( region == backward )
	{                // Backward limit 
	  *re= 1./(2.* sred*sred)*( 2.* sred*sred+(2.*sred+4.*sred*sred)*ReB(-sred)+(-2.*sred-8.*sred*sred)*ReB(sred) +(-1.-2.*sred)*ReT(-sred)+(-1.+2.*sred+4.*sred*sred)*ReT(sred) )  ;
	  *im= 1./(2.* sred*sred)*(               (2.*sred+4.*sred*sred)*ImB(-sred)+(-2.*sred-8.*sred*sred)*ImB(sred) +(-1.-2.*sred)*ImT(-sred)+(-1.+2.*sred+4.*sred*sred)*ImT(sred) )  ;
	}
      else if( region == high ) 
	{            // high energy limit
	  *re = 1. + (sred-ured)/tred * log(-sred/ured)+(sred*sred+ured*ured) / (2.*tred*tred)*pow( log(-sred/ured) , 2 );
	  *im = -PI*( (sred-ured)/tred + (sred*sred+ured*ured) / (tred*tred)*log(-sred/ured));
	}
      else{  Mxxxx_fermion(tred,sred,re,im); }     
    }
  return;
 
};

void Mpmpm_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
// M+-+- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits(sred,tred,ured);     
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *ured*ured ;
      *im=0;  }
    else if( region == forward )
      {                // Forward limit 
	*re= 1./(2.* sred*sred)*( 2.* sred*sred+(2.*sred+4.*sred*sred)*ReB(-sred)+(-2.*sred-8.*sred*sred)*ReB(sred) +(-1.-2.*sred)*ReT(-sred)+(-1.+2.*sred+4.*sred*sred)*ReT(sred) )  ;
	*im= 1./(2.* sred*sred)*(               (2.*sred+4.*sred*sred)*ImB(-sred)+(-2.*sred-8.*sred*sred)*ImB(sred) +(-1.-2.*sred)*ImT(-sred)+(-1.+2.*sred+4.*sred*sred)*ImT(sred) )  ;    
      }
    else if( region == backward )
      {                // Backward limit 
	*re=0.; *im=0.;
      }
    else if( region == high ) 
      {            // high energy limit
	*re = 1. + (tred-sred)/ured * log(-tred/sred)+(sred*sred+tred*tred) / (2.*ured*ured)*pow( log(-tred/sred) , 2 );
	*im = PI*( (tred-sred)/ured + (sred*sred+tred*tred) / (ured*ured)*log(-tred/sred));
      }   
    else{   
      Mxxxx_fermion(ured,tred,re,im);     
    }
    
  }
  return;
  
};

void Mpppm_fermion(double sred, double tred, double * re, double * im, int exclude_loops){
  // M+--- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured=-sred-tred;


  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits (sred, tred, ured);    
    if( region == low ){ // EFT limit
      *re= 0.; *im=0.;}
    else if( region == forward || region == backward )
      {                // Forward and backward limit 
	*re= 0.; *im=0.;}
      else if (region == high ) 
      {            // high energy limit
	*re=-1;
	*im=0;
      }
    else{
      
      *re=-1;
      *im=0;
      
      temp=-1/sred-1/tred-1/ured;
      *re += temp*( ReT(sred) + ReT(tred) + ReT(ured) );
      *im += temp*( ImT(sred) + ImT(tred) + ImT(ured) );
      
      temp = 1/ured + 1/( 2 * sred * tred );
      *re += temp*ReI(sred,tred);
      *im += temp*ImI(sred,tred);
      
      temp = 1/tred + 1/( 2 * sred * ured );
      *re += temp*ReI(sred,ured);
      *im += temp*ImI(sred,ured);
      
      temp = 1/sred + 1/( 2 * tred * ured );
      *re += temp*ReI(tred,ured);
      *im += temp*ImI(tred,ured);
    }
    
  }


  return ;

};

void Mppmm_fermion(double sred, double tred, double * re, double * im, int exclude_loops){
  // M++-- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured = -sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred, tred, ured);
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +(7./90.) )*(sred*sred+tred*tred+ured*ured); *im=0.;}
    else if( region == forward || region == backward )
      {                // Forward and backward limit 
	*re=1./(2.*sred*sred)*( -2.*sred*sred-2.*sred*ReB(sred)+2.*sred*ReB(-sred)-ReT(sred)-ReT(-sred)  ) ; 
	*im=1./(2.*sred*sred)*(              -2.*sred*ImB(sred)+2.*sred*ImB(-sred)-ImT(sred)-ImT(-sred)  );
      }
    else if ( region == high ) 
      {            // high energy limit
	*re=-1;
	*im=0;
      }
    else{ 
      
      *re=-1;
      *im=0;
      
      temp = 1/( 2 * sred * tred );
      *re += temp*ReI(sred,tred);
      *im += temp*ImI(sred,tred);
      
      temp = 1/( 2 * sred * ured );
      *re += temp*ReI(sred,ured);
      *im += temp*ImI(sred,ured);
      
      temp = 1/( 2 * tred * ured );
      *re += temp*ReI(tred,ured);
      *im += temp*ImI(tred,ured);
      
    }
    
  }


  return;

};



void Mxxxx_vector(double x, double y, double * re, double * im){

  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=-1.5;
  *im=0;

  double z = - x - y;
  
  double temp;



  temp = -3* (y-z)/x;
  *re += temp * ( ReB(y) - ReB(z) );
  *im += temp * ( ImB(y) - ImB(z) );

 
  temp = -1/x*(8*x-3-6*y*z/x);
  *re += temp * ( ReT(y) + ReT(z) );
  *im += temp * ( ImT(y) + ImT(z) );

  temp =  1/x*(8*x-6-6*y*z/x)-4*(x-0.25)*(x-0.75)/(y*z);
  *re += temp * ReI(y,z); 
  *im += temp * ImI(y,z);


  temp = -4*(x-0.25)*(x-0.75)/(x*y); ; 
  *re += temp * ReI(x,y); 
  *im += temp * ImI(x,y);

  temp = -4*(x-0.25)*(x-0.75)/(x*z); ; 
  *re += temp * ReI(x,z); 
  *im += temp * ImI(x,z);
    

  return;


  
};

void Mpppp_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-sred-tred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred, tred, ured);
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *sred*sred ;
      *im=0;  }
    else if ( region == forward || region == backward )
      {           // Forward and backward limit 
	*re=-3./2.+8.*(sred-0.25)*(sred-0.75)/sred*ReB(sred)+
	  (-8.*(sred-0.25)*(sred-0.75)/sred+3.)*ReB(-sred)+
	  4.*(sred-0.25)*(sred-0.75)/(sred*sred)*ReT(sred)+
	  (4.*(sred-0.25)*(sred-0.75)/(sred*sred)-(8.*sred-3.)/sred)*ReT(-sred);
	*im=8.*(sred-0.25)*(sred-0.75)/sred*ImB(sred)+
	  (-8.*(sred-0.25)*(sred-0.75)/sred+3.)*ImB(-sred)+
	  4.*(sred-0.25)*(sred-0.75)/(sred*sred)*ImT(sred)+
	  (4.*(sred-0.25)*(sred-0.75)/(sred*sred)-(8.*sred-3.)/sred)*ImT(-sred);
      }
    else if (region == high)
      {          // high energy limit
	*re = - 1.*( 1.5 +  
		  1.5* (ured-tred)/sred * log(ured/tred) +
		  2. * ( 1. - 0.75 * tred*ured / (sred*sred) ) * ( pow( log(ured/tred),2) + PI*PI ) +
		  2. * sred * sred * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
				       log(4.*sred)*log(-4.*ured)/(sred*ured)+
				       log(-4.*ured)*log(-4.*tred)/(ured*tred) )
		  );
	*im =  ( 2. * PI * sred*sred *( log(-4.*ured)  / (sred * ured) +
					log(-4.*tred)  / (sred * tred ) ) 
		  )  ;
      }
    else{   Mxxxx_vector(sred,tred,re,im);  
    }
  }
  return;

};

void Mpmmp_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-sred-tred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits (sred, tred, ured);	
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *tred*tred ;
      *im=0;  }
    else if( region == forward )
      {                // Forward limit 
        *re=0.; *im=0.;
      }
    else if( region == backward )
      {                // Backward limit 
	*re=-3./2.-8.*(-sred-0.25)*(-sred-0.75)/sred*ReB(-sred)+
	  (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ReB(sred)+
	  4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ReT(-sred)+
	  (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ReT(sred);
	*im=-8.*(-sred-0.25)*(-sred-0.75)/sred*ImB(-sred)+
	  (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ImB(sred)+
	  4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ImT(-sred)+
	  (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ImT(sred);       
      }
    else if (region == high)
      {  // high energy limit
	*re = - ( 1.5 +  
		  1.5* (ured-sred)/tred * log(-ured/sred) +
		  2. * ( 1. - 0.75 * sred * ured / (tred*tred) ) *  pow( log( - ured/sred),2)   +
		  2. * tred * tred * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
				       log(4.*sred)*log(-4.*ured)/(sred*ured)+
				       log(-4.*ured)*log(-4.*tred)/(ured*tred) )  
		  );
	
	*im = - ( 1.5 * (sred-ured)/tred * (- PI)+ 
		  2. * ( 1. - 0.75 * sred*ured / (tred*tred) ) * PI * 2. * log(-ured/sred) + 
		  2. * (-PI) * tred*tred*  ( log(-4.*ured) / (ured*sred) +
					     log(-4.*tred) / (tred*sred)  )
		  )  ;
      }
    else{    Mxxxx_vector(tred,sred,re,im); }
  }
  return;

};

void Mpmpm_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-tred-sred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred,tred,ured);
    if( region == low )
      {     // EFT limit 
	*re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *ured*ured ;
	*im=0;  
      }
    else if( region == forward)
      {                // Forward limit 
 	*re=-3./2.-8.*(-sred-0.25)*(-sred-0.75)/sred*ReB(-sred)+
	  (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ReB(sred)+
	  4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ReT(-sred)+
	  (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ReT(sred);
	*im=-8.*(-sred-0.25)*(-sred-0.75)/sred*ImB(-sred)+
	  (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ImB(sred)+
	  4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ImT(-sred)+
	  (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ImT(sred);       
      }
    else if( region == backward )
      {                // Backward limit 
	*re=0.; *im=0.;
      }
    else if ( region == high )
      {         // high energy limit
	*re = - ( 1.5 +  
		  1.5* (tred-sred)/ured * log(-tred/sred) +
		  2. * ( 1. - 0.75 * sred * tred / (ured*ured) ) *  pow( log( - tred/sred),2)   +
		  2. * ured * ured * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
					      log(4.*sred)*log(-4.*ured)/(sred*ured)+
					      log(-4.*ured)*log(-4.*tred)/(ured*tred) )
		  );
	*im = - ( 1.5 * (sred-tred)/ured * (- PI)+ 
		  2. * ( 1. - 0.75 * sred*tred / (ured*ured) ) * PI * 2. * log(-tred/sred) + 
		  2. * (-PI) * ured*ured*  ( log(-4.*ured) / (ured*sred) +
					     log(-4.*tred) / (tred*sred)  )
		  )  ;	
      }
    else{   Mxxxx_vector(ured,tred,re,im); }
    
  }
  return;

};

void Mpppm_vector(double sred, double tred, double * re, double * im, int exclude_loops){

double ured=-tred-sred;


if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
else{ 

  //if(sred<0.001){ // EFT limit
  //		*re= 0.; *im=0.;}
  //else if(sred<10000. && sred>0.001 && (-tred<0.0001*sred||-ured<0.0001*sred ))
  //{                // Forward and backward limit 
  //            *re= 0.; *im=0.;}
  //else{

   Mpppm_fermion(sred,tred,re,im,exclude_loops);
  *re *= -1.5;
  *im *= -1.5;
  //}

    }

  return;
};

void Mppmm_vector(double sred, double tred, double * re, double * im, int exclude_loops){

  double ured = -sred-tred;

if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
else{   

  //if(sred<0.001){ // EFT limit 
  //		*re= -4.*(4.*(-5./32.)  +(27./40.) )*(sred*sred+tred*tred+ured*ured); *im=0.;}
  //else if(sred<10000. && sred>0.001 && (-tred<0.0001*sred||-ured<0.0001*sred ))
  //{                // Forward and backward limit 
  //*re=1./(2.*sred*sred)*( -2.*sred*sred-2.*sred*ReB(sred)+2.*sred*ReB(-sred)-ReT(sred)-ReT(-sred)  ) ; 
  //*im=1./(2.*sred*sred)*(              -2.*sred*ImB(sred)+2.*sred*ImB(-sred)-ImT(sred)-ImT(-sred)  );
  //*re *= -1.5;
  //*im *= -1.5;
  //}
  //else{ 

  Mppmm_fermion(sred,tred,re,im,exclude_loops);
  *re *= -1.5;
  *im *= -1.5;
  // }

   }

  return;

};


void Mpppp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*s*s   ;
*im=0;
  return;
};

void Mpmmp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*t*t   ;
*im=0;
  return;
};

void Mpmpm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*u*u   ;
*im=0;
  return;
};

void Mpppm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= 0 ;
*im=0;
  return;
};

void Mppmm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+zeta2)*(s*s+t*t+u*u)   ;
*im=0;
  return;
};
