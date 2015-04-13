/************************************************
*    CompHEP version 4.4p3      *
*------------------------------------------------
* Copyright (C) 2001-2003, CompHEP Collaboration*
************************************************/
#include<math.h>
#include"out_int.h"
#include"out_ext.h"

namespace anom_aaww {

int gwidth=0;
 double va[11] ={3.162278E+02
,7.757350E-03
,4.717500E-01
,9.118760E+01
,0.000000E+00
,0.000000E+00
,0.000000E+00
,0.000000E+00
,0.000000E+00
,8.817324E-01
,8.040306E+01};
const int nin_ = 2;

const int nout_ = 2;

const int nprc_ = 1;

int pinf_(int nsub,int nprtcl,char* pname, double* pmass)
{
int n;
char names[1][4][7]={
 {"A","A","W+","W-"}
};
int nvalue[1][4]={
 {0,0,10,10},
};
if  (nsub<0 ||nsub>1||nprtcl<0||nprtcl>4) return 1;
if(pname) strcpy(pname,names[nsub-1][nprtcl-1]);
if(pmass)
{
  n=nvalue[nsub-1][nprtcl-1];
if (n>7) if (calcFunc()) return FUCTION_ERROR;
if (n==0) *pmass=0; else *pmass=va[n];
if (*pmass<0) (*pmass)=-(*pmass);
}
return 0;
}
const int nvar_ = 7;

const int nfunc_ = 3;

int vinf_(int numvar,char * name, double * val)
{
char names[11][10]={"Sqrt(S)"
,"alfEMZ"
,"SW"
,"MZ"
,"DkappA"
,"LambA"
,"a0W"
,"aCW"
,"EE"
,"CW"
,"MW"};
   if(numvar<0||numvar>10  ) return 1;
   if(name) strcpy(name,names[numvar]);
   if(val) *val=va[numvar];
   return 0;
}

int asgn_(int numvar,double newval)
{
  if(numvar < 0|| numvar>7) return 1;
   va[numvar]=newval;
   return 0;
}

void cStrings(int nsub,int *nC, int * power, int **  chains)
{
   switch(nsub)
   {
   case 1 :    *nC=0; *power=0; *chains=NULL; break;
   default: *nC=0; *power=0; *chains=NULL;
   }
}

#include"extern.h"
int calcFunc(void)
{
int err=0;
   va[8]=sqrt((double)16*atan((double)1.)*va[1]);
   if(!finite(va[8])) return FUCTION_ERROR;
   va[9]=sqrt((double)1-pow(va[2],(double)2));
   if(!finite(va[9])) return FUCTION_ERROR;
   va[10]=va[3]*va[9];
   if(!finite(va[10])) return FUCTION_ERROR;
if(err) return 1; else return 0;
}
} //namespace anom_aaww
