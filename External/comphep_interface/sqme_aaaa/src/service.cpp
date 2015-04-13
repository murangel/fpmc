/************************************************
*    CompHEP version 4.5.1   *
*------------------------------------------------
* Copyright (C) 2001-2009, CompHEP Collaboration*
************************************************/
#include<math.h>
#include"out_int.h"
#include"out_ext.h"

namespace anom_aaaa {

int gwidth = 0;
 double va[3] ={1.400000E+04
,0.000000E+00
,0.000000E+00};
const int nin_ = 2;

const int nout_ = 2;

const int nprc_ = 1;

int pinf_(int nsub,int nprtcl,char* pname, double* pmass)
{
int n;
char names[1][4][7]={
 {"A","A","A","A"}
};
int nvalue[1][4]={
 {0,0,0,0},
};
if  (nsub<0 ||nsub>1||nprtcl<0||nprtcl>4) return 1;
if(pname) strcpy(pname,names[nsub-1][nprtcl-1]);
if(pmass)
{
  n=nvalue[nsub-1][nprtcl-1];
if (n>2) if (calcFunc()) return FUCTION_ERROR;
if (n==0) *pmass=0; else *pmass=va[n];
if (*pmass<0) (*pmass)=-(*pmass);
}
return 0;
}
const int nvar_ = 2;

const int nfunc_ = 0;

int vinf_ (int numvar, char * name, double * val)
{
char names[3][10]={
"Sqrt(S)"
,"aA1"
,"aA2"};
   if (numvar < 0 || numvar > 2) return 1;
   if (name) strcpy (name, names[numvar]);
   if (val) *val=va[numvar];
   return 0;
}

int asgn_ (int numvar, double newval)
{
  if (numvar < 0|| numvar>2) return 1;
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
if(err) return 1; else return 0;
}
} //namespace anom_aaaa
