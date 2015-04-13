/************************************************
*    CompHEP version 4.4p3      *
*------------------------------------------------
* Copyright (C) 2001-2003, CompHEP Collaboration*
************************************************/
#include<math.h>


#define real double
#include"out_int.h"
#include"out_ext.h"

namespace anom_aaww {

extern real *Q0,*Q1,*Q2;
extern double va[11];
DNN  d_1;
int d_1(real * momenta)
{int I,err=0;
real s0max=0;
 for(I=0;I<nin_;I++) s0max+=momenta[4*I];
s0max=computer_eps*s0max*s0max;
 if(Q1!=NULL) free(Q1);
 if(Q2!=NULL) free(Q2);
 Q1=(real*)malloc(sizeof(real)*3);
 Q2=(real*)malloc(sizeof(real)*3);
Q1[2]=va[10]*va[10]- sqrMom("\1\3",momenta);
Q1[1]=va[10]*va[10]- sqrMom("\1\4",momenta);
for ( I=1; I<=2;I++)
  {
  if((Q1[I]>0? Q1[I]:-Q1[I]) < 10*s0max) err=DENOMINATOR_ERROR;
  if(!Q1[I]) Q1[I]=s0max;
  Q1[I]=1/Q1[I];
  Q2[I]=Q1[I]*Q1[I];
}
return err;
}

} //namespace anom_aaww
