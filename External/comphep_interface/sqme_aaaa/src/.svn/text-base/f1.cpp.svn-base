/************************************************
*    CompHEP version 4.5.1   *
*------------------------------------------------
* Copyright (C) 2001-2009, CompHEP Collaboration*
************************************************/
/*
                      A          !        A                          
                    -1---\       !     /---1-                        
                      P1 |       !     |  P1                         
                         |       !     |                             
                      A  |  A    !  A  |  A                          
                    -2---@-4-----!---4-@---2-                        
                      P2 |  P4   !  P4 |  P2                         
                         |       !     |                             
                         |  A    !  A  |                             
                         \-3-----!---3-/                             
                            P3   !  P3                               
*/
#include<math.h>
#include"out_ext.h"
#include"out_int.h"

namespace anom_aaaa {

extern double *Q0, *Q1, *Q2;
extern double va[3];

FNN F1;
double F1(void)
{
double TOTNUM,TOTDEN,RNUM,result;
static double C[2];double S[2];                                             
     
if(calcCoef[1])
{
//Error in R.S. Gupta (2011) lagrangian (a_2 term)
//S[0]=va[2]*va[2];
//C[0]=+va[1]*(64*va[1]-32*va[2])+4*S[0];
//C[1]=+va[1]*(80*va[1]-56*va[2])+13*S[0];

//S.Fichet G.Gersdorff Lagrangian
S[0]=va[2]*va[2];
C[0]=+va[1]*(32*va[2]+64*va[1])+4*S[0];
C[1]=+va[1]*(56*va[2]+80*va[1])+13*S[0];
}
TOTNUM=+16;
TOTDEN=+1;
S[0]=DP[3]*DP[3];
S[1]=DP[2]*DP[2];
RNUM=+DP[0]*(DP[5]*(C[0]*(-DP[3]*DP[2]-DP[4]*DP[1])+C[1]*DP[5]*DP[0]))+
 DP[1]*(DP[4]*(C[1]*DP[4]*DP[1]-C[0]*DP[3]*DP[2]))+C[1]*S[0]*S[1];
result=RNUM*(TOTNUM/TOTDEN);
 if(result>Fmax) Fmax=result; else if(result<-Fmax) Fmax=-result;
 return result;
}

} //namespace anom_aaaa
