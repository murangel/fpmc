/************************************************
*    CompHEP version 4.4p3      *
*------------------------------------------------
* Copyright (C) 2001-2003, CompHEP Collaboration*
************************************************/
/*
                      A          !        A                          
                    -1---\       !     /---1-                        
                      P1 |       !     |  P1                         
                         |       !     |                             
                      A  |  Z    !  Z  |  A                          
                    -2---@-4-----!---4-@---2-                        
                      P2 |  P4   !  P4 |  P2                         
                         |       !     |                             
                         |  Z    !  Z  |                             
                         \-3-----!---3-/                             
                            P3   !  P3                               
*/
#include<math.h>
#include"out_ext.h"
#include"out_int.h"

namespace anom_aazz {

extern double *Q0, *Q1, *Q2;
extern double va[8];

FNN F1;
double F1(void)
{
double TOTNUM,TOTDEN,RNUM,result;
static double C[6];double S[5];                                             
     
if(calcCoef[1])
{
S[0]=va[5]*va[5];
S[1]=va[3]*va[3]*va[3]*va[3];
C[0]=+S[1]*(va[4]*(16*va[5]+32*va[4])+3*S[0]);
S[2]=va[3]*va[3];
C[1]=+2*S[0]*S[2];
C[2]=+4*S[0];
C[3]=+va[4]*(8*va[5]+16*va[4])+S[0];
S[3]=va[7]*va[7]*va[7]*va[7];
C[4]=+256*S[3]*S[1];
S[4]=va[6]*va[6]*va[6]*va[6];
C[5]=+S[4];
}
TOTNUM=+C[5];
TOTDEN=+C[4];
S[0]=DP[5]*DP[5];
RNUM=+DP[0]*(DP[0]*(C[0]+C[3]*S[0])+C[1]*(DP[2]*DP[1]+DP[4]*DP[3]))+C[2]*
 DP[4]*DP[3]*DP[2]*DP[1];
result=RNUM*(TOTNUM/TOTDEN);
 if(result>Fmax) Fmax=result; else if(result<-Fmax) Fmax=-result;
 return result;
}

} //namespace anom_aazz
