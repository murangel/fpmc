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
                      A  |  W-   !  W- |  A                          
                    -2---@-4<----!--<4-@---2-                        
                      P2 |  P4   !  P4 |  P2                         
                         |       !     |                             
                         |  W+   !  W+ |                             
                         \-3>----!-->3-/                             
                            P3   !  P3                               
*/
#include<math.h>
#include"out_ext.h"
#include"out_int.h"

namespace anom_aaww {

extern double *Q0, *Q1, *Q2;
extern double va[11];

FNN F6;
double F6(void)
{
double TOTNUM,TOTDEN,RNUM,result;
static double C[11];double S[4];                                            
     
if(calcCoef[1])
{
S[0]=va[10]*va[10]*va[10]*va[10];
C[0]=+1216*S[0];
C[1]=+S[0]*(48*va[7]+320*va[6]);
S[1]=va[7]*va[7];
C[2]=+S[0]*(va[6]*(16*va[7]+32*va[6])+3*S[1]);
S[2]=va[10]*va[10];
C[3]=+S[2]*(112*va[7]+64*va[6]);
C[4]=+2*S[2]*S[1];
C[5]=+4*S[1];
C[6]=+32*(va[7]+va[6]);
C[7]=+16*va[7]+128*va[6];
C[8]=+va[6]*(8*va[7]+16*va[6])+S[1];
C[9]=+128*S[0];
S[3]=va[8]*va[8]*va[8]*va[8];
C[10]=+S[3];
}
TOTNUM=+C[10];
TOTDEN=+C[9];
S[0]=DP[5]*DP[5];
RNUM=+DP[0]*(S[0]*(C[8]*DP[0]-C[7])+C[4]*(DP[2]*DP[1]+DP[4]*DP[3])+C[2]*
 DP[0]-C[1])+DP[1]*(DP[2]*(C[3]+C[5]*DP[4]*DP[3])-C[6]*DP[5]*DP[4])+DP[3]*(
 C[3]*DP[4]-C[6]*DP[5]*DP[2])+C[0]+320*S[0];
result=RNUM*(TOTNUM/TOTDEN);
 if(result>Fmax) Fmax=result; else if(result<-Fmax) Fmax=-result;
 return result;
}

} //namespace anom_aaww
