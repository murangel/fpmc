/************************************************
*    CompHEP version 4.5.1   *
*------------------------------------------------
* Copyright (C) 2001-2009, CompHEP Collaboration*
************************************************/
#include <math.h>
#include <stdio.h>
#include"out_int.h"
#include"out_ext.h"

//#include"../src/num/include/alphas2.h"
//#include"../src/num/include/alphas_menu.h"
//#include"../src/num/include/strfun.h"
//#include"../src/num/include/q_kin.h"

namespace anom_aaaa {

extern double va[3];
char  processch[] = "A,A -> A,A";
double DP[45];
//static double xstr0,xstr1;
static int strfun_calc;
extern DNN d_1;
extern FNN F1;
static void sprod_(double*);

static double smpl(int nsub, double * momenta,int * err)
{
// double q;
 double ans, ans0=0.0, ans1=0.0;
switch(nsub)
{
  case 1:
    *err=*err|d_1(momenta);
    sprod_(momenta);
    ans0=0;
    ans1=0+F1();
    break;
}
/*
  if (!strfun_calc) {
    strfun_calc = 1;
    q = qcd_Scale_chep();
    if (ans0) xstr0 = strfun_(0, xbjo[0], xbjo[1], q);
    if (ans1) xstr1 = strfun_(1, xbjo[0], xbjo[1], q);
  }
*/
  ans = ans0 + ans1;
  if(!(*err) && 10000 * Fmax * computer_eps > (ans > 0 ? ans : -ans)) * 
 err = 1;
//  return ans0 * xstr0 + ans1 * xstr1;
  return ans0 + ans1;
}
#include"sqme0.c"

} //namespace anom_aaaa
