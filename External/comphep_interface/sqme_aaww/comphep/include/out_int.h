/*$Log */
#ifndef __ANOM_AAWW_OUT_INT__
#define __ANOM_AAWW_OUT_INT__

#include <stdlib.h>
#include <string.h>
#include <math.h>

namespace anom_aaww {

typedef int (DNN) (double *);
typedef double (FNN) (void);
extern double sqrMom (char *, double *);
extern double computer_eps;
extern double Fmax;
extern int *calcCoef;
extern double DP[];

} // namespace anom_aaww

#endif
