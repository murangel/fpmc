/*********************************************************
 Copyright (C), 200-2008, CompHEP Collaboration
----------------------------------------------------------
*
* $Id$
*
* $Log$
*
**********************************************************/
#ifndef __ANOM_AAAA_OUT_EXT__
#define __ANOM_AAAA_OUT_EXT__

#include<stdlib.h>
#include<string.h>
#include<math.h>

namespace anom_aaaa {

#ifdef INTERPRET

extern int nin_;
extern int nout_;
extern int nprc_;
extern int nvar_;
extern int nfunc_;

#else

extern const int nin_;
extern const int nout_;
extern const int nprc_;
extern const int nvar_;
extern const int nfunc_;

#endif

extern int pinf_ (int nsub, int nprtcl, char *pname, double *pmass);
extern int vinf_ (int numvar, char *name, double *val);
extern int asgn_ (int numvar, double valnew);
extern double sqme_ (int nsub, double *momenta, int *err);
extern int calcFunc (void);
extern int gwidth;
extern char processch[];

extern double *color_weights;


extern void cStrings (int nsub, int *nC, int *power, int **chains);


#define  FUCTION_ERROR       1
#define  DENOMINATOR_ERROR   2
#define  COMPHEP_ZERO        (1.E-20)

} // namespace anom_aaaa

#endif
