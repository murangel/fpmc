/*********************************************************
 Copyright (C), 200-2003, CompHEP Collaboration
----------------------------------------------------------
*
* $Id: out_ext.h,v 1.1 2009/04/22 14:10:25 okepka Exp $
*
* $Log: out_ext.h,v $
* Revision 1.1  2009/04/22 14:10:25  okepka
* Adding fergotten files
*
* Revision 1.3  2009/04/22 13:54:59  okepka
* Namespaces in comphep, anomalous ZZ, change in energy cut to select survival p.
*
* Revision 1.1  2008/12/15 15:58:39  okepka
* comphep
*
* Revision 1.1.1.1  2008-12-13 19:25:18  olda
* Start
*
* Revision 1.3  2003/07/20 14:34:37  kryukov
* Change width menu. Remove rwidth global vars.
* gwidth=0(fix), 1(overall), 2(running).
*
**********************************************************/
#ifndef __ANOM_AAZZ_OUT_EXT__
#define __ANOM_AAZZ_OUT_EXT__

#include<stdlib.h>
#include<string.h>
#include<math.h>

namespace anom_aazz {

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

} // namespace anom_aazz

#endif
