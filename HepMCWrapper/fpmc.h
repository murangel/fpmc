#ifndef FPMC_INC
#define FPMC_INC

//------------------------- FPMC ------------------------------
//
extern "C" {
   void hwmodini_(void);
   void hwchek_(void);
   void hwfxer_(int*);
}
#define hwmodini hwmodini_
#define hwchek hwchek_
#define hwfxer hwfxer_

//      DOUBLE PRECISION XPQNRM(-6:6),PI,FMCONV,AION,ZION,XMSCUT
//      DOUBLE PRECISION RBMIN 
//      DOUBLE PRECISION AMASS,PMASS,GAPSPR,PROSPR,CDFFAC
//      CHARACTER*3 TYPEPR, TYPINT
//      INTEGER ISOFTM, NFLUX
      
//c---Normalizations of partonic PDFs inside the pomeron:
//      COMMON /PDFNRM/ XPQNRM
      
//c---Atomic and proton numbers of the ion
//      COMMON /ION/ AION,ZION,RBMIN
extern "C" {
  extern struct {
     double AION,ZION,RBMIN;
  } ion_;
}
#define ion ion_
      
//c---Flags carrying info whether the process is inclusive or exclusive:
//      COMMON /PRTYPE/ TYPEPR, TYPINT
//      PARAMETER(PI=3.141592654D0,FMCONV=0.197,XMSCUT=50.)
//      PARAMETER(AMASS=931.D-3,PMASS=938D-3)
extern "C" {
  extern struct {
     char TYPEPR[3], TYPINT[3];
  } prtype_;
}
#define prtype prtype_
      
//c---Model parameters
//      COMMON /XSECT/ GAPSPR,PROSPR,CDFFAC,ISOFTM,NFLUX
extern "C" {
  extern struct {
     double GAPSPR,PROSPR,CDFFAC;
     int    ISOFTM,NFLUX;
  } xsect_;
}
#define xsect xsect_

//c---PDF parameter
//      INTEGER IFITPDF
//      COMMON /PDFS/ IFITPDF
extern "C" {
  extern struct {
     int IFITPDF;
  } pdfs_;
}
#define pdfs pdfs_

//c---FLUX - POMERON + REGEON PARAMETERS
//      DOUBLE PRECISION alphaP,Bpom,alphaPp, zh1
//      DOUBLE PRECISION alphaR,alphaRp,Breg,Cr
//      COMMON /PRPARAMS/ alphaP,Bpom,alphaPp,alphaR,alphaRp,Breg,Cr,zh1

//c--- Anomalous coupling variables for anomalous AA->WW, AA-> ZZ, AA-> AA
//      INTEGER AAANOM
//      DOUBLE PRECISION D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
//     &        A1A, A2A
//      COMMON /AAANOMAL/ D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
//     &        A1A, A2A, AAANOM
extern "C" {
  extern struct {
     double D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
            A1A, A2A;
     int    AAANOM;
  } aaanomal_;
}
#define aaanomal aaanomal_
 
//c--- Variables for exotic  AA->AA
//      INTEGER AAEXOTIC
//      DOUBLE PRECISION AAM, AAQ, AAN
//      COMMON /AAEXOTICAL/ AAM, AAQ, AAN, AAEXOTIC
extern "C" {
  extern struct {
     double AAM, AAQ, AAN;
     int    AAEXOTIC;
  } aaexotical_;
}
#define aaexotical aaexotical_
      
//C     CHIDe
//      integer CHIDeIGLU
//      double precision CHIDeX, CHIDeXP, CHIDeS2, CHIDeS
//      common /CHIDeFPMC/ CHIDeX, CHIDeXP, CHIDeS2, CHIDeS, CHIDeIGLU
extern "C" {
  extern struct {
     double CHIDeX, CHIDeXP, CHIDeS2, CHIDeS;
     int    CHIDeIGLU;
  } chidefpmc_;
}
#define chidefpmc chidefpmc_

//C     KMR2
//      integer KMR2DELTA
//      double precision KMR2Q2CUT, KMR2SURV, KMR2SCALE
//      common /KMR2FPMC/ KMR2Q2CUT, KMR2SURV, KMR2SCALE, KMR2DELTA
extern "C" {
  extern struct {
     double KMR2Q2CUT, KMR2SURV, KMR2SCALE;
     int    KMR2DELTA;
  } kmr2fpmc_;
}
#define kmr2fpmc kmr2fpmc_

//c--- temporal variables
//      DOUBLE PRECISION TMAXMAX, TMINMIN
//      DOUBLE PRECISION SMAXMAX, SMINMIN
//      INTEGER ISSET, MYDEBUG
//      COMMON /OKTEMP/ ICOUNT
//      INTEGER ICOUNT

//-------------------------------------------------------------
//------------------------ FFCARD -----------------------------
//
extern "C" {
   void fpmc_var_ini_(int*);
}
#define fpmc_var_ini fpmc_var_ini_

extern "C" {
  extern struct {
     float  URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
            UYJMIN, UYJMAX,UPTMIN, UPTMAX, UEMMIN, UEMMAX, UDKAPPA,
            UACW, UA0W, UA0Z, UACZ, UA1A, UA2A,
            UAAM, UAAQ, UAAN,
            UCHIDeX, UCHIDeXP, UCHIDeS2,
            UXI1Min, UXI1Max, UXI2Min, UXI2Max,
            UCHIDeGapMin, UCHIDeGapMax,
            UKMR2Q2CUT, UKMR2SURV, UKMR2SCALE;
  } myffread1_;
}
#define myffread1 myffread1_
 
extern "C" {
  extern struct {
     float  UDLAMBDA, UANOMCUTOFF, UYWWMIN, UYWWMAX,
            UQ2WWMN, UQ2WWMX;
  } myffread2_;
}
#define myffread2 myffread2_

extern "C" {
  extern struct {
     int UOUTPUT,UMAXEV,UIPROC,UNFLUX,UNRN1,UNRN2,UIFIT, 
         UISOFTM, UZION, UAION, UBMIN, UAAANOM, UAAEXOTIC,
         UCHIDeIGLU, UKMR2DELTA;
  } myffread3_;
}
#define myffread3 myffread3_

extern "C" {
  extern struct {
     char UHADR[1];
  } cc0_;
}
#define cc0 cc0_

extern "C" {
  extern struct {
     char UTYPEPR[3];
  } cc1_;
}
#define cc1 cc1_

extern "C" {
  extern struct {
     char UTYPINT[3];
  } cc2_;
}
#define cc2 cc2_

extern "C" {
  extern struct {
     char UPART1[4];
  } cc3_;
}
#define cc3 cc3_

extern "C" {
  extern struct {
     char UPART2[4];
  } cc4_;
}
#define cc4 cc4_

extern "C" {
  extern struct {
     int UMODPDF1, UMODPDF2;
  } cc5_;
}
#define cc5 cc5_

extern "C" {
  extern struct {
     char UNTNAME[32], UCHIDePATH[32];
  } cyfflong1_;
}
#define cyfflong1 cyfflong1_

#endif
