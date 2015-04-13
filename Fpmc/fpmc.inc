C-----------------------------------------------------------------------
C * 08/08/2003 Maarten Boonekamp and Tibor Kucs
C * Common blocks and declarations needed for the modifications made to
C   HERWIG and POMWIG by Maarten Boonekamp and Tibor Kucs.
C * FMCONV = conversion factor from GeV to Pb
C * AION = atomic number of a heavy ion
C * ZION = proton number of a heavy ion
C * XMSCUT = cut on the invariant mass of the produced system
C            (dijet, dilepton, diphoton)
C * TYPEPR = type of the process 
C            'INC' = inclusive 'EXC' = exclusive
C * TYPINT = type of the interaction
C            'QED' or 'QCD'
C----------------------------------------------------------------------- 
      DOUBLE PRECISION XPQNRM(-6:6),PI,FMCONV,AION,ZION,XMSCUT
      DOUBLE PRECISION RBMIN 
      DOUBLE PRECISION AMASS,PMASS,GAPSPR,PROSPR,CDFFAC
      CHARACTER*3 TYPEPR, TYPINT
      INTEGER ISOFTM, NFLUX
      
c---Normalizations of partonic PDFs inside the pomeron:
      COMMON /PDFNRM/ XPQNRM
      
c---Atomic and proton numbers of the ion
      COMMON /ION/ AION,ZION,RBMIN
      
c---Flags carrying info whether the process is inclusive or exclusive:
      COMMON /PRTYPE/ TYPEPR, TYPINT
      PARAMETER(PI=3.141592654D0,FMCONV=0.197,XMSCUT=50.)
      PARAMETER(AMASS=931.D-3,PMASS=938D-3)
      
c---Model parameters
      COMMON /XSECT/ GAPSPR,PROSPR,CDFFAC,ISOFTM,NFLUX


c---PDF parametr
      INTEGER IFITPDF
      COMMON /PDFS/IFITPDF

c---FLUX - POMERON + REGEON PARAMETERS
      DOUBLE PRECISION alphaP,Bpom,alphaPp, zh1
      DOUBLE PRECISION alphaR,alphaRp,Breg,Cr
      COMMON /PRPARAMS/ alphaP,Bpom,alphaPp,alphaR,alphaRp,Breg,Cr,zh1

c--- Anomalous coupling variables for anomalous AA->WW, AA-> ZZ, AA-> AA
      INTEGER AAANOM
      DOUBLE PRECISION D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
     &        A1A, A2A
      COMMON /AAANOMAL/ D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
     &        A1A, A2A, AAANOM
 
c--- Variables for exotic  AA->AA
      INTEGER AAEXOTIC
      DOUBLE PRECISION AAM, AAQ, AAN
      COMMON /AAEXOTICAL/ AAM, AAQ, AAN, AAEXOTIC
      
C     CHIDe
      integer CHIDeIGLU
      double precision CHIDeX, CHIDeXP, CHIDeS2, CHIDeS
      common /CHIDeFPMC/ CHIDeX, CHIDeXP, CHIDeS2, CHIDeS, CHIDeIGLU
    
C     KMR2
      integer KMR2DELTA
      double precision KMR2Q2CUT, KMR2SURV, KMR2SCALE
      common /KMR2FPMC/ KMR2Q2CUT, KMR2SURV, KMR2SCALE, KMR2DELTA

c--- temporal variables
      DOUBLE PRECISION TMAXMAX, TMINMIN
      DOUBLE PRECISION SMAXMAX, SMINMIN
      INTEGER ISSET, MYDEBUG
      COMMON /OKTEMP/ ICOUNT
      INTEGER ICOUNT

      

