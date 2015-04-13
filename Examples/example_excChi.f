C-----------------------------------------------------------------------
C * Oct 2004 M.Boonekamp
C * Exclusive Chi production
C----------------------------------------------------------------------- 
      PROGRAM DPECS
c---Common block
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'
c---User's declarations
      DOUBLE PRECISION CMSENR
      INTEGER N
      CHARACTER ANSWER
c---External
      EXTERNAL HWUDAT

c---Events in the run
      MAXEV=10000

c---Beam particles (In DPE, E+ stands for initial state proton - this is corrected afterwards)
      PART1='E+'
      PART2='E+'

c---Set cms energy
      CMSENR=2D3

c---Beam momenta
      PBEAM1=CMSENR/2d0
      PBEAM2=PBEAM1

c---Process type : EXClusive or INClusive, QCD or QED when relevant
c   It is important to set those right
      TYPEPR='EXC'
      TYPINT='QCD'

c---Setting the hard subprocess (see manual)
      IPROC=19701

c---Option to include hadronization and showering effects
      ANSWER='N' ! say 'N' to skip this part

c---Initialization of other common blocks
c...The default parameters has to be changed after this call
      CALL HWIGIN
      SUSYIN=.TRUE.

C ... increase inefficiency tolerance
      EFFMIN = 1d-6

c---User's default kinematic parameters
c... beam momentum transfer range ( Q2 = |t| )
      Q2WWMN=0d0
      Q2WWMX=4d0
c... beam momentum loss range
      YWWMIN=0.d0
      YWWMAX=0.1D0

c---Choosing the flux : 9,10 is Cox-Forshaw (DPE)
c                       11 is Bialas-Landshoff (exc. DPE) or Boonekamp et al (inc. DPE)
c                       12 is Cahn, Jackson (heavy ions, QED)
c                       13 is Drees, Ellis, Zeppenfled (heavy ions, QED)
c                       14 is Papageorgiu (proton, QED)
c                       15 is Budnev flux (RECOMMENDED for proton, QED)
c                       16 is KMR flux (tables from L.Lonnblad or ExHume)
      NFLUX=11

c---Choosing pdf to use: 
c                       according to hep/ph 0609291
c                       10:  H1      
c                       20:  Zeus 
c                       30:  H1Zeus combined 
c                       older versions for compat., see h1qcd.f   
c                        2:  NLO fit as in H1
c                        5:  LO  fit as in H1
c                        8:  extended version of H1
      IFITPDF = 10
      
c---Set normalisations for Exclusive processes
c... option for soft corrections : ISOFTM = 0 : no correction
c                                           1 : simple factor, see below
c                                           2 : KMR low mass diffractive
c                                           3 : effective opacity model
      ISOFTM = 1

c---Events printed
      MAXPR=3

c---Initialize model/pdf dependant parameters
      CALL HWMODINI

c---User's initial calculations
      CALL HWABEG

c---Compute parameter dependent constants
      CALL HWUINC

c---Check POMWIG Settings + Initialisations for consistency
      CALL HWCHEK

c---Call HWUSTA to make any particle stable
      CALL HWUSTA('PI0     ')      

c---Initialize elementary process
      CALL HWEINI

c---Initialize event record fixing : this will replace the beam 
c   electrons by protons, radiated photons by pomerons/reggeons etc
      CALL HWFXER(.TRUE.)

c---Loop over events
      DO 100 N=1,MAXEV
c...Initialize event
         CALL HWUINE
c...Generate hard subprocesses
         CALL HWEPRO
c...Include showering and hadronization
         IF (ANSWER.EQ.'Y') THEN
            CALL HWBGEN
            CALL HWDHOB
            CALL HWCFOR
            CALL HWCDEC
            CALL HWDHAD
            CALL HWDHVY
            CALL HWMEVT
         END IF
c...Finish event
         CALL HWUFNE
c...Fix event record (i.e. restore correct intermediate states); print result
         CALL HWFXER(.FALSE.)
         IF(N.LE.MAXPR) THEN
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, 'AFTER EVENT RECORD FIXING:'
           CALL HWUEPR
         ENDIF
c...User's event analysis
         CALL HWANAL
 100  CONTINUE

c---Terminate elementary process
      CALL HWEFIN

c---User's terminal calculations
      CALL HWAEND      
      STOP
      END
C-----------------------------------------------------------------------
C * 07/03/2003, Tibor Kucs                                             
C * User's routine for initialization                                
C----------------------------------------------------------------------- 
      SUBROUTINE HWABEG
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'
      INTEGER NWPAWC
      REAL*4 HMEMOR
      PARAMETER(NWPAWC = 200000)
      COMMON /PAWC/HMEMOR(NWPAWC)
  
c---HBOOK initialization
      CALL HLIMIT(NWPAWC)
      CALL HROPEN(33,'HWIG','fpmc.hbook','NC',1024,ISTAT)
c---Histograms
      CALL HBOOK1(101,'Proton v1sq     ',50,0.,2.,0.)
      CALL HBOOK1(102,'Proton v2sq     ',50,0.,2.,0.)
      CALL HBOOK1(103,'Proton 1-xi1    ',50,0.,0.5,0.)
      CALL HBOOK1(104,'Proton 1-xi2    ',50,0.,0.5,0.)
      CALL HBOOK1(105,'Gluon xg1       ',51,0.,1.02,0.)
      CALL HBOOK1(106,'Gluon xg2       ',51,0.,1.02,0.)
      CALL HBOOK1(107,'Central Mass    ',200,0.,20.,0.)
      CALL HBOOK1(108,'Dijet mass      ', 10,50.,100.,0.)
      CALL HBOOK1(109,'Mass fraction   ',51,0.,1.02,0.)
      CALL HBOOK1(110,'Rapidity y1     ',100,-5.,5.,0.)
      CALL HBOOK1(111,'Rapidity y2     ',100,-5.,5.,0.)
      CALL HBOOK1(112,'Trans. mom pt1  ',50,0.,200.,0.)
      CALL HBOOK1(113,'Trans. mom pt2  ',50,0.,200.,0.)
      CALL HBOOK2(114,'Proton (1-xi1)(1-xi2)',50,0.9,1.,50,0.9,1.,0)
      CALL HBOOK2(115,'Gluons (xg1)(xg2)',100,0.,1.,100,0.,1.,0)
      CALL HBOOK2(116,'Partons (xg1)(y1)',100,0.,1.,100,-10.,10.,0)
      CALL HBOOK2(117,'Partons (xg2)(y2)',100,0.,1.,100,-10.,10.,0)
      CALL HBOOK1(118,'Rapidity y0     ',100,-5.,5.,0.)
      CALL HBOOK1(119,'Trans. mom pt0  ',50,0.,5.,0.)
      CALL HBOOK1(301,'Proton phi1',100,-3.2,3.2,0.)
      CALL HBOOK1(302,'Proton phi2',100,-3.2,3.2,0.)
      CALL HBOOK1(303,'Proton dphi',90,0.,180.,0.)

      RETURN
      END
C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.                                             
C * User's routine to analyze data from event.
C * Distributions relevant for DPE scattering. Assumptions:
C    - initial protons are index 1 and 2
C    - scattered protons are index 5 and 7 (from 1 and 2 resp.)
C    - pomerons are 4 and 6 (radiated from index 1 and 2 resp.)
C    - gluons are 8 and 9 (from index 4 and 6 resp.)
C    - Higgs/Z/... is 10
C    - di-parton is 11 and 12           
C----------------------------------------------------------------------- 
      SUBROUTINE HWANAL
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'

      REAL V1SQ,V2SQ,XI1,XI2,OXI1,OXI2,WEIGHT,XG1,XG2
      REAL ROOTS,CENTM,HARDM,MFRAC
      REAL PT0,Y0,PT1,MT1SQ,Y1,PT2,MT2SQ,Y2,PTH,MTHSQ,YH
      REAL PHI1,PHI2,DPHI

      IF(IERROR.NE.0) RETURN
      
      WEIGHT=1d0

      V1SQ=PHEP(1,5)**2+PHEP(2,5)**2
      V2SQ=PHEP(1,7)**2+PHEP(2,7)**2
      XI1=1.-ABS(PHEP(3,5)/PHEP(3,1))
      XI2=1.-ABS(PHEP(3,7)/PHEP(3,2))
      OXI1=1.-XI1
      OXI2=1.-XI2

      PHI1 = DACOS(PHEP(1,5)/SQRT(V1SQ))
      IF(PHEP(2,5).LT.0d0) PHI1 = -PHI1
      PHI2 = DACOS(PHEP(1,7)/SQRT(V2SQ))
      IF(PHEP(2,7).LT.0d0) PHI2 = -PHI2
      DPHI=PHI2-PHI1
      IF(DPHI.LT.-PI) DPHI=DPHI+2*PI
      IF(DPHI.GT.+PI) DPHI=DPHI-2*PI
      DPHI=ABS(DPHI)*180/3.14159265359

      XG1=ABS(PHEP(4,8)/PHEP(4,4))
      XG2=ABS(PHEP(4,9)/PHEP(4,6))

      ROOTS=PHEP(4,1)+PHEP(4,2)
      CENTM=ROOTS*SQRT(XI1*XI2)
      HARDM=CENTM*SQRT(XG1*XG2)
      MFRAC=HARDM/CENTM

      PT0=SQRT(PHEP(1,10)**2+PHEP(2,10)**2)
      Y0=0.5*LOG((PHEP(4,10)+PHEP(3,10))/(PHEP(4,10)-PHEP(3,10)))

      PT1=SQRT(PHEP(1,11)**2+PHEP(2,11)**2)
      Y1=0.5*LOG((PHEP(4,11)+PHEP(3,11))/(PHEP(4,11)-PHEP(3,11)))

      PT2=SQRT(PHEP(1,12)**2+PHEP(2,12)**2)
      Y2=0.5*LOG((PHEP(4,12)+PHEP(3,12))/(PHEP(4,12)-PHEP(3,12)))

      CALL HF1(101,V1SQ,WEIGHT)
      CALL HF1(102,V2SQ,WEIGHT)
      CALL HF1(103,XI1,WEIGHT)
      CALL HF1(104,XI2,WEIGHT)
      CALL HF1(105,XG1,WEIGHT)
      CALL HF1(106,XG2,WEIGHT)
      CALL HF1(107,CENTM,WEIGHT)
      CALL HF1(108,HARDM,WEIGHT)
      CALL HF1(109,MFRAC,WEIGHT)
      CALL HF1(110,Y1,WEIGHT)
      CALL HF1(111,Y2,WEIGHT)
      CALL HF1(112,PT1,WEIGHT)
      CALL HF1(113,PT2,WEIGHT)
      CALL HF2(114,OXI1,OXI2,WEIGHT)
      CALL HF2(115,XG1,XG2,WEIGHT)
      CALL HF2(116,XG1,Y1,WEIGHT)
      CALL HF2(117,XG2,Y2,WEIGHT)
      CALL HF1(118,Y0,WEIGHT)
      CALL HF1(119,PT0,WEIGHT)
      CALL HF1(301,PHI1,WEIGHT)
      CALL HF1(302,PHI2,WEIGHT)
      CALL HF1(303,DPHI,WEIGHT)

      RETURN
      END
C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.                                             
C * User's routine for terminal calculations, histogram output, etc.           
C----------------------------------------------------------------------- 
      SUBROUTINE HWAEND
c---Finish HBOOK
      CALL HROUT(0,ICYCLE,' ')
      CALL HRENDC('HWIG')
      CLOSE(33)
      RETURN
      END
C-----------------------------------------------------------------------
