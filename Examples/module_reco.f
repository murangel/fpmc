C-----------------------------------------------------------------------
C * 15/11/2006 O. Kepka
C * Multifunctional module with cms fast simulation
C * Set parameters via FF cards      
C----------------------------------------------------------------------- 
      PROGRAM FPMC
c---Common block
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'
      INCLUDE 'CHIDe.inc'
c---User's declarations
      DOUBLE PRECISION CMSENR
      INTEGER N
      CHARACTER ANSWER
      LOGICAL READCARD
c---External
      INCLUDE 'ffcard.inc'

c---Initialize FPMC setup parameters and read ffcards if READCARD = true
      READCARD = .TRUE.
      CALL FPMC_VAR_INI(READCARD) 

c---Apply the setup parameters - either default/ffcard

c---Events in the run
      MAXEV = UMAXEV

c---Beam particles (In DPE as implemented by POMWIG E+ stands for initial state proton)
      PART1=UPART1
      PART2=UPART2

c---Set cms energy
      CMSENR = UECMS    !set by a card

c---Beam momenta
      PBEAM1=CMSENR/2d0
      PBEAM2=PBEAM1

c---Process type : EXClusive or INClusive, QCD or QED when relevant
c   It is important to set those right
c      TYPEPR='INC'
c      TYPINT='QCD'
      TYPEPR = UTYPEPR
      TYPINT = UTYPINT

c---Setting the hard subprocess (see manual)
      IPROC = UIPROC

c---Option to include hadronization and showering effects
c      ANSWER='Y' ! say 'N' to skip this part
      ANSWER=UHADR 

c---Initialization of other Herwig common blocks
c...The default parameters has to be changed after this call
      CALL HWIGIN
cc      SUSYIN=.TRUE.

c---Random number genrator initializaiton      
      NRN(1) = UNRN1 ! set again later in the code
      NRN(2) = UNRN2 ! set again later in the code

C ... increase inefficiency tolerance
      EFFMIN = 1d-6

c---User's default kinematic parameters
c... beam momentum transfer range ( Q2 = |t| )
      RMASS(201) = UHMASS
      RMASS(6)   = UTMASS ! top mass
      RMASS(198) = UWMASS ! W mass
      RMASS(406) = UMST1 ! stop1 mass
      RMASS(405) = UMSB1 ! stop2 mass

c stop mass
c      RMASS(449)=10000.D0
c      RMASS(401)=10000.D0
c      RMASS(402)=10000.D0
c      RMASS(403)=10000.D0
c      RMASS(404)=10000.D0   
c      RMASS(405)=10000.D0
c      RMASS(407)=10000.D0
c      RMASS(408)=10000.D0
c      RMASS(409)=10000.D0  
c      RMASS(410)=10000.D0
c      RMASS(411)=10000.D0

c      RMASS(413)=10000.D0
c      RMASS(414)=10000.D0
c      RMASS(415)=10000.D0  
c      RMASS(416)=10000.D0
c      RMASS(417)=10000.D0
c      RMASS(418)=10000.D0
c      RMASS(419)=10000.D0
c      RMASS(420)=10000.D0
c      RMASS(421)=10000.D0
c      RMASS(422)=10000.D0
c      RMASS(423)=10000.D0
c      RMASS(424)=10000.D0
c      RMASS(406)=393.D0
c      RMASS(412)=393.D0

      Q2WWMN=UQ2WWMN
      Q2WWMX=UQ2WWMX

c      Q2WWMN=0
c      Q2WWMX=4
      
c... beam momentum loss range
c      YWWMIN=0.d0
c      YWWMAX=1.D0
      YWWMIN=UYWWMIN
      YWWMAX=UYWWMAX
      
c--- central products : rapidity, pT or mass cuts
c      YJMAX = 5d0
c---JC      YJMAX = 4d0
      YJMAX = UYJMAX
      YJMIN = UYJMIN

c      PTMIN = 25d0
c---JC      
      PTMIN=UPTMIN
      PTMAX=UPTMAX

c masse des leptons / pas photons, pas jets
c---JC      EMMIN=0d0
      EMMIN=UEMMIN

c---Choosing the flux : 9,10 is Cox-Forshaw (DPE)
c                       11 is Bialas-Landshoff (exc. DPE) or Boonekamp et al (inc. DPE)
c                       12 is Cahn, Jackson (heavy ions, QED)
c                       13 is Drees, Ellis, Zeppenfled (heavy ions, QED)
c                       14 is Papageorgiu (proton, QED)
c                       15 is Budnev flux (RECOMMENDED for proton, QED)
c                       16 is KMR flux (tables from L.Lonnblad or ExHume)
      NFLUX = UNFLUX

c---Choosing pdf to use: 
c                       according to hep/ph 0609291
c                       10:  H1      
c                       20:  Zeus 
c                       30:  H1Zeus combined 
c                       older versions for compat., see h1qcd.f   
c                        2:  NLO fit as in H1
c                        5:  LO  fit as in H1
c                        8:  extended version of H1
      IFITPDF = UIFIT

c ---Anomalous coupling settings      
      AAANOM = UAAANOM
      D_KAPPA = UDKAPPA
      LAMBDA = UDLAMBDA
      A0W = UA0W
      ACW = UACW
      A0Z = UA0Z
      ACZ = UACZ
      A1A = UA1A
      A2A = UA2A
      ANOMCUTOFF = UANOMCUTOFF

c ---Exotic excl AA settings
      AAEXOTIC = UAAEXOTIC
      AAM = UAAM
      AAQ = UAAQ
      AAN = UAAN

C ... begin R.S.
C     CHIDe Model
      CHIDeIGLU = UCHIDeIGLU
      CHIDeX   =  UCHIDeX
      CHIDeXP  =  UCHIDeXP
      CHIDeS2  =  UCHIDeS2
      CHIDeS   =  UECMS*UECMS
      XI1MIN = UXI1MIN
      XI1MAX = UXI1MAX
      XI2MIN = UXI2MIN
      XI2MAX = UXI2MAX
      CHIDeGapMin = UCHIDeGapMin
      CHIDeGapMax = UCHIDeGapMax
      CHIDePATH = UCHIDePATH

C ... end R.S.


c---Set normalisations for Exclusive processes
c... option for soft corrections : ISOFTM = 0 : no correction
c                                           1 : simple factor, see below
c                                           2 : KMR low mass diffractive
c                                           3 : effective opacity model
      ISOFTM = UISOFTM

c---Events printed
      MAXPR=10

c---Set the number of nucleons - for QED photon flux
      ZION = UZION
      AION = UAION
      RBMIN - URBMIN

      write(*,*) ''
      write(*,*) 'USER SETTINGS'
      write(*,*) '-------------'
      write(*,*) 'NTNAME   = ',UNTNAME
      write(*,*) 'MAXEV    = ',MAXEV
      write(*,*) 'TYPEPR   = ',TYPEPR
      write(*,*) 'TYPEINT  = ',TYPINT
      write(*,*) 'HADR     = ', ANSWER
      write(*,*) 'PART1    = ', PART1
      write(*,*) 'PART2    = ', PART2
      write(*,*) 'NRN(1)   = ',NRN(1)
      write(*,*) 'NRN(2)   = ',NRN(2)
      write(*,*) 'WMASS    = ',UWMASS
      write(*,*) 'TMASS    = ',UTMASS
      write(*,*) 'HMASS    = ',UHMASS
      write(*,*) 'MSTOP1   = ',UMST1
      write(*,*) 'MSBOY1   = ',UMSB1
      write(*,*) 'YJMIN    = ', YJMIN
      write(*,*) 'YJMAX    = ',YJMAX
      write(*,*) 'PTMIN    = ',PTMIN
      write(*,*) 'EMMIN    = ',EMMIN
      write(*,*) 'IFIT     = ',IFITPDF
      write(*,*) 'ISOFTM   = ',ISOFTM
      write(*,*) 'Q2WWMN   = ',Q2WWMN
      write(*,*) 'Q2WWMX   = ',Q2WWMX
      write(*,*) 'YWWMIN   = ',YWWMIN
      write(*,*) 'YWWMAX   = ',YWWMAX
      write(*,*) 'ZION     = ',ZION
      write(*,*) 'AION     = ',AION
      write(*,*) 'BMIN     = ',RBMIN
      write(*,*) 'AAANOM = ', AAANOM
      write(*,*) 'DKAPPA = ', D_KAPPA
      write(*,*) 'DLAMBDA = ', LAMBDA
      write(*,*) 'ANOMCUTOFF= ', ANOMCUTOFF
      write(*,*) 'AAEXOTIC = ', AAEXOTIC
      write(*,*)

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
      CALL HWFXER(.TRUE.,IPROC)

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
         CALL HWFXER(.FALSE.,IPROC)
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
      INCLUDE 'ffcard.inc'
!      COMMON/ QUEST / IQUEST(10)   !pridano
!      COMMON/ QUEST / IQUEST(10)   !pridano
      INTEGER NWPAWC
      REAL*4 HMEMOR
      PARAMETER(NWPAWC = 10000000)
      COMMON /PAWC/HMEMOR(NWPAWC)
      INTEGER IQUEST
      COMMON /QUEST/ IQUEST(100)


c---HBOOK initialization
      CALL HLIMIT(NWPAWC)
      CALL HCDIR('//PAWC',' ')

      ! for high capacity hbook
      IQUEST(10) = 64000  

c     UNTNAME is adopted from the key
      CALL HROPEN(33,'HWIG',UNTNAME,'NQ',8192,ISTAT)
      ! for low capacitiy hbook
c      CALL HROPEN(33,'HWIG',UNTNAME,'NC',4096,ISTAT)
      
c---Histograms
      CALL HBOOK1(101,'Proton v1sq     ',50,0.,2.,0.)
      CALL HBOOK1(102,'Proton v2sq     ',50,0.,2.,0.)
      CALL HBOOK1(103,'Proton 1-xi1    ',50,0.,0.8,0.)
      CALL HBOOK1(104,'Proton 1-xi2    ',50,0.,0.8,0.)
      CALL HBOOK1(105,'Gluon xg1       ',51,0.,1.02,0.)
      CALL HBOOK1(106,'Gluon xg2       ',51,0.,1.02,0.)
      CALL HBOOK1(107,'Central Mass    ',200,0.,2000.,0.)
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
      CALL HBOOK1(901,'N events        ',1,0.,1.,0.)
      CALL HBOOK1(902,'E cms           ',1,0.,1.,0.)
      CALL HBOOK1(903,'IPROC           ',1,0.,1.,0.)
      CALL HBOOK1(904,'NFLUX           ',1,0.,1.,0.)
      CALL HBOOK1(905,'YJMAX           ',1,0.,1.,0.)
      CALL HBOOK1(906,'PTMIN           ',1,0.,1.,0.)
      CALL HBOOK1(907,'AAANOM        ',1,0.,1.,0.)
      CALL HBOOK1(908,'D_KAPPA         ',1,0.,1.,0.)
      CALL HBOOK1(909,'LAMBDA        ',1,0.,1.,0.)
      CALL HBOOK1(910,'ANOMCUTOFF      ',1,0.,1.,0.)
      CALL HBOOK1(911,'AAEXOTIC       ',1,0.,1.,0.)
      CALL HBOOK1(999,'XSECT[pb]       ',1,0.,1.,0.)

      ! save some settings into the ntuple
      CALL HF1(901, 0.5, REAL(UMAXEV))
      CALL HF1(902, 0.5, REAL(UECMS))
      CALL HF1(903, 0.5, REAL(UIPROC))
      CALL HF1(904, 0.5, REAL(UNFLUX))
      CALL HF1(905, 0.5, REAL(UYJMAX))
      CALL HF1(906, 0.5, REAL(UPTMIN))
      CALL HF1(907, 0.5, REAL(UAAANOM))
      CALL HF1(908, 0.5, REAL(UDKAPPA))
      CALL HF1(909, 0.5, REAL(UDLAMBDA))
      CALL HF1(911, 0.5, REAL(UAAEXOTIC))

c  simulation ntuple initialisation
      CALL NTINIT

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

      integer IERR
      REAL V1SQ,V2SQ,XI1,XI2,OXI1,OXI2,WEIGHT,XG1,XG2
      REAL ROOTS,CENTM,HARDM,MFRAC
      REAL PT0,Y0,PT1,MT1SQ,Y1,PT2,MT2SQ,Y2,PTH,MTHSQ,YH
      REAL PHI1,PHI2,DPHI
      real xg1b, xg2b
      common /remnant/ xg1b,xg2b

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

      xg1b=xg1
      xg2b=xg2

      ROOTS=2*SQRT(PHEP(4,1)*PHEP(4,2))
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

c      print *,'central mass', sqrt(xi1*xi2)*ROOTS

c --- my reconstruction is called
      CALL MYRECO(IERR)
      IF(IERR.NE.0) return

      call hfnt(777)

      RETURN
      END
C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.                                             
C * User's routine for terminal calculations, histogram output, etc.           
C----------------------------------------------------------------------- 
      SUBROUTINE HWAEND
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'
      INCLUDE 'ffcard.inc'
      INTEGER ICYCLE

      ICYCLE = 2

      ! save the cross section
      CALL HF1(999, 0.5, REAL(1000.*AVWGT) )

      PRINT *, ''
      PRINT *, '========Final summary========='
      PRINT *, 'Cross section[pb]=', 1000.*AVWGT , ' NEVENTS=', UMAXEV

c---Finish HBOOK
      CALL HROUT(0,ICYCLE,' ')
      CALL HRENDC('HWIG')
      CLOSE(33)
      RETURN
      END
C-----------------------------------------------------------------------
