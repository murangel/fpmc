C-----------------------------------------------------------------------
C * 09/13/2003 Tibor Kucs                                             
C * This example : inclusive DPE dijet cross section (all quarks and gluons).
C----------------------------------------------------------------------- 
      PROGRAM DPECS
c---Common block
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'FPMC.INC'
c---User's declarations
      INTEGER N
      CHARACTER ANSWER
c---External
      EXTERNAL HWUDAT

c---Events in the run
      MAXEV=10000

c---Beam particles (In DPE as implemented by POMWIG E+ stands for initial state proton)
      PART1='E+'    ! e --> e+gam will mimic p --> p+Pom
      PART2='E+'    ! e --> e+gam will mimic p --> p+Pom

c---Beam momenta
      PBEAM1=7d3
      PBEAM2=7d3

c---Process type : EXClusive or INClusive, QCD or QED when relevant
c   It's important to set those right
      TYPEPR='INC'
      TYPINT='QCD'

c---Setting the hard subprocess (see manual)
*      IPROC=11500 ! dijets
      IPROC=12800 ! W-pairs

c---Option to include hadronization and showering effects
      ANSWER='Y' ! say 'N' to skip this part

c---Initialization of other common blocks
c...The default parameters has to be changed after this call
      CALL HWIGIN

c---User's default kinematic parameters
c... beam momentum transfer range ( Q2 = |t| )
      Q2WWMN=0d0
      Q2WWMX=4d00
c... beam momentum loss range
      YWWMIN=0d0
      YWWMAX=0.9D0
c... central products : rapidity, pT or mass cuts
      YJMAX = 5d0
      PTMIN = 5d0

c---Choosing the flux : 9,10 is Cox-Forshaw (DPE)
c                       11 is Bialas-Landshoff (exc. DPE) or Boonekamp et al (inc. DPE)
c                       12 is Cahn, Jackson (heavy ions, QED)
c                       13 is Drees, Ellis, Zeppenfled (heavy ions, QED)
c                       14 is Papageorgiu (proton, QED)
c                       15 is Budnev flux (RECOMMENDED for proton, QED)
c                       16 is KMR flux (tables from L.Lonnblad or ExHume)
      NFLUX=9

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


c---Set normalisations for DPE processes - Irrelevant for SD, DDIS
c... option for soft corrections : ISOFTM = 0 : no correction
c                                           1 : simple factor, see below
c                                           2 : KKMR low mass diffractive
c                                           3 : effective opacity model
      ISOFTM = 1

c---Events printed
      MAXPR=1
      
c---Initialize model/pdf dependant parameters
      CALL HWMODINI

c---User's initial calculations
      CALL HWABEG

c---Compute parameter dependent constants
      CALL HWUINC

c---Check FPEMC Settings + Initialisations for consistency
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
         print *, N
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
      INCLUDE 'FPMC.INC'

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

      RETURN
      END
C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.                                             
C * User's routine for terminal calculations, histogram output, etc.           
C----------------------------------------------------------------------- 
      SUBROUTINE HWAEND

      RETURN
      END
C-----------------------------------------------------------------------
