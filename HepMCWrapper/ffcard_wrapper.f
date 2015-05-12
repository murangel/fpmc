C-----------------------------------------------------------------------
C * 10/06/2008 O. Kepka, oldrich.kepka@cern.ch
C * Routine initialize FPMC parameters to default
C * if REACARD = .TRUE. it allows to read setup from external format
C * free card (FFC) with
C * ./module < datacard.ffcard
C----------------------------------------------------------------------- 

      SUBROUTINE FPMC_VAR_INI(READCARD)
      IMPLICIT NONE
      INCLUDE 'ffcard_wrapper.inc'

      INTEGER SPACE(1000)
      COMMON/CFREAD/SPACE
      LOGICAL READCARD

      CALL FPMC_WELCOME

c---FFC default initialiszation
      UOUTPUT    =    1
      UNTNAME    =    'tmpntuple.ntp'
      UMAXEV     =    1000
      UTYPEPR    =    'EXC'
      UTYPINT    =    'QED'
      UHADR      =    'Y'
      UPART1     =    'E+'
      UPART2     =    'E+'
      UIPROC     =    16010
      UHMASS     =    120.
      UNFLUX     =    15
      UNRN1      =    33799
      UNRN2      =    11799
      UWMASS     =    80.425
      UTMASS     =    174.3
      UMST1      =    250.
      UMSB1      =    250.
      UECMS      =    14.E3
      UPTMIN     =    0.
      UPTMAX     =    1.E8
      UYJMIN     =    -6.
      UYJMAX     =    6.
      UEMMIN     =    10.
      UEMMAX     =    1.E8
      UIFIT      =    10
      UISOFTM    =    1
      UZION      =    1
      UAION      =    1
      UBMIN      =    1.
      UAAANOM  =    0       ! original HERWIG AAWW formula
      UDKAPPA    =    0.
      UDLAMBDA   =    0.
      UA0W       =    0.
      UACW       =    0.
      UA0Z       =    0.
      UACZ       =    0.
      UA1A       =    0.
      UA2A       =    0.
      UANOMCUTOFF =  -1.   ! if < 0 turns formfactor off
      UAAEXOTIC  =    0    ! SM Excl AAAA
      UAAM       =    0.
      UAAQ       =    0.
      UAAN       =    0.
      UQ2WWMN    =    0.
      UQ2WWMX    =    4.
      UYWWMIN    =    0.
      UYWWMAX    =    0.1
      UCHIDeIGLU = -1
      UCHIDeX    = -1.0
      UCHIDeXP   = -1.0
      UCHIDeS2   = -1.0
      UXI1MIN= -1.0 
      UXI1MAX= -1.0 
      UXI2MIN= -1.0 
      UXI2MAX= -1.0 
      UKMR2Q2CUT=2.0
      UKMR2SURV=0.03
      UKMR2SCALE=0.618!1.0
      UKMR2DELTA=1!2
      UCHIDePATH='External/CHIDe/Data/'
      UMODPDF1   = -1
      UMODPDF2   = -1


      IF(READCARD) THEN
        !---Key reading      
        PRINT *, "Reading datacard ..."
        call FFINIT(1000)
        !---JC read config file
        CALL FFSET('SIZE', 32)
        call FFKEY('OUTPUT',UOUTPUT,1,'integer')
        call FFKEY('MAXEV',UMAXEV,1,'integer')
        call FFKEY('TYPEPR',UTYPEPR,3,'mixed')
        call FFKEY('TYPINT',UTYPINT,3,'mixed')
        call FFKEY('HADR',UHADR,1,'mixed')
        call FFKEY('PART1',UPART1,1,'mixed')
        call FFKEY('PART2',UPART2,1,'mixed')
        call FFKEY('ECMS',UECMS,1,'real')
        call FFKEY('IPROC',UIPROC,1,'integer')
        call FFKEY('NFLUX',UNFLUX,1,'integer')
        call FFKEY('NRN1',UNRN1,1,'integer')
        call FFKEY('NRN2',UNRN2,1,'integer')
        call FFKEY('IFIT',UIFIT,1,'integer')
        call FFKEY('IZION',UZION,1,'integer')
        call FFKEY('BMIN',UBMIN,1,'real')
        call FFKEY('IAION',UAION,1,'integer')
        call FFKEY('WMASS',UWMASS,1,'real')
        call FFKEY('TMASS',UTMASS,1,'real')
        call FFKEY('HMASS',UHMASS,1,'real')
        call FFKEY('MSTOP1',UMST1,1,'real')
        call FFKEY('MSBOT1',UMSB1,1,'real')
        call FFKEY('YJMAX',UYJMAX,1,'real')
        call FFKEY('YJMIN',UYJMIN,1,'real')
        call FFKEY('PTMIN',UPTMIN,1,'real')
        call FFKEY('PTMAX',UPTMAX,1,'real')
        call FFKEY('EMMIN',UEMMIN,1,'real')
        call FFKEY('EMMAX',UEMMAX,1,'real')
        call FFKEY('Q2WWMN',UQ2WWMN,1,'real')
        call FFKEY('Q2WWMX',UQ2WWMX,1,'real')
        call FFKEY('YWWMIN',UYWWMIN,1,'real')
        call FFKEY('YWWMAX',UYWWMAX,1,'real')
        call FFKEY('ISOFTM',UISOFTM,1,'integer')
        call FFKEY('AAANOM', UAAANOM,1,'integer')
        call FFKEY('DKAPPA', UDKAPPA,1,'real')
        call FFKEY('LAMBDA', UDLAMBDA,1,'real')
        call FFKEY('A0W', UA0W,1,'real')
        call FFKEY('ACW', UACW,1,'real')
        call FFKEY('A0Z', UA0Z,1,'real')
        call FFKEY('ACZ', UACZ,1,'real')
        call FFKEY('A1A', UA1A,1,'real')
        call FFKEY('A2A', UA2A,1,'real')
        call FFKEY('ANOMCUTOFF', UANOMCUTOFF,1,'real')
        call FFKEY('AAEXOTIC', UAAEXOTIC,1,'integer')
        call FFKEY('AAM', UAAM,1,'real')
        call FFKEY('AAQ', UAAQ,1,'real')
        call FFKEY('AAN', UAAN,1,'real')
        UNTNAME = ''
        call FFKEY('NTNAME',UNTNAME,32,'mixed')
        call FFKEY('CHIDePATH', UCHIDePATH, 32, 'mixed')
        call FFKEY('IGLU', UCHIDeIGLU, 1, 'integer')
        call FFKEY('USCALE',  UCHIDeX, 1, 'real')
        call FFKEY('LSCALE', UCHIDeXp, 1, 'real')
        call FFKEY('SURV', UCHIDeS2, 1, 'real')
        call FFKEY('XI1Min', UXI1Min, 1, 'real')
        call FFKEY('XI1Max', UXI1Max, 1, 'real')
        call FFKEY('XI2Min', UXI2Min, 1, 'real')
        call FFKEY('XI2Max', UXI2Max, 1, 'real')
        call FFKEY('CHIDeGapMin', UCHIDeGapMin, 1, 'real')
        call FFKEY('CHIDeGapMax', UCHIDeGapMax, 1, 'real')
        call FFKEY('Q2CUT', UKMR2Q2CUT, 1, 'real')
        call FFKEY('SCALE', UKMR2SCALE, 1, 'real')
        call FFKEY('DELTA', UKMR2DELTA, 1, 'integer')
        call FFKEY('MODPDF1', UMODPDF1, 1, 'integer')
        call FFKEY('MODPDF2', UMODPDF2, 1, 'integer')

        call FFGO
       
C    Xi cuts
          IF(UXI1Min.LT.0.0) UXI1Min = UYWWMIN
          IF(UXI1Max.LT.0.0) UXI1Max = UYWWMAX
          IF(UXI2Min.LT.0.0) UXI2Min = UYWWMIN
          IF(UXI2Max.LT.0.0) UXI2Max = UYWWMAX

C     CHIDe Model
        IF(UNFLUX.EQ.18) THEN          
C   Impact factor parameterisation
          IF(UCHIDeIGLU.LT.0.0) UCHIDeIGLU = 4
C   Scaling lower limit of Sudakov factor integration
          IF(UCHIDeXP.LT.0.0.AND.(.NOT.(UIPROC.EQ.16059)))UCHIDeXP=0.5
          IF(UCHIDeXP.LT.0.0 .AND. UIPROC.EQ.16059) UCHIDeXP = 1.0
C   Gap Survival Probability
          IF(UCHIDeS2.LT.0.0) THEN 
            IF(abs(UECMS-14000.).lt.1.) THEN
              UCHIDeS2 = 0.03 ! ~LHC
            ELSEIF(abs(UECMS-1960.).lt.1) THEN
              UCHIDeS2 = 0.10  ! ~Tevatron
            ELSE
              WRITE (*,*) "ERROR: CHIDeS2 not specified!",
     &  "  For nonstandard ECMS (different than 1960 and 14000)",
     &  "  CHIDeS2 must be specified in data card."
              WRITE (*,*) "STOP"
              STOP
            ENDIF
          ENDIF
          IF(UIPROC.EQ.16012) THEN ! Scaling upper limit of Sudakov factor integration:
            IF(UCHIDeX.LT.0.0) UCHIDeX = 0.5      
          ELSEIF(UIPROC.EQ.16059) THEN
            IF(UCHIDeX.LT.0.0) UCHIDeX = 0.5
          ELSEIF(UIPROC.EQ.19999) THEN
            IF(UCHIDeX.LT.0.0) UCHIDeX = 1    
          ENDIF
       ENDIF
       
       IF(UNFLUX.EQ.16) THEN 
         IF(UCHIDeS2.GT.0) UKMR2SURV = UCHIDeS2
         IF(UCHIDeS2.LT.0) UKMR2SURV = 0.03  
       ENDIF

      ELSE
        PRINT *, "Reading datacard not requested, using default",
     & " parameters."
      ENDIF    

      END SUBROUTINE

      SUBROUTINE FPMC_WELCOME
      WRITE (*,*) ""
      WRITE (*,*) ""
      WRITE (*,*) ""
      WRITE (*,*) "     _|_|_|_|  _|_|_|    _|      _|    _|_|_| " 
      WRITE (*,*) "     _|        _|    _|  _|_|  _|_|  _|       " 
      WRITE (*,*) "     _|_|_|    _|_|_|    _|  _|  _|  _|       " 
      WRITE (*,*) "     _|        _|        _|      _|  _|       " 
      WRITE (*,*) "     _|        _|        _|      _|    _|_|_| " 
      WRITE (*,*) ""
      WRITE (*,*) ""
      WRITE (*,*) "     FPMC - Forward Physics Monte Carlo v1.0",
     &   " 12 Jun 2012"
      WRITE (*,*) ""
      WRITE (*,*) "     www.cern.ch/fpmc"
      WRITE (*,*) ""
      WRITE (*,*) ""

      END SUBROUTINE
