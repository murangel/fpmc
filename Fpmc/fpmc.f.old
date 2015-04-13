
CDECK  ID>, HWDHAD.
*CMZ :-        -26/04/91  11.11.54  by  Peter Richardson
*-- Author :    Ian Knowles, Bryan Webber & Mike Seymour
C-----------------------------------------------------------------------
      SUBROUTINE HWDHAD
C-----------------------------------------------------------------------
C     GENERATES DECAYS OF UNSTABLE HADRONS AND LEPTONS
C     Modified for TAUOLA interface 16/10/01 PR
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      COMMON/FFS/TB,BT
      COMMON/SFF/IT1,IB1,IT2,IB2
      DOUBLE PRECISION TB,BT
      INTEGER IT1,IB1,IT2,IB2
      DOUBLE PRECISION HWRGEN,HWULDO,RN,BF,COSANG,RSUM,DIST(4),VERTX(4),
     & PMIX,WTMX,WTMX2,XS,DOT1,DOT2,HWDPWT,HWDWWT,HWDHWT,XXX,YYY
      INTEGER IHEP,ID,MHEP,IDM,I,IDS,IM,MO,IPDG
      LOGICAL STABLE
      EXTERNAL HWRGEN,HWDPWT,HWDWWT,HWDHWT,HWULDO
c O.K. INEG counts number of cycles. Kills event if > 50       
      INTEGER INEG
      IF (IERROR.NE.0) RETURN
      DO 100 IHEP=1,NMXHEP
      IF (IHEP.GT.NHEP) THEN
        ISTAT=90
        RETURN
      ELSEIF (ISTHEP(IHEP).EQ.120 .AND.
     &  JDAHEP(1,IHEP).EQ.IHEP.AND.JDAHEP(2,IHEP).EQ.IHEP) THEN
C---COPY COLOUR SINGLET CMF
        NHEP=NHEP+1
        IF (NHEP.GT.NMXHEP) CALL HWWARN('HWDHAD',100,*999)
        CALL HWVEQU(5,PHEP(1,IHEP),PHEP(1,NHEP))
        CALL HWVEQU(4,VHEP(1,IHEP),VHEP(1,NHEP))
        IDHW(NHEP)=IDHW(IHEP)
        IDHEP(NHEP)=IDHEP(IHEP)
        ISTHEP(NHEP)=190
        JMOHEP(1,NHEP)=IHEP
        JMOHEP(2,NHEP)=NHEP
        JDAHEP(2,NHEP)=NHEP
        JDAHEP(1,IHEP)=NHEP
        JDAHEP(2,IHEP)=NHEP
      ELSEIF (ISTHEP(IHEP).GE.190.AND.ISTHEP(IHEP).LE.193) THEN
C---FIRST CHECK FOR STABILITY
        ID=IDHW(IHEP)
        IF (RSTAB(ID)) THEN
          ISTHEP(IHEP)=1
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C---SPECIAL FOR GAUGE BOSON DECAY
          IF (ID.GE.198.AND.ID.LE.200) CALL HWDBOS(IHEP)
C---SPECIAL FOR HIGGS BOSON DECAY
          IF (ID.EQ.201) CALL HWDHIG(ZERO)
        ELSE
C---UNSTABLE.
C Calculate position of decay vertex
          IF (DKLTM(ID).EQ.ZERO) THEN
            CALL HWVEQU(4,VHEP(1,IHEP),VERTX)
            MHEP=IHEP
            IDM=ID
          ELSE
            CALL HWUDKL(ID,PHEP(1,IHEP),DIST)
            CALL HWVSUM(4,VHEP(1,IHEP),DIST,VERTX)
            IF (MAXDKL) THEN
              CALL HWDXLM(VERTX,STABLE)
              IF (STABLE) THEN
                ISTHEP(IHEP)=1
                JDAHEP(1,IHEP)=0
                JDAHEP(2,IHEP)=0
                GOTO 100
              ENDIF
            ENDIF
            IF (MIXING.AND.(ID.EQ.221.OR.ID.EQ.223.OR.
     &                      ID.EQ.245.OR.ID.EQ.247)) THEN
C Select flavour of decaying b-meson allowing for flavour oscillation
              IDS=MOD(ID,3)
              XXX=XMRCT(IDS)*DIST(4)/PHEP(4,IHEP)
              YYY=YMRCT(IDS)*DIST(4)/PHEP(4,IHEP)
              IF (ABS(YYY).LT.10) THEN
                PMIX=HALF*(ONE-COS(XXX)/COSH(YYY))
              ELSE
                PMIX=HALF
              ENDIF
              IF (HWRGEN(1).LE.PMIX) THEN
                IF (ID.LE.223) THEN
                  IDM=ID+24
                ELSE
                  IDM=ID-24
                ENDIF
              ELSE
                IDM=ID
              ENDIF
C Introduce a decaying neutral b-meson
              IF (NHEP+1.GT.NMXHEP) CALL HWWARN('HWDHAD',101,*999)
              MHEP=NHEP+1
              ISTHEP(MHEP)=ISTHEP(IHEP)
              ISTHEP(IHEP)=200
              JDAHEP(1,IHEP)=MHEP
              JDAHEP(2,IHEP)=MHEP
              IDHW(MHEP)=IDM
              IDHEP(MHEP)=IDPDG(IDM)
              JMOHEP(1,MHEP)=IHEP
              JMOHEP(2,MHEP)=JMOHEP(2,IHEP)
              CALL HWVEQU(5,PHEP(1,IHEP),PHEP(1,MHEP))
              CALL HWVEQU(4,VERTX,VHEP(1,MHEP))
              NHEP=NHEP+1
            ELSE
              MHEP=IHEP
              IDM=ID
            ENDIF
          ENDIF
C Use CLEO/EURODEC packages for b-hadrons if requested
          IF ((IDM.GE.221.AND.IDM.LE.231).OR.
     &        (IDM.GE.245.AND.IDM.LE.254)) THEN
            IF (BDECAY.EQ.'CLEO') THEN
              CALL HWDCLE(MHEP)
              GOTO 100
            ELSEIF (BDECAY.EQ.'EURO') THEN
              CALL HWDEUR(MHEP)
              GOTO 100
            ENDIF
          ENDIF
C Use TAUOLA package for tau decays if requested
          IF((IDM.EQ.125.OR.IDM.EQ.131).AND.TAUDEC.EQ.'TAUOLA') THEN
            CALL HWDTAU(1,MHEP,0.0D0)
            GOTO 100
          ENDIF
C Choose decay mode
          ISTHEP(MHEP)=ISTHEP(MHEP)+5
          RN=HWRGEN(2)
          BF=0.
          IM=LSTRT(IDM)
          DO 10 I=1,NMODES(IDM)
          BF=BF+BRFRAC(IM)
          IF (BF.GE.RN) GOTO 20
  10      IM=LNEXT(IM)
          CALL HWWARN('HWDHAD',50,*20)
  20      IF ((IDKPRD(1,IM).GE.1.AND.IDKPRD(1,IM).LE.13).OR.
     &        (IDKPRD(3,IM).GE.1.AND.IDKPRD(3,IM).LE.13)) THEN
C Partonic decay of a heavy-(b,c)-hadron, store details
            NQDK=NQDK+1
            IF (NQDK.GT.NMXQDK) CALL HWWARN('HWDHAD',102,*999)
            LOCQ(NQDK)=MHEP
            IMQDK(NQDK)=IM
            CALL HWVEQU(4,VERTX,VTXQDK(1,NQDK))
            GOTO 100
          ELSE
C Exclusive decay, add decay products to event record
            IF (NHEP+NPRODS(IM).GT.NMXHEP)
     &        CALL HWWARN('HWDHAD',103,*999)
            JDAHEP(1,MHEP)=NHEP+1
            DO 30 I=1,NPRODS(IM)
            NHEP=NHEP+1
            IDHW(NHEP)=IDKPRD(I,IM)
            IDHEP(NHEP)=IDPDG(IDKPRD(I,IM))
            ISTHEP(NHEP)=193
            JMOHEP(1,NHEP)=MHEP
            JMOHEP(2,NHEP)=JMOHEP(2,MHEP)
            PHEP(5,NHEP)=RMASS(IDKPRD(I,IM))
  30        CALL HWVEQU(4,VERTX,VHEP(1,NHEP))
            JDAHEP(2,MHEP)=NHEP
          ENDIF
C Next choose momenta:
          IF (NPRODS(IM).EQ.1) THEN
C 1-body decay: K0(BR) --> K0S,K0L
            CALL HWVEQU(4,PHEP(1,MHEP),PHEP(1,NHEP))
          ELSEIF (NPRODS(IM).EQ.2) THEN
C 2-body decay
C---SPECIAL TREATMENT OF POLARIZED MESONS
            COSANG=TWO
            IF (ID.EQ.IDHW(JMOHEP(1,MHEP))) THEN
              MO=JMOHEP(1,MHEP)
              RSUM=0
              DO 40 I=1,3
  40          RSUM=RSUM+RHOHEP(I,MO)
              IF (RSUM.GT.ZERO) THEN
                RSUM=RSUM*HWRGEN(3)
                IF (RSUM.LT.RHOHEP(1,MO)) THEN
C---(1+COSANG)**2
                  COSANG=MAX(HWRGEN(4),HWRGEN(5),HWRGEN(6))*TWO-ONE
                ELSEIF (RSUM.LT.RHOHEP(1,MO)+RHOHEP(2,MO)) THEN
C---1-COSANG**2
                  COSANG=2*COS((ACOS(HWRGEN(7)*TWO-ONE)+PIFAC)/THREE)
                ELSE
C---(1-COSANG)**2
                  COSANG=MIN(HWRGEN(8),HWRGEN(9),HWRGEN(10))*TWO-ONE
                ENDIF
              ENDIF
            ENDIF
            CALL HWDTWO(PHEP(1,MHEP),PHEP(1,NHEP-1),
     &                  PHEP(1,NHEP),CMMOM(IM),COSANG,.FALSE.)
          ELSEIF (NPRODS(IM).EQ.3) THEN
C 3-body decay
            IF (NME(IM).EQ.100) THEN
C  Use free massless (V-A)*(V-A) Matrix Element
              CALL HWDTHR(PHEP(1,MHEP),PHEP(1,NHEP-1),PHEP(1,NHEP-2),
     &                    PHEP(1,NHEP),HWDWWT)
            ELSEIF (NME(IM).EQ.101) THEN
C  Use bound massless (V-A)*(V-A) Matrix Element
              WTMX=((PHEP(5,MHEP)-PHEP(5,NHEP))
     &             *(PHEP(5,MHEP)+PHEP(5,NHEP))
     &             +(PHEP(5,NHEP-1)-PHEP(5,NHEP-2))
     &             *(PHEP(5,NHEP-1)+PHEP(5,NHEP-2)))/TWO
              WTMX2=WTMX**2
              IPDG=ABS(IDHEP(MHEP))
              XS=ONE-MAX(RMASS(MOD(IPDG/1000,10)),
     &                   RMASS(MOD(IPDG/100,10)),RMASS(MOD(IPDG/10,10)))
     &              /(RMASS(MOD(IPDG/1000,10))+RMASS(MOD(IPDG/100,10))
     &               +RMASS(MOD(IPDG/10,10)))
  50          CALL HWDTHR(PHEP(1,MHEP),PHEP(1,NHEP-1),PHEP(1,NHEP-2),
     &                    PHEP(1,NHEP),HWDWWT)
              DOT1=HWULDO(PHEP(1,MHEP),PHEP(1,NHEP-1))
              DOT2=HWULDO(PHEP(1,MHEP),PHEP(1,NHEP-2))

C O.K. trick to prevent infinite loop              
              IF (DOT1*(WTMX-DOT1-XS*DOT2).LT.ZERO) THEN
c               print *, (DOT1*(WTMX-DOT1-XS*DOT2))," ", HWRGEN(11)*WTMX2
                 IF(INEG.LT.50) THEN
                     INEG=INEG+1
                 ELSE    
                     PRINT *, "O.K.: Killing the event, >50 negative"
                     CALL HWWARN('HWBGEN',100,*999)
                 ENDIF
              END IF
              IF (DOT1*(WTMX-DOT1-XS*DOT2).LT.HWRGEN(11)*WTMX2) GOTO 50
            ELSE IF (NME(IM).EQ.200) THEN
C Use free massless ((V-A)*TB1+(V+A)*CT1)*((V-A)*TB2+(V+A)*CT2)) Matrix Element
C sort tan(beta)
              IF((IDK(IM).EQ.  2).OR.(IDK(IM).EQ.  4).OR.
     &           (IDK(IM).EQ.  6).OR.(IDK(IM).EQ.  8).OR.
     &           (IDK(IM).EQ. 10).OR.(IDK(IM).EQ. 12).OR.
     &           (IDK(IM).EQ.122).OR.(IDK(IM).EQ.124).OR.
     &           (IDK(IM).EQ.126).OR.(IDK(IM).EQ.128).OR.
     &           (IDK(IM).EQ.130).OR.(IDK(IM).EQ.132))THEN
                TB=TANB
              ELSE
                TB=1./TANB
              END IF
              IF((IDKPRD(1,IM).EQ.  2).OR.(IDKPRD(1,IM).EQ.  4).OR.
     &           (IDKPRD(1,IM).EQ.  6).OR.(IDKPRD(1,IM).EQ.  8).OR.
     &           (IDKPRD(1,IM).EQ. 10).OR.(IDKPRD(1,IM).EQ. 12).OR.
     &           (IDKPRD(1,IM).EQ.122).OR.(IDKPRD(1,IM).EQ.124).OR.
     &           (IDKPRD(1,IM).EQ.126).OR.(IDKPRD(1,IM).EQ.128).OR.
     &           (IDKPRD(1,IM).EQ.130).OR.(IDKPRD(1,IM).EQ.132))THEN
                BT=TANB
              ELSE
                BT=1./TANB
              END IF
              IT1=IDK(IM)
              IB1=IDKPRD(3,IM)
              IT2=IDKPRD(1,IM)
              IB2=IDKPRD(2,IM)
              CALL HWDTHR(PHEP(1,MHEP),PHEP(1,NHEP),PHEP(1,NHEP-2),
     &                    PHEP(1,NHEP-1),HWDHWT)
            ELSE
              CALL HWDTHR(PHEP(1,MHEP),PHEP(1,NHEP-2),PHEP(1,NHEP-1),
     &                    PHEP(1,NHEP),HWDPWT)
            ENDIF
          ELSEIF (NPRODS(IM).EQ.4) THEN
C 4-body decay
            CALL HWDFOR(PHEP(1,MHEP  ),PHEP(1,NHEP-3),PHEP(1,NHEP-2),
     &                  PHEP(1,NHEP-1),PHEP(1,NHEP))
          ELSEIF (NPRODS(IM).EQ.5) THEN
C 5-body decay
            CALL HWDFIV(PHEP(1,MHEP  ),PHEP(1,NHEP-4),PHEP(1,NHEP-3),
     &                  PHEP(1,NHEP-2),PHEP(1,NHEP-1),PHEP(1,NHEP))
          ELSE
            CALL HWWARN('HWDHAD',104,*999)
          ENDIF
        ENDIF
      ENDIF
  100 CONTINUE
C---MAY HAVE OVERFLOWED /HEPEVT/
      CALL HWWARN('HWDHAD',105,*999)
  999 END


C-----------------------------------------------------------------------
      SUBROUTINE HWCHEK
C
C     M.Boonekamp, Nov 2003 : 
C
C     Performs some checks on the user's FPMC settings.
C     Prints out what it finds
C
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INTEGER IPR
      INTEGER PDFID, ACTID

c------------------------------------------------------------------- 

C ... Check NFLUX / TYPINT compatibility
      IPR=MOD(IPROC/100,100)
      PRINT*, ' '
      PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
      PRINT*, ' '

      IF(TYPINT.NE.'QED') THEN
        IF(TYPINT.NE.'QCD') THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested ilegal type of interaction '
        PRINT*, '         TYPINT = ', TYPINT
        PRINT*, ' TYPINT should be set to QCD/QED '
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        STOP
        ENDIF   
      ENDIF

      IF(TYPEPR.NE.'INC') THEN
         IF(TYPEPR.NE.'EXC') THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested ilegal type of production '
        PRINT*, '         TYPEPR = ', TYPEPR
        PRINT*, ' TYPEPR should be set to INC/EXC '
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        STOP
        ENDIF   
      ENDIF

      IF(NFLUX.GE.9.AND.NFLUX.LE.11.AND.TYPINT.EQ.'QED') THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested a Pomeron/Reggeon flux'
        PRINT*, '          with a photon initiated hard process:'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          TYPINT = ', TYPINT 
        PRINT*, ' NFLUX=9,10,11 should be used with TYPINT=QCD'
        PRINT*, ' NFLUX=12,13,14 should be used with TYPINT=QED'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
      ELSEIF(NFLUX.EQ.11.AND.TYPEPR.EQ.'INC'.AND.CDFFAC.EQ.0d0) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested an Inclusive Pomeron flux'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          TYPEPR = ', TYPEPR 
        PRINT*, ' CDFFAC not set ; now set to 1d0'
        PRINT*, ' (suggested value : 3.8d0)'
        CDFFAC = 1d0
        PRINT*, ' '
      ELSEIF(NFLUX.EQ.11.AND.TYPEPR.EQ.'EXC'.AND.GAPSPR.EQ.0d0) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested an Exclusive Pomeron flux'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          TYPEPR = ', TYPEPR 
        PRINT*, ' GAPSPR not set ; now set to 1d0'
        PRINT*, ' (suggested values : 0.03 at LHC, 0.1 at Tevatron)'
        GAPSPR = 1d0
        PRINT*, ' '
      ELSEIF(NFLUX.GE.12.AND.NFLUX.LE.15.AND.TYPINT.EQ.'QCD') THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested a Photon flux'
        PRINT*, '          with a QCD hard process:'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          TYPINT = ', TYPINT 
        PRINT*, ' NFLUX=12,13,14,15 should be used with TYPINT=QED,',
     & ' TYPEPR=EXC'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
C ... begin R.S.
      ELSEIF(NFLUX.EQ.18 .AND. TYPINT.EQ.'QED') THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested a CHIDe model'
        PRINT*, '          with a photon initiated hard process:'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          TYPINT = ', TYPINT 
        PRINT*, ' NFLUX=18 should be used with TYPINT=QCD'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
        ELSEIF(NFLUX.EQ.18 .AND. .NOT.(IPROC.EQ.19999 
     &  .OR. IPROC.EQ.16012.OR.IPROC.EQ.16059)) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested a CHIDe model'
        PRINT*, '          with wrong IPROC:'
        PRINT*, '          NFLUX  = ', NFLUX
        PRINT*, '          IPROC = ', IPROC 
        PRINT*, ' Availiable IPROC: 19999, 16012, 16059'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP

C ... end R.S.
      ELSEIF(TYPEPR.EQ.'EXC'.AND.IPR.NE.60.AND.IPR.NE.99.AND.IPR.NE.98
     &       .AND.IPR.NE.97.AND.IPR.NE.96) THEN
C CHR to do t tbar higgs, does not exist in gamma gamma herwig
c     &       .AND.IPR.NE.97.AND.IPR.NE.96.and.IPR.NE.25) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested an exclusive cross-section'
        PRINT*, '          with a bad IPROC:'
        PRINT*, '          IPROC  = ', IPROC
        PRINT*, '          TYPEPR = ', TYPEPR
        PRINT*, ' TYPEPR=EXC needs IPROC=(1)60**, (1)96**, (1)97**, 
     & (1)98** or (1)99**'
        PRINT*,'IPR :', IPR
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
      ELSEIF(TYPEPR.EQ.'INC'.AND.(IPR.EQ.99.OR.IPR.EQ.98)) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested an inclusive cross-section'
        PRINT*, '         with IPROC set for exclusive production :'
        PRINT*, '          IPROC  = ', IPROC
        PRINT*, '          TYPEPR = ', TYPEPR
        PRINT*, ' IPROC=(1)98**, (1)99** needs  TYPEPR=EXC'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
      ELSEIF(TYPEPR.EQ.'EXC'.AND.NFLUX.EQ.9) THEN
        PRINT*, ' '
        PRINT*, ' FPMC - You requested an exclusive production'
        PRINT*, '        with NFLUX=9 which corresponds to inclusive'
        PRINT*, '        factorized model only'
        PRINT*, ' '
        PRINT*, ' - STOP'
        PRINT*, ' '
        PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
        PRINT*, ' '
        STOP
      ENDIF
C ... Compute Higgs BR's : not done in HWUINC when IPRO=99...
      IF(IPR.EQ.99) THEN
        GAMH=RMASS(201+IHIGGS)
        CALL HWDHIG(GAMH)
      ENDIF
C ... Check the pdf choice
c      PDFID=MOD(IFITPDF, 100)
c      ACTID=IFITPDF-PDFID
       PDFID=IFITPDF
      IF ((PDFID.EQ.2).OR.(PDFID.EQ.5).OR.(PDFID.EQ.8).OR.
     .   (PDFID.EQ.10).OR.(PDFID.EQ.20).OR.(PDFID.EQ.30).OR.
     .   (PDFID.EQ.100).OR.(PDFID.EQ.101)) THEN
         !nothing
      ELSE
         PRINT*, ' '
         print*, ' FPMC - You requested ilegal PDF parameter'
         print*, '         (see h1qcd.f for details)'
         PRINT*, ' '
         print*, ' - STOP'
         PRINT*, ' '
         PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
         PRINT*, ' '
         STOP
      ENDIF


      IF((NFLUX.GE.12).AND.(NFLUX.LE.15)) THEN

         IF( (IPROC.EQ.16015).AND.AAANOM.NE.2) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Anomalous AAZZ coupling available '
            PRINT*, '        only with AAANOM = 2'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
         ENDIF


         IF( (IPROC.EQ.16016).AND.AAANOM.NE.2) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Anomalous AAAA coupling available '
            PRINT*, '        only with AAANOM = 2'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
         ENDIF


         IF(AAANOM.EQ.0)THEN
         ELSEIF(AAANOM.EQ.1)THEN
            PRINT *, 'Full SM formula for AAWW used'
         ELSEIF(AAANOM.EQ.2)THEN
           IF( IPROC.NE.16010.AND.IPROC.NE.16015
     & .AND.IPROC.NE.16016) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Anomalous AAWW coupling available '
            PRINT*, '        only for IPROC = 16010 -- WW'
            PRINT*, '                 IPROC = 16015 -- ZZ'
            PRINT*, '                 IPROC = 16016 -- AA'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
           ELSEIF(IPROC.EQ.16010) THEN
            PRINT *, 'Comphep ME used for TGC+QGC AA->WW '
            PRINT *, 'Vertex couplings set to:'
            PRINT *, '   D_KAPPA  = ', D_KAPPA 
            PRINT *, '   LAMBDA = ', LAMBDA
            PRINT *, '   A0W = ', A0W
            PRINT *, '   ACW = ', ACW
            PRINT *, '   ANOMCUTOFF = ', ANOMCUTOFF
            PRINT *, ' For ANOMCUTOFF < 0 GeV, dipole form of the
     & coupling formfactor turned off (see comphep_wrapper.cpp)'
           ELSEIF(IPROC.EQ.16015) THEN
            PRINT *, 'Comphep ME used for QGC AA->ZZ '
            PRINT *, 'Vertex couplings set to:'
            PRINT *, '   A0Z = ', A0Z
            PRINT *, '   ACZ = ', ACZ
            PRINT *, '   ANOMCUTOFF = ', ANOMCUTOFF
            PRINT *, ' For ANOMCUTOFF < 0 GeV, dipole form of the
     & coupling formfactor turned off (see comphep_wrapper.cpp)'
           ELSEIF(IPROC.EQ.16016) THEN
            PRINT *, 'Comphep ME used for QGC AA->AA '
            PRINT *, 'Vertex couplings set to:'
            PRINT *, '   A1A = ', A1A
            PRINT *, '   A2A = ', A2A
            PRINT *, '   ANOMCUTOFF = ', ANOMCUTOFF
            PRINT *, ' For ANOMCUTOFF < 0 GeV, dipole form of the
     & coupling formfactor turned off (see comphep_wrapper.cpp)'
           ENDIF
         ELSE
               PRINT *, 'Unknown AAANOM = ', AAANOM
               PRINT *, 'Should be =0, 1, 2'
               STOP
         ENDIF

C--M.S exclusive diphoton production
         IF( (IPROC.EQ.16063).AND.AAEXOTIC.NE.1) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Exotic AAAA coupling available '
            PRINT*, '        only with AAEXOTIC = 1'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
         ENDIF

         IF( (IPROC.EQ.16064).AND.AAEXOTIC.NE.1) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Exotic AAAA coupling available '
            PRINT*, '        only with AAEXOTIC = 1'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
         ENDIF

         IF(AAEXOTIC.EQ.0)THEN
         ELSEIF(AAEXOTIC.EQ.1)THEN
            PRINT *, 'EXOTICS FOR EXCL AAAA'
           IF( IPROC.NE.16063.AND.IPROC.NE.16064) THEN
            PRINT*, ' '
            PRINT*, ' FPMC - Exotic AAAA coupling available '
            PRINT*, '        only for IPROC = 16063 -- Bosons'
            PRINT*, '                 IPROC = 16064 -- Fermions'
            PRINT*, ' - STOP'
            PRINT*, ' '
            PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
            PRINT*, ' '
            STOP
           ELSEIF(IPROC.EQ.16063) THEN
            PRINT *, 'G.von Gersdorff ME used for AA->AA '
            PRINT *, 'With extra bosons         '
            PRINT *, 'Exotic bosons parameters set to:'
            PRINT *, '   AAM  = ', AAM 
            PRINT *, '   AAQ = ', AAQ
            PRINT *, '   AAN = ', AAN
           ELSEIF(IPROC.EQ.16064) THEN
            PRINT *, 'G.von Gersdorff ME used for AA->AA '
            PRINT *, 'With extra fermions         '
            PRINT *, 'Exotic fermions parameters set to:'
            PRINT *, '   AAM  = ', AAM 
            PRINT *, '   AAQ = ', AAQ
            PRINT *, '   AAN = ', AAN
           ENDIF
         ELSE
               PRINT *, 'Unknown AAEXOTIC = ', AAEXOTIC
               PRINT *, 'Should be =0, 1'
               STOP
         ENDIF
      ENDIF   
               
C ... Print out settings
      CALL PRINTSETTING()
      PRINT*, ' '
      PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
      PRINT*, ' '
      END
C-----------------------------------------------------------------------


c------------------------------------------------------------------- 
c * 15/11/2006 OK
c * Routin for setting model dependant parameters...
c * ... CDF factorization according to IFITPDF
c * Previously located in examples     
c------------------------------------------------------------------- 
      SUBROUTINE HWMODINI

      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'CHIDe.inc'
c---To pass the relative GAP Srv. prob factor for BL inclusive modified dist.
      DOUBLE PRECISION GAPSPRREL
      COMMON /BLMODIF/  GAPSPRREL

      INTEGER IFIT,NFLAVR
      DOUBLE PRECISION XPQ(-6:6),XP,Q2,SCORR,X,Y,Z
      COMMON /NORM/ Q2,NFLAVR

      CHARACTER(80) STABLE


c ... begin R.S.
      character*500  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
c ... end R.S.
       

c--- Initialization of the QCD subroutine
      IF(TYPEPR.NE.'EXC') THEN
         Q2=75D0        !dummy - just initialization
          IFIT = IFITPDF

         XP=1D-1
         CALL QCDFIT(XP,Q2,XPQ,IFIT)
      ENDIF
      
c... Inc.DPE : cross-section normalization to the CDF measurement (NFLUX=11, 'INC')
c      old factor CDFFAC = 3.8d0
c OK 15/11/06  tuning to CDF Phys.Rev.Lett 85, 4215 (ptmin = 7GeV)
c CHR/YURA
      IF(NFLUX.LE.11.OR.NFLUX.EQ.19.OR.NFLUX.EQ.21) THEN
        IF(IFITPDF.EQ.10) THEN
            CDFFAC = 0.074d0/8d0
        ELSEIF(IFITPDF.EQ.20) THEN
            CDFFAC = 0.085d0/8d0
        ELSEIF(IFITPDF.EQ.30) THEN
            CDFFAC = 0.072d0/8d0
        ELSEIF(IFITPDF.EQ.100) THEN
            CDFFAC = 0.072d0/8d0
        ELSEIF(IFITPDF.EQ.101) THEN
            CDFFAC = 0.072d0/8d0
        ELSEIF(IFITPDF.EQ.2) THEN
            CDFFAC = 0.111d0/8d0
        ELSEIF(IFITPDF.EQ.5) THEN
            CDFFAC = 0.111d0/8d0
        ELSEIF(IFITPDF.EQ.8) THEN  
            CDFFAC = 0.111d0/8d0
        ELSE 
            print *, 'Non-standard choice of IFITPDF(2/5/8/10/20/30)' 
            STOP
        ENDIF  
      ENDIF

c... option for soft corrections : ISOFTM = 0 : no correction
c                                           1 : simple factor, see below
c                                           2 : KMR low mass diffractive
c                                           3 : effective opacity model
      IF(ISOFTM.EQ.1) THEN
c        Proton rescattering factor : at the LHC, put 0.03 for NFLUX=11, QCD/EXC
c                                                 0.90 for NFLUX=14,15, QED
         PROSPR = 0.90d0 ! used when NFLUX=14,15, with QED and ISOFTM=1

         IF(PBEAM1.LE.4.5d3) THEN
                        ! used when NFLUX=11, with QCD/EXC and ISOFTM=1
            GAPSPR = 0.10d0 
            GAPSPRREL = 1
         ELSE
            GAPSPRREL = 0.03/0.1
            GAPSPR = 0.03d0   
         ENDIF   

         
c--- Initialize soft correction code (from A.Kupco)
      ELSEIF(ISOFTM.EQ.2) THEN
         IF(PBEAM1.LE.1d3) THEN
            != './Tables/soft.2TeV.tab'
            CALL GETSOFTC(SCORR,X,Y,Z,1,1)
         ELSEIF(PBEAM1.GE.4.5d3) THEN
            != './Tables/soft.14TeV.tab'
            CALL GETSOFTC(SCORR,X,Y,Z,1,2)
         ELSE
            PRINT*, 'Soft corr : wrong ECMS - change ISOFTM'
            STOP
         ENDIF
      ELSEIF(ISOFTM.EQ.3) THEN
         IF(PBEAM1.LE.1d3) THEN
            !='./Tables/soft.2TeV.effopa.tab'
            CALL GETSOFTC(SCORR,X,Y,Z,1,3)
         ELSEIF(PBEAM1.GE.4.5d3) THEN
            !='./Tables/soft.14TeV.effopa.tab'
            CALL GETSOFTC(SCORR,X,Y,Z,1,4)
         ELSE
            PRINT*, 'Soft corr : wrong ECMS - change ISOFTM'
            STOP
         ENDIF
      ENDIF

c--- Initialize KMR calculation

      IF(NFLUX.EQ.16) THEN
         CALL KMRINI(99,CHIDeS,KMR2SURV,KMR2SCALE,KMR2DELTA,KMR2Q2CUT)
c      IF(PBEAM1.LE.1d3) THEN
c         CALL KMRINI(1)
c      ELSEIF(PBEAM1.GE.6.5d3) THEN
c         CALL KMRINI(2)
c      ENDIF
      ENDIF

C ... begin R.S.
C ... Initialize CHIDe model
      IF(NFLUX.EQ.18 .AND. IPROC.EQ.19999) THEN
        dgdtab1=TRIM(CHIDePATH) // "dgdtab1.d"
        dgdtab2=TRIM(CHIDePATH) // "dgdtab2.d"
        dgdtab3=TRIM(CHIDePATH) // "dgdtab3.d"
        dgdtab4=TRIM(CHIDePATH) // "dgdtab4.d"
        sudatab=TRIM(CHIDePATH) // "Higgssudatab.d"
C     Parameters: Higgs mass, top mass, cms energy^2, impact factor
C     parameterisation, upper and lower limit scales of sudakov factor,
C     Rapidity Gap Survival
        CALL CHIDeHiggsInit(RMASS(201),RMASS(6),CHIDeS,CHIDeIGLU,
     &  CHIDeX, CHIDeXP, CHIDeS2)
      ELSEIF(NFLUX.EQ.18 .AND. IPROC.EQ.16012) THEN
        dgdtab1=TRIM(CHIDePATH) // "dgdtab1.d"
        dgdtab2=TRIM(CHIDePATH) // "dgdtab2.d"
        dgdtab3=TRIM(CHIDePATH) // "dgdtab3.d"
        dgdtab4=TRIM(CHIDePATH) // "dgdtab4.d"
        sudatab=TRIM(CHIDePATH) // "ggsudatab.d"
        CALL CHIDeGGInit(CHIDeS,DBLE(PTMIN),CHIDeS2,
     &       CHIDeX,CHIDeXP,CHIDeIGLU,
     &       XI1Min, XI1Max, XI2Min, XI2Max,
     &       DBLE(YJMIN), DBLE(YJMAX))
        XI1MIN=XI1MIN/1d3! CHIDeZ1 = CHIDeB1 + CHIDeB2 
c        XI1MIN=1d-10 ! CHIDeZ1 = CHIDeB1 + CHIDeB2 
c       XI1MAX=1d0 ! CHIDeZ1 = CHIDeB1 + CHIDeB2 
      ELSEIF(NFLUX.EQ.18 .AND. IPROC.EQ.16059) THEN
        dgdtab1=TRIM(CHIDePATH) // "dgdtab1.d"
        dgdtab2=TRIM(CHIDePATH) // "dgdtab2.d"
        dgdtab3=TRIM(CHIDePATH) // "dgdtab3.d"
        dgdtab4=TRIM(CHIDePATH) // "dgdtab4.d"
        sudatab=TRIM(CHIDePATH) // "ggsudatab.d"
        CALL CHIDeDiphotonInit(CHIDeS,DBLE(PTMIN),CHIDeS2,
     &       CHIDeX,CHIDeXP,CHIDeIGLU,
     &       XI1Min, XI1Max, XI2Min, XI2Max,
     &       DBLE(YJMIN), DBLE(YJMAX))
        XI1MIN=XI1MIN/1d3! CHIDeZ1 = CHIDeB1 + CHIDeB2 
c        XI1MIN=1d-10 ! CHIDeZ1 = CHIDeB1 + CHIDeB2 
c       XI1MAX=1d0 ! CHIDeZ1 = CHIDeB1 + CHIDeB2 
      ENDIF
C ... end R.S.


c---Initialize the pomeron/regeon parameters
c   OK 29/11/06
c CHR/Yura add Reggeon...
C nflux 19: Pomeron Reggeon
c nflux 20: Photon Pomeron
c nflux 21: Reggeon Pomeron
c nflux 22: Pomeron Photon

      IF ((NFLUX.EQ.9).OR.(NFLUX.EQ.10).OR.(NFLUX.EQ.19).OR.
     +     (NFLUX.EQ.20).OR.(NFLUX.EQ.21).OR.(NFLUX.EQ.22).OR.
     +     (NFLUX.EQ.25).OR.(NFLUX.EQ.26)) THEN

               !1994 tables FitB: hep/ph 9708016
         IF((IFITPDF.EQ.2).OR.(IFITPDF.EQ.5)) THEN
            alphaP=1.203
            alphaPp=0.26
            Bpom=4.6                  ! pomeron
            alphaR=0.50
            alphaRp=0.90
            Breg=2.0
            Cr=16.0             ! reggeon only normalization from POMWIG
            zh1=0.003
c           FITB            
c            alphaP=1.200
c            alphaPp=0.26
c            Bpom=4.6                  ! pomeron
c            alphaR=0.57
c            alphaRp=0.90
c            Breg=2.0
c            Cr=48.0             ! reggeon only normalization from POMWIG
c            zh1=0.003
         ENDIF

               !1997 tables hep/ph 0602228
         IF(IFITPDF.EQ.8) THEN
            alphaP=1.19
            alphaPp=0.26
            Bpom=4.6                  ! pomeron
            alphaR=0.62
            alphaRp=0.90
            Breg=2.6
            Cr=48.0             ! reggeon only normalization from POMWIG
            zh1=0.003
          !H1 best fits: hep/ph 0609291
          ELSEIF((IFITPDF.EQ.10).OR.(IFITPDF.EQ.20).OR.
     *           (IFITPDF.EQ.30)) THEN
            IF(IFITPDF.EQ.10) THEN              ! ... H1 fits
                  alphaP=1.12
               ELSEIF(IFITPDF.EQ.20) THEN       ! ... Zeus fits         
                  alphaP=1.104
               ELSEIF(IFITPDF.EQ.30) THEN       ! ... H1Zeus
                  alphaP=1.118
            ENDIF
            alphaPp=0.06
            Bpom=5.5                  ! pomeron
            alphaR=0.50
            alphaRp=0.30
            Breg=1.6
            Cr=48.0             ! reggeon only normalization from POMWIG
            zh1=0.003
         ELSEIF((IFITPDF.EQ.100).OR.(IFITPDF.EQ.101)) THEN
           IF(IFITPDF.EQ.100) THEN              ! ... H1 fit A
                  alphaP=1.118
C CHR YURA
                  Cr=0.0017
               ELSEIF(IFITPDF.EQ.101) THEN       ! ... H1 fit B
                  alphaP=1.111
C CHR YURA
		  Cr=0.0014
            ENDIF
            
            alphaPp=0.06
            Bpom=5.5                  ! pomeron
c CHR/YURA modify values for reggeon	    
            alphaR=0.5
            alphaRp=0.3
            Breg=1.6
            Cr=0.0014             ! reggeon only normalization from POMWIG
            zh1=0.003
         ENDIF   
c___________________________________________________________________          
c         old version
c...H1 best fits: 
c         alpha=1.203
c         alphaPp=0.26
c         B=4.6                  ! pomeron
c         alphar=0.50
c         alpharp=0.90
c         Br=2.0
c         Cr=16.0                ! reggeon
c...  H1 parameters with no interference (best fit to H1 F2D3 using POMWIG):
*     PARAMETER (alpha=1.200,alphap=0.26,B=4.6)
*     PARAMETER (alphar=0.57,alpharp=0.90,Br=2.0,Cr=48.0)
c------------------------------------------------------------------- 
      ELSEIF(NFLUX.EQ.11)THEN
c---Determination of the normalizations for PDFs inside pomeron:
      NFLAVR=4 ! number of flavors considered + gluon
      CALL NRMPDF(XPQNRM)
      CALL NRMPRT(XPQNRM)

c...  Bialas-Landshoff model:
         alphaP=1.08
         alphaPp=0.25
         Bpom=6.0                  !pomeron
         alphaR=0.57
         alphaRp=0.90
         Breg=2.0
         Cr=48.0                !reggeon here dummy
c         zh1=0                 ! dummy for B.L.
      END IF

      RETURN
      END


      
C----------------------------------------------------------------------
C * 08/19/2003 Maarten Boonekamp and Tibor Kucs:
C * Modified to account for different normalizations 
C   in the cross section and gluon PDF inside the pomeron in 
C   Cox-Forshaw and Bialas-Landshoff models. We also incorporated the 
C   fact that the functional dependence of the flux factors
C   is different in the above mentioned models.
C * NFLUX = 11   Bialas model
C * NFLUX = 12   heavy ion-gamgam (valid for all Z,A)
C * NFLUX = 13   heavy ion-gamgam (heavy ions only)
C * NFLUX = 14   Photon flux (E.Papageorgiu)
C * NFLUX = 15   Correct Budnev flux (from K.Piotrzkowski)
C * NFLUX = 16   KMR-like flux (from L.Lonnblad)
C nflux 19: Pomeron Reggeon
c nflux 20: Photon Pomeron
c nflux 21: Reggeon Pomeron
c nflux 22: Pomeron Photon
c nflux 23: proton ion - photon photon
c nflux 24: ion proton - photon photon
c nflux 25: proton ion - pomeron photon
c nflux 26: ion proton - photon pomeron
C-----------------------------------------------------------------------
      SUBROUTINE FLUX(F,Z,TMIN,TMAX,IPRO,IND)
      IMPLICIT NONE
c---Returns H1 pomeron flux if NFLUX=9
c...Returns H1 reggeon flux if NFLUX=10
c...Returns flux for user defined structure function if NFLUX>10
c...B.Cox and J. Forshaw 11/05/00
      DOUBLE PRECISION F,Z,TMIN,TMAX,TTMIN,TTMAX
      INTEGER IPRO,IND
c      DOUBLE PRECISION alpha,B,alphap
c      DOUBLE PRECISION alphar,alpharp,Br,Cr
      DOUBLE PRECISION V,W,X,FAC
      DOUBLE PRECISION QZERO,EXPARG,ALPHAE,EI,QMIN2
      DOUBLE PRECISION RZERO,R,BMIN,XM,ZZERO,ARG,F1,F2
      DOUBLE PRECISION HWUALF,DGAGNC,DBESK0,DBESK1
      INCLUDE 'fpmc.inc'
      EXTERNAL HWUALF,DGAGNC,DBESK0,DBESK1
      PARAMETER(ALPHAE=1.D0/137.D0)
c... parameters for the Budnev flux (thanks K.Piotrzkowski)
      DOUBLE PRECISION QSCALE, BUD, Q2MIN, Q2MAX, PHI
      PARAMETER (QSCALE=0.71d0,BUD=1d0)!5.4d-6)
      EXTERNAL PHI
c... parameters for the KMR model
      DOUBLE PRECISION Z1KMR,Z2KMR,FLKMR
      SAVE Z1KMR,Z2KMR

      integer IIND
      
      COMMON /YURA/ IIND
c index of protons Yura/CHR      
      IIND = IND


c---Choice of flux parameters:
c   OK 29/11/06 moved to MODINI  

c        print *,'**** NFLUX:',nflux
c---Fluxes for different models:
c...Cox-Forshaw pomeron flux:
      IF(NFLUX.EQ.9) THEN
c                  print *, 'alpha:', alphaP,  ' alphap:', alphaPp, 
c     .             'TMIN:', TMIN, ' TMAX:', TMAX, ' B:', Bpom, ' Z:', Z

         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
         F = X*W*V
c...Cox-Forshaw reggeon flux:
      ELSEIF (NFLUX.EQ.10) THEN  
         V = DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Breg+2.D0*alphaRp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaR-1.D0))
         F = Cr*X*W*V         
      ELSEIF (NFLUX.EQ.11) THEN
C---Bialas-Landshoff pomeron flux:
c ... T.K. : In order to match the functional dependence of Bialas-Landshoff 
c ... and Cox-Forshaw models we need to rescale the integration boundaries:
         TTMIN = TMIN*(1-Z)
         TTMAX = TMAX*(1-Z)
         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TTMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TTMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
c ... Factors accounting for:
c ...... 1) different normalization of gluonic PDF
c ...... 2) jacobian of going from t to v^2
c ...... 3) numerical factor due to the nonperturbative couplings in the 
c ......    Bialas-Landshoff model
c ...... 4) depending on whether we produce a color singlet or not, we must
c ......    divide a factor root(2*[Nc^2-1])
         FAC = 1.D0/(1.D0-Z)*DSQRT(2304.D0/(PI)**5)
         IF(TYPEPR.EQ.'INC'.AND.TYPINT.EQ.'QCD') THEN
            FAC = FAC / XPQNRM(0)
         ELSEIF(TYPEPR.EQ.'EXC'.AND.TYPINT.EQ.'QCD') THEN
            FAC = FAC / DSQRT(16.D0)
         ENDIF
         F = FAC*X*W*V
      ELSEIF (NFLUX.EQ.12) THEN
C---Photon flux from ions : should be valid for all Z
c ... T.K. : Implemented factorized flux (11) in Cahn, Jackson; PR D42 (1990) 3690
c ... M.B. : coherency conditions are adapted for ions
c corrections heavy ions CHR 10/2014
c         goto 1516
c         RZERO=1.2/FMCONV
c         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R ! f0actorized flux is a good approximation for larger radius
c         XM=AION*AMASS
c         ZZERO=1d0/XM/BMIN
c         ARG=Z/ZZERO
c         FAC=2.*(ZION**2)*ALPHAE/PI/ZZERO
c         F1=DBESK0(ARG)*DBESK1(ARG)
c         F2=0.5*ARG*(DBESK1(ARG)-DBESK0(ARG))*(DBESK1(ARG)+DBESK0(ARG))
c         F=FAC*(F1-F2)
c1516     continue
c         print *,'con',FMCONV,AMASS
         RZERO=1.2/FMCONV
         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R
CHR 10/2014 The cutoff should be transmitted in the cards (depends on ion)
         BMIN=RBMIN*R
c         print *,'BMIN = ',RBMIN
         ARG=Z*AMASS*BMIN
         F=2.*ARG*DBESK0(ARG)*DBESK1(ARG)
         F=F-ARG**2.D0*(DBESK1(ARG)**2.D0-DBESK0(ARG)**2.D0)
         F=F*ZION**2.D0*ALPHAE/PI/Z


      ELSEIF (NFLUX.EQ.23) THEN
C---Photon flux from ions : should be valid for all Z
c ... T.K. : Implemented factorized flux (11) in Cahn, Jackson; PR D42 (1990) 3690
c ... M.B. : coherency conditions are adapted for ions
c CHR ions 10/2014
c         print *,'con',FMCONV,AMASS
         IF(IND.EQ.1) THEN
c         Q2MIN = Z*Z*0.88d0/(1d0-Z)/QSCALE
c         Q2MAX = TMAX/QSCALE
c         F = ALPHAE/PI*(1d0-Z)/Z*(PHI(Q2MAX,Z)-PHI(Q2MIN,Z))
c ... Rangel and Goncalves: PRD 39 (1989) 2536
         QMIN2 = (Z*AMASS)**2.D0
         ARG   = 1.D0+0.71D0/QMIN2
         F = DLOG(ARG)-11.D0/6.D0+3.D0/ARG
         F = F-3.D0/(2.D0*ARG*ARG)+1.D0/(3.D0*ARG*ARG*ARG)
         F = F*ALPHAE/PI
         F = F*(1.D0-Z+(1.D0/2.D0)*Z*Z)
         F = F/Z
         ELSEIF(IND.EQ.2) THEN
         RZERO=1.2/FMCONV
         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R
c         BMIN=1.1*R
         BMIN=RBMIN*R
         ARG=Z*AMASS*BMIN
         F=2.*ARG*DBESK0(ARG)*DBESK1(ARG)
         F=F-ARG**2.D0*(DBESK1(ARG)**2.D0-DBESK0(ARG)**2.D0)
         F=F*ZION**2.D0*ALPHAE/PI/Z
         ENDIF
      ELSEIF (NFLUX.EQ.24) THEN
C---Photon flux from ions : should be valid for all Z
c ... T.K. : Implemented factorized flux (11) in Cahn, Jackson; PR D42 (1990) 3690
c ... M.B. : coherency conditions are adapted for ions
c CHR ions 10/2014
c         print *,'con',FMCONV,AMASS
         IF(IND.EQ.2) THEN
         Q2MIN = Z*Z*0.88d0/(1d0-Z)/QSCALE
         Q2MAX = TMAX/QSCALE
         F = ALPHAE/PI*(1d0-Z)/Z*(PHI(Q2MAX,Z)-PHI(Q2MIN,Z))
         ELSEIF(IND.EQ.1) THEN
         RZERO=1.2/FMCONV
         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R
c         BMIN=1.1*R
         BMIN=RBMIN*R
         ARG=Z*AMASS*BMIN
         F=2.*ARG*DBESK0(ARG)*DBESK1(ARG)
         F=F-ARG**2.D0*(DBESK1(ARG)**2.D0-DBESK0(ARG)**2.D0)
         F=F*ZION**2.D0*ALPHAE/PI/Z
         ENDIF
      ELSEIF (NFLUX.EQ.25) THEN
C---Photon flux from ions : should be valid for all Z
c ... T.K. : Implemented factorized flux (11) in Cahn, Jackson; PR D42 (1990) 3690
c ... M.B. : coherency conditions are adapted for ions
c CHR ions 10/2014
c         print *,'con',FMCONV,AMASS
         IF(IND.EQ.1) THEN
         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
         F = X*W*V
         ELSEIF(IND.EQ.2) THEN
         RZERO=1.2/FMCONV
         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R
c         BMIN=1.1*R
         BMIN=RBMIN*R
         ARG=Z*AMASS*BMIN
         F=2.*ARG*DBESK0(ARG)*DBESK1(ARG)
         F=F-ARG**2.D0*(DBESK1(ARG)**2.D0-DBESK0(ARG)**2.D0)
         F=F*ZION**2.D0*ALPHAE/PI/Z
         ENDIF
      ELSEIF (NFLUX.EQ.26) THEN
C---Photon flux from ions : should be valid for all Z
c ... T.K. : Implemented factorized flux (11) in Cahn, Jackson; PR D42 (1990) 3690
c ... M.B. : coherency conditions are adapted for ions
c CHR ions 10/2014
c         print *,'con',FMCONV,AMASS
         IF(IND.EQ.2) THEN
         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
         F = X*W*V
         ELSEIF(IND.EQ.1) THEN
         RZERO=1.2/FMCONV
         R=RZERO*(AION**(1./3.))
c         BMIN=1.2*R
c         BMIN=1.1*R
         BMIN=RBMIN*R
         ARG=Z*AMASS*BMIN
         F=2.*ARG*DBESK0(ARG)*DBESK1(ARG)
         F=F-ARG**2.D0*(DBESK1(ARG)**2.D0-DBESK0(ARG)**2.D0)
         F=F*ZION**2.D0*ALPHAE/PI/Z
         ENDIF
      ELSEIF (NFLUX.EQ.13) THEN
C---Photon flux for protons
c ... Rangel and Goncalves: PRD 39 (1989) 2536
         QMIN2 = (Z*AMASS)**2.D0
         ARG   = 1.D0+0.71D0/QMIN2
         F = DLOG(ARG)-11.D0/6.D0+3.D0/ARG
         F = F-3.D0/(2.D0*ARG*ARG)+1.D0/(3.D0*ARG*ARG*ARG)
         F = F*ALPHAE/PI
         F = F*(1.D0-Z+(1.D0/2.D0)*Z*Z)
         F = F/Z
c         print *,'*** data:',z,exparg,ei,f
      ELSEIF (NFLUX.EQ.17) THEN
C---Photon flux from heavy ions (Ca, Pb):
c ... T.K. : Implemented (6) in Drees, Ellis, Zeppenfeld; PL B223 (1989) 455
         QZERO=6.D-2 ! best fit
         EXPARG=(Z*AION*AMASS/QZERO)**2
         EI=DGAGNC(0.D0,EXPARG) / DEXP(EXPARG)
         F=(ALPHAE/PI/Z)*(-DEXP(-EXPARG)+(1.D0+EXPARG)*EI)
         F=F*(ZION**2)
c         print *,'*** data:',z,exparg,ei,f
      ELSEIF (NFLUX.EQ.14) THEN
C---Photon flux in pp (use with ZION=AION=1 only)
c ... M.B. : implement Papageorgiu; PL B250 (1995) 394
c ...        This is equivalent to Cahn-Jackson, but taking R0=0.2 fm 
c0 ...        and relaxing the p-p overlap cut a bit
         RZERO=0.2/FMCONV
         BMIN=1.*RZERO
         XM=PMASS
         ZZERO=1d0/XM/BMIN
         F=2.*ALPHAE/PI/Z*DLOG(ZZERO/Z)
      ELSEIF (NFLUX.EQ.15) THEN
C---Budnev photon flux : in principle this is most precise for pp
         Q2MIN = Z*Z*0.88d0/(1d0-Z)/QSCALE
         Q2MAX = TMAX/QSCALE
         F = ALPHAE/PI*(1d0-Z)/Z*(PHI(Q2MAX,Z)-PHI(Q2MIN,Z))
c nflux=20,22 CHR/Yura	 
      ELSEIF(NFLUX.EQ.20) THEN
         IF(IND.EQ.1) THEN
         Q2MIN = Z*Z*0.88d0/(1d0-Z)/QSCALE
         Q2MAX = TMAX/QSCALE
         F = ALPHAE/PI*(1d0-Z)/Z*(PHI(Q2MAX,Z)-PHI(Q2MIN,Z))
         ELSEIF(IND.EQ.2) THEN
         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
         F = X*W*V
	 ENDIF
      ELSEIF(NFLUX.EQ.22) THEN
         IF(IND.EQ.2) THEN
         Q2MIN = Z*Z*0.88d0/(1d0-Z)/QSCALE
         Q2MAX = TMAX/QSCALE
         F = ALPHAE/PI*(1d0-Z)/Z*(PHI(Q2MAX,Z)-PHI(Q2MIN,Z))
         ELSEIF(IND.EQ.1) THEN
         V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     +        DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
         W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
         X = 1.D0/(Z**(2.D0*alphaP-1.D0))
         F = X*W*V
	 ENDIF
      ELSEIF (NFLUX.EQ.16) THEN
C---KMR flux
         IF(IND.EQ.1) THEN
            Z1KMR = Z
            Z2KMR = 0d0
            F = 1d0
         ELSEIF(IND.EQ.2) THEN
            Z2KMR = Z
            CALL KMRINT(Z1KMR,Z2KMR,FLKMR)
            F=FLKMR/Z1KMR/Z2KMR
         ENDIF      
C ... begin R.S.
      ELSEIF(NFLUX.EQ.18) THEN
         F = 1d0
C ... end R.S.
C CR YURA POm Reg
      ELSEIF(NFLUX.EQ.19) THEN
         IF(IND.EQ.1) THEN
	 	V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     + 		  DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
     		W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
		X = 1.D0/(Z**(2.D0*alphaP-1.D0))
		F = X*W*V
c          print *,'F ind 1 :',Z,F
c...Cox-Forshaw reggeon flux:
         ELSEIF (IND.EQ.2) THEN  
                V = DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMIN)-
     +            DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMAX)
                W = 1.D0/(Breg+2.D0*alphaRp*DLOG(1.D0/Z))
                X = 1.D0/(Z**(2.D0*alphaR-1.D0))
                F = Cr*X*W*V  
c                print *,'F IND 2 :',Z,F
	 ENDIF       
      ELSEIF(NFLUX.EQ.21) THEN
         IF(IND.EQ.2) THEN
	 	V = DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMIN)-
     + 		  DEXP(-(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))*TMAX)
     		W = 1.D0/(Bpom+2.D0*alphaPp*DLOG(1.D0/Z))
		X = 1.D0/(Z**(2.D0*alphaP-1.D0))
		F = X*W*V
c...Cox-Forshaw reggeon flux:
         ELSEIF (IND.EQ.1) THEN  
                V = DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMIN)-
     +            DEXP(-(Breg+2.D0*alphaRp*DLOG(1.D0/Z))*TMAX)
                W = 1.D0/(Breg+2.D0*alphaRp*DLOG(1.D0/Z))
                X = 1.D0/(Z**(2.D0*alphaR-1.D0))
                F = Cr*X*W*V  
	 ENDIF       
      ELSE
         WRITE(*,*) 'In FLUX: NFLUX must be 9-16,18-26 in FPMC!'
         STOP
      ENDIF

      RETURN
      END
C------------------------------------------------------------------------

C Code for the Budnev flux (thanks K.Piotrzkowski)
      double precision function phi(q2,z)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      data a/7.16/,b/-3.96/,c/0.028/
      
      y = z*z/(1-z)
      sum1 = 0.
      sum2 = 0.
      do k=1,3
         sum1=sum1+1./(k*(1.+q2)**k)
         sum2=sum2+b**k/(k*(1.+q2)**k)
      enddo
      
      phi = (1.+a*y)*(-dlog(1+1/q2)+sum1)+(1-b)*y/(4*q2*(1+q2)**3)
      phi = phi +c*(1+0.25*y)*(dlog((1+q2-b)/(1+q2))+sum2)
      
      end
C------------------------------------------------------------------------
C Non-integrated Budnev flux (for proton Q2 dependence)
C------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GESQ(Q2) 
C     Electric proton formfactor
C------------------------------------------------------------------------
      DOUBLE PRECISION Q2, Q0SCALE
      PARAMETER (Q0SCALE=0.71d0)
      GESQ =1.D0/(1.D0+Q2/Q0SCALE)**4
      RETURN 
      END

C------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GMSQ(Q2) 
C     Magnetic proton formfactor
C------------------------------------------------------------------------
      DOUBLE PRECISION Q2, MUSQ, GESQ
      EXTERNAL GESQ
      PARAMETER (MUSQ=7.78D0)
      GMSQ =  MUSQ*GESQ(Q2)
      RETURN
      END

C------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BUDNEVFLUX(Z, Q2)
C     Non-integrated Budnev flux (Oldrich Kepka)
C     Is unitless, shoud be multiplied with 1/ebeam when used
C------------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION Z, Q2, FE, FM, GESQ, GMSQ, NONINTFLUX
      DOUBLE PRECISION Q2MINPROT
      EXTERNAL GESQ, GMSQ
      FE=(4D0*0.88D0*GESQ(Q2)+Q2*GMSQ(Q2))/(4D0*0.88D0+Q2)
      FM=GMSQ(Q2);
      Q2MINPROT = Z*Z*0.88d0/(1d0-Z)
      NONINTFLUX = 1.D0/Z*ALPHEM/PIFAC/Q2*
     &    ((1.D0-Z)*(1.D0-Q2MINPROT/Q2)*FE + 0.5D0*Z**2*FM);
      IF(NONINTFLUX.GT.0) THEN
         BUDNEVFLUX = NONINTFLUX
      ELSE
         BUDNEVFLUX = 0D0
      ENDIF
      RETURN
      END

CDECK  ID>, HWEGAM.
*CMZ :-        -26/04/91  11.11.55  by  Bryan Webber
*-- Author :    Bryan Webber & Luca Stanco
*Modified to substitue pomeron (reggeon) flux for photon flux
* B. Cox and J. Forshaw 11/07/00
* Bug fix for double pomeron 12/06/01 B. Cox
* Modified for herwig65 21/11/02 B. Cox
* Modified for exclusive production 8/8/3 M.Boonekamp T. Kucs
C-----------------------------------------------------------------------
      SUBROUTINE HWEGAM(IHEP,ZMI,ZMA,WWA)
C-----------------------------------------------------------------------
C     GENERATES A PHOTON IN WEIZSACKER-WILLIAMS (WWA=.TRUE.) OR
C     ELSE EQUIVALENT PHOTON APPROX FROM INCOMING E+, E-, MU+ OR MU-
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION GAPSPRREL
      COMMON /BLMODIF/  GAPSPRREL

      DOUBLE PRECISION HWRGEN,HWRUNI,EGMIN,ZMIN,ZMAX,ZGAM,SS,ZMI,ZMA,
     & PPL,PMI,QT2,Q2,QQMIN,QQMAX,S0,A,RPM(2),HWRGAU
      INTEGER IHEP,IHADIS,HQ,I
      LOGICAL WWA
      DOUBLE PRECISION F,C,FN      
      DATA EGMIN/5.D0/
C ... M.B. T.K. 
      CHARACTER(80) STABLE
c      DOUBLE PRECISION Bp,ALPHAPp,Br,ALPHAPr,Bbl,ALPHAPbl,zh1,Bkmr
      DOUBLE PRECISION Bbl,ALPHAPbl,Bkmr
      DOUBLE PRECISION Q0SQ,QPSQ,Q2PRIM,RNGEN,PXGAM,PYGAM,DPHI,SOFTC
c      old        
c      PARAMETER (ALPHAPp =0.26,Bp =4.6,zh1=0.003)
c      PARAMETER (ALPHAPr =0.90,Br =2.0)
c     new
c      PARAMETER (ALPHAPp =0.06,Bp =5.5,zh1=0.009)
c      PARAMETER (ALPHAPr =0.30,Br =1.6)
c      PARAMETER (ALPHAPbl=0.25,Bbl=6.0)
      PARAMETER (Bkmr=4.0)
      PARAMETER (Q0SQ=0.06d0**2)
      PARAMETER (QPSQ=0.027d0) ! = (0.197/1.2)**2       
      DOUBLE PRECISION ZGAM1,ZGAM2,TGAM1,TGAM2,PHIGAM1,PHIGAM2
      SAVE ZGAM1,ZGAM2,TGAM1,TGAM2,PHIGAM1,PHIGAM2
      EXTERNAL HWRGEN,HWRUNI,HWRGAU
      INTEGER ICHI
      LOGICAL EXCKINE
      DOUBLE PRECISION Q2MINPROT, RNDFLUX, BUDNEVFLUX, BFLUX, MAXFLUX
      EXTERNAL BUDNEVFLUX
      LOGICAL LOOP
c M.R. ExcChi
      DOUBLE PRECISION THETA,PHIGAM,MASS2,THETA2,THETA1,K0,K1,OMEGA,Q2WT

C ... begin R.S.

      include 'CHIDe.inc' 
      DOUBLE PRECISION CHIDe_C, CHIDe_B
      DOUBLE PRECISION S, AA, BB, CC, DELTA

      LOGICAL IsCHIDeHiggs, IsCHIDeGG, IsCHIDe, IsCHIDeDiphoton
      DOUBLE PRECISION CHIDeKmax, CHIDeK2max
      DOUBLE PRECISION CHIDedotdiff
      EXTERNAL CHIDedotdiff

      DOUBLE PRECISION CHIDeQ, CHIDeQp,CHIDeSigma,SIGMA
      DOUBLE PRECISION CHIDePhi, CHIDePhip
      DOUBLE PRECISION CHIDeK(2), CHIDeKp(2)
      integer II,N
      integer IND
      COMMON /YURA/ IND
      
      IsCHIDe = (NFLUX.EQ.18)
      IsCHIDeGG = (IPROC.EQ.16012 .AND. NFLUX.EQ.18)
      IsCHIDeDiphoton = (IPROC.EQ.16059 .AND. NFLUX.EQ.18)
      IsCHIDeHiggs = (IPROC.EQ.19999 .AND. NFLUX.EQ.18)
      
C ... end R.S.


c ... pomeron/regeon structure from the common block
      ALPHAPbl = alphaPp
      Bbl = Bpom


C YURA/CHR IND=IHEP!!!!
      IND = IHEP

      
C ... end M.B. T.K.

      IF (IERROR.NE.0)  RETURN 
      IF (IHEP.LT.1.OR.IHEP.GT.2) CALL HWWARN('HWEGAM',500,*999)
      SS=PHEP(5,3)
      IF (IHEP.EQ.1) THEN
        IHADIS=2
      ELSE
        IHADIS=1
        IF (JDAHEP(1,IHADIS).NE.0) IHADIS=JDAHEP(1,IHADIS)
      ENDIF

C *** Murilo Rangel ExcChi: added IPRO=96,97 (exclusive chi production)
      IF(IPRO.EQ.96.OR.IPRO.EQ.97) ICHI=MOD(IPROC,100)
      EXCKINE=(TYPEPR.EQ.'EXC'.AND.(IPRO.EQ.97))

C Initialize Soft correction
      SOFTC=1d0

C---Exclusive 2->1 : second call, recall ZGAM2, skip computations
C---Works only for gluglu/gamgam -> Higgs, Chi_c/Chi_b
C---additions from Murilo Rangel for correct kinematics
C ... M.B.
      IF(TYPEPR.EQ.'EXC'.AND.IHEP.EQ.2.AND.IPRO.EQ.99) THEN
         ZGAM = ZGAM2
         GOTO 998
      ENDIF
      IF(EXCKINE.AND.IHEP.EQ.2) GOTO 998 !M.R. ExcChi
C ... end M.B.


C---DEFINE LIMITS FOR GAMMA MOMENTUM FRACTION
      IF (ZMI.LE.0D0 .OR. ZMA.GT.1D0) THEN
        IF (IPRO.EQ.13.OR.IPRO.EQ.14) THEN
          S0 = EMMIN**2
C ... M.B. : added IPRO=98 below (exclusive photon pair production)
        ELSEIF(IPRO.EQ.15.OR.IPRO.EQ.18.OR.IPRO.EQ.22.OR.IPRO.EQ.24.OR.
     &         IPRO.EQ.50.OR.IPRO.EQ.53.OR.IPRO.EQ.55.OR.IPRO.EQ.98)THEN
          S0 = 4.D0*PTMIN**2
C ... end M.B.
        ELSEIF (IPRO.EQ.17.OR.IPRO.EQ.51) THEN
          HQ = MOD(IPROC,100)
          S0 = 4.D0*(PTMIN**2+RMASS(HQ)**2)
C ... M.B. : added IPRO=99 below (exclusive Higgs boson production)
        ELSEIF (IPRO.EQ.16.OR.IPRO.EQ.19.OR.IPRO.EQ.95
     &                                  .OR.IPRO.EQ.99) THEN
          S0 = MAX(2*RMASS(1),RMASS(201)-GAMMAX*GAMH)**2
C ... end M.B.
        ELSEIF (IPRO.EQ.23) THEN
          S0 = MAX(2*RMASS(1),RMASS(201)-GAMMAX*GAMH)**2
          S0 = (PTMIN+SQRT(PTMIN**2+S0))**2
        ELSEIF (IPRO.EQ.20) THEN
          S0 = RMASS(201)**2
        ELSEIF (IPRO.EQ.21) THEN
          S0 = (PTMIN+SQRT(PTMIN**2+RMASS(198)**2))**2
C--PR MOD 7/7/99
        ELSEIF(IPRO.EQ.30) THEN
          S0 = 4.0D0*(PTMIN**2+RMMNSS**2)
        ELSEIF(IPRO.EQ.40.OR.IPRO.EQ.41) THEN
          HQ = IPROC-100*IPRO
          RPM(1) = RMMNSS
          RPM(2) = ZERO
          IF(HQ.GE.10.AND.HQ.LT.20) THEN
            RPM(1) = ABS(RMASS(450))
            IF(HQ.GT.10) RPM(1) = ABS(RMASS(449+MOD(HQ,10)))
          ELSEIF(HQ.GE.20.AND.HQ.LT.30) THEN
            RPM(1) = ABS(RMASS(454))
            IF(HQ.GT.20) RPM(1) = ABS(RMASS(453+MOD(HQ,20)))
          ELSEIF(HQ.EQ.30) THEN
            RPM(1) = RMASS(449)
          ELSEIF(HQ.EQ.40) THEN
            IF(IPRO.EQ.40) THEN
              RPM(1) = RMASS(425)
              DO I=1,5
                RPM(1) = MIN(RPM(1),RMASS(425+I))
              ENDDO
            ELSE
              RPM(1) = MIN(RMASS(405),RMASS(406))
            ENDIF
            RPM(2) = RMASS(198)
          ELSEIF(HQ.EQ.50) THEN
            IF(IPRO.EQ.40) THEN
              RPM(1) = RMASS(425)
              DO I=1,5
                RPM(1) = MIN(RPM(1),RMASS(425+I))
              ENDDO
              DO I=1,3
                RPM(2) = MIN(RPM(1),RMASS(433+2*I))
              ENDDO
              RPM(1) = MIN(RPM(1),RPM(2))
              RPM(2) = RMASS(203)
              DO I=1,2
                RPM(2) = MIN(RPM(2),RMASS(204+I))
              ENDDO
            ELSE
              RPM(1) = RMASS(401)
              RPM(2) = RMASS(413)
              DO I=1,5
                RPM(1) = MIN(RPM(1),RMASS(401+I))
                RPM(2) = MIN(RPM(2),RMASS(413+I))
              ENDDO
              RPM(1) = MIN(RPM(1),RPM(2))
              RPM(2) = RMASS(203)
              DO I=1,2
                RPM(2) = MIN(RPM(2),RMASS(204+I))
              ENDDO
            ENDIF
            RPM(2) = RMASS(203)
            DO I=1,2
              RPM(2) = MIN(RPM(2),RMASS(204+I))
            ENDDO
          ELSEIF(HQ.GE.60) THEN
            RPM(1) = ZERO
          ENDIF
          RPM(1) = RPM(1)**2
          RPM(2) = RPM(2)**2
          S0 = RPM(1)+RPM(2)+TWO*(PTMIN**2+
     &         SQRT(RPM(1)*RPM(2)+PTMIN**2*(RPM(1)+RPM(2)+PTMIN**2)))
C--end of mod
        ELSEIF (IPRO.EQ.52) THEN
          HQ = MOD(IPROC,100)
          S0 = (PTMIN+SQRT(PTMIN**2+RMASS(HQ)**2))**2
C ... M.B. : add IPRO=60 below
        ELSEIF (IPRO.EQ.60) THEN
          S0 = 4d0*PTMIN**2
C ... end M.B.
        ELSEIF (IPRO.EQ.80) THEN
          S0 = WHMIN**2
        ELSEIF (IPRO.EQ.90) THEN
          S0 = Q2MIN
        ELSEIF (IPRO.EQ.91.OR.IPRO.EQ.92) THEN
          S0 = Q2MIN+4.D0*PTMIN**2
          HQ = MOD(IPROC,100)
          IF (HQ.GT.0) S0 = S0+4.D0*RMASS(HQ)**2
          IF (IPRO.EQ.91) S0 = MAX(S0,EMMIN**2)
        ELSEIF (IPRO.EQ.96.OR.IPRO.EQ.97) THEN
          IF(ICHI.EQ.1) THEN
            S0=RMASS(165)**2
          ELSEIF(ICHI.EQ.2) THEN
            S0=RMASS(306)**2
          ENDIF          
        ELSE
          S0 = 0
        ENDIF
        IF (S0.GT.0) THEN
          S0 = (SQRT(S0)+ABS(PHEP(5,IHADIS)))**2-PHEP(5,IHADIS)**2
          S0 = MAX(S0,WHMIN**2)
          ZMIN = S0 / (SS**2 - PHEP(5,IHEP)**2 - PHEP(5,IHADIS)**2)
          ZMAX = ONE
C ... M.B. : set boundaries properly for Exclusive case
          IF(TYPEPR.EQ.'EXC'.AND.(IPRO.EQ.16.OR.IPRO.EQ.99)) THEN
             ZMAX = YWWMAX
             ZMIN = RMASS(201)**2/SS**2/ZMAX
             !R.S.
             XI1Min = max(XI1Min,RMASS(201)**2/SS**2/XI2Max)
             XI2Min = min(XI2Min,RMASS(201)**2/SS**2/XI1Max)
          ENDIF
          IF(EXCKINE) THEN  ! M.R. ExcChi
            ZMAX = YWWMAX
            IF(ICHI.EQ.1) ZMIN = RMASS(165)**2/SS**2/ZMAX
            IF(ICHI.EQ.2) ZMIN = RMASS(306)**2/SS**2/ZMAX
          ENDIF
C ... end M.B.
        ELSE
C---UNKNOWN PROCESS: USE ENERGY CUTOFF, AND WARN USER
          IF (FSTWGT) CALL HWWARN('HWEGAM',1,*999)
          ZMIN = EGMIN / PHEP(4,IHEP)
          ZMAX = ONE
        ENDIF
      ELSE
        ZMIN=ZMI
        ZMAX=ZMA
      ENDIF

C---APPLY USER DEFINED CUTS YWWMIN,YWWMAX AND INDIRECT LIMITS ON Z
      IF (.NOT.WWA) THEN
        ZMIN=MAX(ZMIN,YWWMIN,SQRT(Q2WWMN)/ABS(PHEP(3,IHEP)))

C ... begin R.S.
       IF(IHEP.EQ.1 .AND. (NFLUX.EQ.18 .OR. NFLUX.EQ.16)) THEN
         ZMIN=XI1Min
         ZMax=XI1Max
       ELSEIF(IHEP.EQ.2 .AND. (NFLUX.EQ.18 .OR. NFLUX.EQ.16)) THEN
         ZMIN=XI2Min
         ZMax=XI2Max
       ENDIF 
C ... end R.S.        
        ZMAX=MIN(ZMAX,YWWMAX)
        IF (ZMIN.GT.ZMAX) THEN
          GAMWT=ZERO
          RETURN
        ENDIF
      ENDIF
C---Q2WWMN AND Q2WWMX ARE USER-DEFINED LIMITS IN THE Q**2 INTEGRATION
      QQMAX=Q2WWMX
      QQMIN=Q2WWMN
C     IF (QQMIN.GT.QQMAX) CALL HWWARN('HWEGAM',50,*10)
C---GENERATE GAMMA MOMENTUM FRACTION
      ZGAM=(ZMIN/ZMAX)**HWRGEN(1)*ZMAX
c tempo debug CHR
c      ZGAM=0.01D0      
      C=1.D0/DLOG(ZMAX/ZMIN)
C---POMERON (REGGEON) FLUX ; calculate GAMWT     
C CHR comment this line and put it in each flux to be simpler/clearer
cc      CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
C changes CHR/YURA
      IF ((NFLUX.EQ.9).OR.(NFLUX.EQ.10)) THEN
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP) 
         CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/(C*FN*zh1)   
C CHR NFLUX=9 should be also true for nflux=10
      ELSEIF ((NFLUX.GE.11).AND.(NFLUX.LE.16)) THEN 
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/C
C CHR Yura modif 10/2013 Pom Reg
      ELSEIF (NFLUX.EQ.19) THEN
c         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
c         IF(IND.EQ.1) THEN
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1)   
c         ELSEIF(IND.EQ.2) THEN
c             GAMWT = GAMWT*F*ZGAM/C
c	 ENDIF 
C CHR Yura modif Reg Pom              
      ELSEIF (NFLUX.EQ.21) THEN
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1)   
c         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
c         IF(IND.EQ.2) THEN
c             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
c             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1)   
c         ELSEIF(IND.EQ.1) THEN
c             GAMWT = GAMWT*F*ZGAM/C
c	 ENDIF               
C CHR Yura modif 10/2013 gamma+P
      ELSEIF (NFLUX.EQ.20) THEN
c         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         IF(IND.EQ.2) THEN
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1) 
c	     print *,'f fn gamwt :',f,fn,gamwt
	 ELSEIF(IND.EQ.1) THEN      
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/C
c	     print *,'f fn gamwt ind 1:',f,fn,gamwt
	 ENDIF    
C CHR Yura P + gamma           
      ELSEIF (NFLUX.EQ.22) THEN
         IF(IND.EQ.1) THEN
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1) 
c	     print *,'f fn gamwt :',f,fn,gamwt
	 ELSEIF(IND.EQ.2) THEN      
             CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
             GAMWT = GAMWT*F*ZGAM/C
c	     print *,'f fn gamwt ind 2:',f,fn,gamwt
	 ENDIF    
C CHR ion proton-ion gamma gamma
      ELSEIF (NFLUX.EQ.23) THEN 
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/C
C CHR ion proton-ion gamma gamma
      ELSEIF (NFLUX.EQ.24) THEN 
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/C
C CHR ion proton-ion pomeron gamma
      ELSEIF (NFLUX.EQ.25) THEN 
         IF(IND.EQ.1) THEN
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/(C*FN*zh1) 
	 ELSEIF(IND.EQ.2) THEN      
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/C
	 ENDIF    
C CHR ion proton-ion pomeron gamma
      ELSEIF (NFLUX.EQ.26) THEN 
         IF(IND.EQ.2) THEN
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/(C*FN*zh1) 
	 ELSEIF(IND.EQ.1) THEN      
         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
         GAMWT = GAMWT*F*ZGAM/C
	 ENDIF    


c         CALL FLUX(F,ZGAM,QQMIN,QQMAX,IPRO,IHEP)
c         IF(IND.EQ.1) THEN
c             CALL FLUX(FN,ZH1,QQMIN,QQMAX,IPRO,IHEP)
c             GAMWT = GAMWT*F*ZGAM/(C*FN*zh1)   
c	 ENDIF               
C ... begin R.S.
C     In CHIDe model the integrand is not factorised into flux and hard
C     process, therefore there is no flux here.
      ELSEIF(IsCHIDeHiggs .AND. IHEP.EQ.1) THEN
        CHIDe_B = 4d0! Only matters for generation efficiency
        CHIDe_C = CHIDe_B/(DEXP(-CHIDe_B*QQMIN)-DEXP(-CHIDe_B*QQMAX))
        CHIDeQQ1 = -1.0/CHIDe_B*DLOG(DEXP(-CHIDe_B*QQMIN) 
     &            - CHIDe_B/CHIDe_C*HWRGEN(1))
        CHIDeQQ3 = -1.0/CHIDe_B*DLOG(DEXP(-CHIDe_B*QQMIN) 
     &            - CHIDe_B/CHIDe_C*HWRGEN(2))
        GAMWT = pi*GAMWT/CHIDe_C/DEXP(-CHIDe_B*CHIDeQQ1)
        GAMWT = pi*GAMWT/CHIDe_C/DEXP(-CHIDe_B*CHIDeQQ3)
        CHIDePhi1 = 2*pi*HWRGEN(1)
        CHIDePhi3 = 2*pi*HWRGEN(2)
        
        CHIDeQ1 = sqrt(CHIDeQQ1)
        CHIDeQ3 = sqrt(CHIDeQQ3)
                
        CHIDeK1(1) =  CHIDeQ1*cos(CHIDePhi1)
        CHIDeK1(2) =  CHIDeQ1*sin(CHIDePhi1)
        CHIDeK3(1) =  CHIDeQ3*cos(CHIDePhi3)
        CHIDeK3(2) =  CHIDeQ3*sin(CHIDePhi3)
        
        CHIDeZ1 = (ZMIN/ZMAX)**HWRGEN(1)*ZMAX
        GAMWT = GAMWT*CHIDeZ1*DLOG(ZMAX/ZMIN)

C     Exact kinematics (from A.Dechambre)           
c       S = (PBEAM1+PBEAM2)**2
c       AA = S*CHIDeZ1
c       BB = S*CHIDeZ1 - CHIDeQQ3 + CHIDeZ1/(1.-CHIDeZ1)*CHIDeQQ1
c    &      - (CHIDeK1(1)-CHIDeK1(2))**2 - (CHIDeK2(1)-CHIDeK2(2))**2 
c    &      + RMASS(201)**2
c       CC = CHIDeZ1/(1.-CHIDeZ1)*CHIDeQQ1
c    &      - (CHIDeK2(1)-CHIDeK2(2))**2 + RMASS(201)**2
c       DELTA = BB**2-(4.*AA*CC)
c
c       IF(DELTA.LE.ZERO) GAMWT = ZERO
c       CHIDeZ2 = (-BB + sqrt(DELTA))/2./AA
        
        CHIDeZ2 = (RMASS(201)/SS)**2 / CHIDeZ1
c       print*, RMASS(201), sqrt(CHIDeZ1*CHIDeZ2)*SS
c       stop
        ZGAM = CHIDeZ1
      ELSEIF( (IsCHIDeGG.OR.IsCHIDeDiphoton) .AND. IHEP.EQ.1) THEN
        CHIDeK2max =PTMIN+100d0
        CHIDe_B = 0.2d0! Only matters for generation efficiency
        CHIDe_C=CHIDe_B/(DEXP(-CHIDe_B*PTMIN)-DEXP(-CHIDe_B*CHIDeK2max))
        CHIDeQ2 = -1.0/CHIDe_B*DLOG(DEXP(-CHIDe_B*PTMIN) 
     &            - CHIDe_B/CHIDe_C*HWRGEN(1))
        GAMWT = 2*pi*GAMWT*CHIDeQ2/CHIDe_C/DEXP(-CHIDe_B*CHIDeQ2)
        CHIDePhi2 = 2*pi*HWRGEN(2)
        CHIDeQQ2 = CHIDeQ2**2
       
        CHIDe_B = 5d0! Only matters for generation efficiency
        CHIDe_C = CHIDe_B/(DEXP(-CHIDe_B*QQMIN)-DEXP(-CHIDe_B*QQMAX))
        CHIDeQQ1 = -1.0/CHIDe_B*DLOG(DEXP(-CHIDe_B*QQMIN) 
     &            - CHIDe_B/CHIDe_C*HWRGEN(1))
        CHIDeQQ3 = -1.0/CHIDe_B*DLOG(DEXP(-CHIDe_B*QQMIN) 
     &            - CHIDe_B/CHIDe_C*HWRGEN(2))
        GAMWT = pi*GAMWT/CHIDe_C/DEXP(-CHIDe_B*CHIDeQQ1)
        GAMWT = pi*GAMWT/CHIDe_C/DEXP(-CHIDe_B*CHIDeQQ3)
        CHIDePhi1 = 2*pi*HWRGEN(1)
        CHIDePhi3 = 2*pi*HWRGEN(2)
        CHIDeQ1 = sqrt(CHIDeQQ1)
        CHIDeQ3 = sqrt(CHIDeQQ3)
        
        CHIDeK1(1) =  CHIDeQ1*cos(CHIDePhi1)
        CHIDeK1(2) =  CHIDeQ1*sin(CHIDePhi1)
        CHIDeK2(1) =  CHIDeQ2*cos(CHIDePhi2)
        CHIDeK2(2) =  CHIDeQ2*sin(CHIDePhi2)
        CHIDeK3(1) =  CHIDeQ3*cos(CHIDePhi3)
        CHIDeK3(2) =  CHIDeQ3*sin(CHIDePhi3)
       
        CHIDeB1 = (XI1MIN/XI1MAX)**HWRGEN(1)*XI1MAX
        CHIDeB2 = (XI1MIN/XI1MAX)**HWRGEN(2)*XI1MAX
        GAMWT = GAMWT*DLOG(XI1MAX/(XI1MIN))*CHIDeB1
        GAMWT = GAMWT*DLOG(XI1MAX/(XI1MIN))*CHIDeB2

        CHIDeA1=CHIDedotdiff(CHIDeK1,CHIDeK2,CHIDeK1,CHIDeK2)
     &           /CHIDeS/CHIDeB1
        CHIDeA2=CHIDedotdiff(CHIDeK3,CHIDeK2,CHIDeK3,CHIDeK2)
     &           /CHIDeS/CHIDeB2
        
        CHIDeZ1 = CHIDeB1 + CHIDeB2
        CHIDeZ2 = CHIDeA1 + CHIDeA2
        ZGAM = CHIDeZ1
      ELSEIF(IsCHIDe .AND. IHEP.EQ.2) THEN
        ZGAM =  CHIDeZ2
c ... end R.S.
      ELSE
         WRITE(*,*) 'In HWEGAM : NFLUX MUST BE 9-16,18-26 IN FPMC'
         STOP
      ENDIF

C---Exclusive 2->1 process : first call, save ZGAM2, compute total weight      
C---Works only for gluglu/gamgam -> Higgs. Other use cases?
C ... M.B.
      IF(TYPEPR.EQ.'EXC'.AND.IHEP.EQ.1.AND.IPRO.EQ.99) THEN
C ... begin R.S
         IF(IsCHIDeHiggs) THEN
           ZGAM1 = CHIDeZ1
           ZGAM2 = CHIDeZ2
         ELSE
           ZGAM2 = (RMASS(201)/SS)**2 / ZGAM
           ZGAM1 = ZGAM
           CALL FLUX(F,ZGAM2,QQMIN,QQMAX,IPRO,2)
           GAMWT = GAMWT * F * ZGAM2
         ENDIF
        IF(ZGAM.GT.XI1Max.OR.ZGAM2.GT.XI2Max) THEN 
C ... end R.S.
c       IF(ZGAM.GT.YWWMAX.OR.ZGAM2.GT.YWWMAX) THEN 
          PRINT*, '!!! HWEGAM : ZGAM OUT OF RANGE !!!'
        ENDIF
      ENDIF
C ... end M.B.
C *** M.B. M.R. ExcChi: Generate THETA for exclusive production
 998  IF (EXCKINE) THEN
         C=4.d0*PHEP(4,IHEP)**2
         THETA=DSQRT((1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+
     &   HWRGEN(1)*(DEXP(-C*QQMIN)-DEXP(-C*QQMAX)))))
         PHIGAM=HWRGEN(1)*2*PI ! Generate phi randomly
C *** M.B. M.R. Compute ZGAM in second call for exclusive production
        IF(IHEP.EQ.2) THEN
          IF(IPRO.EQ.99) THEN      ! Higgs
            MASS2 = RMASS(201)**2
          ELSEIF(IPRO.EQ.97) THEN  ! Chi_c, Chi_b
            IF(ICHI.EQ.1) MASS2 = RMASS(165)**2
            IF(ICHI.EQ.2) MASS2 = RMASS(306)**2
          ENDIF
          IF (NFLUX.EQ.11) THEN
            THETA2=PI-THETA
            K0=PHEP(4,IHEP)
            K1=K0*(ONE-ZGAM1)/DCOS(THETA1)
            OMEGA=ONE-DCOS(THETA1)*DCOS(THETA2)-DSIN(THETA1)*
     &      DSIN(THETA2)*(DCOS(PHIGAM1)*DCOS(PHIGAM)+
     &      DSIN(PHIGAM1)*DSIN(PHIGAM))

            ZGAM=ONE-DCOS(THETA)*(MASS2-4*K0**2+4*K0*K1)/
     &      (K0*(-4*K0+2*K1*OMEGA))
     
          ELSE
            ZGAM=(MASS2/SS**2)/ZGAM1
            PRINT*,'Exclusive CHI non-BL not yet tested'
            PRINT*,'Use flux 11'
            STOP
          ENDIF
          ZGAM2=ZGAM 

          CALL FLUX(F,ZGAM,QQMIN,QQMAX,NFLUX,IPRO)
          IF(.NOT. IsCHIDe) GAMWT=GAMWT*F*ZGAM
          IF(ZGAM.LT.YWWMIN) GAMWT = 0d0
          IF(ZGAM.GT.YWWMAX.OR.ZGAM.LT.YWWMIN) THEN
            PRINT*, ''
            PRINT*,' !!! HWEGAM : ZGAM OUT OF RANGE !!!'
            CALL HWWARN('HWEGAM',100,*999)
          ENDIF      
        ENDIF
      ENDIF
C *** M.B. M.R. ExcChi end

C ... MB 03/05 : pick and save t (=Q2)        
C ... M.B. T.K.: choice depends on model considered
      IF(NFLUX.EQ.9) THEN        
C ...... H1 Pomeron
        C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
        Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
      ELSEIF(NFLUX.EQ.10) THEN   
C ...... H1 Reggeon
        C=Breg+2.0*alphaRp*DLOG(1.D0/ZGAM)
        Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
      ELSEIF(NFLUX.EQ.11.AND.EXCKINE) THEN !M.R. ExcChi
C...... Bialas-Landshoff Pomeron
        Q2=THETA**2*(ONE-ZGAM)*(SS**2)/4
      ELSEIF(NFLUX.EQ.11) THEN   
C...... Bialas-Landshoff Pomeron
         C=Bbl+2.0*ALPHAPbl*DLOG(1.D0/ZGAM)
        Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
      ELSEIF((NFLUX.EQ.12.OR.NFLUX.EQ.13).AND.ZION.GE.10) THEN   
C...... QED Heavy Ion : exponential Form Factor, (60 MeV)**2 slope 
C ..... (what about Z,A dependence?)
         Q2=1.D-3

C CHR proton ion
      ELSEIF((NFLUX.EQ.23).AND.ZION.GE.10) THEN   
C...... QED Heavy Ion : exponential Form Factor, (60 MeV)**2 slope 
C ..... (what about Z,A dependence?)
       IF(IND.EQ.1) THEN   
C O.K.... QED Proton : "Dirac Form Factor"
c         NFLUX = 15 according to Budnev flux, Rejection Technique
         Q2MINPROT = ZGAM*ZGAM*0.88d0/(1d0-ZGAM)
        IF(QQMIN.LT.Q2MINPROT) QQMIN = Q2MINPROT
        IF(QQMIN.GE.QQMAX) THEN
            GAMWT=ZERO
            RETURN
        ENDIF    
            
        LOOP=.TRUE.
        DO WHILE(LOOP)
         RNGEN=HWRGEN(1) 
         Q2 = (QQMIN**RNGEN)*(QQMAX**(1-RNGEN));
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
         MAXFLUX = HWRUNI(1, ZERO, 
     &                     1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2)
         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
            IF(RNDFLUX.LT.BFLUX) THEN
               LOOP=.FALSE.
            ENDIF   
        END DO
        ELSEIF(IND.EQ.2) THEN
         Q2=1.D-3
	ENDIF


C CHR  ion proton
      ELSEIF((NFLUX.EQ.24).AND.ZION.GE.10) THEN   
C...... QED Heavy Ion : exponential Form Factor, (60 MeV)**2 slope 
C ..... (what about Z,A dependence?)
       IF(IND.EQ.2) THEN   
C O.K.... QED Proton : "Dirac Form Factor"
c         NFLUX = 15 according to Budnev flux, Rejection Technique
         Q2MINPROT = ZGAM*ZGAM*0.88d0/(1d0-ZGAM)
        IF(QQMIN.LT.Q2MINPROT) QQMIN = Q2MINPROT
        IF(QQMIN.GE.QQMAX) THEN
            GAMWT=ZERO
            RETURN
        ENDIF    
            
        LOOP=.TRUE.
        DO WHILE(LOOP)
         RNGEN=HWRGEN(1) 
         Q2 = (QQMIN**RNGEN)*(QQMAX**(1-RNGEN));
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
         MAXFLUX = HWRUNI(1, ZERO, 
     &                     1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2)
         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
            IF(RNDFLUX.LT.BFLUX) THEN
               LOOP=.FALSE.
            ENDIF   
        END DO
        ELSEIF(IND.EQ.1) THEN
         Q2=1.D-3
	ENDIF

C CHR proton ion Pomeron gamma
      ELSEIF((NFLUX.EQ.25).AND.ZION.GE.10) THEN   
C...... QED Heavy Ion : exponential Form Factor, (60 MeV)**2 slope 
C ..... (what about Z,A dependence?)
       IF(IND.EQ.1) THEN   
C ......   H1 Pomeron
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
       ELSEIF(IND.EQ.2) THEN
         Q2=1.D-3
	ENDIF



C CHR ion proton Pomeron gamma
      ELSEIF((NFLUX.EQ.26).AND.ZION.GE.10) THEN   
C...... QED Heavy Ion : exponential Form Factor, (60 MeV)**2 slope 
C ..... (what about Z,A dependence?)
       IF(IND.EQ.2) THEN   
C ......   H1 Pomeron
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
       ELSEIF(IND.EQ.1) THEN
         Q2=1.D-3
	ENDIF


      ELSEIF((NFLUX.EQ.12.AND.ZION.EQ.1)
     &        .OR.NFLUX.EQ.14) THEN
C .... QED Proton : "Dirac Form Factor"
        RNGEN=HWRGEN(1) 
        Q2 = (QPSQ**RNGEN)*((ZGAM*ZGAM/(1-ZGAM)*PMASS*PMASS)
     & **(1.-RNGEN))
C Yura/CHR Pom+gamma
      ELSEIF(NFLUX.EQ.20) THEN
       IF(IND.EQ.1) THEN   
C O.K.... QED Proton : "Dirac Form Factor"
c         NFLUX = 15 according to Budnev flux, Rejection Technique
         Q2MINPROT = ZGAM*ZGAM*0.88d0/(1d0-ZGAM)
        IF(QQMIN.LT.Q2MINPROT) QQMIN = Q2MINPROT
        IF(QQMIN.GE.QQMAX) THEN
            GAMWT=ZERO
            RETURN
        ENDIF    
            
        LOOP=.TRUE.
        DO WHILE(LOOP)
         RNGEN=HWRGEN(1) 
         Q2 = (QQMIN**RNGEN)*(QQMAX**(1-RNGEN));
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
         MAXFLUX = HWRUNI(1, ZERO, 
     &                     1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2)
         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
            IF(RNDFLUX.LT.BFLUX) THEN
               LOOP=.FALSE.
            ENDIF   
        END DO
       ELSEIF(IND.EQ.2) THEN
C ......   H1 Pomeron
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
       ENDIF
c       print *,'nflux 20 gamwt :',gamwt
C CHR/Yura Pom photon
      ELSEIF(NFLUX.EQ.22) THEN
       IF(IND.EQ.2) THEN   
C O.K.... QED Proton : "Dirac Form Factor"
c         NFLUX = 15 according to Budnev flux, Rejection Technique
         Q2MINPROT = ZGAM*ZGAM*0.88d0/(1d0-ZGAM)
        IF(QQMIN.LT.Q2MINPROT) QQMIN = Q2MINPROT
        IF(QQMIN.GE.QQMAX) THEN
            GAMWT=ZERO
            RETURN
        ENDIF    
            
        LOOP=.TRUE.
        DO WHILE(LOOP)
         RNGEN=HWRGEN(1) 
         Q2 = (QQMIN**RNGEN)*(QQMAX**(1-RNGEN));
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
         MAXFLUX = HWRUNI(1, ZERO, 
     &                     1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2)
         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
            IF(RNDFLUX.LT.BFLUX) THEN
               LOOP=.FALSE.
            ENDIF   
        END DO
       ELSEIF(IND.EQ.1) THEN
C ......   H1 Pomeron
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
       ENDIF
c       print *,'nflux 22 gamwt :',gamwt



       ELSEIF(NFLUX.EQ.15) THEN   


C O.K.... QED Proton : "Dirac Form Factor"
c      NFLUX = 15 according to Budnev flux, Rejection Technique
        Q2MINPROT = ZGAM*ZGAM*0.88d0/(1d0-ZGAM)
        IF(QQMIN.LT.Q2MINPROT) QQMIN = Q2MINPROT

        IF(QQMIN.GE.QQMAX) THEN
c           print *, QQMIN, QQMAX, ZGAM
            GAMWT=ZERO
            RETURN
        ENDIF    
            
        LOOP=.TRUE.
        DO WHILE(LOOP)
         RNGEN=HWRGEN(1) 
         Q2 = (QQMAX/QQMIN**RNGEN)*QQMIN;
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
         MAXFLUX =1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2
         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
            IF(RNDFLUX.LT.BFLUX) THEN
               LOOP=.FALSE.
            ENDIF   
            IF(MAXFLUX.LT.BFLUX) THEN
               PRINT *, " BFLUX larger than MAXFLUX, PROBLEM"
            ENDIF   
        END DO


c        LOOP=.TRUE.
c        DO WHILE(LOOP)
c         RNGEN=HWRGEN(1) 
c         Q2 = (QQMIN**(1-RNGEN))*(QQMAX**(RNGEN));
c         Q2 = (QQMIN**RNGEN)*(QQMAX**(1-RNGEN));
c         Q2 = (QQMAX/QQMIN**RNGEN)*QQMIN;
c         Q2 = (QQMAX/QQMIN*RNGEN)*(QQMIN)
c         Q2 = (QQMAX-QQMIN)*RNGEN+(QQMIN)
C        function envelope is 1/ebeam/z*alpha/pifac/Q2         
c         MAXFLUX = HWRUNI(1, ZERO, 
c     &                     1D0/PHEP(4, IHEP)/ZGAM*ALPHEM/PIFAC/Q2)
c         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
c         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
c            IF(RNDFLUX.LT.BFLUX) THEN
c               LOOP=.FALSE.
c            ENDIF   
c          LOOP=.FALSE.
c        END DO
c        CALL HF1(1111, REAL(Q2), 1.)
c        CALL HF1(1112, REAL(LOG(Q2)/LOG(10d0)), 1.)

c        LOOP=.TRUE.
c        DO WHILE(LOOP)
c         RNGEN=HWRGEN(1) 
c         Q2 = (QQMAX-QQMIN*RNGEN)+ QQMIN;
cC        function envelope is 1/ebeam/z*alpha/pifac/Q2         
c         MAXFLUX = HWRUNI(1, ZERO, 
c     &                      BUDNEVFLUX(ZGAM, QQMIN) )
c         RNDFLUX = HWRUNI(1, ZERO, MAXFLUX)
c         BFLUX = BUDNEVFLUX(ZGAM, Q2)/PHEP(4, IHEP) 
c            IF(RNDFLUX.LT.BFLUX) THEN
c               LOOP=.FALSE.
c            ENDIF   
c        END DO

      ELSEIF(NFLUX.EQ.16) THEN   
        C=Bkmr
        Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +     (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))

c CHR + Yura Pom Reg
      ELSEIF(NFLUX.EQ.19) THEN        
C ...... H1 Pomeron
        IF(IND.EQ.1) THEN
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +       (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
        ELSEIF(IND.EQ.2) THEN   
C ...... H1 Reggeon
           C=Breg+2.0*alphaRp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +        (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
        ENDIF
C CHR Yura Pom Reg
      ELSEIF(NFLUX.EQ.21) THEN        
C ...... H1 Pomeron
        IF(IND.EQ.2) THEN
           C=Bpom+2.0*alphaPp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +       (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
        ELSEIF(IND.EQ.1) THEN   
C ...... H1 Reggeon
           C=Breg+2.0*alphaRp*DLOG(1.D0/ZGAM)
           Q2=(1.D0/C)*DLOG(1.D0/(DEXP(-C*QQMAX)+HWRGEN(1)*
     +        (DEXP(-C*QQMIN)-DEXP(-C*QQMAX))))
        ENDIF


C ... begin R.S.
      ELSEIF(IsCHIDe) THEN
        IF(IHEP.EQ.1) Q2 = CHIDeQQ1
        IF(IHEP.EQ.2) Q2 = CHIDeQQ3
C ... end R.S.
      ENDIF
      IF(IHEP.EQ.1) THEN
         TGAM1=-Q2
         IF(EXCKINE) THEN
           THETA1=THETA !Save THETA1 for  M.R. ExcChi
           ZGAM1=ZGAM
         ENDIF
      ELSEIF(IHEP.EQ.2) THEN
         TGAM2=-Q2
      ELSE
         PRINT*, 'HWEGAM -- THIS SHOULD NEVER HAPPEN'
         STOP
      ENDIF
C ... end M.B. T.K.

C *** M.R. ExcChi Update GAMWT for exclusive production
      IF(IHEP.EQ.1) THEN
         TGAM1=-Q2
         IF(EXCKINE) THEN
           THETA1=THETA !Save THETA1 for  M.R. ExcChi
           ZGAM1=ZGAM
         ENDIF
      ELSEIF(IHEP.EQ.2) THEN
         TGAM2=-Q2
      ELSE
         PRINT*, 'HWEGAM -- THIS SHOULD NEVER HAPPEN'
         STOP
      ENDIF
C ... end M.B. T.K.

C *** M.R. ExcChi Update GAMWT for exclusive production
      IF (EXCKINE.AND.NFLUX.EQ.11) THEN
         C=PHEP(4,IHEP)**2*(ONE-ZGAM)*
     &   (Bbl+2.0*ALPHAPbl*DLOG(1.D0/ZGAM))
         Q2WT=(DEXP(-C*THETA**2))/
     &        ((1.d0/C)*(DEXP(-C*QQMIN)-DEXP(-C*QQMAX)))
         GAMWT=GAMWT*Q2WT
         C=4.d0*PHEP(4,IHEP)**2
         Q2WT=((1.d0/C)*(DEXP(-C*QQMIN)-
     &        DEXP(-C*QQMAX)))/(DEXP(-C*THETA**2))
         GAMWT=GAMWT*Q2WT
       ENDIF
C *** end M.R.

C ... MB 03/05 : pick and save phi here.

C ... begin R.S.     
      IF(IsCHIDe .AND. IHEP.EQ.1) THEN
        PXGAM = -CHIDeK1(1)
        PYGAM = -CHIDeK1(2)
      ELSEIF(IsCHIDe .AND. IHEP.EQ.2) THEN
        PXGAM = CHIDeK3(1)
        PYGAM = CHIDeK3(2)
      ELSE
        CALL HWRAZM(SQRT(Q2),PXGAM,PYGAM)
      ENDIF
C ... end R.S.
      PHIGAM1=-999d0
      PHIGAM2=-999d0
C ... begin R.S.     
      IF(IsCHIDe) THEN
              IF(IHEP.EQ.1) PHIGAM1 = 2*pi-CHIDePhi1
              IF(IHEP.EQ.2) PHIGAM2 = CHIDePhi3
C ... end R.S.
      ELSEIF(IHEP.EQ.1) THEN
        PHIGAM1=DACOS(PXGAM/SQRT(Q2))
        IF(PYGAM.LT.0d0) PHIGAM1=-PHIGAM1
      ELSEIF(IHEP.EQ.2) THEN
        PHIGAM2=DACOS(PXGAM/SQRT(Q2))
        IF(PYGAM.LT.0d0) PHIGAM2=-PHIGAM2
      ELSE
         PRINT*, 'HWEGAM -- THIS SHOULD NEVER HAPPEN'
         STOP
      ENDIF
      
C ... Apply soft survival probability
      IF(IHEP.EQ.2) THEN
         DPHI=PHIGAM2-PHIGAM1
         IF(DPHI.LT.-PI) DPHI=DPHI+2*PI
         IF(DPHI.GT.+PI) DPHI=DPHI-2*PI
         SOFTC = 1d0
c ...... CF model : survival factor
c CHR YURA change
         IF((NFLUX.GE.9.AND.NFLUX.LE.10).or.(NFLUX.EQ.19).or.
     &       (NFLUX.EQ.21)) THEN
            IF(ISOFTM.EQ.1) THEN
               SOFTC = GAPSPR
            ELSEIF(ISOFTM.GT.1) THEN
               CALL GETSOFTC(SOFTC,TGAM1,TGAM2,DABS(DPHI),2,0)
            ENDIF
            
c ...... BPR model : CDF normalization for inclusive, 
c                    survival factor for exclusive
         ELSEIF(NFLUX.EQ.11) THEN
            IF(TYPEPR.EQ.'INC') THEN
               SOFTC = CDFFAC
c OK 19/04/07 for BL density(beta)=density(beta) new "relative factor"
c             for LHC: 0.03/0.1     
c             for TEV: 1
c             values set in HWMODINI
               IF((IFITPDF.GE.0).AND.(IFITPDF.LT.100)) THEN
                  SOFTC=SOFTC*GAPSPRREL
               ENDIF   
            ELSEIF(TYPEPR.EQ.'EXC') THEN
               IF(ISOFTM.EQ.1) THEN
                  SOFTC = GAPSPR
               ELSEIF(ISOFTM.GT.1) THEN
                  CALL GETSOFTC(SOFTC,TGAM1,TGAM2,DABS(DPHI),2,0)
               ENDIF
            ENDIF
c ...... QED : survival factor
         ELSEIF((NFLUX.GE.12).AND.(NFLUX.LE.15)) THEN
            IF(ISOFTM.EQ.1) THEN
               SOFTC = PROSPR
            ENDIF
         ELSEIF(NFLUX.EQ.16) THEN
c           GAPSPR already applied in the tables 
         ENDIF
         GAMWT = GAMWT * SOFTC
      ENDIF
C... end MB 03/05
      IF (GAMWT.LT.ZERO) GAMWT=ZERO

C---FILL PHOTON
      NHEP=NHEP+1
      IDHW(NHEP)=59
      ISTHEP(NHEP)=3
      IDHEP(NHEP)=22
      JMOHEP(1,NHEP)=IHEP
      JMOHEP(2,NHEP)=0
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
      JDAHEP(1,IHEP)=NHEP
      IF (WWA) THEN
C---FOR COLLINEAR KINEMATICS, ZGAM IS THE ENERGY FRACTION
        PHEP(4,NHEP)=PHEP(4,IHEP)*ZGAM
        PHEP(3,NHEP)=PHEP(3,IHEP)-SIGN(SQRT(
     &     (PHEP(4,IHEP)-PHEP(4,NHEP))**2-PHEP(5,IHEP)**2),PHEP(3,IHEP))
        PHEP(2,NHEP)=0
        PHEP(1,NHEP)=0
        CALL HWUMAS(PHEP(1,NHEP))
      ELSE
C---FOR EXACT KINEMATICS, ZGAM IS TAKEN TO BE FRACTION OF (E+PZ)
        PPL=ZGAM*(ABS(PHEP(3,IHEP))+PHEP(4,IHEP))
        QT2=(ONE-ZGAM)*Q2-(ZGAM*PHEP(5,IHEP))**2
        PMI=(QT2-Q2)/PPL
        PHEP(5,NHEP)=-SQRT(Q2)
        PHEP(4,NHEP)=(PPL+PMI)/TWO
        PHEP(3,NHEP)=SIGN((PPL-PMI)/TWO,PHEP(3,IHEP))
        PHEP(1,NHEP)=PXGAM
        PHEP(2,NHEP)=PYGAM
      ENDIF
C---UPDATE OVERALL CM FRAME
      JMOHEP(IHEP,3)=NHEP
      CALL HWVDIF(4,PHEP(1,3),PHEP(1,IHEP),PHEP(1,3))
      CALL HWVSUM(4,PHEP(1,NHEP),PHEP(1,3),PHEP(1,3))
      CALL HWUMAS(PHEP(1,3))
      
C---FILL OUTGOING LEPTON
      NHEP=NHEP+1
      IDHW(NHEP)=IDHW(IHEP)
      ISTHEP(NHEP)=1
      IDHEP(NHEP)=IDHEP(IHEP)
      JMOHEP(1,NHEP)=IHEP
      JMOHEP(2,NHEP)=0
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
      JDAHEP(2,IHEP)=NHEP
      CALL HWVDIF(4,PHEP(1,IHEP),PHEP(1,NHEP-1),PHEP(1,NHEP))
      PHEP(5,NHEP)=PHEP(5,IHEP)
 999  END
C------------------------------------------------------------------------------
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWEPRO
C-----------------------------------------------------------------------
C     WHEN NEVHEP=0, CHOOSES X VALUES AND FINDS WEIGHT FOR PROCESS IPROC
C     OTHERWISE, CHOOSES AND LOADS ALL VARIABLES FOR HARD PROCESS
C     modifications for Les Houches accord by PR (7/15/02)
C
C M.Boonekamp, Aug-2003: - include a routine to cut on phase space
C                        - include a switch for gamgam to Higgs
C                        - include a switch for Jz=0 photon pairprod.
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
C ... begin R.S.
      include 'CHIDe.inc'
      INCLUDE 'fpmc.inc' ! To be able to check NFLUX
C ... end R.S.
      DOUBLE PRECISION CIRCKP(2)
      COMMON /HWCIR2/CIRCKP
      DOUBLE PRECISION Z1,Z2,C1,C2,B1,B2,CIRCEE,CIRCGG,RS,MISS,ETA,
     $     HWUGAU,HWECIR,QMX1,QMN1,QMX2,QMN2,TEST
      INTEGER IHAD
      SAVE MISS
      DOUBLE PRECISION HWRGEN
      EXTERNAL HWRGEN,HWECIR
C--Les Houches Common Block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)

C ... begin R.S.
C ... end R.S.
      
      IF (IERROR.NE.0)  RETURN
C--pick the type of event to generate if using Les Houches accord
C--first choice according to maxiumum weight
      IF(IPROC.LT.0) THEN
        IF(ABS(IDWTUP).EQ.1) THEN
          IF(ITYPLH.EQ.0) THEN
            TEST = HWRGEN(1)*LHMXSM
            DO ITYPLH=1,NPRUP
              IF(TEST.LE.ABS(LHXMAX(ITYPLH))) GOTO 5
              TEST = TEST-ABS(LHXMAX(ITYPLH))
            ENDDO
 5          WGTMAX = ABS(LHXMAX(ITYPLH))
            WBIGST = ABS(LHXMAX(ITYPLH))
          ENDIF
C--second choice according to cross section
        ELSEIF(ABS(IDWTUP).EQ.2) THEN
          IF(ITYPLH.EQ.0) THEN
            TEST = HWRGEN(1)*LHMXSM
            DO ITYPLH=1,NPRUP
              IF(TEST.LE.ABS(LHXSCT(ITYPLH))) GOTO 6
              TEST = TEST-ABS(LHXSCT(ITYPLH))
            ENDDO
 6          WGTMAX = ABS(LHXMAX(ITYPLH))
            WBIGST = ABS(LHXMAX(ITYPLH))
          ENDIF
        ELSE
          WGTMAX = 1.0D0
          WBIGST = 1.0D0
          ITYPLH = 1
        ENDIF
      ENDIF
C---ROUTINE LOOPS BACK TO HERE IF GENERATED WEIGHT WAS NOT ACCEPTED
   10 GENEV=.FALSE.
C---FSTWGT IS .TRUE. DURING FIRST CALL TO HARD PROCESS ROUTINE
      FSTWGT=NWGTS.EQ.0
C---FSTEVT IS .TRUE. THROUGHOUT THE FIRST EVENT
      FSTEVT=NEVHEP.EQ.1
C---SET COLOUR CORRECTION TO FALSE
      COLUPD = .FALSE.
      HRDCOL(1,1)=0
      HRDCOL(1,3)=0
C---SET UP INITIAL STATE
      NHEP=1
      ISTHEP(NHEP)=101
      PHEP(1,NHEP)=0.
      PHEP(2,NHEP)=0.
      PHEP(3,NHEP)=PBEAM1
      PHEP(4,NHEP)=EBEAM1
      PHEP(5,NHEP)=RMASS(IPART1)
      JMOHEP(1,NHEP)=0
      JMOHEP(2,NHEP)=0
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
      IDHW(NHEP)=IPART1
      IDHEP(NHEP)=IDPDG(IPART1)
      NHEP=NHEP+1
      ISTHEP(NHEP)=102
      PHEP(1,NHEP)=0.
      PHEP(2,NHEP)=0.
      PHEP(3,NHEP)=-PBEAM2
      PHEP(4,NHEP)=EBEAM2
      PHEP(5,NHEP)=RMASS(IPART2)
      JMOHEP(1,NHEP)=0
      JMOHEP(2,NHEP)=0
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
      IDHW(NHEP)=IPART2
      IDHEP(NHEP)=IDPDG(IPART2)
C---NEXT ENTRY IS OVERALL CM FRAME
      NHEP=NHEP+1
      IDHW(NHEP)=14
      IDHEP(NHEP)=0
      ISTHEP(NHEP)=103
      JMOHEP(1,NHEP)=NHEP-2
      JMOHEP(2,NHEP)=NHEP-1
      JDAHEP(1,NHEP)=0
      JDAHEP(2,NHEP)=0
      CALL HWVSUM(4,PHEP(1,NHEP-1),PHEP(1,NHEP-2),PHEP(1,NHEP))
      CALL HWUMAS(PHEP(1,NHEP))
C Select a primary interaction point
      IF (PIPSMR) THEN
        CALL HWRPIP
      ELSE
        CALL HWVZRO(4,VTXPIP)
      ENDIF
      CALL HWVEQU(3,VTXPIP,VHEP(1,NHEP))
      VHEP(4,NHEP)=0.0
C---GENERATE PHOTONS (WEIZSACKER-WILLIAMS APPROX)
C   FOR HADRONIC PROCESSES WITH LEPTON BEAMS
      GAMWT=ONE
C ... M.B. : added IPRO=96-99 below
      IF ((IPRO.GT.12.AND.IPRO.LT.90).OR.
     &    (IPRO.GE.96.AND.IPRO.LE.99)) THEN
        IF (CIRCOP.EQ.0) THEN
           IF (ABS(IDHEP(1)).EQ.11.OR.ABS(IDHEP(1)).EQ.13) then
               CALL HWEGAM(1,ZERO, ONE,.FALSE.)
c     &          CALL HWEGAM(1,ZERO, ONE,.FALSE.)
             endif     
           IF (ABS(IDHEP(2)).EQ.11.OR.ABS(IDHEP(2)).EQ.13) then
               CALL HWEGAM(2,ZERO, ONE,.FALSE.)
             endif     
        ELSE
C---MODIFIED TO USE CIRCE FOR BEAMSTRAHLUNG EFFECTS
          IF (ABS(IDHEP(1)).NE.11.OR.IDHEP(1)+IDHEP(2).NE.0) STOP
     $         'This version only works for e+e- annihilation'
          IF (FSTWGT) THEN
            RS=NINT(PHEP(5,3)*10)/1D1
            CALL CIRCES(ZERO,ZERO,RS,CIRCAC,CIRCVR,CIRCRV,CIRCCH)
          ENDIF
          CALL HWEGAM(1,ZERO, ONE,.TRUE.)
          CALL HWEGAM(2,ZERO, ONE,.TRUE.)
          Z1=PHEP(4,4)/PHEP(4,1)
          Z2=PHEP(4,6)/PHEP(4,2)
C---FACTORIZE THE DISTRIBUTIONS FROM CIRCE
          C1=CIRCGG(Z1,-1D0)/SQRT(CIRCGG(-1D0,-1D0))
          C2=CIRCGG(-1D0,Z2)/SQRT(CIRCGG(-1D0,-1D0))
C---REMOVE SPURIOUS WEIGHT GIVEN IN HWEGAM
          GAMWT=GAMWT/(.5*ALPHEM/PIFAC*(1+(1-Z1)**2)/Z1*
     $         LOG((ONE-Z1)/Z1*4*PHEP(4,1)*PHEP(4,2)/PHEP(5,1)**2))
     $               /(.5*ALPHEM/PIFAC*(1+(1-Z2)**2)/Z2*
     $         LOG((ONE-Z2)/Z2*4*PHEP(4,4)*PHEP(4,2)/PHEP(5,1)**2))
C---REPLACE IT BY THE SUM OF BEAM AND BREM STRAHLUNG
          QMX1=MIN(Q2WWMX,(Z1*PHEP(3,1))**2)
          QMN1=MAX(Q2WWMN,(PHEP(5,1)*Z1)**2/(1-Z1))
          QMX2=MIN(Q2WWMX,(Z2*PHEP(3,2))**2)
          QMN2=MAX(Q2WWMN,(PHEP(5,2)*Z2)**2/(1-Z2))
          B1=.5*ALPHEM/PIFAC*(1+(1-Z1)**2)/Z1*LOG(QMX1/QMN1)
          B2=.5*ALPHEM/PIFAC*(1+(1-Z2)**2)/Z2*LOG(QMX2/QMN2)
          IF (CIRCOP.EQ.1) THEN
            GAMWT=GAMWT*B1*B2
          ELSEIF (CIRCOP.EQ.2) THEN
            GAMWT=GAMWT*C1*C2
          ELSEIF (CIRCOP.EQ.3) THEN
            GAMWT=GAMWT*(C1+B1)*(C2+B2)
          ELSE
            STOP 'Illegal value of circop!'
          ENDIF
        ENDIF
      ELSEIF (IPRO.GE.90) THEN
        IF (CIRCOP.NE.0) STOP 'Circe not interfaced for DIS processes'
        IF (ABS(IDHEP(2)).EQ.11.OR.ABS(IDHEP(2)).EQ.13)
     &       CALL HWEGAM(2,ZERO, ONE,.FALSE.)
      ENDIF
C---GENERATE ISR PHOTONS FOR LEPTONIC PROCESSES
      IF (IPRO.GT.0.AND.IPRO.LE.12) THEN
        IF (CIRCOP.EQ.0) THEN
           CALL HWEISR(1)
           CALL HWEISR(2)
        ELSE
C---MODIFIED TO USE CIRCE FOR BEAMSTRAHLUNG EFFECTS
          IF (ABS(IDHEP(1)).NE.11.OR.IDHEP(1)+IDHEP(2).NE.0) STOP
     $         'This version only works for e+e- annihilation'
          IF (FSTWGT) THEN
            RS=NINT(PHEP(5,3)*10)/1D1
            CALL CIRCES(ZERO,ZERO,RS,CIRCAC,CIRCVR,CIRCRV,CIRCCH)
C---PRECALCULATE THE PART OF THE SPECTRUM MISSED BETWEEN ZMXISR AND 1
            ETA=0.6D0
            MISS=HWUGAU(HWECIR,1D-15**(1-ETA),(1-ZMXISR)**(1-ETA),1D-12)
          ENDIF
          COLISR=.TRUE.
          CALL HWEISR(1)
          CALL HWEISR(2)
          IHAD=1
          IF (JDAHEP(1,IHAD).NE.0) IHAD=JDAHEP(1,IHAD)
          Z1=PHEP(4,IHAD)/PHEP(4,1)
          IHAD=2
          IF (JDAHEP(1,IHAD).NE.0) IHAD=JDAHEP(1,IHAD)
          Z2=PHEP(4,IHAD)/PHEP(4,2)
C---FACTORIZE THE DISTRIBUTIONS FROM CIRCE
          C1=CIRCEE(Z1,-1D0)/SQRT(CIRCEE(-1D0,-1D0))
          C2=CIRCEE(-1D0,Z2)/SQRT(CIRCEE(-1D0,-1D0))
          IF (Z1.EQ.ONE) C1=C1+MISS
          IF (Z2.EQ.ONE) C2=C2+MISS
C---REMOVE WEIGHT GIVEN IN HWEISR
          B1=CIRCKP(1)
          B2=CIRCKP(2)
          GAMWT=GAMWT/(B1*B2)
C---REPLACE IT BY THE SUM OF BEAM AND BREM STRAHLUNG
          IF (CIRCOP.EQ.1) THEN
            GAMWT=GAMWT*B1*B2
          ELSEIF (CIRCOP.EQ.2) THEN
            GAMWT=GAMWT*C1*C2
          ELSEIF (CIRCOP.EQ.3) THEN
C---IN THE APPROXIMATION OF DOMINANCE BY THE DELTA-FUNCTION TERM
            IF (Z1.EQ.ONE) C1=C1-1
            IF (Z2.EQ.ONE) C2=C2-1
C---IF IT DOES NOT DOMINATE, ZMXISR SHOULD BE DECREASED
            IF (B1+C1.LT.ZERO) CALL HWWARN('HWEPRO',501,*999)
            IF (B2+C2.LT.ZERO) CALL HWWARN('HWEPRO',502,*999)
            GAMWT=GAMWT*(C1+B1)*(C2+B2)
          ELSE
            STOP 'Illegal value of circop!'
          ENDIF
        ENDIF
      ENDIF
      
C---IF USER LIMITS WERE TOO TIGHT, MIGHT NOT BE ANY PHASE-SPACE
      IF (GAMWT.LE.ZERO) GOTO 30
C---IF CMF HAS ACQUIRED A TRANSVERSE BOOST, OR USER REQUESTS IT ANYWAY,
C   BOOST EVENT RECORD BACK TO CMF
      IF (PHEP(1,3)**2+PHEP(2,3)**2.GT.ZERO .OR. USECMF) CALL HWUBST(1)
C---ROUTINE LOOPS BACK TO HERE IF GENERATED WEIGHT WAS ACCEPTED
   20 CONTINUE

      IPRO=MOD(IPROC/100,100)

C---PROCESS GENERATED BY LES HOUCHES INTERFACE
      IF(IPRO.LE.0) THEN
        CALL HWHGUP
      ELSEIF (IPRO.EQ.1) THEN
        IF (IPROC.LT.110.OR.IPROC.GE.120) THEN
C--- E+E- -> Q-QBAR OR L-LBAR
          CALL HWHEPA
        ELSE
C--- E+E- -> Q-QBAR-GLUON
          CALL HWHEPG
        ENDIF
      ELSEIF (IPRO.EQ.2) THEN
C--- E+E- -> W+ W-
        CALL HWHEWW
      ELSEIF (IPRO.EQ.3) THEN
C---E+E- -> Z H
        CALL HWHIGZ
      ELSEIF (IPRO.EQ.4) THEN
C---E+E- -> NUEB NUE H
        CALL HWHIGW
      ELSEIF (IPRO.EQ.5 .AND. IPROC.LT.550) THEN
C---EE -> EE GAMGAM -> EE FFBAR/WW
        CALL HWHEGG
      ELSEIF (IPRO.EQ.5) THEN
C---EE -> ENU GAMW -> ENU FF'BAR/WZ
        CALL HWHEGW
      ELSEIF (IPRO.EQ.6) THEN
C---EE -> FOUR JETS
        CALL HWH4JT
      ELSEIF(IPRO.EQ.7) THEN
C--EE -> SUSY PARTICLES(PAIR PRODUCTION)
        CALL HWHESP
      ELSEIF(IPRO.EQ.8) THEN
C--EE -> RPV SUSY PARTICLE PRODUCTION
        CALL HWHREP
      ELSEIF (IPRO.EQ.9) THEN
        IF((MOD(IPROC,10000).EQ.955).OR.
     &     (MOD(IPROC,10000).EQ.965).OR.
     &     (MOD(IPROC,10000).EQ.975))THEN
C---MSSM Higgs pair production in l+l-: H+ H- and A0 Higgs, Higgs=h0,H0.
          CALL HWHIHH
        ELSEIF((MOD(IPROC,10000).EQ.910).OR.
     &         (MOD(IPROC,10000).EQ.920))THEN
C---MSSM scalar Higgs production from vector-vector fusion.
          CALL HWHIGW
        ELSEIF((MOD(IPROC,10000).EQ.960).OR.
     &         (MOD(IPROC,10000).EQ.970))THEN
C---MSSM scalar Higgs production from Higgs-strahlung.
          CALL HWHIGZ
        END IF
      ELSEIF ((IPRO.EQ.10).OR.(IPRO.EQ.11)) THEN
C---SM/MSSM Higgs production with heavy quark flavours via e+e-.
        CALL HWHIGE
      ELSEIF (IPRO.EQ.13) THEN
C---GAMMA/Z0/Z' DRELL-YAN PROCESS
        CALL HWHDYP
      ELSEIF (IPRO.EQ.14) THEN
C---W+/- PRODUCTION VIA DRELL-YAN PROCESS
        CALL HWHWPR
      ELSEIF (IPRO.EQ.15) THEN
C---QCD HARD 2->2 PROCESSES
        CALL HWHQCD
      ELSEIF ((IPRO.EQ.16).OR.(IPRO.EQ.36)) THEN
C---SM/MSSM HIGGS PRODUCTION VIA QUARK/GLUON FUSION
        CALL HWHIGS
      ELSEIF (IPRO.EQ.17) THEN
C---QCD HEAVY FLAVOUR PRODUCTION
        CALL HWHHVY
      ELSEIF (IPRO.EQ.18) THEN
C---QCD DIRECT PHOTON + JET PRODUCTION
        CALL HWHPHO
      ELSEIF ((IPRO.EQ.19).OR.(IPRO.EQ.37)) THEN
C---SM/MSSM HIGGS PRODUCTION VIA W/Z FUSION
        CALL HWHIGW
      ELSEIF (IPRO.EQ.20) THEN
C---TOP PRODUCTION FROM W EXCHANGE
        CALL HWHWEX
      ELSEIF (IPRO.EQ.21) THEN
C---VECTOR BOSON + JET PRODUCTION
        CALL HWHV1J
      ELSEIF (IPRO.EQ.22) THEN
C QCD direct photon pair production
C M.B. : call corrected routine
         CALL HWH2PH
c        CALL HWHPH2
      ELSEIF (IPRO.EQ.23) THEN
C QCD Higgs plus jet production
        CALL HWHIGJ
      ELSEIF (IPRO.EQ.24) THEN
C---COLOUR-SINGLET EXCHANGE
        CALL HWHSNG
      ELSEIF (IPRO.EQ.25) THEN
C---SM Higgs production with heavy quark flavours via qq and gg.
        CALL HWHIGQ
      ELSEIF ((IPRO.EQ.26).OR.(IPRO.EQ.27)) THEN
C---SM Higgs production with heavy gauge bosons via qq(').
        CALL HWHIGV
C---Gauge boson pair in hadron hadron
      ELSEIF (IPRO.EQ.28) THEN
C M.B skip following test since HWHGBP is not there anyway
C        IF (IPROC.LT.2850) THEN
C          CALL HWHGBP
C        ELSE
C          CALL HWHVVJ
C        ENDIF
        CALL HWHGBP
C--Vector boson + two jets
      ELSEIF(IPRO.EQ.29) THEN
        CALL HWHV2J
      ELSEIF (IPRO.EQ.30) THEN
C---HADRON-HADRON SUSY PROCESSES
        CALL HWHSSP
      ELSEIF ((IPRO.EQ.31).OR.(IPRO.EQ.32)) THEN
C---MSSM charged/neutral Higgs production in association with squarks.
        CALL HWHISQ
      ELSEIF (IPRO.EQ.33) THEN
        IF(MOD(IPROC,10000).EQ.3350)THEN
C---MSSM charged Higgs production in association with W: W+H- + W-H+.
          CALL HWHIBK
        ELSEIF((MOD(IPROC,10000).EQ.3310).OR.
     &         (MOD(IPROC,10000).EQ.3320).OR.
     &         (MOD(IPROC,10000).EQ.3360).OR.
     &         (MOD(IPROC,10000).EQ.3370))THEN
C---MSSM Higgs production with heavy gauge bosons via qq(').
          CALL HWHIGV
        ELSE
C---MSSM charged/neutral Higgs pair production.
          CALL HWHIGH
        END IF
      ELSEIF (IPRO.EQ.34) THEN
C---MSSM charged/neutral Higgs production via bg fusion.
        CALL HWHIBG
      ELSEIF (IPRO.EQ.35) THEN
C---MSSM charged Higgs production via bq fusion.
        CALL HWHIBQ
      ELSEIF (IPRO.EQ.38) THEN
C---MSSM charged/neutral Higgs production with heavy quarks via qq and gg.
        CALL HWHIGQ
      ELSEIF(IPRO.EQ.40.OR.IPRO.EQ.41) THEN
C---HADRON-HADRON R-PARITY VIOLATING SUSY PROCESSES
        CALL HWHRSP
      ELSEIF (IPRO.EQ.42) THEN
C---SPIN-TWO RESONANCE
        CALL HWHGRV
      ELSEIF (IPRO.EQ.50) THEN
C Point-like photon two-jet production
        CALL HWHPPT
      ELSEIF (IPRO.EQ.51) THEN
C Point-like photon/QCD heavy flavour pair production
        CALL HWHPPH
      ELSEIF (IPRO.EQ.52) THEN
C Point-like photon/QCD heavy flavour single excitation
        CALL HWHPPE
      ELSEIF (IPRO.EQ.53) THEN
C Compton scattering of point-like photon and (anti)quark
        CALL HWHPQS
      ELSEIF (IPRO.EQ.55) THEN
C Point-like photon/higher twist meson production
        CALL HWHPPM
      ELSEIF (IPRO.EQ.60) THEN
C---QPM GAMMA-GAMMA-->QQBAR
        CALL HWHQPM     
      ELSEIF (IPRO.GE.70.AND.IPRO.LE.79) THEN
C---BARYON-NUMBER VIOLATION, AND OTHER MULTI-W PRODUCTION PROCESSES
        CALL HVHBVI
      ELSEIF (IPRO.EQ.80) THEN
C---MINIMUM-BIAS: NO HARD SUBPROCESS
C   FIND WEIGHT
        CALL HWMWGT
      ELSEIF (IPRO.EQ.90) THEN
C---DEEP INELASTIC
        CALL HWHDIS
      ELSEIF(IPRO.EQ.91) THEN
C---BOSON - GLUON(QUARK) FUSION -->  ANTIQUARK(GLUON) + QUARK
        CALL HWHBGF
      ELSEIF(IPRO.EQ.92) THEN
C---DEEP INELASTIC WITH EXTRA JET: OBSOLETE PROCESS
        WRITE (6,40)
 40     FORMAT (1X,' IPROC=92** is no longer supported.'
     &         /1X,' Please use IPROC=91** instead.')
        CALL HWWARN('HWEPRO',500,*999)
      ELSEIF(IPRO.EQ.95) THEN
C---HIGGS PRODUCTION VIA W FUSION IN E P
        CALL HWHIGW
C ... New processes for DPE : 
C ... M.B.T.K.
      ELSEIF(IPRO.EQ.96) THEN
C---INCLUSIVE CHI PRODUCTION
        CALL HWINCH
      ELSEIF(IPRO.EQ.97) THEN
C---EXCLUSIVE CHI PRODUCTION
        CALL HWEXCH
      ELSEIF(IPRO.EQ.98) THEN
C---EXCLUSIVE PHOTON PAIR PRODUCTION
        CALL HWHQPP
      ELSEIF(IPRO.EQ.99) THEN
C---EXCLUSIVE HIGGS BOSON PRODUCTION
        CALL HWHIGP
C ... end M.B.T.K.
      ELSE
C---UNKNOWN PROCESS
        CALL HWWARN('HWEPRO',102,*999)
      ENDIF
 30   IF (GENEV) THEN
         IF (NOWGT) THEN
            IF (NEGWTS) THEN
               IF (EVWGT.LT.ZERO) THEN
                  EVWGT=-AVABW
               ELSE
                  EVWGT= AVABW
               ENDIF
            ELSE
               EVWGT=AVWGT
            ENDIF
         ENDIF
         ISTAT=10
C--New call spin correlation code if needed
         IF(SYSPIN.AND.(IPRO.EQ. 1.OR.IPRO.EQ.13.OR.IPRO.EQ.14.OR.
     &                  IPRO.EQ.15.OR.IPRO.EQ.17.OR.IPRO.EQ.20.OR.
     &                  IPRO.EQ. 7.OR.IPRO.EQ.30.OR.IPRO.EQ.40.OR.
     &                  IPRO.EQ.41.OR.IPRO.EQ.8)) CALL HWHSPN
         RETURN
      ELSE
C---IF AN EVENT IS CANCELLED BEFORE IT IS GENERATED, GIVE IT ZERO WEIGHT
C ... M.B.T.K. : insert phase space cut
*        CALL HWPCUT  !!! At least temporarily disabled
C end M.B.T.K.
        IF (IERROR.NE.0) THEN
          EVWGT=ZERO
          IERROR=0
        ENDIF
        EVWGT=EVWGT*GAMWT
        if(EVWGT-EVWGT.NE.ZERO) EVWGT=ZERO ! R.S. IF =NAN or =INF
        NWGTS=NWGTS+1
        ABWGT=ABS(EVWGT)
        IF (EVWGT.LT.ZERO) THEN
           IF (NEGWTS) THEN
              NNEGWT=NNEGWT+1
           ELSE
              IF (EVWGT.LT.-1.D-9) CALL HWWARN('HWEPRO',3,*999)
              EVWGT=ZERO
              ABWGT=ZERO
           ENDIF
        ENDIF
        WGTSUM=WGTSUM+EVWGT
        WSQSUM=WSQSUM+EVWGT**2
        ABWSUM=ABWSUM+ABWGT
C--weight addition for Les Houches accord
        IF(IPROC.LE.0) THEN
          IF(ABS(IDWTUP).EQ.1) THEN
             LHWGT (ITYPLH) = LHWGT (ITYPLH)+EVWGT
             LHWGTS(ITYPLH) = LHWGTS(ITYPLH)+EVWGT**2
             LHIWGT(ITYPLH) = LHIWGT(ITYPLH)+1
          ENDIF
        ENDIF
        IF (ABWGT.GT.WBIGST) THEN
           WBIGST=ABWGT
           IF (NOWGT.AND.WBIGST.GT.WGTMAX) THEN
              IF (NEVHEP.NE.0) CALL HWWARN('HWEPRO',1,*999)
              WGTMAX=WBIGST*1.1
              WRITE (6,99) WGTMAX
C--additional for Les Houche accord
              IF(IPROC.LE.0) THEN
                IF(ABS(IDWTUP).EQ.1) 
     &                LHMXSM = LHMXSM-LHXMAX(ITYPLH)+ABWGT
                LHXMAX(ITYPLH) = EVWGT
              ENDIF
           ENDIF
        ENDIF
        IF (NEVHEP.NE.0) THEN
C---LOW EFFICIENCY WARNINGS:
C   WARN AT 10*EFFMIN, STOP AT EFFMIN
          IF (10*EFFMIN*NWGTS.GT.NEVHEP) THEN
            IF (EFFMIN*NWGTS.GT.NEVHEP) CALL HWWARN('HWEPRO',200,*999)
            IF (EFFMIN.GT.ZERO) THEN
              IF (MOD(NWGTS,INT(10/EFFMIN)).EQ.0) THEN
                CALL HWWARN('HWEPRO',2,*999)
                WRITE (6,98) WGTMAX
              ENDIF
            ENDIF
          ENDIF
          IF (NOWGT) THEN
            GENEV=ABWGT.GT.WGTMAX*HWRGEN(0)
          ELSE
            GENEV=ABWGT.NE.ZERO
          ENDIF
          IF (GENEV)  GOTO 20
c          PRINT*, 'begin7'
          GOTO 10
c        PRINT*, 'begin8'
        ENDIF
      ENDIF
 98   FORMAT(10X,'    MAXIMUM WEIGHT =',1PG24.16)
 99   FORMAT(10X,'NEW MAXIMUM WEIGHT =',1PG24.16)
  999 END
C-----------------------------------------------------------------------
      SUBROUTINE HWEXCH
C-----------------------------------------------------------------------
C     CHI_C PRODUCTION VIA GLUON FUSION
C     MEAN EVWGT = CHI_C PRODN C-S
C     DERIVED FROM HWHIGS.
C
C     M.Boonekamp, Jun 2004
C
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION HWUALF,HWRGEN,EMH,CSFAC,EMFAC
      DOUBLE PRECISION TAUT,TAUW,EMW,TAUTR,TAUTI,TAUWR,TAUWI,SUMR,
     &     SUMI,HWIDTH,EPSM,P2OLD,P2NEW
      INTEGER ICHI,I,J,ID1,ID2,INDEX
      EXTERNAL HWUALF,HWRGEN
      SAVE CSFAC
      EQUIVALENCE (EMW,RMASS(198))
C ... Decide on which chi to produce
      ICHI=MOD(IPROC,100)
      IF(ICHI.EQ.1) THEN
         INDEX=165
         HWIDTH=0.0106 !! Measured width : 10.6 MeV (PDG 2004)
      ELSEIF(ICHI.EQ.2) THEN
         INDEX=306
         HWIDTH=0.0016 !! Predicted width, central value
                       !! Eichten, Quigg, PR D52 1995
      ELSE
         PRINT*, 'In HWEXCH : WRONG IPROC ', IPROC
         STOP
      ENDIF
      IF (GENEV) THEN
C ... fill as gamma+gamma -> Chi, whatever TYPINT is
C ... An additional factor 16 = 2(Nc^2-1) is included in the cross-section
C ... (gluon color&spin coherence)
        IDN(1)=59
        IDN(2)=59
        IDCMF=INDEX
        XX(1)=1.
        XX(2)=1.
        CALL HWEONE
      ELSE
        EVWGT=0.
        EMH=RMASS(INDEX)
        EMFAC=1.D0
        PHEP(5,3)=EMH
        EMSCA=EMH
        IF (EMSCA.NE.EMLST) THEN
          EMLST=EMH
          IF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC') THEN
             CSFAC=GEV2NB*2.D0*(PI**2)*HWIDTH/(RMASS(INDEX)**3)*16
c---End modif by Tibor Kucs
          ELSE
            PRINT*, 'In HWEXCH : wrong TYP* settings, ', TYPINT, TYPEPR, 
     &              ' - STOP'
            STOP
          ENDIF
         ENDIF
         EVWGT=EVWGT + CSFAC
        ENDIF
  999 END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE HWFXER(INIFX)
C
C     M.Boonekamp, Oct 2003 : 
C
C     Edits the event record after event finalization. The following
C     replacements are done, if needed:
C     - primary beams : E- -> P   E+ -> PBAR
C     - secondary beams : GAMMA -> POMERON/REGGEON
C     - hard initial parton : GAMMA -> GLUON
C     Masses are set correctly, and momenta are resaled.
C
C     A LOT IS HARDCODED AND A POTENTIAL SOURCE OF PROBLEMS
C
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION PMOM,PRESC
      INTEGER I,J,ISEC,IPOM,IREG
      LOGICAL INIFX
      integer ind
      DATA IPOM/67/,IREG/68/

c CHR/YURA common      
      COMMON /YURA/ IND
      
C Initialise :
      IF(INIFX) THEN
        RNAME(IPOM)='POMERON '
        RNAME(IREG)='REGGEON '
        IDPDG(IPOM)=990
        IDPDG(IREG)=110
        IFLAV(IPOM)=0
        IFLAV(IREG)=0
        ICHRG(IPOM)=0
        ICHRG(IREG)=0
        RMASS(IPOM)=0d0
        RMASS(IREG)=0d0
        RLTIM(IPOM)=1d30
        RLTIM(IREG)=1d30
        RSPIN(IPOM)=0d0
        RSPIN(IREG)=0d0
        RETURN
C DIS case : don't touch Beam1 (always e+/e-), modify only Beam2 (p,pbar)
      ELSEIF((PART1.EQ.'E+'.OR.PART1.EQ.'E-').AND.
     &       (PART2.EQ.'E+'.OR.PART2.EQ.'E-').AND.
     &       NFLUX.GE.9.AND.NFLUX.LE.11.AND.IPRO.EQ.90) THEN
        I=2
        J=5
        IDHW(I)=73
        IDHEP(I)=IDPDG(73)
        PHEP(5,I)=RMASS(73)
        IDHW(J)=73
        IDHEP(J)=IDPDG(73)
        PHEP(5,J)=RMASS(73)
        IF (JDAHEP(1,I).NE.0) ISEC=JDAHEP(1,I)
        IF(ISEC.GT.0) THEN
          IF(NFLUX.EQ.10) THEN
            IDHW(ISEC)=IREG
            IDHEP(ISEC)=IDPDG(IREG)
C CHR/Yura Pom Reg
          ELSEIF(NFLUX.EQ.19) THEN
	    IF(I.EQ.1) THEN
            IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
            ELSEIF(I.EQ.2) THEN	    
            IDHW(ISEC)=IREG
            IDHEP(ISEC)=IDPDG(IREG)
	    ENDIF
C Yura/CHR Reg Pom
          ELSEIF(NFLUX.EQ.21) THEN
	    IF(I.EQ.2) THEN
            IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
            ELSEIF(I.EQ.1) THEN	    
            IDHW(ISEC)=IREG
            IDHEP(ISEC)=IDPDG(IREG)
	    ENDIF
          ELSEIF(NFLUX.LE.11.OR.NFLUX.EQ.16) THEN
            IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
C Yura/CHR Photon Pom
          ELSEIF(NFLUX.EQ.20) THEN
            IF(I.EQ.2) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)	    
            ENDIF
C Yura/CHR Pom Photon
          ELSEIF(NFLUX.EQ.22) THEN
            IF(I.EQ.1) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)	    
            ENDIF
          ENDIF
        ENDIF
C Single Diffraction case
C CHR Rafal Should it be 11 or 16 ???????
      ELSEIF((PART1.EQ.'E+'.OR.PART1.EQ.'E-').AND.
     &       (PART2.EQ.'P'.OR.PART2.EQ.'PBAR').AND.
     &       NFLUX.GE.9.AND.NFLUX.LE.11) THEN
         I=1
         J=5
         IDHW(I)=73
         IDHEP(I)=IDPDG(73)
         PHEP(5,I)=RMASS(73)
         IDHW(J)=73
         IDHEP(J)=IDPDG(73)
         PHEP(5,J)=RMASS(73)
         IF (JDAHEP(1,I).NE.0) ISEC=JDAHEP(1,I)
         IF(ISEC.GT.0) THEN
            IF(NFLUX.EQ.10) THEN
               IDHW(ISEC)=IREG
               IDHEP(ISEC)=IDPDG(IREG)
C CHR/Yura Pom Reg
            ELSEIF(NFLUX.EQ.19) THEN
	     IF(I.EQ.1) THEN
              IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
              ELSEIF(I.EQ.2) THEN	    
               IDHW(ISEC)=IREG
               IDHEP(ISEC)=IDPDG(IREG)
	     ENDIF
C CHR/Yura Reg Pom (non sense here)
            ELSEIF(NFLUX.EQ.21) THEN
	     IF(I.EQ.2) THEN
              IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
              ELSEIF(I.EQ.1) THEN	    
               IDHW(ISEC)=IREG
               IDHEP(ISEC)=IDPDG(IREG)
	     ENDIF
C CHR/Yura Photon Pom (non sense here)
            ELSEIF(NFLUX.EQ.20) THEN
             IF(I.EQ.2) THEN
	      IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
	     ENDIF
C CHR/Yura Pom Photon (non sense here)
            ELSEIF(NFLUX.EQ.22) THEN
             IF(I.EQ.1) THEN
	      IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
	     ENDIF
            ELSEIF(NFLUX.LE.11.OR.NFLUX.EQ.16) THEN
               IDHW(ISEC)=IPOM
               IDHEP(ISEC)=IDPDG(IPOM)
            ENDIF
         ENDIF


C CHR/Rafal is it 16 or 11??????
      ELSEIF((PART1.EQ.'P'.OR.PART1.EQ.'PBAR').AND.
     &       (PART2.EQ.'E+'.OR.PART2.EQ.'E-').AND.
     &       NFLUX.GE.9.AND.NFLUX.LE.11) THEN
         I=2
         J=5
         IDHW(I)=73
         IDHEP(I)=IDPDG(73)
         PHEP(5,I)=RMASS(73)
         IDHW(J)=73
         IDHEP(J)=IDPDG(73)
         PHEP(5,J)=RMASS(73)
         IF (JDAHEP(1,I).NE.0) ISEC=JDAHEP(1,I)
         IF(ISEC.GT.0) THEN
            IF(NFLUX.EQ.10) THEN
               IDHW(ISEC)=IREG
               IDHEP(ISEC)=IDPDG(IREG)
C CHR/Yura Pom Reg
            ELSEIF(NFLUX.EQ.19) THEN
	     IF(I.EQ.1) THEN
              IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
             ELSEIF(I.EQ.2) THEN	    
              IDHW(ISEC)=IREG
              IDHEP(ISEC)=IDPDG(IREG)
	     ENDIF
C CHR/Yura Reg Pom (non sense here)
            ELSEIF(NFLUX.EQ.21) THEN
	     IF(I.EQ.2) THEN
              IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
             ELSEIF(I.EQ.1) THEN	    
              IDHW(ISEC)=IREG
              IDHEP(ISEC)=IDPDG(IREG)
	     ENDIF
C CHR/Yura Pom Photon (non sense here)
            ELSEIF(NFLUX.EQ.20) THEN
             IF(I.EQ.2) THEN
	      IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
	     ENDIF
C CHR/Yura Photon Pom (non sense here)
            ELSEIF(NFLUX.EQ.22) THEN
             IF(I.EQ.1) THEN
	      IDHW(ISEC)=IPOM
              IDHEP(ISEC)=IDPDG(IPOM)
	     ENDIF
            ELSEIF(NFLUX.LE.11.OR.NFLUX.EQ.16) THEN
               IDHW(ISEC)=IPOM
               IDHEP(ISEC)=IDPDG(IPOM)
            ENDIF
         ENDIF

c CHR YURA modifs
C DPE case
      ELSEIF((PART1.EQ.'E+'.OR.PART1.EQ.'E-').AND.
     &       (PART2.EQ.'E+'.OR.PART2.EQ.'E-').AND.
     &       ((NFLUX.GE.9.AND.NFLUX.LE.16).OR.NFLUX.EQ.18.OR.
     &      NFLUX.GE.19)) THEN
C ... Modify primary and secondary beams :
        DO I=1,2
           J=3+2*I
           IDHW(I)=73
           IDHEP(I)=IDPDG(73)
           PHEP(5,I)=RMASS(73)
           IDHW(J)=73
           IDHEP(J)=IDPDG(73)
           PHEP(5,J)=RMASS(73)
        ENDDO
C ... Modify secondary beams :
        DO I=1,2
           IF (JDAHEP(1,I).NE.0) ISEC=JDAHEP(1,I)
           IF(ISEC.GT.0) THEN
              IF(NFLUX.EQ.10) THEN
                 IDHW(ISEC)=IREG
                 IDHEP(ISEC)=IDPDG(IREG)
C CHR/Yura Pom Reg
          ELSEIF(NFLUX.EQ.19) THEN
	    IF(I.EQ.1) THEN
            IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
            ELSEIF(I.EQ.2) THEN	    
            IDHW(ISEC)=IREG
            IDHEP(ISEC)=IDPDG(IREG)
	    ENDIF
C CHR/Yura  Reg Pom - finally, it makes sense
          ELSEIF(NFLUX.EQ.21) THEN
	    IF(I.EQ.2) THEN
            IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
            ELSEIF(I.EQ.1) THEN	    
            IDHW(ISEC)=IREG
            IDHEP(ISEC)=IDPDG(IREG)
	    ENDIF
C CHR/Yura  Photon Pom - finally, it makes sense
          ELSEIF(NFLUX.EQ.20) THEN
            IF(I.EQ.2) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
	    ENDIF
C CHR/Yura  Pom Photon - finally, it makes sense
          ELSEIF(NFLUX.EQ.22) THEN
            IF(I.EQ.1) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
	    ENDIF
C CHR 2014 Photon Pom heavy ion - finally, it makes sense
          ELSEIF(NFLUX.EQ.26) THEN
            IF(I.EQ.2) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
	    ENDIF
C CHR 2014 Pom Photon heavy ion - finally, it makes sense
          ELSEIF(NFLUX.EQ.25) THEN
            IF(I.EQ.1) THEN
	    IDHW(ISEC)=IPOM
            IDHEP(ISEC)=IDPDG(IPOM)
	    ENDIF
              ELSEIF(NFLUX.LE.11.OR.NFLUX.EQ.16) THEN
                 IDHW(ISEC)=IPOM
                 IDHEP(ISEC)=IDPDG(IPOM)
              ENDIF
           ENDIF
        ENDDO
C ... Modify hard initial parton if necessary
        IF((TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC').OR.
     &       IPRO.EQ.96.OR.IPRO.EQ.97) THEN
          DO I=8,9
            IF(IDHEP(I).EQ.22) THEN
              IDHW(I)=13
              IDHEP(I)=IDPDG(13)
            ENDIF
            IF(IDHEP(I+5).EQ.22) THEN
              IDHW(I+5)=13
              IDHEP(I+5)=IDPDG(13)
            ENDIF
          ENDDO
        ENDIF
      ELSE
        PRINT*, 'HWFXER : not in diffraction mode - nothing done'
      ENDIF
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE HWHIGP
C-----------------------------------------------------------------------
C     HIGGS PRODUCTION VIA PHOTON OR GLUON FUSION
C     MEAN EVWGT = HIGGS PRODN C-S * BRANCHING FRACTION IN NB.
C     DERIVED FROM HWHIGS.
C
C     M.Boonekamp, T. Kucs Aug 2003 : 
C
C     If TYPINT = 'QED' then gamma + gamma -> Higgs cross section is 
C     calculated.
C     If TYPINT = 'QCD' then g + g -> Higgs cross-section needed for 
C     exclusive DPE is calculated.
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION HWUALF,HWHIGT,HWRGEN,HWUSQR,HWUAEM,BRHIGQ,EMH,
     &     CSFAC,EVSUM(13),EMFAC,CV,CA,BR,RWGT,E1,E2,EMQ,GFACTR
      DOUBLE PRECISION TAUT,TAUW,EMW,TAUTR,TAUTI,TAUWR,TAUWI,SUMR,
     &     SUMI,HWIDTH, NLOConst
      INTEGER IDEC,I,J,ID1,ID2
      EXTERNAL HWUALF,HWHIGT,HWRGEN,HWUSQR,HWUAEM
      SAVE CSFAC,BR
      EQUIVALENCE (EMW,RMASS(198))

C ... begin R.S.
      include 'CHIDe.inc'
      INTEGER CHIDeN, CHIDeI
      LOGICAL IsCHIDeHiggs
      DOUBLE PRECISION CHIDeSigma, CHIDeJac

      DOUBLE PRECISION CHIDeQ, CHIDeQp
      DOUBLE PRECISION CHIDePhi, CHIDePhip
      DOUBLE PRECISION CHIDeK(2), CHIDeKp(2)

      DOUBLE PRECISION CHIDeKmax
      DOUBLE PRECISION CHIDedotdiff
      EXTERNAL CHIDedotdiff

      IsCHIDeHiggs = (NFLUX.EQ.18 .AND. IPROC.EQ.19999)
C ... end R.S.
      

      IF (GENEV) THEN
C ... fill as gamma+gamma -> Higgs, whatever TYPINT is
        IDN(1)=59
        IDN(2)=59
        IDCMF=201+IHIGGS
        XX(1)=1.
        XX(2)=1.
       CALL HWEONE
      ELSE
        EVWGT=0.
        EMH=RMASS(201+IHIGGS)
        EMFAC=1.D0
C ... M.B. : commented below; should be done in HWEGAM
!!!!        IF(IMSSM.EQ.0)CALL HWHIGM(EMH,EMFAC)
        PHEP(5,3)=EMH
        EMSCA=EMH

*        IF (EMSCA.NE.EMLST) THEN

          EMLST=EMH
          IF(TYPINT.EQ.'QED'.AND.TYPEPR.EQ.'EXC') THEN
c---Begin modif by Tibor Kucs 08/13/2003
c...Gamma + gamma -> Higgs cross-section
c...Taken from SUBROUTINE HWDHIG
             TAUT=(2*RMASS(6)/EMH)**2
             TAUW=(2*EMW/EMH)**2
             CALL HWDHGC(TAUT,TAUTR,TAUTI)
             CALL HWDHGC(TAUW,TAUWR,TAUWI)
             SUMR=4./3*(  - 2*TAUT*( 1 + (1-TAUT)*TAUTR ) ) * ENHANC(6)
     &            +(2 + 3*TAUW*( 1 + (2-TAUW)*TAUWR ) ) * ENHANC(10)
             SUMI=4./3*(  - 2*TAUT*(     (1-TAUT)*TAUTI ) ) * ENHANC(6)
     &            +(    3*TAUW*(     (2-TAUW)*TAUWI ) ) * ENHANC(10)
             GFACTR=ALPHEM/(8.*SWEIN*EMW**2)
             HWIDTH=GFACTR*.03125*(ALPHEM/PIFAC)**2
     &            *EMH**3 * (SUMR**2 + SUMI**2)
             CSFAC=GEV2NB*8.D0*(PI**2)*HWIDTH/(RMASS(201)**3)
c---End modif by Tibor Kucs
          ELSEIF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC') THEN
C ... here the gluon-gluon to Higgs cross-section
C ... w.r.t. HWHIGS, the XLMIN factor is not there, since it comes
C     from the pdf's.
C     An additional factor 16 = 2(Nc^2-1) is included
C ... (gluon color&spin coherence)
C ... O.K. NLO k-factor NLOConst         
             GFACTR=GEV2NB*HWUAEM(EMH**2)/(576.*SWEIN*RMASS(198)**2)
             EMQ=RMASS(13)
             CSFAC=GFACTR*HWHIGT(EMH)*HWUALF(1,EMH)**2*EMFAC*16
             NLOConst = PI + 5.5/PI;
             CSFAC=CSFAC*(1+HWUALF(1,EMH)*NLOConst)
C ... begin R.S.
             IF(IsCHIDeHiggs) THEN
C     Average over unphysical variables to increase the efficiency
                CSFAC=ZERO
                CHIDeN = 100
                DO CHIDeI = 1, CHIDeN
                   CHIDeKmax=20.0 ! Maximal momentum in gluon loop
                   CHIDeJac=1.

                   CHIDeQ=CHIDeKmax*HWRGEN(1)
                   CHIDePhi=2.*PIFAC*HWRGEN(2)
                   CHIDeJac=CHIDeJac*2.*PIFAC*CHIDeQ*CHIDeKmax
                   CHIDeK(1)=CHIDeQ*cos(CHIDePhi)
                   CHIDeK(2)=CHIDeQ*sin(CHIDePhi)
                   
                   CHIDeQp=CHIDeKmax*HWRGEN(1)
                   CHIDePhip=2.*PIFAC*HWRGEN(2)                   
                   CHIDeJac=CHIDeJac*2.*PIFAC*CHIDeQp*CHIDeKmax
                   CHIDeKp(1)=CHIDeQp*cos(CHIDePhip)
                   CHIDeKp(2)=CHIDeQp*sin(CHIDePhip) 
                   
                   call CHIDeHiggs(CHIDeSigma,CHIDeZ1,CHIDeZ2,
     &                             CHIDeK1,CHIDeK3,CHIDeK,CHIDeKp)
                   CHIDeSigma=CHIDeJac*CHIDeSigma
                   CSFAC=CSFAC+CHIDeSigma
                ENDDO
                CSFAC=CSFAC/CHIDeN
                IF(CSFAC.LE.0) CSFAC=ZERO
             ENDIF        
C     ... end R.S.

          ELSE
            PRINT*, 'In HWHIGP : wrong TYP* setting ', TYPINT, TYPEPR, 
     &              ' - STOP'
            STOP
          ENDIF
C INCLUDE BRANCHING RATIO OF HIGGS
          IDEC=MOD(IPROC,100)
          BR=1
          IF(IMSSM.EQ.0)THEN
C SM case
            IF(IDEC.EQ.0) THEN
              BRHIGQ=0
              DO 30 I=1,6
 30             BRHIGQ=BRHIGQ+BRHIG(I)
              BR=BRHIGQ
            ELSEIF (IDEC.EQ.10) THEN
              CALL HWDBOZ(198,ID1,ID2,CV,CA,BR,1)
              CALL HWDBOZ(199,ID1,ID2,CV,CA,BR,1)
              BR=BR*BRHIG(IDEC)
            ELSEIF (IDEC.EQ.11) THEN
              CALL HWDBOZ(200,ID1,ID2,CV,CA,BR,1)
              CALL HWDBOZ(200,ID1,ID2,CV,CA,BR,1)
              BR=BR*BRHIG(IDEC)
            ELSEIF (IDEC.LE.12) THEN
              BR=BRHIG(IDEC)
            ENDIF
          ENDIF
         ENDIF
         EVWGT=EVWGT + CSFAC*BR
*      ENDIF
  999 END
C-----------------------------------------------------------------------
CDECK  ID>, HWHPH2.
*CMZ :-        -12/01/93  10.12.43  by  Bryan Webber
*-- Author :    Ian Knowles
C-----------------------------------------------------------------------
      SUBROUTINE HWH2PH
C-----------------------------------------------------------------------
C     QQD direct photon pair production: mean EVWGT = sigma in nb
C
C     M.B., Nov 2003 : small bug correction
C
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWRGEN,HWRUNI,HWUALF,HWHPPB,EPS,RCS,ET,EJ,KK,KK2,
     & YJ1INF,YJ1SUP,Z1,YJ2INF,YJ2SUP,Z2,FACT,FACTR,RS,S,T,U,CSTU,TQSQ,
     & DSTU,HCS
      INTEGER ID,ID1,ID2
      EXTERNAL HWRGEN,HWRUNI,HWUALF,HWHPPB
      SAVE HCS,CSTU,DSTU,FACT
      PARAMETER (EPS=1.D-9)
      IF (GENEV) THEN
        RCS=HCS*HWRGEN(0)
      ELSE
        EVWGT=0.
        CALL HWRPOW(ET,EJ)
        KK=ET/PHEP(5,3)
        KK2=KK**2
        IF (KK.GE.ONE) RETURN
        YJ1INF=MAX( YJMIN , LOG((1.-SQRT(1.-KK2))/KK) )
        YJ1SUP=MIN( YJMAX , LOG((1.+SQRT(1.-KK2))/KK) )
        IF (YJ1INF.GE.YJ1SUP) RETURN
        Z1=EXP(HWRUNI(1,YJ1INF,YJ1SUP))
        YJ2INF=MAX( YJMIN , -LOG(2./KK-1./Z1) )
        YJ2SUP=MIN( YJMAX ,  LOG(2./KK-Z1) )
        IF (YJ2INF.GE.YJ2SUP) RETURN
        Z2=EXP(HWRUNI(2,YJ2INF,YJ2SUP))
        XX(1)=0.5*(Z1+Z2)*KK
        IF (XX(1).GE.ONE) RETURN
        XX(2)=XX(1)/(Z1*Z2)
        IF (XX(2).GE.ONE) RETURN
        COSTH=(Z1-Z2)/(Z1+Z2)
        S=XX(1)*XX(2)*PHEP(5,3)**2
        RS=0.5*SQRT(S)
        T=-0.5*S*(1.-COSTH)
        U=-S-T
        EMSCA=SQRT(2.*S*T*U/(S*S+T*T+U*U))
        FACT=GEV2NB*PIFAC*0.5*ET*EJ*(YJ1SUP-YJ1INF)*(YJ2SUP-YJ2INF)
     &      *(ALPHEM/S)**2
        CALL HWSGEN(.FALSE.)
        CSTU=2.*(U/T+T/U)/CAFAC
        IF (DISF(13,1).GT.EPS.AND.DISF(13,2).GT.EPS) THEN
           TQSQ=0.
           DO 10 ID=1,6
  10       IF (RMASS(ID).LT.RS) TQSQ=TQSQ+QFCH(ID)**2
           DSTU=DISF(13,1)*DISF(13,2)*FACT*HWHPPB(S,T,U)
     &         /64.*(HWUALF(1,EMSCA)*TQSQ/PIFAC)**2
        ENDIF
      ENDIF
      HCS=0.
      DO 30 ID=1,6
      FACTR=FACT*CSTU*QFCH(ID)**4
C q+qbar ---> gamma+gamma
      ID1=ID
      ID2=ID+6
      IF (DISF(ID1,1).LT.EPS.OR.DISF(ID2,2).LT.EPS) GOTO 20
      HCS=HCS+FACTR*DISF(ID1,1)*DISF(ID2,2)
      IF (GENEV.AND.HCS.GT.RCS) CALL HWHQCP(59,59,2134,61,*99)
C qbar+q ---> gamma+gamma
  20  ID1=ID+6
      ID2=ID
      IF (DISF(ID1,1).LT.EPS.OR.DISF(ID2,2).LT.EPS) GOTO 30
      HCS=HCS+FACTR*DISF(ID1,1)*DISF(ID2,2)
      IF (GENEV.AND.HCS.GT.RCS) CALL HWHQCP(59,59,2134,62,*99)
  30  CONTINUE
C g+g ---> gamma+gamma
      ID1=13
      ID2=13
C M.B. Protection added below
      IF(DISF(ID1,1).GT.EPS.AND.DISF(ID2,2).GT.EPS) HCS=HCS+DSTU
      IF (GENEV.AND.HCS.GT.RCS) CALL HWHQCP(59,59,2134,63,*99)
      EVWGT=HCS
      RETURN
C Generate event
  99  IDN(1)=ID1
      IDN(2)=ID2
      IDCMF=15
      CALL HWETWO
  999 END
C-----------------------------------------------------------------------
      SUBROUTINE HWHQPM
C     HARD PROCESS: GAMGAM --> QQBAR/LLBAR/W+W-/ZZ/GAMGAM
C     MEAN EVENT WEIGHT = CROSS-SECTION IN NB AFTER CUTS ON PT
C
C     M.Boonekamp T.Kucs, Aug 2003 : 
C        Modifications made to allow for Jz=0 dijet production
C        in a gluon gluon collision.
C     O.Kepka, April, 2009 : 
C        Anomalous ZZ
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION RCS,HCS,RS,S,EMSQ,BE,TMIN,TMAX,T,U,FACTR,Q,CFAC,
     $     HWRGEN
      INTEGER HQ,ID3,ID4,I1,I2
C --- Begin modif by Tibor Kucs 08/14/2003 
C     Updated by Maarten Boonekamp 11/03/2003
      INTEGER IFLAVR,NFLAVR, IFIRST, ILAST
      DATA NFLAVR / 6 /
      DOUBLE PRECISION XMS,ALPHAS,TFACT,SGGQQ,XMSQ,SINTH,ETSQ,SGGGG,
     $     HWUALF,BESQ,XSQQ(500),SFLAV,SINSQ,RR,PTSQ
      EXTERNAL HWUALF
      SAVE SGGGG, SGGQQ, XSQQ
c---End modif by Tibor Kucs 08/14/2003
      SAVE HCS,FACTR,HQ,RS
c --- edited by Oldrich Kepka      
      DOUBLE PRECISION BETA, AMP2
      DOUBLE PRECISION EMSQ1, RGEN1, RGEN2
      EXTERNAL sqme_aaww_c
C ... begin R.S.
      include 'CHIDe.inc'
      INTEGER CHIDeN, CHIDeI
      LOGICAL IsCHIDeGG,IsCHIDeDiphoton
      DOUBLE PRECISION CHIDeSigma, CHIDeJac

      DOUBLE PRECISION CHIDeQ, CHIDeQp
      DOUBLE PRECISION CHIDePhi, CHIDePhip
      DOUBLE PRECISION CHIDeK(2), CHIDeKp(2)

      DOUBLE PRECISION CHIDeKmax
      DOUBLE PRECISION CHIDedotdiff,CHIDedot
      EXTERNAL CHIDedotdiff,CHIDedot

      IsCHIDeGG = (NFLUX.EQ.18 .AND. IPROC.EQ.16012)
      IsCHIDeDiphoton = (NFLUX.EQ.18 .AND. IPROC.EQ.16059)
C ... end R.S.      
      
      IF (GENEV) THEN
        RCS=HCS*HWRGEN(0)
      ELSE
        IFIRST = 1
        ILAST = NFLAVR
        EVWGT=0.
        RS=PHEP(5,3)
        S=RS**2
        HQ=MOD(IPROC,100)
C Checks on HQ
        IF (HQ.EQ.0.OR.HQ.EQ.1) THEN 
          EMSQ=0
          BESQ=1
          BE=1
          CFAC=3
        ELSE
c ... M.B. : if HQ=13, do not change it
c ... O.K. : HQ=15 ZZ
c ... M.S./O.K. : HQ=16 AA
c ... M.S.      : HQ=60,61,62,63,64 SM AA + AAANOM=3 def
          IF (HQ.GT.6.AND.HQ.LE.10) HQ=2*HQ+107
          IF (HQ.EQ.15) HQ=200 ! ZZ
          IF (HQ.EQ.16) HQ=59 ! AA
          IF (HQ.GE.60.AND.HQ.LE.64) THEN
            HQ=59 ! AA
            AAANOM=3
          ENDIF
          IF (HQ.EQ.127) HQ=198 
          IF (HQ.GE.21.AND.HQ.LE.26) THEN 
            HQ=HQ+400-20
            IFIRST = 401
            ILAST = 406
          ENDIF
          EMSQ=RMASS(HQ)**2
          IF(HQ.EQ.12) EMSQ=RMASS(13)**2 ! R.S.
          BESQ=1-4*EMSQ/S
          IF (BESQ.LT.ZERO) RETURN
          BE=SQRT(BESQ)
          CFAC=1
          IF (HQ.LE.6) CFAC=3
        ENDIF
C Kinematics 
           TMIN=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMIN**2)/S,ZERO)))
           TMAX=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMAX**2)/S,ZERO)))
           IF (TMIN.GE.TMAX) RETURN

           RGEN1=HWRGEN(1)
           RGEN2=HWRGEN(2)
           IF(IsCHIDeGG.OR.IsCHIDeDiphoton) THEN
             T = -CHIDeQQ2*(CHIDeB1+CHIDeB2)/CHIDeB1
           ELSE
             T=-(TMAX/TMIN)**RGEN1*TMIN
             IF(RGEN2.GT.HALF) T=-S-T
           ENDIF
           U=-S-T
           COSTH=(T-U)/(BE*S)
           EMSCA=SQRT(2.*S*T*U/(S*S+T*T+U*U)) 
C Cross-sections
        IF (HQ.NE.198.AND.HQ.NE.200.AND.HQ.NE.59) THEN
C --- Begin modif by Tibor Kucs 08/14/2003
C     Updated by Maarten Boonekamp 11/03/2003
c ... If QED use the original implementation :
c     this is for lepton - antilepton
          IF(TYPINT.EQ.'QED') THEN
            FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     $            *2*PIFAC*CFAC*ALPHEM**2/S**2
     $            *((U-4*EMSQ)/T+(T-4*EMSQ)/U-4*(EMSQ/T+EMSQ/U)**2)
c ... If QCD / EXC take the Jz=0, color singlet QCD cross-section
          ELSEIF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC') THEN
          if(.NOT.(IsCHIDeGG.OR.IsCHIDeDiphoton)) 
     &        TFACT=PIFAC*HWUALF(1,EMSCA)**2/3.D0
            SGGQQ=0d0 
            SGGGG=0d0
            SINSQ=(1.D0-COSTH)*(1.D0+COSTH)
c ...... gg -> qq
            IF((HQ.GE.1.AND.HQ.LE.6).OR.(HQ.GE.401.AND.HQ.LE.406)
     &           .OR.HQ.EQ.11) THEN
            DO IFLAVR=IFIRST,ILAST
            XMSQ=RMASS(IFLAVR)**2
            IF(S.GT.4*XMSQ.AND.(IFLAVR.EQ.HQ.OR.HQ.EQ.11)) THEN
              PTSQ=(S/4d0-XMSQ)*SINSQ
              ETSQ=XMSQ+PTSQ
              XSQQ(IFLAVR)=0d0
              IF(IFLAVR.GE.1.AND.IFLAVR.LE.6) THEN
                BESQ=1d0-4*XMSQ/S
                XSQQ(IFLAVR) = TFACT*XMSQ/S/ETSQ**2*BESQ
                SGGQQ=SGGQQ+XSQQ(IFLAVR)
                ELSEIF(IFLAVR.GE.401.AND.IFLAVR.LE.406) THEN
                XSQQ(IFLAVR) = 2*TFACT*XMSQ**2/S**2/ETSQ**2
                SGGQQ=SGGQQ+XSQQ(IFLAVR)
              ENDIF
            ENDIF
            ENDDO 
            ELSEIF(IsCHIDeGG.OR.IsCHIDeDiphoton) THEN
            CHIDeN = 50
            FACTR = 0d0
            DO CHIDeI = 1, CHIDeN
            CHIDeKmax=20.0 ! Maximal momentum in gluon loop
            CHIDeJac=1.

            CHIDeQ=CHIDeKmax*HWRGEN(1)
            CHIDePhi=2.*PIFAC*HWRGEN(2)
            CHIDeJac=CHIDeJac*2.*PIFAC*CHIDeQ*CHIDeKmax
            CHIDeK(1)=CHIDeQ*cos(CHIDePhi)
            CHIDeK(2)=CHIDeQ*sin(CHIDePhi)

            CHIDeQp=CHIDeKmax*HWRGEN(1)
            CHIDePhip=2.*PIFAC*HWRGEN(2)                   
            CHIDeJac=CHIDeJac*2.*PIFAC*CHIDeQp*CHIDeKmax
            CHIDeKp(1)=CHIDeQp*cos(CHIDePhip)
            CHIDeKp(2)=CHIDeQp*sin(CHIDePhip) 

            if(IsCHIDeGG) call CHIDeGG(CHIDeSigma,
     &             CHIDeK,CHIDeKp,CHIDeK1,CHIDeK2,CHIDeK3, 
     &             CHIDeB1,CHIDeB2,CHIDeA1, CHIDeA2)

            if(IsCHIDeDiphoton) call CHIDeDiphoton(CHIDeSigma,
     &             CHIDeK,CHIDeKp,CHIDeK1,CHIDeK2,CHIDeK3, 
     &             CHIDeB1,CHIDeB2,CHIDeA1, CHIDeA2)

            FACTR=FACTR+CHIDeSigma*CHIDeJac
            ENDDO

            FACTR=FACTR/CHIDeN
            IF(FACTR.LE.0) FACTR=ZERO

c ...... gg -> gg : include a factor 1/2! to compensate for double PhSp integration
            ELSEIF(HQ.EQ.13) THEN
            XMSQ=RMASS(HQ)**2
            PTSQ=(S/4d0-XMSQ)*SINSQ
            ETSQ=XMSQ+PTSQ
            SGGGG=TFACT*(27.D0/2.D0)/ETSQ**2 / 2d0
          ENDIF
c ...... Total 
          if(.not.(IsCHIDeGG.OR.IsCHIDeDiphoton)) 
     &        FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)*(SGGQQ+SGGGG)
        ELSE
          PRINT*, 'In HWHQPM : settings not consistent',TYPINT,TYPEPR,
     &              ' - STOP'
          STOP
        ENDIF
      ELSE

c O.K. 01/11/2007 modified to include call to O'Mega matrix
c element for AA->WW
c generation of full kinematics added because neded by O'Mega

c       original formula
c               FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
c     $         *6*PIFAC*CFAC*ALPHEM**2/S**2
c     $         *(1-S/(T*U)*(4D0/3*S+2*EMSQ)
c     $         +(S/(T*U))**2*(2D0/3*S**2+2*EMSQ**2))

          IF(AAANOM.EQ.0)THEN
          ! original herwig formula
            FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     $             *6*PIFAC*CFAC*ALPHEM**2/S**2
     $             *(1-S/(T*U)*(4D0/3*S+2*EMSQ)
     $             +(S/(T*U))**2*(2D0/3*S**2+2*EMSQ**2))

          ELSEIF(AAANOM.EQ.1)THEN
          ! AA->WW standard model formula but no approximation on Mw
          ! assumed. Have to generate kinematics properly
          ! complete formula arXiv:hep-ph/9601355

              TMIN=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMIN**2)/S,ZERO)))-EMSQ 
              TMAX=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMAX**2)/S,ZERO)))-EMSQ
              IF (TMIN.GE.TMAX) RETURN
c             T=-(TMAX/TMIN)**WRGEN(1)*TMIN
              T=-(TMAX/TMIN)**RGEN1*TMIN

cc            IF (HWRGEN(2).GT.HALF) then
              IF (RGEN2.GT.HALF) then
                  T=-S-T+2*EMSQ
              ENDIF   
              U=-S-T+2*EMSQ

              FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     $             *6*pifac*cfac*alphem**2/s**2
     $             *(1-2*s*(2*s+3*emsq)/(3*(emsq-t)*(emsq-u))
     $           +2*s**2*(s**2+3*emsq*emsq)/(3*(emsq-t)**2*(emsq-u)**2))

          ELSEIF(AAANOM.EQ.2)THEN
c          PTMIN>0 only for AAAA
              IF(HQ.EQ.59.AND.PTMIN.EQ.0) THEN
                  PTMIN=1
              ENDIF
              TMIN=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMIN**2)/S,ZERO)))-EMSQ 
              TMAX=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMAX**2)/S,ZERO)))-EMSQ
              IF (TMIN.GE.TMAX) RETURN
c           T=-(TMAX/TMIN)**WRGEN(1)*TMIN
              T=-(TMAX/TMIN)**RGEN1*TMIN

c           IF (HWRGEN(2).GT.HALF) then
              IF (RGEN2.GT.HALF) then
                     T=-S-T+2*EMSQ
              ENDIF   
              U=-S-T+2*EMSQ
             

              COSTH=(2*T+S-2*EMSQ)/SQRT(S**2-4*EMSQ*S)
              EMSCA=SQRT(2.*S*T*U/(S*S+T*T+U*U)) 

C ... O.K./M.S. Calling anomalous aaww or aazz coupling or aaaa coupling
              IF(HQ.EQ.198) call sqme_aaww_c(AMP2, S, T, alphem,
     $        SQRT(EMSQ), SWEIN, D_KAPPA, LAMBDA, A0W, ACW, ANOMCUTOFF)
              IF(HQ.EQ.200) call sqme_aazz_c(AMP2, S, T, alphem,
     $        SQRT(EMSQ), SWEIN, A0Z, ACZ, ANOMCUTOFF)
              IF(HQ.EQ.59) THEN 
C ... M.S. Former routine from compHEP
c              call sqme_aaaa_c(AMP2, S, T, alphem,
c     $        SQRT(EMSQ), SWEIN, D_KAPPA, LAMBDA, A1A, A2A, ANOMCUTOFF)
              call eft_sqme_aaaa_c(AMP2, S, T, 1, A1A, A2A, ANOMCUTOFF)
              ENDIF
              FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     $         *2*PIFAC/(64.*PIFAC**2)/S**2*2d0*AMP2

c              FACTR1=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
c     $             *6*pifac*cfac*alphem**2/s**2
c     $             *(1-2*s*(2*s+3*emsq)/(3*(emsq-t)*(emsq-u))
c     $           +2*s**2*(s**2+3*emsq*emsq)/(3*(emsq-t)**2*(emsq-u)**2))


          ELSEIF(AAANOM.EQ.3)THEN
c          PTMIN>0 only for AAAA
              IF(HQ.EQ.59.AND.PTMIN.EQ.0) THEN
                  PTMIN=1
              ENDIF
              TMIN=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMIN**2)/S,ZERO)))-EMSQ 
              TMAX=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMAX**2)/S,ZERO)))-EMSQ
              IF (TMIN.GE.TMAX) RETURN
c           T=-(TMAX/TMIN)**WRGEN(1)*TMIN
              T=-(TMAX/TMIN)**RGEN1*TMIN

c           IF (HWRGEN(2).GT.HALF) then
              IF (RGEN2.GT.HALF) then
                     T=-S-T+2*EMSQ
              ENDIF   
              U=-S-T+2*EMSQ
             

              COSTH=(2*T+S-2*EMSQ)/SQRT(S**2-4*EMSQ*S)
              EMSCA=SQRT(2.*S*T*U/(S*S+T*T+U*U)) 

C ... M.S. Calling SM and EXOTIC exclusive photon pair production
              IF(HQ.EQ.59.AND.IPROC.EQ.16060) THEN
              call sm_sqme_aaaa_c(AMP2, S, T, 0)
              ENDIF
              IF(HQ.EQ.59.AND.IPROC.EQ.16061) THEN
              call sm_sqme_aaaa_c(AMP2, S, T, 1)
              ENDIF       
              IF(HQ.EQ.59.AND.IPROC.EQ.16062) THEN
              call sm_sqme_aaaa_c(AMP2, S, T, 2)
              ENDIF               
              IF(HQ.EQ.59.AND.IPROC.EQ.16063) THEN
              call bsmv_sqme_aaaa_c(AMP2, S, T, 1, 0,
     $        AAM, AAQ, AAN)
              ENDIF
              IF(HQ.EQ.59.AND.IPROC.EQ.16064) THEN
              call bsmf_sqme_aaaa_c(AMP2, S, T, 1, 0,
     $        AAM, AAQ, AAN)
              ENDIF
             FACTR=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     $         *2*PIFAC/(64.*PIFAC**2)/S**2*2d0*AMP2

          ELSE
             print *, 'Nonstandard choice of AAANOM=', AAANOM
             stop
          ENDIF   
        ENDIF
      ENDIF
C ... For 'QED' include EM charge
      IF(TYPINT.EQ.'QED') THEN 
        HCS=0.
        XX(1)=1.
        XX(2)=1.
        IF (HQ.EQ.0) THEN
          I1=1
          I2=6
        ELSE
          I1=HQ
          I2=HQ
        ENDIF
        DO 10 ID3=I1,I2
          IF (RS.GT.2*RMASS(ID3)) THEN
            Q=ICHRG(ID3)
            IF (HQ.LE.6) Q=Q/THREE
            ID4=ID3+6
            IF (HQ.EQ.198) ID4=199
c O.K ZZ                    
            IF (HQ.EQ.200) THEN 
               ID4=200
               Q=1
            ENDIF
c O.K AA                    
            IF (HQ.EQ.59) THEN 
               ID4=59
               Q=1
            ENDIF
            HCS=HCS+Q**4
            IF (GENEV.AND.HCS.GT.RCS) CALL HWHQCP(ID3,ID4,1243,61,*99)
            ENDIF
 10      CONTINUE
C ... For 'QCD' no EM charge, it is set to ONE
      ELSE
        HCS=1.
        XX(1)=1.
        XX(2)=1.
        IF(GENEV) THEN
           IF(IsCHIDeDiphoton) THEN
            ID3=59           ! photons in the final state 
            ID4=59
            CALL HWHQCP(ID3,ID4,1243,61,*99)
            GOTO 99
          ELSEIF(IsCHIDeGG .OR. (HWRGEN(0)*(SGGQQ+SGGGG).LE.SGGGG)) THEN
            ID3=13           ! gluons in the final state 
            ID4=13
            CALL HWHQCP(ID3,ID4,1243,61,*99)
            GOTO 99
          ELSE
            RR=HWRGEN(0)*SGGQQ
            SFLAV=0d0
            DO IFLAVR=IFIRST,ILAST
              SFLAV=SFLAV+XSQQ(IFLAVR)
              IF(RR.LT.SFLAV) THEN
                ID3=IFLAVR
                ID4=ID3+6
                CALL HWHQCP(ID3,ID4,1243,61,*99)
                GOTO 99
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      EVWGT=FACTR*HCS
      RETURN
c---End modif by Tibor Kucs 08/14/2003
 99   IDN(1)=59
      IDN(2)=59      
      IDCMF=15
      IF(IsCHIDeGG.OR.IsCHIDeDiphoton) THEN
        CALL HWETWO_MOD(CHIDePhi2) ! not sure if this is
        ! the right angle
      ELSE 
        CALL HWETWO
      ENDIF
      END

CDECK  ID>, HWETWO.
*CMZ :-        -26/04/91  11.11.55  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWETWO_MOD(phi)
C-----------------------------------------------------------------------
C     SETS UP 2->2 HARD SUBPROCESS
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
C ... begin R.S.      
      double precision phi
C ... end R.S.
      DOUBLE PRECISION HWUMBW,HWUPCM,PA,PCM
      INTEGER ICMF,IBM,I,J,K,IHEP,NTRY
      EXTERNAL HWUPCM
C---INCOMING LINES
      ICMF=NHEP+3
      DO 15 I=1,2
      IBM=I
C---FIND BEAM AND TARGET
      IF (JDAHEP(1,I).NE.0) IBM=JDAHEP(1,I)
      IHEP=NHEP+I
      IDHW(IHEP)=IDN(I)
      IDHEP(IHEP)=IDPDG(IDN(I))
      ISTHEP(IHEP)=110+I
      JMOHEP(1,IHEP)=ICMF
      JMOHEP(I,ICMF)=IHEP
      JDAHEP(1,IHEP)=ICMF
C---SPECIAL - IF INCOMING PARTON IS INCOMING BEAM THEN COPY IT
      IF (XX(I).EQ.ONE.AND.IDHW(IBM).EQ.IDN(I)) THEN
        CALL HWVEQU(5,PHEP(1,IBM),PHEP(1,IHEP))
        IF (I.EQ.2) PHEP(3,IHEP)=-PHEP(3,IHEP)
      ELSE
        PHEP(1,IHEP)=0.
        PHEP(2,IHEP)=0.
        PHEP(5,IHEP)=RMASS(IDN(I))
        PA=XX(I)*(PHEP(4,IBM)+ABS(PHEP(3,IBM)))
        PHEP(4,IHEP)=0.5*(PA+PHEP(5,IHEP)**2/PA)
        PHEP(3,IHEP)=PA-PHEP(4,IHEP)
      ENDIF
 15   CONTINUE
      PHEP(3,NHEP+2)=-PHEP(3,NHEP+2)
C---HARD CENTRE OF MASS
      IDHW(ICMF)=IDCMF
      IDHEP(ICMF)=IDPDG(IDCMF)
      ISTHEP(ICMF)=110
      CALL HWVSUM(4,PHEP(1,NHEP+1),PHEP(1,NHEP+2),PHEP(1,ICMF))
      CALL HWUMAS(PHEP(1,ICMF))
C---OUTGOING LINES
      NTRY=0
 19   CONTINUE
      DO 20 I=3,4
      IHEP=NHEP+I+1
      IDHW(IHEP)=IDN(I)
      IDHEP(IHEP)=IDPDG(IDN(I))
      ISTHEP(IHEP)=110+I
      JMOHEP(1,IHEP)=ICMF
      JDAHEP(I-2,ICMF)=IHEP
 20   PHEP(5,IHEP)=HWUMBW(IDN(I))
      PCM=HWUPCM(PHEP(5,NHEP+3),PHEP(5,NHEP+4),PHEP(5,NHEP+5))
      IF (PCM.LT.ZERO) THEN
        NTRY=NTRY+1
        IF (NTRY.LE.NETRY) GO TO 19
        CALL HWWARN('HWETWO',103,*999)
      ENDIF
      IHEP=NHEP+4
      PHEP(4,IHEP)=SQRT(PCM**2+PHEP(5,IHEP)**2)
      PHEP(3,IHEP)=PCM*COSTH
      PHEP(1,IHEP)=SQRT((PCM+PHEP(3,IHEP))*(PCM-PHEP(3,IHEP)))
C ... begin R.S
C      CALL HWRAZM(PHEP(1,IHEP),PHEP(1,IHEP),PHEP(2,IHEP)) 
       PHEP(2,IHEP) = PHEP(1,IHEP)*sin(phi)
       PHEP(1,IHEP) = PHEP(1,IHEP)*cos(phi)
C ... end R.S      
      CALL HWULOB(PHEP(1,NHEP+3),PHEP(1,IHEP),PHEP(1,IHEP))
      CALL HWVDIF(4,PHEP(1,NHEP+3),PHEP(1,IHEP),PHEP(1,NHEP+5))
C---SET UP COLOUR STRUCTURE LABELS
      DO 30 I=1,4
      J=I
      IF (J.GT.2) J=J+1
      K=ICO(I)
      IF (K.GT.2) K=K+1
      JMOHEP(2,NHEP+J)=NHEP+K
   30 JDAHEP(2,NHEP+K)=NHEP+J
      NHEP=NHEP+5
  999 END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE HWHQPP
C     Hard processes: gamma + gamma --> gamma + gamma
C                     gluon + gluon --> gamma + gamma exclusive
C
C     Mean event weight = cross section in NB after cuts on pT
C
C     NOTE: The implementation of the QED process does not contain
C           the contribution from the box diagram containing W bosons.
C
C     Tibor Kucs, Dec 2003 : 
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION RCS,HCS,RS,S,EMSQ,BE,TMIN,TMAX,T,U,FACTR,Q,CFAC,
     $     HWRGEN,A2,XMSQ,ALPHAS,CHGSUM,FACT,HWUALF,HWHBOX
      COMPLEX*16 CA1,CM1111,CM1122,CM2211,CM1212,CM2222,CM
      INTEGER HQ,ID3,ID4,I1,I2
      EXTERNAL HWUALF,HWRGEN,CA1,A2,HWHBOX
      SAVE HCS,FACTR,HQ,RS            
      IF (GENEV) THEN
        GOTO 99
      ELSE
         EVWGT=0.
         RS=PHEP(5,3)
         S=RS**2
         EMSQ=0
         BE=1.
         TMIN=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMIN**2)/S,ZERO)))
         TMAX=S/2*(1-SQRT(MAX(1-4*(EMSQ+PTMAX**2)/S,ZERO)))
         IF (TMIN.GE.TMAX) RETURN
         T=-(TMAX/TMIN)**HWRGEN(1)*TMIN
         IF (HWRGEN(2).GT.HALF) T=-S-T
         U=-S-T
         COSTH=(T-U)/(BE*S)
         EMSCA=SQRT(2.*S*T*U/(S*S+T*T+U*U))        
c...gamma + gamma --> gamma + gamma
         IF(TYPINT.EQ.'QED'.AND.TYPEPR.EQ.'EXC') THEN
            EVWGT=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     &            *(ALPHEM**4)*HWHBOX(S,T,U)/16./PI/S**2
         ELSEIF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC') THEN
c...gluon + gluon --> gamma + gamma exclusive
            ALPHAS=HWUALF(1,EMSCA)
            CHGSUM=(11./9.)**2  ! charge^2 sum for five flavors
            FACT=(ALPHEM*ALPHAS/2.)**2/16./PI/S**2
            EVWGT=-GEV2NB*2*LOG(TMAX/TMIN)*MAX(T,U)
     &            *FACT*CHGSUM*HWHBOX(S,T,U)
         ELSE
            PRINT*, 'In HWHQPP : not implemented ',TYPINT,TYPEPR, 
     &              ' - STOP'
            STOP
         END IF
      END IF
      XX(1)=1.
      XX(2)=1.
      RETURN
c...the ingoing particles are recorded as photons, to avoid
c   problems with gluons having XX(i)=1  
 99   IF(TYPINT.EQ.'QED'.OR.TYPINT.EQ.'QCD') THEN 
        IDN(1)=59
        IDN(2)=59
        IDN(3)=59   ! photons in the final state as well
        IDN(4)=59
      END IF
      IDCMF=15
      CALL HWETWO
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION HWHBOX(S,T,U)
C-----------------------------------------------------------------------
C     Quark box diagram contribution to photon/gluon scattering
C     Internal quark mass neglected: m_q << U,T,S
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION S,T,U,A2
      COMPLEX*16 CA1,CM1111,CM1122,CM2211,CM1212,CM2222,CM
      EXTERNAL CA1,A2
c...Matrix elements for given polarizations
      CM1111=CA1(S,T)+CA1(S,U)+A2(T,U)
      CM1122=CA1(S,T)+CA1(S,U)-A2(T,U)-32./3.
      CM2211=CM1122
      CM1212=CA1(S,T)-CA1(S,U)+A2(T,U)-32./3.
      CM2222=CM1111-32.
c...gamma + gamma --> gamma + gamma
      IF(TYPINT.EQ.'QED'.AND.TYPEPR.EQ.'EXC') THEN
         HWHBOX=0.25*(CM1111*CONJG(CM1111)+CM1122*CONJG(CM1122)+
     &        CM2211*CONJG(CM2211)+CM1212*CONJG(CM1212)+
     &        CM2222*CONJG(CM2222))
      ELSEIF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'EXC') THEN
c...gluon + gluon --> gamma + gamma exclusive
         CM=CM1111+CM1122+CM2211+CM2222
         HWHBOX=0.25*CM*CONJG(CM)
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
      COMPLEX*16 FUNCTION CA1(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,TEMP1,TEMP2,TEMP3,REA1,IMA1,PI
      PARAMETER(PI=3.141592654D0)
      TEMP1=4.*(X-Y)/(X+Y)
      TEMP2=2.*(X*X+Y*Y)/((X+Y)**2)
      TEMP3=DLOG(-X/Y)
      REA1=8./3.+TEMP1*TEMP3-TEMP2*(TEMP3**2)
      IMA1=-PI*(TEMP1+2.*TEMP2*TEMP3)
      CA1=DCMPLX(REA1,IMA1)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION A2(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,TEMP1,TEMP2,TEMP3,REA2,IMA2,PI
      PARAMETER(PI=3.141592654D0)
      TEMP1=4.*(X-Y)/(X+Y)
      TEMP2=2.*(X*X+Y*Y)/(X+Y)**2
      TEMP3=DLOG(X/Y)
      A2=8./3.+TEMP1*TEMP3-TEMP2*(TEMP3**2+PI**2)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE HWINCH
C-----------------------------------------------------------------------
C     CHI PRODUCTION VIA GLUON OR FUSION
C     MEAN EVWGT = CHI PRODUCTION C-S * BRANCHING FRACTION IN NB
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION HWUALF,HWHIGT,HWRGEN,HWUSQR,HWUAEM,BRHIGQ,EMH,
     & CSFAC,EMFAC,CV,CA,RWGT,E1,E2,EMQ,HWIDTH
      INTEGER IDEC,I,J,ID1,ID2,ICHI,INDEX
      EXTERNAL HWUALF,HWHIGT,HWRGEN,HWUSQR,HWUAEM
      SAVE CSFAC
C ... Decide on which chi to produce
      ICHI=MOD(IPROC,100)
      IF(ICHI.EQ.1) THEN
         INDEX=165
         HWIDTH=0.0106 !! Measured width : 10.6 MeV (PDG 2004)
      ELSEIF(ICHI.EQ.2) THEN
         INDEX=306
         HWIDTH=0.0016 !! Predicted width, central value
                       !! Eichten, Quigg, PR D52 1995
      ELSE
         PRINT*, 'In HWINCH : WRONG IPROC ', IPROC
         STOP
      ENDIF
      IF (GENEV) THEN
C ... M.B. : this is ugly : gluglu -> chi recorded as gamgam -> chi, to prevent
C ...        trouble with the shower algorithm, which lacks a bit of phase space in this process.
C ...        This is corrected later on
        IDN(1)=59
        IDN(2)=59
        IDCMF=INDEX
        CALL HWEONE
      ELSE
        EVWGT=0.
        EMH=RMASS(INDEX)
        EMFAC=1.D0
        IF (EMH.LE.0 .OR. EMH.GE.PHEP(5,3)) RETURN
        EMSCA=EMH
        EMLST=EMH
        XXMIN=(EMH/PHEP(5,3))**2
        XLMIN=LOG(XXMIN)
        EMQ=0d0
        IF(TYPINT.EQ.'QCD'.AND.TYPEPR.EQ.'INC') THEN
           CSFAC=-GEV2NB*2.D0*(PI**2)*HWIDTH/(RMASS(INDEX)**3)*XLMIN
        ELSE
           PRINT*, 'In HWINCH : wrong TYP* settings, ',TYPINT,TYPEPR, 
     &          ' - STOP'
           STOP
        ENDIF
        CALL HWSGEN_CHI(.TRUE.)
        EVWGT=0
        E1=PHEP(4,MAX(1,JDAHEP(1,1)))
        E2=PHEP(4,MAX(2,JDAHEP(1,2)))
        EMQ=0d0
        IF (EMH.GT.2*EMQ) THEN
           IF (XX(1).LT.0.5*(1-EMQ/E1+HWUSQR(1-2*EMQ/E1)) .AND.
     &          XX(2).LT.0.5*(1-EMQ/E2+HWUSQR(1-2*EMQ/E2)))
     &          EVWGT=EVWGT+DISF(13,1)*DISF(13,2)*CSFAC
        ENDIF
 40     CONTINUE
      ENDIF
  999 END
C-----------------------------------------------------------------------
C * 07/03/2003 Tibor Kucs                                             
C * normalizations of PDFs inside the pomeron                                
C----------------------------------------------------------------------- 
      SUBROUTINE NRMPDF(XPQNRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(1),XPQNRM(-6:6)     
      COMMON /NORM/ Q2,NFLAVR
      COMMON /IFLAV/ IFLAVR
      COMMON /RESULT/ S1,S2,S3,S4
      EXTERNAL FUN
c---Initialize the VEGAS parameters:
      ACC=5D-3
      NDIM=1
      NCALL=1000
      ITS=4
c---Loop through PDFs    
      DO 100 IFLAVR=-NFLAVR,NFLAVR 
      CALL VEGAS(FUN,ACC,NDIM,NCALL,ITS,0,0)
c---Norm for parton with flavor IFLAVR
      XPQNRM(IFLAVR)=S1
 100  CONTINUE
      RETURN
      END
C----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUN(Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Z(1),XPQ(-6:6)
      COMMON /NORM/ Q2,NFLAVR
      COMMON /IFLAV/ IFLAVR
c------ CALL TO DIFFERENT QCD FITS (CHR)
      CALL QCDFIT(Z(1),Q2,XPQ,0)
c      CALL QCD_1994(Z(1),Q2,XPQ,0)
c---PDF for a parton of type IFLAVR inside the pomeron
      FUN=XPQ(IFLAVR)
      RETURN
      END
C----------------------------------------------------------------------------- 
      SUBROUTINE NRMPRT(XPQNRM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XPQNRM(-6:6)
      COMMON /NORM/ Q2,NFLAVR
      PRINT*, ' '
      PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
      PRINT*, ' '
      PRINT*, '  Normalization of PDFs inside the pomeron:'
      DO 200 IFLAVR=-NFLAVR,NFLAVR
         PRINT*,'    XPQNRM(',IFLAVR,') = ',XPQNRM(IFLAVR)
 200  CONTINUE
      PRINT*, ' '
      PRINT*, ' - - - - - - - - - FPMC - - - - - - - - - '
      PRINT*, ' '
      RETURN
      END

CDECK  ID>, HWSFUN.
*CMZ :-        -02/05/91  11.30.51  by  Federico Carminati
*-- Author :    Miscellaneous, combined by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWSFUN(X,MUSCALE,IDHAD,NSET,DIST,IBEAM)
C-----------------------------------------------------------------------
C     NUCLEON AND PION STRUCTURE FUNCTIONS DIST=X*QRK(X,Q=SCALE)
C
C     IDHAD = TYPE OF HADRON:
C     73=P  91=PBAR  75=N  93=NBAR  38=PI+  30=PI-  59=PHOTON
C
C     NEW SPECIAL CODES:
C     71=`REMNANT PHOTON' 72=`REMNANT NUCLEON'
C
C     NSET = STRUCTURE FUNCTION SET
C          = 1,2 FOR DUKE+OWENS SETS 1,2 (SOFT/HARD GLUE)
C          = 3,4 FOR EICHTEN ET AL SETS 1,2 (NUCLEON ONLY)
C          = 5   FOR OWENS SET 1.1 (PREPRINT FSU-HEP-910606)
C
C     FOR PHOTON DREES+GRASSIE IS USED
C
C     N.B. IF IBEAM.GT.0.AND.MODPDF(IBEAM).GE.0 THEN NSET IS
C     IGNORED AND CERN PDFLIB WITH AUTHOR GROUP=AUTPDF(IBEAM) AND
C     SET=MODPDF(IBEAM) IS USED.  FOR COMPATABILITY WITH VERSIONS 3
C     AND EARLIER, AUTPDF SHOULD BE SET TO 'MODE'
C     NOTE THAT NO CONSISTENCY CHECK IS MADE, FOR EXAMPLE THAT THE
C     REQUESTED SET FOR A PHOTON IS ACTUALLY A PHOTON SET
C
C     IF (ISPAC.GT.0) SCALE IS REPLACED BY MAX(SCALE,QSPAC)
C
C     FOR PHOTON, IF (PHOMAS.GT.0) THEN QUARK DISTRIBUTIONS ARE
C     SUPPRESSED BY      LOG((Q**2+PHOMAS**2)/(P**2+PHOMAS**2))
C                    L = -------------------------------------- ,
C                        LOG((Q**2+PHOMAS**2)/(     PHOMAS**2))
C     WHILE GLUON DISTRIBUTIONS ARE SUPPRESSED BY L**2,
C     WHERE Q=SCALE AND P=VIRTUALITY OF THE PHOTON
C
C   DUKE+OWENS = D.W.DUKE AND J.F.OWENS, PHYS. REV. D30 (1984) 49 (P/N)
C              + J.F.OWENS, PHYS. REV. D30 (1984) 943 (PI+/-)
C   WITH EXTRA SIGNIFICANT FIGURES VIA ED BERGER
C   WARNING....MOMENTUM SUM RULE BADLY VIOLATED ABOVE 1 TEV
C   DUKE+OWENS SETS 1,2 OBSOLETE. SET 1 UPDATED TO OWENS 1.1 (1991)
C   PION NOT RELIABLE ABOVE SCALE = 50 GEV
C
C   EICHTEN ET AL = E.EICHTEN,I.HINCHLIFFE,K.LANE AND C.QUIGG,
C                   REV. MOD. PHYS. 56 (1984) 579
C   REVISED AS IN   REV. MOD. PHYS. 58 (1986) 1065
C   RELIABLE RANGE : SQRT(5)GEV < SCALE < 10TEV, 1E-4 < X < 1
C
C   DREES+GRASSIE = M.DREES & K.GRASSIE, ZEIT. PHYS. C28 (1985) 451
C   MODIFIED IN     M.DREES & C.S.KIM, DESY 91-039
C                         AND C.S.KIM, DTP/91/16   FOR HEAVY QUARKS
C
C   FOR CERN PDFLIB DETAILS SEE PDFLIB DOC Q ON CERNVM OR
C   CERN_ROOT:[DOC]PDFLIB.TXT ON VXCERN
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      DOUBLE PRECISION HWSGAM,X,MUSCALE,XOLD,QOLD,XMWN,QSCA,SS,SMIN,S,T,
     & TMIN,TMAX,VX,AA,VT,WT,UPV,DNV,SEA,STR,CHM,BTM,TOP,GLU,WX,XQSUM,
     & DMIN,TPMIN,TPMAX,DIST(13),G(2),Q0(5),QL(5),F(5),A(6,5),
     & B(3,6,5,4),XQ(6),TX(6),TT(6),TB(6),NEHLQ(8,2),CEHLQ(6,6,2,8,2),
     & BB(4,6,5),VAL(20),USEA,DSEA,TOTAL,SCALEF,FAC,TBMIN(2),TTMIN(2)
      REAL HWSDGG,HWSDGQ,XSP,Q2,P2,W2,EMB2,EMC2,ALAM2,XPGA(-6:6),F2GM,
     & XPVMD,XPANL,XPANH,XPBEH,XPDIR
      COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &     XPDIR(-6:6)
      LOGICAL PDFWRX(2,2),PDFWRQ(2,2)
      DOUBLE PRECISION PDFXMN,PDFXMX,PDFQMN,PDFQMX
      COMMON /W50513/PDFXMN,PDFXMX,PDFQMN,PDFQMX
      INTEGER IDHAD,NSET,IBEAM,IOLD,NOLD,IP,I,J,K,NX,IT,IX,IFL,NFL,
     & MPDF,IHAD,ISET,IOP1,IOP2,IP2
      CHARACTER*20 PARM(20)
      EXTERNAL HWSGAM,HWSDGG,HWSDGQ


* B.C. Pomwig
      DOUBLE PRECISION XPQ(-6:6),BCQ, XXX
      INTEGER ifit,iloop
      PARAMETER(ifit=0)
      CHARACTER*20 PARMA(20)
      DOUBLE PRECISION valu(20)

c init pdflib
      logical first
      data first/.true./

* transmit proton Id in case of Pom+Reg
      integer ind
      
      common /yura/ind

* B.C. Pomwig

      SAVE QOLD,IOLD,NOLD,XOLD,SS,S,T,TMIN,TMAX,G,A,TX,TT,TB,IP,NX
      DATA PDFWRX,PDFWRQ/8*.TRUE./
      DATA (((B(I,J,K,1),I=1,3),J=1,6),K=1,5)/
     &3.D0,0.D0,0.D0,.419D0,.004383D0,-.007412D0,
     &3.46D0,.72432D0,-.065998D0,4.4D0,-4.8644D0,1.3274D0,
     &6*0.D0,1.D0,
     &0.D0,0.D0,.763D0,-.23696D0,.025836D0,4.D0,.62664D0,-.019163D0,
     &0.D0,-.42068D0,.032809D0,6*0.D0,1.265D0,-1.1323D0,.29268D0,
     &0.D0,-.37162D0,-.028977D0,8.05D0,1.5877D0,-.15291D0,
     &0.D0,6.3059D0,-.27342D0,0.D0,-10.543D0,-3.1674D0,
     &0.D0,14.698D0,9.798D0,0.D0,.13479D0,-.074693D0,
     &-.0355D0,-.22237D0,-.057685D0,6.3494D0,3.2649D0,-.90945D0,
     &0.D0,-3.0331D0,1.5042D0,0.D0,17.431D0,-11.255D0,
     &0.D0,-17.861D0,15.571D0,1.564D0,-1.7112D0,.63751D0,
     &0.D0,-.94892D0,.32505D0,6.D0,1.4345D0,-1.0485D0,
     &9.D0,-7.1858D0,.25494D0,0.D0,-16.457D0,10.947D0,
     &0.D0,15.261D0,-10.085D0/
      DATA (((B(I,J,K,2),I=1,3),J=1,6),K=1,5)/
     &3.D0,0.D0,0.D0,.3743D0,.013946D0,-.00031695D0,
     &3.329D0,.75343D0,-.076125D0,6.032D0,-6.2153D0,1.5561D0,
     &6*0.D0,1.D0,0.D0,
     &0.D0,.7608D0,-.2317D0,.023232D0,3.83D0,.62746D0,-.019155D0,
     &0.D0,-.41843D0,.035972D0,6*0.D0,1.6714D0,-1.9168D0,.58175D0,
     &0.D0,-.27307D0,-.16392D0,9.145D0,.53045D0,-.76271D0,
     &0.D0,15.665D0,-2.8341D0,0.D0,-100.63D0,44.658D0,
     &0.D0,223.24D0,-116.76D0,0.D0,.067368D0,-.030574D0,
     &-.11989D0,-.23293D0,-.023273D0,3.5087D0,3.6554D0,-.45313D0,
     &0.D0,-.47369D0,.35793D0,0.D0,9.5041D0,-5.4303D0,
     &0.D0,-16.563D0,15.524D0,.8789D0,-.97093D0,.43388D0,
     &0.D0,-1.1612D0,.4759D0,4.D0,1.2271D0,-.25369D0,
     &9.D0,-5.6354D0,-.81747D0,0.D0,-7.5438D0,5.5034D0,
     &0.D0,-.59649D0,.12611D0/
      DATA (((B(I,J,K,3),I=1,3),J=1,6),K=1,5)/
     &1.D0,0.D0,0.D0,0.4D0,-0.06212D0,-0.007109D0,0.7D0,0.6478D0,
     &0.01335D0,27*0.D0,0.9D0,-0.2428D0,0.1386D0,0.D0,-0.2120D0,
     &0.003671D0,5.0D0,0.8673D0,0.04747D0,
     &0.D0,1.266D0,-2.215D0,0.D0,2.382D0,0.3482D0,3*0.D0,
     &0.D0,0.07928D0,-0.06134D0,-0.02212D0,-0.3785D0,-0.1088D0,2.894D0,
     &9.433D0,
     &-10.852D0,0.D0,5.248D0,-7.187D0,0.D0,8.388D0,-11.61D0,3*0.D0,
     &0.888D0,-1.802D0,1.812D0,0.D0,-1.576D0,1.20D0,3.11D0,-0.1317D0,
     &0.5068D0,6.0D0,2.801D0,-12.16D0,0.D0,-17.28D0,20.49D0,3*0.D0/
      DATA (((B(I,J,K,4),I=1,3),J=1,6),K=1,5)/
     &1.D0,0.D0,0.D0,0.4D0,-0.05909D0,-0.006524D0,0.628D0,0.6436D0,
     &0.01451D0,27*0.D0,
     &0.90D0,-0.1417D0,-0.1740D0,0.D0,-0.1697D0,-0.09623D0,5.0D0,
     &-2.474D0,1.575D0,
     &0.D0,-2.534D0,1.378D0,0.D0,0.5621D0,-0.2701D0,3*0.D0,
     &0.D0,0.06229D0,-0.04099D0,-0.0882D0,-0.2892D0,-0.1082D0,1.924D0,
     &0.2424D0,
     &2.036D0,0.D0,-4.463D0,5.209D0,0.D0,-0.8367D0,-0.04840D0,3*0.D0,
     &0.794D0,-0.9144D0,0.5966D0,0.D0,-1.237D0,0.6582D0,2.89D0,0.5966D0,
     &-0.2550D0,
     &6.0D0,-3.671D0,-2.304D0,0.D0,-8.191D0,7.758D0,3*0.D0/
C---COEFFTS FOR NEW OWENS 1.1 SET
      DATA BB/3.D0,3*0.D0,.665D0,-.1097D0,-.002442D0,0.D0,
     &3.614D0,.8395D0,-.02186D0,0.D0,.8673D0,-1.6637D0,.342D0,0.D0,
     &0.D0,1.1049D0,-.2369D0,5*0.D0,1.D0,3*0.D0,
     &.8388D0,-.2092D0,.02657D0,0.D0,4.667D0,.7951D0,.1081D0,0.D0,
     &0.D0,-1.0232D0,.05799D0,0.D0,0.D0,.8616D0,.153D0,5*0.D0,
     &.909D0,-.4023D0,.006305D0,0.D0,
     &0.D0,-.3823D0,.02766D0,0.D0,7.278D0,-.7904D0,.8108D0,0.D0,
     &0.D0,-1.6629D0,.5719D0,0.D0,0.D0,-.01333D0,.5299D0,0.D0,
     &0.D0,.1211D0,-.1739D0,0.D0,0.D0,.09469D0,-.07066D0,.01236D0,
     &-.1447D0,-.402D0,.1533D0,-.06479D0,6.7599D0,1.6596D0,.6798D0,
     &-.8525D0,0.D0,-4.4559D0,3.3756D0,-.9468D0,
     &0.D0,7.862D0,-3.6591D0,.03672D0,0.D0,-.2472D0,-.751D0,.0487D0,
     &3.017D0,-4.7347D0,3.3594D0,-.9443D0,0.D0,-.9342D0,.5454D0,
     &-.1668D0,
     &5.304D0,1.4654D0,-1.4292D0,.7569D0,0.D0,-3.9141D0,2.8445D0,
     &-.8411D0,
     &0.D0,9.0176D0,-10.426D0,4.0983D0,0.D0,-5.9602D0,7.515D0,-2.7329D0/
C...THE FOLLOWING DATA LINES ARE COEFFICIENTS NEEDED IN THE
C...EICHTEN, HINCHLIFFE, LANE, QUIGG PROTON STRUCTURE FUNCTION
C...POWERS OF 1-X IN DIFFERENT CASES
      DATA NEHLQ/3,4,7,5,7,7,7,7,3,4,7,6,7,7,7,7/
C...EXPANSION COEFFICIENTS FOR UP VALENCE QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,1,1),IX=1,6),IT=1,6),NX=1,2)/
     1 7.677D-01,-2.087D-01,-3.303D-01,-2.517D-02,-1.570D-02,-1.000D-04,
     2-5.326D-01,-2.661D-01, 3.201D-01, 1.192D-01, 2.434D-02, 7.620D-03,
     3 2.162D-01, 1.881D-01,-8.375D-02,-6.515D-02,-1.743D-02,-5.040D-03,
     4-9.211D-02,-9.952D-02, 1.373D-02, 2.506D-02, 8.770D-03, 2.550D-03,
     5 3.670D-02, 4.409D-02, 9.600D-04,-7.960D-03,-3.420D-03,-1.050D-03,
     6-1.549D-02,-2.026D-02,-3.060D-03, 2.220D-03, 1.240D-03, 4.100D-04,
     1 2.395D-01, 2.905D-01, 9.778D-02, 2.149D-02, 3.440D-03, 5.000D-04,
     2 1.751D-02,-6.090D-03,-2.687D-02,-1.916D-02,-7.970D-03,-2.750D-03,
     3-5.760D-03,-5.040D-03, 1.080D-03, 2.490D-03, 1.530D-03, 7.500D-04,
     4 1.740D-03, 1.960D-03, 3.000D-04,-3.400D-04,-2.900D-04,-1.800D-04,
     5-5.300D-04,-6.400D-04,-1.700D-04, 4.000D-05, 6.000D-05, 4.000D-05,
     6 1.700D-04, 2.200D-04, 8.000D-05, 1.000D-05,-1.000D-05,-1.000D-05/
      DATA (((CEHLQ(IX,IT,NX,1,2),IX=1,6),IT=1,6),NX=1,2)/
     1 7.237D-01,-2.189D-01,-2.995D-01,-1.909D-02,-1.477D-02, 2.500D-04,
     2-5.314D-01,-2.425D-01, 3.283D-01, 1.119D-01, 2.223D-02, 7.070D-03,
     3 2.289D-01, 1.890D-01,-9.859D-02,-6.900D-02,-1.747D-02,-5.080D-03,
     4-1.041D-01,-1.084D-01, 2.108D-02, 2.975D-02, 9.830D-03, 2.830D-03,
     5 4.394D-02, 5.116D-02,-1.410D-03,-1.055D-02,-4.230D-03,-1.270D-03,
     6-1.991D-02,-2.539D-02,-2.780D-03, 3.430D-03, 1.720D-03, 5.500D-04,
     1 2.410D-01, 2.884D-01, 9.369D-02, 1.900D-02, 2.530D-03, 2.400D-04,
     2 1.765D-02,-9.220D-03,-3.037D-02,-2.085D-02,-8.440D-03,-2.810D-03,
     3-6.450D-03,-5.260D-03, 1.720D-03, 3.110D-03, 1.830D-03, 8.700D-04,
     4 2.120D-03, 2.320D-03, 2.600D-04,-4.900D-04,-3.900D-04,-2.300D-04,
     5-6.900D-04,-8.200D-04,-2.000D-04, 7.000D-05, 9.000D-05, 6.000D-05,
     6 2.400D-04, 3.100D-04, 1.100D-04, 0.000D+00,-2.000D-05,-2.000D-05/
C...EXPANSION COEFFICIENTS FOR DOWN VALENCE QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,2,1),IX=1,6),IT=1,6),NX=1,2)/
     1 3.813D-01,-8.090D-02,-1.634D-01,-2.185D-02,-8.430D-03,-6.200D-04,
     2-2.948D-01,-1.435D-01, 1.665D-01, 6.638D-02, 1.473D-02, 4.080D-03,
     3 1.252D-01, 1.042D-01,-4.722D-02,-3.683D-02,-1.038D-02,-2.860D-03,
     4-5.478D-02,-5.678D-02, 8.900D-03, 1.484D-02, 5.340D-03, 1.520D-03,
     5 2.220D-02, 2.567D-02,-3.000D-05,-4.970D-03,-2.160D-03,-6.500D-04,
     6-9.530D-03,-1.204D-02,-1.510D-03, 1.510D-03, 8.300D-04, 2.700D-04,
     1 1.261D-01, 1.354D-01, 3.958D-02, 8.240D-03, 1.660D-03, 4.500D-04,
     2 3.890D-03,-1.159D-02,-1.625D-02,-9.610D-03,-3.710D-03,-1.260D-03,
     3-1.910D-03,-5.600D-04, 1.590D-03, 1.590D-03, 8.400D-04, 3.900D-04,
     4 6.400D-04, 4.900D-04,-1.500D-04,-2.900D-04,-1.800D-04,-1.000D-04,
     5-2.000D-04,-1.900D-04, 0.000D+00, 6.000D-05, 4.000D-05, 3.000D-05,
     6 7.000D-05, 8.000D-05, 2.000D-05,-1.000D-05,-1.000D-05,-1.000D-05/
      DATA (((CEHLQ(IX,IT,NX,2,2),IX=1,6),IT=1,6),NX=1,2)/
     1 3.578D-01,-8.622D-02,-1.480D-01,-1.840D-02,-7.820D-03,-4.500D-04,
     2-2.925D-01,-1.304D-01, 1.696D-01, 6.243D-02, 1.353D-02, 3.750D-03,
     3 1.318D-01, 1.041D-01,-5.486D-02,-3.872D-02,-1.038D-02,-2.850D-03,
     4-6.162D-02,-6.143D-02, 1.303D-02, 1.740D-02, 5.940D-03, 1.670D-03,
     5 2.643D-02, 2.957D-02,-1.490D-03,-6.450D-03,-2.630D-03,-7.700D-04,
     6-1.218D-02,-1.497D-02,-1.260D-03, 2.240D-03, 1.120D-03, 3.500D-04,
     1 1.263D-01, 1.334D-01, 3.732D-02, 7.070D-03, 1.260D-03, 3.400D-04,
     2 3.660D-03,-1.357D-02,-1.795D-02,-1.031D-02,-3.880D-03,-1.280D-03,
     3-2.100D-03,-3.600D-04, 2.050D-03, 1.920D-03, 9.800D-04, 4.400D-04,
     4 7.700D-04, 5.400D-04,-2.400D-04,-3.900D-04,-2.400D-04,-1.300D-04,
     5-2.600D-04,-2.300D-04, 2.000D-05, 9.000D-05, 6.000D-05, 4.000D-05,
     6 9.000D-05, 1.000D-04, 2.000D-05,-2.000D-05,-2.000D-05,-1.000D-05/
C...EXPANSION COEFFICIENTS FOR UP AND DOWN SEA QUARK DISTRIBUTIONS
      DATA (((CEHLQ(IX,IT,NX,3,1),IX=1,6),IT=1,6),NX=1,2)/
     1 6.870D-02,-6.861D-02, 2.973D-02,-5.400D-03, 3.780D-03,-9.700D-04,
     2-1.802D-02, 1.400D-04, 6.490D-03,-8.540D-03, 1.220D-03,-1.750D-03,
     3-4.650D-03, 1.480D-03,-5.930D-03, 6.000D-04,-1.030D-03,-8.000D-05,
     4 6.440D-03, 2.570D-03, 2.830D-03, 1.150D-03, 7.100D-04, 3.300D-04,
     5-3.930D-03,-2.540D-03,-1.160D-03,-7.700D-04,-3.600D-04,-1.900D-04,
     6 2.340D-03, 1.930D-03, 5.300D-04, 3.700D-04, 1.600D-04, 9.000D-05,
     1 1.014D+00,-1.106D+00, 3.374D-01,-7.444D-02, 8.850D-03,-8.700D-04,
     2 9.233D-01,-1.285D+00, 4.475D-01,-9.786D-02, 1.419D-02,-1.120D-03,
     3 4.888D-02,-1.271D-01, 8.606D-02,-2.608D-02, 4.780D-03,-6.000D-04,
     4-2.691D-02, 4.887D-02,-1.771D-02, 1.620D-03, 2.500D-04,-6.000D-05,
     5 7.040D-03,-1.113D-02, 1.590D-03, 7.000D-04,-2.000D-04, 0.000D+00,
     6-1.710D-03, 2.290D-03, 3.800D-04,-3.500D-04, 4.000D-05, 1.000D-05/
      DATA (((CEHLQ(IX,IT,NX,3,2),IX=1,6),IT=1,6),NX=1,2)/
     1 1.008D-01,-7.100D-02, 1.973D-02,-5.710D-03, 2.930D-03,-9.900D-04,
     2-5.271D-02,-1.823D-02, 1.792D-02,-6.580D-03, 1.750D-03,-1.550D-03,
     3 1.220D-02, 1.763D-02,-8.690D-03,-8.800D-04,-1.160D-03,-2.100D-04,
     4-1.190D-03,-7.180D-03, 2.360D-03, 1.890D-03, 7.700D-04, 4.100D-04,
     5-9.100D-04, 2.040D-03,-3.100D-04,-1.050D-03,-4.000D-04,-2.400D-04,
     6 1.190D-03,-1.700D-04,-2.000D-04, 4.200D-04, 1.700D-04, 1.000D-04,
     1 1.081D+00,-1.189D+00, 3.868D-01,-8.617D-02, 1.115D-02,-1.180D-03,
     2 9.917D-01,-1.396D+00, 4.998D-01,-1.159D-01, 1.674D-02,-1.720D-03,
     3 5.099D-02,-1.338D-01, 9.173D-02,-2.885D-02, 5.890D-03,-6.500D-04,
     4-3.178D-02, 5.703D-02,-2.070D-02, 2.440D-03, 1.100D-04,-9.000D-05,
     5 8.970D-03,-1.392D-02, 2.050D-03, 6.500D-04,-2.300D-04, 2.000D-05,
     6-2.340D-03, 3.010D-03, 5.000D-04,-3.900D-04, 6.000D-05, 1.000D-05/
C...EXPANSION COEFFICIENTS FOR GLUON DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,4,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.482D-01,-9.578D-01, 1.009D-01,-1.051D-01, 3.456D-02,-3.054D-02,
     2-9.627D-01, 5.379D-01, 3.368D-01,-9.525D-02, 1.488D-02,-2.051D-02,
     3 4.300D-01,-8.306D-02,-3.372D-01, 4.902D-02,-9.160D-03, 1.041D-02,
     4-1.925D-01,-1.790D-02, 2.183D-01, 7.490D-03, 4.140D-03,-1.860D-03,
     5 8.183D-02, 1.926D-02,-1.072D-01,-1.944D-02,-2.770D-03,-5.200D-04,
     6-3.884D-02,-1.234D-02, 5.410D-02, 1.879D-02, 3.350D-03, 1.040D-03,
     1 2.948D+01,-3.902D+01, 1.464D+01,-3.335D+00, 5.054D-01,-5.915D-02,
     2 2.559D+01,-3.955D+01, 1.661D+01,-4.299D+00, 6.904D-01,-8.243D-02,
     3-1.663D+00, 1.176D+00, 1.118D+00,-7.099D-01, 1.948D-01,-2.404D-02,
     4-2.168D-01, 8.170D-01,-7.169D-01, 1.851D-01,-1.924D-02,-3.250D-03,
     5 2.088D-01,-4.355D-01, 2.239D-01,-2.446D-02,-3.620D-03, 1.910D-03,
     6-9.097D-02, 1.601D-01,-5.681D-02,-2.500D-03, 2.580D-03,-4.700D-04/
      DATA (((CEHLQ(IX,IT,NX,4,2),IX=1,6),IT=1,6),NX=1,2)/
     1 2.367D+00, 4.453D-01, 3.660D-01, 9.467D-02, 1.341D-01, 1.661D-02,
     2-3.170D+00,-1.795D+00, 3.313D-02,-2.874D-01,-9.827D-02,-7.119D-02,
     3 1.823D+00, 1.457D+00,-2.465D-01, 3.739D-02, 6.090D-03, 1.814D-02,
     4-1.033D+00,-9.827D-01, 2.136D-01, 1.169D-01, 5.001D-02, 1.684D-02,
     5 5.133D-01, 5.259D-01,-1.173D-01,-1.139D-01,-4.988D-02,-2.021D-02,
     6-2.881D-01,-3.145D-01, 5.667D-02, 9.161D-02, 4.568D-02, 1.951D-02,
     1 3.036D+01,-4.062D+01, 1.578D+01,-3.699D+00, 6.020D-01,-7.031D-02,
     2 2.700D+01,-4.167D+01, 1.770D+01,-4.804D+00, 7.862D-01,-1.060D-01,
     3-1.909D+00, 1.357D+00, 1.127D+00,-7.181D-01, 2.232D-01,-2.481D-02,
     4-2.488D-01, 9.781D-01,-8.127D-01, 2.094D-01,-2.997D-02,-4.710D-03,
     5 2.506D-01,-5.427D-01, 2.672D-01,-3.103D-02,-1.800D-03, 2.870D-03,
     6-1.128D-01, 2.087D-01,-6.972D-02,-2.480D-03, 2.630D-03,-8.400D-04/
C...EXPANSION COEFFICIENTS FOR STRANGE SEA QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,5,1),IX=1,6),IT=1,6),NX=1,2)/
     1 4.968D-02,-4.173D-02, 2.102D-02,-3.270D-03, 3.240D-03,-6.700D-04,
     2-6.150D-03,-1.294D-02, 6.740D-03,-6.890D-03, 9.000D-04,-1.510D-03,
     3-8.580D-03, 5.050D-03,-4.900D-03,-1.600D-04,-9.400D-04,-1.500D-04,
     4 7.840D-03, 1.510D-03, 2.220D-03, 1.400D-03, 7.000D-04, 3.500D-04,
     5-4.410D-03,-2.220D-03,-8.900D-04,-8.500D-04,-3.600D-04,-2.000D-04,
     6 2.520D-03, 1.840D-03, 4.100D-04, 3.900D-04, 1.600D-04, 9.000D-05,
     1 9.235D-01,-1.085D+00, 3.464D-01,-7.210D-02, 9.140D-03,-9.100D-04,
     2 9.315D-01,-1.274D+00, 4.512D-01,-9.775D-02, 1.380D-02,-1.310D-03,
     3 4.739D-02,-1.296D-01, 8.482D-02,-2.642D-02, 4.760D-03,-5.700D-04,
     4-2.653D-02, 4.953D-02,-1.735D-02, 1.750D-03, 2.800D-04,-6.000D-05,
     5 6.940D-03,-1.132D-02, 1.480D-03, 6.500D-04,-2.100D-04, 0.000D+00,
     6-1.680D-03, 2.340D-03, 4.200D-04,-3.400D-04, 5.000D-05, 1.000D-05/
      DATA (((CEHLQ(IX,IT,NX,5,2),IX=1,6),IT=1,6),NX=1,2)/
     1 6.478D-02,-4.537D-02, 1.643D-02,-3.490D-03, 2.710D-03,-6.700D-04,
     2-2.223D-02,-2.126D-02, 1.247D-02,-6.290D-03, 1.120D-03,-1.440D-03,
     3-1.340D-03, 1.362D-02,-6.130D-03,-7.900D-04,-9.000D-04,-2.000D-04,
     4 5.080D-03,-3.610D-03, 1.700D-03, 1.830D-03, 6.800D-04, 4.000D-04,
     5-3.580D-03, 6.000D-05,-2.600D-04,-1.050D-03,-3.800D-04,-2.300D-04,
     6 2.420D-03, 9.300D-04,-1.000D-04, 4.500D-04, 1.700D-04, 1.100D-04,
     1 9.868D-01,-1.171D+00, 3.940D-01,-8.459D-02, 1.124D-02,-1.250D-03,
     2 1.001D+00,-1.383D+00, 5.044D-01,-1.152D-01, 1.658D-02,-1.830D-03,
     3 4.928D-02,-1.368D-01, 9.021D-02,-2.935D-02, 5.800D-03,-6.600D-04,
     4-3.133D-02, 5.785D-02,-2.023D-02, 2.630D-03, 1.600D-04,-8.000D-05,
     5 8.840D-03,-1.416D-02, 1.900D-03, 5.800D-04,-2.500D-04, 1.000D-05,
     6-2.300D-03, 3.080D-03, 5.500D-04,-3.700D-04, 7.000D-05, 1.000D-05/
C...EXPANSION COEFFICIENTS FOR CHARM SEA QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,6,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.270D-03,-1.817D-02, 9.590D-03,-6.390D-03, 1.690D-03,-1.540D-03,
     2 5.710D-03,-1.188D-02, 6.090D-03,-4.650D-03, 1.240D-03,-1.310D-03,
     3-3.960D-03, 7.100D-03,-3.590D-03, 1.840D-03,-3.900D-04, 3.400D-04,
     4 1.120D-03,-1.960D-03, 1.120D-03,-4.800D-04, 1.000D-04,-4.000D-05,
     5 4.000D-05,-3.000D-05,-1.800D-04, 9.000D-05,-5.000D-05,-2.000D-05,
     6-4.200D-04, 7.300D-04,-1.600D-04, 5.000D-05, 5.000D-05, 5.000D-05,
     1 8.098D-01,-1.042D+00, 3.398D-01,-6.824D-02, 8.760D-03,-9.000D-04,
     2 8.961D-01,-1.217D+00, 4.339D-01,-9.287D-02, 1.304D-02,-1.290D-03,
     3 3.058D-02,-1.040D-01, 7.604D-02,-2.415D-02, 4.600D-03,-5.000D-04,
     4-2.451D-02, 4.432D-02,-1.651D-02, 1.430D-03, 1.200D-04,-1.000D-04,
     5 1.122D-02,-1.457D-02, 2.680D-03, 5.800D-04,-1.200D-04, 3.000D-05,
     6-7.730D-03, 7.330D-03,-7.600D-04,-2.400D-04, 1.000D-05, 0.000D+00/
      DATA (((CEHLQ(IX,IT,NX,6,2),IX=1,6),IT=1,6),NX=1,2)/
     1 9.980D-03,-1.945D-02, 1.055D-02,-6.870D-03, 1.860D-03,-1.560D-03,
     2 5.700D-03,-1.203D-02, 6.250D-03,-4.860D-03, 1.310D-03,-1.370D-03,
     3-4.490D-03, 7.990D-03,-4.170D-03, 2.050D-03,-4.400D-04, 3.300D-04,
     4 1.470D-03,-2.480D-03, 1.460D-03,-5.700D-04, 1.200D-04,-1.000D-05,
     5-9.000D-05, 1.500D-04,-3.200D-04, 1.200D-04,-6.000D-05,-4.000D-05,
     6-4.200D-04, 7.600D-04,-1.400D-04, 4.000D-05, 7.000D-05, 5.000D-05,
     1 8.698D-01,-1.131D+00, 3.836D-01,-8.111D-02, 1.048D-02,-1.300D-03,
     2 9.626D-01,-1.321D+00, 4.854D-01,-1.091D-01, 1.583D-02,-1.700D-03,
     3 3.057D-02,-1.088D-01, 8.022D-02,-2.676D-02, 5.590D-03,-5.600D-04,
     4-2.845D-02, 5.164D-02,-1.918D-02, 2.210D-03,-4.000D-05,-1.500D-04,
     5 1.311D-02,-1.751D-02, 3.310D-03, 5.100D-04,-1.200D-04, 5.000D-05,
     6-8.590D-03, 8.380D-03,-9.200D-04,-2.600D-04, 1.000D-05,-1.000D-05/
C...EXPANSION COEFFICIENTS FOR BOTTOM SEA QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,7,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.010D-03,-1.401D-02, 7.150D-03,-4.130D-03, 1.260D-03,-1.040D-03,
     2 6.280D-03,-9.320D-03, 4.780D-03,-2.890D-03, 9.100D-04,-8.200D-04,
     3-2.930D-03, 4.090D-03,-1.890D-03, 7.600D-04,-2.300D-04, 1.400D-04,
     4 3.900D-04,-1.200D-03, 4.400D-04,-2.500D-04, 2.000D-05,-2.000D-05,
     5 2.600D-04, 1.400D-04,-8.000D-05, 1.000D-04, 1.000D-05, 1.000D-05,
     6-2.600D-04, 3.200D-04, 1.000D-05,-1.000D-05, 1.000D-05,-1.000D-05,
     1 8.029D-01,-1.075D+00, 3.792D-01,-7.843D-02, 1.007D-02,-1.090D-03,
     2 7.903D-01,-1.099D+00, 4.153D-01,-9.301D-02, 1.317D-02,-1.410D-03,
     3-1.704D-02,-1.130D-02, 2.882D-02,-1.341D-02, 3.040D-03,-3.600D-04,
     4-7.200D-04, 7.230D-03,-5.160D-03, 1.080D-03,-5.000D-05,-4.000D-05,
     5 3.050D-03,-4.610D-03, 1.660D-03,-1.300D-04,-1.000D-05, 1.000D-05,
     6-4.360D-03, 5.230D-03,-1.610D-03, 2.000D-04,-2.000D-05, 0.000D+00/
      DATA (((CEHLQ(IX,IT,NX,7,2),IX=1,6),IT=1,6),NX=1,2)/
     1 8.980D-03,-1.459D-02, 7.510D-03,-4.410D-03, 1.310D-03,-1.070D-03,
     2 5.970D-03,-9.440D-03, 4.800D-03,-3.020D-03, 9.100D-04,-8.500D-04,
     3-3.050D-03, 4.440D-03,-2.100D-03, 8.500D-04,-2.400D-04, 1.400D-04,
     4 5.300D-04,-1.300D-03, 5.600D-04,-2.700D-04, 3.000D-05,-2.000D-05,
     5 2.000D-04, 1.400D-04,-1.100D-04, 1.000D-04, 0.000D+00, 0.000D+00,
     6-2.600D-04, 3.200D-04, 0.000D+00,-3.000D-05, 1.000D-05,-1.000D-05,
     1 8.672D-01,-1.174D+00, 4.265D-01,-9.252D-02, 1.244D-02,-1.460D-03,
     2 8.500D-01,-1.194D+00, 4.630D-01,-1.083D-01, 1.614D-02,-1.830D-03,
     3-2.241D-02,-5.630D-03, 2.815D-02,-1.425D-02, 3.520D-03,-4.300D-04,
     4-7.300D-04, 8.030D-03,-5.780D-03, 1.380D-03,-1.300D-04,-4.000D-05,
     5 3.460D-03,-5.380D-03, 1.960D-03,-2.100D-04, 1.000D-05, 1.000D-05,
     6-4.850D-03, 5.950D-03,-1.890D-03, 2.600D-04,-3.000D-05, 0.000D+00/
C...EXPANSION COEFFICIENTS FOR TOP SEA QUARK DISTRIBUTION
      DATA (((CEHLQ(IX,IT,NX,8,1),IX=1,6),IT=1,6),NX=1,2)/
     1 4.410D-03,-7.480D-03, 3.770D-03,-2.580D-03, 7.300D-04,-7.100D-04,
     2 3.840D-03,-6.050D-03, 3.030D-03,-2.030D-03, 5.800D-04,-5.900D-04,
     3-8.800D-04, 1.660D-03,-7.500D-04, 4.700D-04,-1.000D-04, 1.000D-04,
     4-8.000D-05,-1.500D-04, 1.200D-04,-9.000D-05, 3.000D-05, 0.000D+00,
     5 1.300D-04,-2.200D-04,-2.000D-05,-2.000D-05,-2.000D-05,-2.000D-05,
     6-7.000D-05, 1.900D-04,-4.000D-05, 2.000D-05, 0.000D+00, 0.000D+00,
     1 6.623D-01,-9.248D-01, 3.519D-01,-7.930D-02, 1.110D-02,-1.180D-03,
     2 6.380D-01,-9.062D-01, 3.582D-01,-8.479D-02, 1.265D-02,-1.390D-03,
     3-2.581D-02, 2.125D-02, 4.190D-03,-4.980D-03, 1.490D-03,-2.100D-04,
     4 7.100D-04, 5.300D-04,-1.270D-03, 3.900D-04,-5.000D-05,-1.000D-05,
     5 3.850D-03,-5.060D-03, 1.860D-03,-3.500D-04, 4.000D-05, 0.000D+00,
     6-3.530D-03, 4.460D-03,-1.500D-03, 2.700D-04,-3.000D-05, 0.000D+00/
      DATA (((CEHLQ(IX,IT,NX,8,2),IX=1,6),IT=1,6),NX=1,2)/
     1 4.260D-03,-7.530D-03, 3.830D-03,-2.680D-03, 7.600D-04,-7.300D-04,
     2 3.640D-03,-6.050D-03, 3.030D-03,-2.090D-03, 5.900D-04,-6.000D-04,
     3-9.200D-04, 1.710D-03,-8.200D-04, 5.000D-04,-1.200D-04, 1.000D-04,
     4-5.000D-05,-1.600D-04, 1.300D-04,-9.000D-05, 3.000D-05, 0.000D+00,
     5 1.300D-04,-2.100D-04,-1.000D-05,-2.000D-05,-2.000D-05,-1.000D-05,
     6-8.000D-05, 1.800D-04,-5.000D-05, 2.000D-05, 0.000D+00, 0.000D+00,
     1 7.146D-01,-1.007D+00, 3.932D-01,-9.246D-02, 1.366D-02,-1.540D-03,
     2 6.856D-01,-9.828D-01, 3.977D-01,-9.795D-02, 1.540D-02,-1.790D-03,
     3-3.053D-02, 2.758D-02, 2.150D-03,-4.880D-03, 1.640D-03,-2.500D-04,
     4 9.200D-04, 4.200D-04,-1.340D-03, 4.600D-04,-8.000D-05,-1.000D-05,
     5 4.230D-03,-5.660D-03, 2.140D-03,-4.300D-04, 6.000D-05, 0.000D+00,
     6-3.890D-03, 5.000D-03,-1.740D-03, 3.300D-04,-4.000D-05, 0.000D+00/
      DATA TBMIN,TTMIN/8.1905D0,7.4474D0,11.5528D0,10.8097D0/
      DATA XOLD,QOLD,IOLD,NOLD/-1.D0,0.D0,0,0/
      DATA DMIN,Q0,QL/0.D0,2*2.D0,2*2.236D0,2.D0,.2D0,
     &                .4D0,.2D0,.29D0,.177D0/
      IF (X.LE.ZERO) CALL HWWARN('HWSFUN',100,*999)
      XMWN=ONE-X
      IF (XMWN.LE.ZERO) THEN
        DO 1 I=1,13
          DIST(I)=0
 1      CONTINUE
        RETURN
      ENDIF
C---FREEZE THE SCALE IF REQUIRED
      SCALEF=MUSCALE
      IF (ISPAC.GT.0) SCALEF=MAX(SCALEF,QSPAC)
C---CHECK IF PDFLIB REQUESTED
      IF (IBEAM.EQ.1.OR.IBEAM.EQ.2) THEN
        MPDF=MODPDF(IBEAM)
      ELSE
        MPDF=-1
      ENDIF
      QSCA=ABS(SCALEF)

* B.C. Pomwig Start of mod
C--- pomeron for photon
      IF (IDHAD.EQ.59) THEN
         IF (NFLUX.EQ.9) THEN
            BCQ=QSCA*QSCA
C Initialise xpq
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,ifit)
c CHR call to different QCD fits
            CALL QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)=XPQ(1)
            DIST(2)=XPQ(2)
            DIST(3)=XPQ(3)
            DIST(4)=XPQ(4)
c            DIST(4)=0
            DIST(5)=0
            DIST(6)=0
            DIST(7)=XPQ(-1)
            DIST(8)=XPQ(-2)
            DIST(9)=XPQ(-3)
            DIST(10)=XPQ(-4)
c            DIST(10)=0
            DIST(11)=0
            DIST(12)=0
            DIST(13)=XPQ(0)
            GOTO 999
         ELSEIF (NFLUX.EQ.10) THEN
C CHR YURA
            if(first) then
                PARMA(1)='init0'
                VALU(1)=0.D0
                call pdfset(PARMA, VALU)
                first=.false.
             endif
C     Reggeon 
            parma(1)='NPTYPE'
            parma(2)='NGROUP'
            parma(3)='NSET'
            valu(1)=2
*            valu(2)=5
*            valu(3)=2            ! GRV-P LO
            valu(2)=1
            valu(3)=1            ! OWENS

            CALL PDFSET(PARMA,VALU)
            CALL STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
            DIST(1)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(2)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(3)=Cr*STR
            DIST(4)=Cr*CHM
            DIST(5)=0
            DIST(6)=0
            DIST(7)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(8)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(9)=Cr*STR
            DIST(10)=Cr*CHM
            DIST(11)=0
            DIST(12)=0
            DIST(13)=Cr*GLU
            GOTO 999
c Yura CHR Pom Reg	    
	   ELSEIF (NFLUX.EQ.19) THEN
c	   print *,'index in 19 :',ind,ibeam
	   IF(IBEAM.EQ.1) THEN 
           BCQ=QSCA*QSCA
C Initialise xpq
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,ifit)
c CHR call to different QCD fits
            CALL QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)=XPQ(1)
            DIST(2)=XPQ(2)
            DIST(3)=XPQ(3)
            DIST(4)=XPQ(4)
c            DIST(4)=0
            DIST(5)=0
            DIST(6)=0
            DIST(7)=XPQ(-1)
            DIST(8)=XPQ(-2)
            DIST(9)=XPQ(-3)
            DIST(10)=XPQ(-4)
c            DIST(10)=0
            DIST(11)=0
            DIST(12)=0
            DIST(13)=XPQ(0)
c            print *,'ind 1 glu :',xxx,qsca,dist(13)
            GOTO 999
            ELSEIF (IBEAM.EQ.2) THEN
            if(first) then
                PARMA(1)='init0'
                VALU(1)=0.D0
                call pdfset(PARMA, VALU)
                first=.false.
             endif
C     Reggeon 
            parma(1)='NPTYPE'
            parma(2)='NGROUP'
            parma(3)='NSET'
            valu(1)=2
c            valu(2)=5
c            valu(3)=2            ! GRV-P LO
            valu(2)=1
            valu(3)=1            ! OWENS

            CALL PDFSET(PARMA,VALU)
            CALL STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
            DIST(1)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(2)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(3)=Cr*STR
            DIST(4)=Cr*CHM
            DIST(5)=0
            DIST(6)=0
            DIST(7)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(8)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(9)=Cr*STR
            DIST(10)=Cr*CHM
            DIST(11)=0
            DIST(12)=0
            DIST(13)=Cr*GLU
c            print *,'ind 2 glu :',x,qsca,dist(13)
c	    print *,'enter glu pion :',glu
            GOTO 999
	    ENDIF

c Yura CHR Pom Reg	    
	   ELSEIF (NFLUX.EQ.21) THEN
c	   print *,'index in 21 :',ind,ibeam
	   IF(IBEAM.EQ.2) THEN 
           BCQ=QSCA*QSCA
C Initialise xpq
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,ifit)
c CHR call to different QCD fits
            CALL QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)=XPQ(1)
            DIST(2)=XPQ(2)
            DIST(3)=XPQ(3)
            DIST(4)=XPQ(4)
c            DIST(4)=0
            DIST(5)=0
            DIST(6)=0
            DIST(7)=XPQ(-1)
            DIST(8)=XPQ(-2)
            DIST(9)=XPQ(-3)
            DIST(10)=XPQ(-4)
c            DIST(10)=0
            DIST(11)=0
            DIST(12)=0
            DIST(13)=XPQ(0)
            GOTO 999
            ELSEIF (IBEAM.EQ.1) THEN
            if(first) then
                PARMA(1)='init0'
                VALU(1)=0.D0
                call pdfset(PARMA, VALU)
                first=.false.
             endif
C     Reggeon 
            parma(1)='NPTYPE'
            parma(2)='NGROUP'
            parma(3)='NSET'
            valu(1)=2
*            valu(2)=5
*            valu(3)=2            ! GRV-P LO
            valu(2)=1
            valu(3)=1            ! OWENS

            CALL PDFSET(PARMA,VALU)
            CALL STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
            DIST(1)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(2)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(3)=Cr*STR
            DIST(4)=Cr*CHM
            DIST(5)=0
            DIST(6)=0
            DIST(7)=Cr*(0.5D0*DNV)+Cr*DSEA
            DIST(8)=Cr*(0.5D0*UPV)+Cr*USEA
            DIST(9)=Cr*STR
            DIST(10)=Cr*CHM
            DIST(11)=0
            DIST(12)=0
            DIST(13)=Cr*GLU
c	    print *,'enter glu pion :',glu
            GOTO 999
	    ENDIF



c Yura/CHR Photon pom
         ELSEIF (NFLUX.EQ.20) THEN
	 IF(IBEAM.EQ.2) THEN
            BCQ=QSCA*QSCA
C Initialise xpq
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,ifit)
c CHR call to different QCD fits
            CALL QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)=XPQ(1)
            DIST(2)=XPQ(2)
            DIST(3)=XPQ(3)
            DIST(4)=XPQ(4)
c            DIST(4)=0
            DIST(5)=0
            DIST(6)=0
            DIST(7)=XPQ(-1)
            DIST(8)=XPQ(-2)
            DIST(9)=XPQ(-3)
            DIST(10)=XPQ(-4)
c            DIST(10)=0
            DIST(11)=0
            DIST(12)=0
            DIST(13)=XPQ(0)
            GOTO 999
            ENDIF

c Yura/CHR Pom Photon 
         ELSEIF (NFLUX.EQ.22) THEN
	 IF(IBEAM.EQ.1) THEN
            BCQ=QSCA*QSCA
C Initialise xpq
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,ifit)
c CHR call to different QCD fits
            CALL QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)=XPQ(1)
            DIST(2)=XPQ(2)
            DIST(3)=XPQ(3)
            DIST(4)=XPQ(4)
c            DIST(4)=0
            DIST(5)=0
            DIST(6)=0
            DIST(7)=XPQ(-1)
            DIST(8)=XPQ(-2)
            DIST(9)=XPQ(-3)
            DIST(10)=XPQ(-4)
c            DIST(10)=0
            DIST(11)=0
            DIST(12)=0
            DIST(13)=XPQ(0)
            GOTO 999
            ENDIF

         ELSEIF (NFLUX.EQ.11) THEN
C     User defined pomeron structure
C ... M.B.
            DO ILOOP=-6,6
               XPQ(ILOOP)=0.D0
            ENDDO
            BCQ = 75d0
            XXX = X
c            call QCD_1994(XXX,BCQ,XPQ,0)
c CHR call to different QCD fits
            call QCDFIT(XXX,BCQ,XPQ,ifit)
            DIST(1)  = XPQ(1)  
            DIST(2)  = XPQ(2) 
            DIST(3)  = XPQ(3) 
            DIST(4)  = XPQ(4)
c            DIST(4)  = 0
            DIST(5)  = 0 
            DIST(6)  = 0
            DIST(7)  = XPQ(-1) 
            DIST(8)  = XPQ(-2) 
            DIST(9)  = XPQ(-3) 
            DIST(10) = XPQ(-4) 
c            DIST(10) = 0
            DIST(11) = 0
            DIST(12) = 0
            DIST(13) = XPQ(0) 
            GOTO 999
C ... end M.B.
         ELSE
               write(*,*) 'FPMC : NFLUX is not set'
               STOP
         ENDIF
      ENDIF
* B. C. Pomwig End of mod.

      IF (IDHAD.EQ.59.OR.IDHAD.EQ.71) THEN
        IF (MPDF.GE.0) THEN
C---USE PDFLIB PHOTON STRUCTURE FUNCTIONS
          PARM(1)=AUTPDF(IBEAM)
          VAL(1)=FLOAT(MPDF)
C---FIX TO CALL SCHULER-SJOSTRAND CODE
          IF (AUTPDF(IBEAM).EQ.'SaSph') THEN
            XSP=X
            IF (    XSP.LE.ZERO) CALL HWWARN('HWSFUN',102,*999)
            IF (ONE-XSP.LE.ZERO) CALL HWWARN('HWSFUN',103,*999)
            Q2=QSCA**2
            ISET=MOD(MODPDF(IBEAM),10)
            IOP1=MOD(MODPDF(IBEAM)/10,2)
            IOP2=MOD(MODPDF(IBEAM)/20,2)
            IP2=MODPDF(IBEAM)/100
            IF (IOP2.EQ.0) THEN
              P2=0.
            ELSE
              IHAD=IBEAM
              IF (JDAHEP(1,IHAD).NE.0) IHAD=JDAHEP(1,IHAD)
              P2=PHEP(5,IHAD)**2
            ENDIF
            CALL SASGAM(ISET,XSP,Q2,P2,IP2,F2GM,XPGA)
            IF (IOP1.EQ.1 .AND. ISTAT.LT.10) THEN
              DO 5 I=-6,6
 5            XPGA(I)=XPVMD(I)+XPANL(I)+XPBEH(I)+XPDIR(I)
            ENDIF
            UPV=XPGA(2)
            DNV=XPGA(1)
            USEA=XPGA(2)
            DSEA=XPGA(1)
            STR=XPGA(3)
            CHM=XPGA(4)
            BTM=XPGA(5)
            TOP=XPGA(6)
            GLU=XPGA(0)
          ELSE
            CALL PDFSET(PARM,VAL)
            IF (X.LT.PDFXMN.AND.PDFWRX(IBEAM,1) .OR.
     &          X.GT.PDFXMX.AND.PDFWRX(IBEAM,2)) THEN
              CALL HWWARN('HWSFUN',2,*999)
              WRITE (6,'(2A)') ' WARNING: PDFLIB CALLED WITH X',
     &             ' OUTSIDE ALLOWED RANGE!'
              WRITE (6,'(1P,3(A,E9.3))') ' X VALUE=',X,
     &             ', MINIMUM=',PDFXMN,', MAXIMUM=',PDFXMX
              WRITE (6,'(A)') ' NO FURTHER WARNINGS WILL BE ISSUED'
              IF (X.LT.PDFXMN) PDFWRX(IBEAM,1)=.FALSE.
              IF (X.GT.PDFXMX) PDFWRX(IBEAM,2)=.FALSE.
            ENDIF
            IF (QSCA**2.LT.PDFQMN.AND.PDFWRQ(IBEAM,1) .OR.
     &          QSCA**2.GT.PDFQMX.AND.PDFWRQ(IBEAM,2)) THEN
              CALL HWWARN('HWSFUN',3,*999)
              WRITE (6,'(2A)') ' WARNING: PDFLIB CALLED WITH Q',
     &             ' OUTSIDE ALLOWED RANGE!'
              WRITE (6,'(1P,3(A,E9.3))') ' Q VALUE=',QSCA,
     &             ', MINIMUM=',SQRT(PDFQMN),', MAXIMUM=',SQRT(PDFQMX)
              WRITE (6,'(A)') ' NO FURTHER WARNINGS WILL BE ISSUED'
              IF (QSCA**2.LT.PDFQMN) PDFWRQ(IBEAM,1)=.FALSE.
              IF (QSCA**2.GT.PDFQMN) PDFWRQ(IBEAM,2)=.FALSE.
            ENDIF
            CALL STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
          ENDIF
          DIST(1)=DSEA
          DIST(2)=USEA
          DIST(7)=DSEA
          DIST(8)=USEA
        ELSE
          XSP=X
          IF (    XSP.LE.ZERO) CALL HWWARN('HWSFUN',102,*999)
          IF (ONE-XSP.LE.ZERO) CALL HWWARN('HWSFUN',103,*999)
          Q2=SCALEF**2
          W2=Q2*(1-X)/X
          EMC2=4*RMASS(4)**2
          EMB2=4*RMASS(5)**2
          ALAM2=0.160
          NFL=3
          IF (Q2.GT.50.) NFL=4
          IF (Q2.GT.500.) NFL=5
          STR=HWSDGQ(XSP,Q2,NFL,1)
          CHM=HWSDGQ(XSP,Q2,NFL,2)
          GLU=HWSDGG(XSP,Q2,NFL)
          DIST(1)=STR
          DIST(2)=CHM
          DIST(7)=STR
          DIST(8)=CHM
          IF (W2.GT.EMB2) THEN
            BTM=STR
            IF (W2*ALAM2.LT.Q2*EMB2)
     &          BTM=BTM*LOG(W2/EMB2)/LOG(Q2/ALAM2)
          ELSE
            BTM=0.
          ENDIF
          IF (W2.GT.EMC2) THEN
            IF (W2*ALAM2.LT.Q2*EMC2)
     &          CHM=CHM*LOG(W2/EMC2)/LOG(Q2/ALAM2)
          ELSE
            CHM=0.
          ENDIF
          TOP=0.
        ENDIF
C---INCLUDE SUPPRESSION FROM PHOTON VIRTUALITY IF NECESSARY
        IF (PHOMAS.GT.ZERO.AND.(IBEAM.EQ.1.OR.IBEAM.EQ.2)) THEN
          IHAD=IBEAM
          IF (JDAHEP(1,IHAD).NE.0) IHAD=JDAHEP(1,IHAD)
          IF (IDHW(IHAD).EQ.59) THEN
            FAC=LOG((QSCA**2+PHOMAS**2)/(PHEP(5,IHAD)**2+PHOMAS**2))/
     $          LOG((QSCA**2+PHOMAS**2)/(                PHOMAS**2))
            IF (FAC.LT.ZERO) FAC=ZERO
            DIST(1)=DIST(1)*FAC
            DIST(2)=DIST(2)*FAC
            DIST(7)=DIST(7)*FAC
            DIST(8)=DIST(8)*FAC
            STR=STR*FAC
            CHM=CHM*FAC
            BTM=BTM*FAC
            TOP=TOP*FAC
            GLU=GLU*FAC**2
          ELSE
            CALL HWWARN('HWSFUN',1,*999)
          ENDIF
        ENDIF
        GOTO 900
      ENDIF
      IF (MPDF.GE.0) THEN
C---USE PDFLIB NUCLEON STRUCTURE FUNCTIONS
        PARM(1)=AUTPDF(IBEAM)
        VAL(1)=FLOAT(MPDF)
        CALL PDFSET(PARM,VAL)
        IF (X.LT.PDFXMN.AND.PDFWRX(IBEAM,1) .OR.
     &      X.GT.PDFXMX.AND.PDFWRX(IBEAM,2)) THEN
          CALL HWWARN('HWSFUN',4,*999)
          WRITE (6,'(2A)') ' WARNING: PDFLIB CALLED WITH X',
     &         ' OUTSIDE ALLOWED RANGE!'
          WRITE (6,'(1P,3(A,E9.3))') ' X VALUE=',X,
     &         ', MINIMUM=',PDFXMN,', MAXIMUM=',PDFXMX
          WRITE (6,'(A)') ' NO FURTHER WARNINGS WILL BE ISSUED'
          IF (X.LT.PDFXMN) PDFWRX(IBEAM,1)=.FALSE.
          IF (X.GT.PDFXMX) PDFWRX(IBEAM,2)=.FALSE.
        ENDIF
        IF (QSCA**2.LT.PDFQMN.AND.PDFWRQ(IBEAM,1) .OR.
     &      QSCA**2.GT.PDFQMX.AND.PDFWRQ(IBEAM,2)) THEN
          CALL HWWARN('HWSFUN',5,*999)
          WRITE (6,'(2A)') ' WARNING: PDFLIB CALLED WITH Q',
     &         ' OUTSIDE ALLOWED RANGE!'
          WRITE (6,'(1P,3(A,E9.3))') ' Q VALUE=',QSCA,
     &         ', MINIMUM=',SQRT(PDFQMN),', MAXIMUM=',SQRT(PDFQMX)
          WRITE (6,'(A)') ' NO FURTHER WARNINGS WILL BE ISSUED'
          IF (QSCA**2.LT.PDFQMN) PDFWRQ(IBEAM,1)=.FALSE.
          IF (QSCA**2.GT.PDFQMN) PDFWRQ(IBEAM,2)=.FALSE.
        ENDIF
        CALL STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
C--new MRST98 LO PDF's
      ELSEIF(NSET.GE.6.AND.NSET.LE.8) THEN
        CALL HWSMRS(X,MUSCALE,NSET-5,UPV,DNV,USEA,DSEA,STR,CHM,BTM,GLU)
        TOP=ZERO
      ELSE
        IF (NSET.LT.1.OR.NSET.GT.5) CALL HWWARN('HWSFUN',400,*999)
        IF (QSCA.LT.Q0(NSET)) QSCA=Q0(NSET)
        IF (QSCA.NE.QOLD.OR.IDHAD.NE.IOLD.OR.NSET.NE.NOLD) THEN
C---INITIALIZE
          QOLD=QSCA
          IOLD=IDHAD
          NOLD=NSET
          SS=LOG(QSCA/QL(NSET))
          SMIN=LOG(Q0(NSET)/QL(NSET))
          IF (NSET.LT.3.OR.NSET.EQ.5) THEN
            S=LOG(SS/SMIN)
          ELSE
            T=2.*SS
            TMIN=2.*SMIN
            TMAX=2.*LOG(1.E4/QL(NSET))
          ENDIF
          IF (IDHAD.GE.72) THEN
            IF (NSET.LT.3) THEN
              IP=NSET
              DO 10 I=1,5
              DO 10 J=1,6
   10         A(J,I)=B(1,J,I,IP)+S*(B(2,J,I,IP)+S*B(3,J,I,IP))
              DO 20 K=1,2
              AA=ONE+A(2,K)+A(3,K)
   20         G(K)=HWSGAM(AA)/((ONE+A(2,K)*A(4,K)/AA)*HWSGAM(A(2,K))
     &            *HWSGAM(ONE+A(3,K)))
            ELSEIF (NSET.EQ.5) THEN
              DO 21 I=1,5
              DO 21 J=1,6
   21         A(J,I)=BB(1,J,I)+S*(BB(2,J,I)+S*(BB(3,J,I)+S*BB(4,J,I)))
              DO 22 K=1,2
              AA=ONE+A(2,K)+A(3,K)
   22         G(K)=HWSGAM(AA)/((ONE+A(2,K)/AA*(A(4,K)+
     &            (ONE+A(2,K))/(ONE+AA)*A(5,K)))*HWSGAM(A(2,K))
     &            *HWSGAM(ONE+A(3,K)))
            ELSE
              IP=NSET-2
              VT=MAX(-ONE,MIN(ONE,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
              WT=VT*VT
C...CHEBYSHEV POLYNOMIALS FOR T EXPANSION
              TT(1)=1.
              TT(2)=VT
              TT(3)=   2.*WT- 1.
              TT(4)=  (4.*WT- 3.)*VT
              TT(5)=  (8.*WT- 8.)*WT+1.
              TT(6)=((16.*WT-20.)*WT+5.)*VT
            ENDIF
          ELSEIF (NSET.LT.3) THEN
              IP=NSET+2
              DO 30 I=1,5
              DO 30 J=1,6
   30         A(J,I)=B(1,J,I,IP)+S*(B(2,J,I,IP)+S*B(3,J,I,IP))
              AA=ONE+A(2,1)+A(3,1)
              G(1)=HWSGAM(AA)/(HWSGAM(A(2,1))*HWSGAM(ONE+A(3,1)))
              G(2)=0.
           ENDIF
        ENDIF
C
        IF (NSET.LT.3.OR.NSET.EQ.5) THEN
          DO 50 I=1,5
   50     F(I)=A(1,I)*X**A(2,I)*XMWN**A(3,I)*(ONE+X*
     &        (A(4,I)+X*(A(5,I)  +  X*A(6,I))))
          F(1)=F(1)*G(1)
          F(2)=F(2)*G(2)
          UPV=F(1)-F(2)
          DNV=F(2)
          SEA=F(3)/6.
          STR=SEA
          CHM=F(4)
          BTM=ZERO
          TOP=ZERO
          GLU=F(5)
        ELSE
          IF (X.NE.XOLD) THEN
            XOLD=X
            IF (X.GT.0.1) THEN
              NX=1
              VX=(2.*X-1.1)/0.9
            ELSE
              NX=2
              VX=MAX(-ONE,(2.*LOG(X)+11.51293)/6.90776)
            ENDIF
            WX=VX*VX
            TX(1)=1.
            TX(2)=VX
            TX(3)=   2.*WX- 1.
            TX(4)=  (4.*WX- 3.)*VX
            TX(5)=  (8.*WX- 8.)*WX+1.
            TX(6)=((16.*WX-20.)*WX+5.)*VX
          ENDIF
C...CALCULATE STRUCTURE FUNCTIONS
          DO 120 IFL=1,6
          XQSUM=0.
          DO 110 IT=1,6
          DO 110 IX=1,6
  110     XQSUM=XQSUM+CEHLQ(IX,IT,NX,IFL,IP)*TX(IX)*TT(IT)
  120     XQ(IFL)=XQSUM*XMWN**NEHLQ(IFL,IP)
          UPV=XQ(1)
          DNV=XQ(2)
          STR=XQ(5)
          CHM=XQ(6)
          SEA=XQ(3)
          GLU=XQ(4)
C...SPECIAL EXPANSION FOR BOTTOM (THRESHOLD EFFECTS)
          IF (NFLAV.LT.5.OR.T.LE.TBMIN(IP)) THEN
            BTM=0.
          ELSE
            VT=MAX(-ONE,MIN(ONE,(2.*T-TMAX-TBMIN(IP))/(TMAX-TBMIN(IP))))
            WT=VT*VT
            TB(1)=1.
            TB(2)=VT
            TB(3)=   2.*WT- 1.
            TB(4)=  (4.*WT- 3.)*VT
            TB(5)=  (8.*WT- 8.)*WT+1.
            TB(6)=((16.*WT-20.)*WT+5.)*VT
            XQSUM=0.
            DO 130 IT=1,6
            DO 130 IX=1,6
  130       XQSUM=XQSUM+CEHLQ(IX,IT,NX,7,IP)*TX(IX)*TB(IT)
            BTM=XQSUM*XMWN**NEHLQ(7,IP)
          ENDIF
C...SPECIAL EXPANSION FOR TOP (THRESHOLD EFFECTS)
          TPMIN=TTMIN(IP)+TMTOP
C---TMTOP=2.*LOG(TOPMAS/30.)
          TPMAX=TMAX+TMTOP
          IF (NFLAV.LT.6.OR.T.LE.TPMIN) THEN
            TOP=0.
          ELSE
            VT=MAX(-ONE,MIN(ONE,(2.*T-TPMAX-TPMIN)/(TPMAX-TPMIN)))
            WT=VT*VT
            TB(1)=1.
            TB(2)=VT
            TB(3)=   2.*WT- 1.
            TB(4)=  (4.*WT- 3.)*VT
            TB(5)=  (8.*WT- 8.)*WT+1.
            TB(6)=((16.*WT-20.)*WT+5.)*VT
            XQSUM=0.
            DO 150 IT=1,6
            DO 150 IX=1,6
  150       XQSUM=XQSUM+CEHLQ(IX,IT,NX,8,IP)*TX(IX)*TB(IT)
            TOP=XQSUM*XMWN**NEHLQ(8,IP)
          ENDIF
        ENDIF
      ENDIF
      IF (MPDF.LT.0.AND.NSET.LE.5) THEN
        USEA=SEA
        DSEA=USEA
      ENDIF
      IF(MPDF.LT.0.AND.NSET.GT.2.AND.(IDHAD.EQ.38.OR.IDHAD.EQ.30)) THEN
        WRITE(6,*) '     THIS SET OF PDFS DOES NOT SUPPORT PIONS'
        WRITE(6,*) 'EITHER USE SET NSTRU=1,2 OR A PION SET FROM PDFLIB'
        STOP
      ENDIF
      IF (IDHAD.EQ.73.OR.IDHAD.EQ.72) THEN
         DIST(1)=DSEA+DNV
         DIST(2)=USEA+UPV
         DIST(7)=DSEA
         DIST(8)=USEA
      ELSEIF (IDHAD.EQ.91) THEN
         DIST(1)=DSEA
         DIST(2)=USEA
         DIST(7)=DSEA+DNV
         DIST(8)=USEA+UPV
      ELSEIF (IDHAD.EQ.75) THEN
         DIST(1)=USEA+UPV
         DIST(2)=DSEA+DNV
         DIST(7)=USEA
         DIST(8)=DSEA
      ELSEIF (IDHAD.EQ.93) THEN
         DIST(1)=USEA
         DIST(2)=DSEA
         DIST(7)=USEA+UPV
         DIST(8)=DSEA+DNV
      ELSEIF (IDHAD.EQ.38) THEN
         DIST(1)=USEA
         DIST(2)=USEA+UPV
         DIST(7)=USEA+UPV
         DIST(8)=USEA
      ELSEIF (IDHAD.EQ.30) THEN
         DIST(1)=USEA+UPV
         DIST(2)=USEA
         DIST(7)=USEA
         DIST(8)=USEA+UPV
      ELSE
         PRINT *,' CALLED HWSFUN FOR IDHAD =',IDHAD
         CALL HWWARN('HWSFUN',400,*999)
      ENDIF
  900 DIST(3)=STR
      DIST(4)=CHM
      DIST(5)=BTM
      DIST(6)=TOP
      DIST(9)=STR
      DIST(10)=CHM
      DIST(11)=BTM
      DIST(12)=TOP
      DIST(13)=GLU
      DO 901 I=1,13
      IF (DIST(I).LT.DMIN) DIST(I)=DMIN
  901 CONTINUE
C---FOR REMNANT NUCLEONS SWITCH OFF VALENCE QUARKS,
C   WHILE MAINTAINING MOMENTUM SUM RULE
      IF (IDHAD.EQ.72) THEN
        TOTAL=0
        DO 910 I=1,13
          TOTAL=TOTAL+DIST(I)
 910    CONTINUE
        DIST(1)=DIST(1)-DNV
        DIST(2)=DIST(2)-UPV
        IF (TOTAL.GT.DNV+UPV) THEN
          DO 920 I=1,13
            DIST(I)=DIST(I)*TOTAL/(TOTAL-DNV-UPV)
 920      CONTINUE
        ENDIF
      ENDIF
  999 CONTINUE
      END



* User defined pomeron structure function routine. 
* B. Cox 15/05/2001
*
*  Input : X = x_{i/IP}
*          Q2 = photon virtuality
*
*  Output: XPQ(-6:6): PDG style array of partons, 0=gluon.
*          

      SUBROUTINE POMSTR(X,Q2,XPQ)
      DOUBLE PRECISION XPQ,X,Q2
      DIMENSION XPQ(-6:6)

      XPQ(1)=X*(1.0-X)
      XPQ(2)=X*(1.0-X)
*      XPQ(1)=X**-2
*      XPQ(2)=X**-2
      XPQ(3)=0.
      XPQ(4)=0.
      XPQ(5)=0.
      XPQ(6)=0.
      XPQ(-1)=X*(1.0-X)
      XPQ(-2)=X*(1.0-X)
*      XPQ(-1)=X**-2
*      XPQ(-2)=X**-2
      XPQ(-3)=0.
      XPQ(-4)=0.
      XPQ(-5)=0.
      XPQ(-6)=0.
      XPQ(0)=X*(1.0-X)

      RETURN
      END

CDECK  ID>, HWSGEN.
*CMZ :-        -26/04/91  14.55.45  by  Federico Carminati
*-- Author :    Bryan Webber
C Maarten Boonekamp, 10/03/06 : introduce this to avoid cutoffs on X, 
C                               in inclusive chi production
C-----------------------------------------------------------------------
      SUBROUTINE HWSGEN_CHI(GENEX)
C-----------------------------------------------------------------------
C     GENERATES X VALUES (IF GENEX)
C     EVALUATES STRUCTURE FUNCTIONS AND ENFORCES CUTOFFS ON X
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWBVMC,HWRUNI,X,QL
      INTEGER I,J
      LOGICAL GENEX
      EXTERNAL HWBVMC,HWRUNI
      IF (GENEX) THEN
        XX(1)=EXP(HWRUNI(0,ZERO,XLMIN))
        XX(2)=XXMIN/XX(1)
      ENDIF
      DO 10 I=1,2
        J=I
        IF (JDAHEP(1,I).NE.0) J=JDAHEP(1,I)
        X=XX(I)
        QL=(1.-X)*EMSCA
        if(EMSCA*EMSCA.gt.30000) print*, 'EMSCA ', EMSCA
        CALL HWSFUN(X,EMSCA,IDHW(J),NSTRU,DISF(1,I),I)
* Only modification wrt HWSGEN is below: we do not cut off at the gluon
* consitutent mass, but ensure there remains enough phase space to create
* at least one pion
*      DO 10 J=1,13
*        IF (QL.LT.RMASS(21)) DISF(J,I)=0.
   10 CONTINUE
      END

C-----------------------------------------------------------------------
      SUBROUTINE PRINTSETTING
C-----------------------------------------------------------------------
C-- Author :    Oldrich Kepka
C-- prints important user setting
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'

      write(*,*) ''
      write(*,*) ' FPMC - USER SETTINGS'
      write(*,*) ' --------------------'
      write(*,*) '          NFLUX    = ',NFLUX
      write(*,*) '          TYPEPR   = ',TYPEPR
      write(*,*) '          TYPEINT  = ',TYPINT
      write(*,*) '          GAPSPR   = ',GAPSPR
      write(*,*) '          PROSPR   = ',PROSPR
      write(*,*) '          CDFFAC   = ',CDFFAC
      write(*,*) '          MAXEV    = ',MAXEV
      write(*,*) '          PART1    = ',PART1
      write(*,*) '          PART2    = ',PART2
      write(*,*) '          NRN(1)   = ',NRN(1)
      write(*,*) '          NRN(2)   = ',NRN(2)
      write(*,*) '          WMASS    = ',RMASS(198)
      write(*,*) '          TMASS    = ',RMASS(6)
      write(*,*) '          HMASS    = ',RMASS(201)
      write(*,*) '          YJMIN    = ', YJMIN
      write(*,*) '          YJMAX    = ',YJMAX
      write(*,*) '          PTMIN    = ',PTMIN
      write(*,*) '          EMMIN    = ',EMMIN
      write(*,*) '          IFIT     = ',IFITPDF
      write(*,*) '          ISOFTM   = ',ISOFTM
      write(*,*) '          Q2WWMN   = ',Q2WWMN
      write(*,*) '          Q2WWMX   = ',Q2WWMX
      write(*,*) '          YWWMIN   = ',YWWMIN
      write(*,*) '          YWWMAX   = ',YWWMAX
      write(*,*) '          ZION     = ',ZION
      write(*,*) '          AION     = ',AION
      write(*,*) '          AAANOM = ',AAANOM
      write(*,*) '          DKAPPA   = ',D_KAPPA
      write(*,*) '          LAMBDA   = ',LAMBDA
      write(*,*) '          A0W      = ',A0W
      write(*,*) '          ACW      = ',ACW
      write(*,*) '          A0Z      = ',A0Z
      write(*,*) '          ACZ      = ',ACZ
      write(*,*) '          A1A      = ',A1A
      write(*,*) '          A2A      = ',A2A
      write(*,*) '          ANOMCUTOFF = ',ANOMCUTOFF
      write(*,*) '          AAEXOTIC = ',AAEXOTIC
      write(*,*) '          AAM      = ',AAM
      write(*,*) '          AAQ      = ',AAQ
      write(*,*) '          AAN      = ',AAN
      write(*,*) '----------others-----------------'
      write(*,*) '          PTMAX    = ',PTMAX
      write(*,*)
      END     


CHECK  ID>, FPMCDEBUGPART
C-- Author :    Oldrich Kepka
C-- particle listing, for debugging pourpouses
C-----------------------------------------------------------------------

      SUBROUTINE FPMCDEBUGPART
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'

      INTEGER N, IPART, I, ID
      DOUBLE PRECISION PX, PY, PT, PZ, E, M
      PRINT *, 'Stable particle list'
      PRINT *, 'PX, PY, PZ, E, M, PT, ID, STATE'
      N=NHEP
      IPART=0
      do 1515 I=1,N
          IF(ISTHEP(I).EQ.1) THEN
             IPART=IPART+1
             PX=PHEP(1,I)
             PY=PHEP(2,I)
             PZ=PHEP(3,I)
             E =PHEP(4,I)
             M=PHEP(5,I)
             ID=IDHEP(I)
             PT=DSQRT(PX*PX+PY*PY) 
             print '(6F10.2,2I6)', PX, PY, PZ, E, M, PT, ID, ISTHEP(I)

c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
          ENDIF
1515  continue 
      PRINT *, '# of stables=', IPART

      PRINT *, ''   
      PRINT *, 'All particle list'

      N=NHEP
      IPART=0 
      do 1615 I=1,N
             IPART=IPART+1
             PX=PHEP(1,I)
             PY=PHEP(2,I)
             PZ=PHEP(3,I)
             E =PHEP(4,I)
             M=PHEP(5,I)
             ID=IDHEP(I)
             PT=DSQRT(PX*PX+PY*PY) 
             print '(6F10.2,2I6)', PX, PY, PZ, E, M, PT, ID, ISTHEP(I)

c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
1615  continue
      END


