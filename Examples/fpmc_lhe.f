C-----------------------------------------------------------------------
CDECK  ID>, HWUEMP.
*CMZ :-        -20/10/2013    
*-- Author :   M.E.POL from HWUEPR: Ian Knowles, Bryan Webber & Kosuke Odagiri
C-----------------------------------------------------------------------
      SUBROUTINE HWUEMP
C-----------------------------------------------------------------------
C     Prints out event data in LHE format in unit 45
C     Works for QED
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INTEGER I,J,IST,ILEP
      INTEGER IGA1,IGA2,IDA1,IDA2,IPR1,IPR2
      INTEGER JMOT,ISIS1,ISIS2,IFLEP1,IFLEP2,IFLEP3,IFLEP4
      INTEGER II1,II2,II3,II4,II5
      INTEGER NUP,IDRPUP
      DOUBLE PRECISION AMASSG,VTIM,ASPI,ALFAS,HWUALF,SCALE
      EXTERNAL HWUALF
      II1=0
      II2=0
      II3=0
      II4=0
      II5=0
      ILEP=0
      JMOT=0
      ISIS1 = 0
      ISIS2 = 0
      IFLEP1 = 0
      IFLEP2 = 0
      IFLEP3 = 0
      IFLEP4 = 0
      NUP = 10
C      NUP = 4
C      IDRPUP = 10042
      IDRPUP = 1
      AMASSG=0.0
      VTIM = 0.0
      ASPI = 9.0
      ALFAS = HWUALF(1,EMSCA)
      SCALE = 1.0
C Filter leptonic decays
      DO 300 J=1,NHEP
        JMOT = JMOHEP(1,J)
        IF(ABS(IDHEP(J)).EQ.11.AND.ABS(IDHEP(JMOT)).EQ.23)ILEP = ILEP+1
        IF(ABS(IDHEP(J)).EQ.12.AND.ABS(IDHEP(JMOT)).EQ.23)ILEP = ILEP+1
 300  CONTINUE
      IF (ILEP.NE.2) GOTO 999
C See if there is a FSR photons
C      DO 350 J=1,NHEP
C        JMOT = JMOHEP(1,J)
C        IF(IDHEP(J).EQ.22.AND.ABS(IDHEP(JMOT)).EQ.24)THEN
C          WRITE(45,192)
C        ENDIF
C 350  CONTINUE
C It is a WW -> electrons and/or muons
C Write start of event in the lhe file
      WRITE(45,40)
      WRITE(45,41)NUP,IDRPUP,EVWGT,SCALE,ALPHEM,ALFAS 
      DO 410 I=1,NHEP
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.1)IPR1=I
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.2)IPR2=I
        IF(RNAME(IDHW(I)).EQ.'HARD')THEN
          IGA1 = JMOHEP(1,I)
          IGA2 = JMOHEP(2,I)
          IDA1 = JDAHEP(1,I)
          IDA2 = JDAHEP(2,I)
          ISIS1 = JDAHEP(1,IDA1)
          ISIS2 = JDAHEP(1,IDA2)
          IF(ABS(IDHEP(ISIS1)).EQ.24)THEN
             IFLEP1 = JDAHEP(1,ISIS1)
             IFLEP2 = JDAHEP(2,ISIS1)
          ENDIF
          IF(ABS(IDHEP(ISIS2)).EQ.24)THEN
             IFLEP3 = JDAHEP(1,ISIS2)
             IFLEP4 = JDAHEP(2,ISIS2)
          ENDIF
        ENDIF
  410 CONTINUE
C Print the gammas
      DO 411 I=1,NHEP
        IF (I.EQ.IGA1.OR.I.EQ.IGA2)THEN
          II1 = -1
          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
     &      (PHEP(J,I),J=1,5),VTIM,ASPI
C Print the daughters of the gammas(WW)
        ELSE IF(I.EQ.IDA1.OR.I.EQ.IDA2)THEN
          II1 = 2
          II2 = 1
          II3 = 2
          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
     &      (PHEP(J,I),J=1,5),VTIM,ASPI
        ELSE IF(I.EQ.IFLEP1.OR.I.EQ.IFLEP2)THEN
          II1 = 1
          II2 = 3
          II3 = 3
          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
     &      (PHEP(J,I),J=1,5),VTIM,ASPI
        ELSE IF(I.EQ.IFLEP3.OR.I.EQ.IFLEP4)THEN
          II1 = 1
          II2 = 4
          II3 = 4
          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
     &      (PHEP(J,I),J=1,5),VTIM,ASPI
        ENDIF
  411 CONTINUE
C Now print the outgoing protons
      DO 412 I=1,NHEP
        IF (I.EQ.IPR1.OR.I.EQ.IPR2)THEN
          II1 = 1
          II2 = 1
          II3 = 2
          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
     &      (PHEP(J,I),J=1,5),VTIM,ASPI
C          WRITE(45,191) I,RNAME(IDHW(I)),IDHEP(I),IST,
C     &      JMOHEP(1,I),JMOHEP(2,I),JDAHEP(1,I),JDAHEP(2,I),
C     &      (PHEP(J,I),J=1,5),(VHEP(J,I),J=1,4)
        ENDIF
  412 CONTINUE
      WRITE(45,50) 
   30 FORMAT(///1X,'EVENT ',I7,':',F8.2,' GEV/C ',A8,' ON ',F8.2,
     & ' GEV/C ',A8,' PROCESS:',I6/1X,'SEEDS: ',I11,' & ',I11,
     & '   STATUS: ',I4,' ERROR:',I4,'  WEIGHT: ',1P,E11.4/)
   40 FORMAT('<event>')
   50 FORMAT('</event>')
   41 FORMAT(1X,I2,2X,I5,1X,4(1X,E15.8))
  190 FORMAT(1X,I8,5I4,5(2X,E15.8),2F3.0)
  191 FORMAT('#',I4,1X,A8,I8,5I4,2F8.2,2F7.1,F8.2,1P,4E10.3)
  192 FORMAT(1X,'There is a final state photon')
  999 CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LHEEVT
C-----------------------------------------------------------------------
C     Prints out event data in LHE format in unit 45
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INTEGER J,IST,ILEP
      INTEGER I
      INTEGER IGA1,IGA2,IDA1,IDA2,IPR1,IPR2
      INTEGER JMOT,ISIS1,ISIS2,IFLEP1,IFLEP2,IFLEP3,IFLEP4
      INTEGER II1,II2,II3,II4,II5
      INTEGER NUP,IDRPUP
      DOUBLE PRECISION AMASSG,VTIM,ASPI,ALFAS,HWUALF,SCALE
      
      EXTERNAL HWUALF
      II1=0
      II2=0
      II3=0
      II4=0
      II5=0
      ILEP=0
      JMOT=0
      ISIS1 = 0
      ISIS2 = 0
      IFLEP1 = 0
      IFLEP2 = 0
      IFLEP3 = 0
      IFLEP4 = 0
      NUP = 0
      IDRPUP = 1
      AMASSG=0.0
      VTIM = 0.0
      ASPI = 9.0
      ALFAS = HWUALF(1,EMSCA)
      SCALE = 1.0

C
C Write start of event in the lhe file
C
      WRITE(45,40)

      DO 510 I=1,NHEP
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.1)IPR1=I
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.2)IPR2=I
     
        NUP = NUP + 1
  510 CONTINUE

      WRITE(45,41)NUP,IDRPUP,EVWGT,SCALE,ALPHEM,ALFAS 
      DO 511 I=1,NHEP
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.1)IPR1=I
        IF(IDHEP(I).EQ.2212.AND.JMOHEP(1,I).EQ.2)IPR2=I
    
        WRITE(45,190) IDHEP(I),ISTHEP(I),
     &    JMOHEP(1,I),JMOHEP(2,I),JDAHEP(1,I),
     &    JDAHEP(2,I),(PHEP(J,I),J=1,5),VTIM,ASPI
      
  511 CONTINUE

C Now print the outgoing protons
C      DO 412 I=1,NHEP
C        IF (I.EQ.IPR1.OR.I.EQ.IPR2)THEN
C          II1 = 1
C          II2 = 1
C          II3 = 2
C          WRITE(45,190)IDHEP(I),II1,II2,II3,II4,II5,
C     &      (PHEP(J,I),J=1,5),VTIM,ASPI
C        ENDIF
C  412 CONTINUE

      WRITE(45,50)
 
   30 FORMAT(///1X,'EVENT ',I7,':',F8.2,' GEV/C ',A8,' ON ',F8.2,
     & ' GEV/C ',A8,' PROCESS:',I6/1X,'SEEDS: ',I11,' & ',I11,
     & '   STATUS: ',I4,' ERROR:',I4,'  WEIGHT: ',1P,E11.4/)
   40 FORMAT('<event>')
   50 FORMAT('</event>')
   41 FORMAT(1X,I5,2X,I5,1X,4(1X,E15.8))
  190 FORMAT(1X,I8,5I5,5(2X,E15.8),2F3.0)
  191 FORMAT('#',I4,1X,A8,I8,5I4,2F8.2,2F7.1,F8.2,1P,4E10.3)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LHEINI
C-----------------------------------------------------------------------
C     Prints out initialization in LHE format in unit 45
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'ffcard.inc'
      INTEGER IPR1,IPR2,I1,I2,I3,I4,III5,IPROC2
      DOUBLE PRECISION ERWGT

C FIXME
      IPR1 = 2212
      IPR2 = 2212
C      I1 = 0
      I1 = 1
      I2 = 0
C      I3 = 3
      I3 = 4
      I4 = 1
      III5 = -1
      IPROC2 = -1

      OPEN(UNIT=45,FILE=ULHEFILE)

      WRITE(45,30)
      WRITE(45,40)
C      WRITE(45,41)IPR1,IPR2,PBEAM1,PBEAM2,III5,III5,III5,III5,I3,I4
C      WRITE(45,42)1000.*AVWGT,1000.*ERWGT,AVWGT,IPROC2
      WRITE(45,41)IPR1,IPR2,PBEAM1,PBEAM2,I2,I2,IPROC2,IPROC2,I3,I4
      WRITE(45,42)1000.*AVWGT,1000.*ERWGT,AVWGT,I1
      WRITE(45,50)
   30 FORMAT('<LesHouchesEvents version=' '"' '1.0' '"' '>')
   40 FORMAT('<header>'/'</header>'/'<init>')
C   41 FORMAT(2(5X,I4),2(2X,E15.8),2(2X,I2),2(2x,I2),2(2X,I1))
   41 FORMAT(2(5X,I4),2(2X,E15.8),2(2X,I1),2(2x,I5),2(2X,I1))
   42 FORMAT(3(2x,E15.8),2X,I5)   
   50 FORMAT('</init>')    
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LHEEND
C-----------------------------------------------------------------------
      WRITE(45,151)
 151  FORMAT('</LesHouchesEvents>')

      CLOSE(45)

      END
