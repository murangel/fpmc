CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Interface to the KMR luminosity functions,
C Reading tables produced with code from L.Lonnblad
C . tables are in ./Tables/kmrLum_*.txt
C . code is in ./Ext/*cc
C
C     MB, Dec 2004
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE KMRINI(ILUM,S,GAPSURV,USCALE,DELTA,Q2CUT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL GET_LUMI_LONNBLAD_TEV, GET_LUMI_LONNBLAD_LHC,
     &  GET_LUMI_EXHUME_TEV, GET_LUMI_EXHUME_LHC
C ... begin R.S.
      INTEGER LUMINOSITY,DELTA
      COMMON /KMR/ LUMINOSITY
      LUMINOSITY = ILUM
C ... end R.S
      print*, AAAAA


      IF(ILUM.EQ.1) THEN
      PRINT *, 'Initializing KMR luminosity from : LONNBLAD_TEV'
         CALL GET_LUMI_LONNBLAD_TEV()
      ELSEIF(ILUM.EQ.2) THEN
      PRINT *, 'Initializing KMR luminosity from : LONNBLAD_LHC'
         CALL GET_LUMI_LONNBLAD_LHC()
      ELSEIF(ILUM.EQ.3) THEN
      PRINT *, 'Initializing KMR luminosity from : EXHUME_TEV'
         CALL GET_LUMI_EXHUME_TEV()
      ELSEIF(ILUM.EQ.4) THEN
      PRINT *, 'Initializing KMR luminosity from : EXHUME_LHC'
         CALL GET_LUMI_EXHUME_LHC()
C ... begin R.S.
      ELSEIF(ILUM.EQ.99) THEN
      PRINT *, 'Initializing KMR luminosity'
         CALL KMR2_INIT(S,GAPSURV,USCALE,DELTA,Q2CUT) 
C ... end R.S.         
      ELSE
         PRINT*, 'KMRINI -- Change ILUM'
         STOP
      ENDIF

c      CHARACTER*80 CHFILE
c      IF(ILUM.EQ.1) THEN
c         CHFILE = './Tables/kmrLumi_Lonnblad_Tev.txt'
c      ELSEIF(ILUM.EQ.2) THEN
c         CHFILE = './Tables/kmrLumi_Lonnblad_LHC.txt'
c      ELSEIF(ILUM.EQ.3) THEN
c         CHFILE = './Tables/kmrLumi_ExHume_Tev.txt'
c      ELSEIF(ILUM.EQ.4) THEN
c         CHFILE = './Tables/kmrLumi_ExHume_LHC.txt'
c      ELSE
c         PRINT*, 'KMRINI -- Change ILUM'
c         STOP
c      ENDIF
c      PRINT*, 'Initializing KMR luminosity from :', CHFILE
c      OPEN(55,FILE=CHFILE,STATUS='OLD')
c      READ(55,*) NLINES, XMIN, XMAX
c      IMAX = NINT(NLINES*(1.*NLINES+1.)/2.)
c      DO I = 1, IMAX
c         READ(55,*) X1(I), X2(I), XMAS(I), XLUM(I)
c      ENDDO
c      CLOSE(55)


      RETURN
      END

      SUBROUTINE KMRINT(X1TEST,X2TEST,XKMRLUM)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X1(25000), X2(25000), XLUM(25000), XMAS(25000)
      COMMON / KMRTAB / XMINLOG, XMAXLOG, STEP, 
     &                  X1, X2, XLUM, XMAS, NLINES
      DOUBLE PRECISION X1OK(4), X2OK(4), XLUMOK(4)

C ... begin R.S.
      INTEGER LUMINOSITY
      COMMON /KMR/ LUMINOSITY
      DOUBLE PRECISION KMR2_GET_LUMI
      EXTERNAL KMR2_GET_LUMI
      IF(LUMINOSITY.EQ.99) THEN
        XKMRLUM = KMR2_GET_LUMI(X1TEST,X2TEST)
        RETURN
      ENDIF
C ... end R.S
      XTESTMIN = MIN(X1TEST,X2TEST)
      XTESTMAX = MAX(X1TEST,X2TEST)
      XTLOG1 = DLOG10(XTESTMIN)
      XTLOG2 = DLOG10(XTESTMAX)

      ibin1 = int((xtlog1-xminlog)/step) + 1
      ibin2 = int((xtlog2-xminlog)/step) + 1
      ioffset = int( (1.*ibin1-1.) * (1.*(nlines+1)-1.*ibin1/2.) )
      iline = ioffset + ibin2 - ibin1 + 1
      ioffs2 = int( (1.*ibin1) * (1.*(nlines+1)-1.*(ibin1+1)/2.) )
      iline2 = ioffs2 + ibin2 - ibin1

      if(ibin1.eq.ibin2) then
         x1ok(1) = x1(iline)
         x1ok(2) = x1(iline+1) 
         x1ok(3) = x1(iline2+1)  
         x1ok(4) = x1(iline2+1)
         x2ok(1) = x2(iline)
         x2ok(2) = x2(iline+1) 
         x2ok(3) = x1(iline2)
         x2ok(4) = x2(iline2+1)
         xlumok(1) = xlum(iline)
         xlumok(2) = xlum(iline+1) 
         xlumok(3) = xlum(iline+1)  
         xlumok(4) = xlum(iline2+1)
      else
         x1ok(1) = x1(iline)
         x1ok(2) = x1(iline+1) 
         x1ok(3) = x1(iline2)  
         x1ok(4) = x1(iline2+1)
         x2ok(1) = x2(iline)
         x2ok(2) = x2(iline+1) 
         x2ok(3) = x2(iline2)  
         x2ok(4) = x2(iline2+1)
         xlumok(1) = xlum(iline)
         xlumok(2) = xlum(iline+1) 
         xlumok(3) = xlum(iline2)  
         xlumok(4) = xlum(iline2+1)
      endif

*      print*, 'kmr ', xtestmin, xtestmax, ibin1, ibin2, iline, iline2
*      print*, 'kmr ', x1ok
*      print*, 'kmr ', x2ok
*      print*, 'kmr ', xlumok

      xlum1 = xlumok(1)+(xlumok(2)-xlumok(1))*
     &     (xtestmax-x2ok(1))/(x2ok(2)-x2ok(1))
      xlum2 = xlumok(3)+(xlumok(4)-xlumok(3))*
     &     (xtestmax-x2ok(3))/(x2ok(4)-x2ok(3))
      xlum3 = xlum1+(xlum2-xlum1)*
     &     (xtestmin-x1ok(1))/(x1ok(3)-x1ok(1))

      XKMRLUM = xlum3

      RETURN
      END



