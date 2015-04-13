      SUBROUTINE I_NLOTABLES(IFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NQ2MAX=100,NZMAX=100)
      DIMENSION UDS_GRID(NZMAX,NQ2MAX),C_GRID(NZMAX,NQ2MAX),
     &          G_GRID(NZMAX,NQ2MAX)
      dimension beta_grid(nzmax),q2_grid(nq2max)
      DIMENSION UDS_GRDL(NZMAX,NQ2MAX),C_GRDL(NZMAX,NQ2MAX),
     &          G_GRDL(NZMAX,NQ2MAX)
      COMMON /CPARAMB/ BETA_GRID, Q2_GRID,UDS_GRID,C_GRID,G_GRID
      INTEGER ILINE
      INTEGER I,J, ROW, COL
      PARAMETER(ROW=10000,COL=5)
      INCLUDE 'pdftable_h194ext.inc'
      INCLUDE 'pdftable_h197.inc'
      INCLUDE 'pdftable_h1zeus.inc'
      INCLUDE 'pdftable_zeus.inc'

c read input tables and put numbers in grid
      ILINE = 1
      IF(IFIT.EQ.8) THEN
         DO IQ2_TEST = 1,NQ2MAX
            DO IZ_TEST  = 1,NZMAX
             BETA_GRID(IZ_TEST)         = TABLE_H194EXT(ILINE, 1) 
             Q2_GRID(IQ2_TEST)          = TABLE_H194EXT(ILINE, 2)
             UDS_GRID(IZ_TEST,IQ2_TEST) = TABLE_H194EXT(ILINE, 3)
             G_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H194EXT(ILINE, 4)
             C_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H194EXT(ILINE, 5)

             UDS_GRID(IZ_TEST,IQ2_TEST)=UDS_GRID(IZ_TEST,IQ2_TEST)/6.
             ILINE=ILINE+1

C      WRITE(91,1004) BETA_GRID(IZ_TEST),Q2_GRID(IQ2_TEST),
C     + UDS_GRID(IZ_TEST,IQ2_TEST),
C     + G_GRID(IZ_TEST,IQ2_TEST),
C     + C_GRID(IZ_TEST,IQ2_TEST)
C1004 format(2(F14.7,1X),3(F10.5,1X))
            ENDDO
         ENDDO
c     table 10
      ELSEIF(IFIT.EQ.10) THEN
         DO IQ2_TEST = 1,NQ2MAX
            DO IZ_TEST  = 1,NZMAX
             BETA_GRID(IZ_TEST)         = TABLE_H197(ILINE, 1) 
             Q2_GRID(IQ2_TEST)          = TABLE_H197(ILINE, 2)
             UDS_GRID(IZ_TEST,IQ2_TEST) = TABLE_H197(ILINE, 3)
             G_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H197(ILINE, 4)
             C_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H197(ILINE, 5)

             UDS_GRID(IZ_TEST,IQ2_TEST)=UDS_GRID(IZ_TEST,IQ2_TEST)/6.
             ILINE=ILINE+1
            ENDDO
         ENDDO
c     table 20
      ELSEIF(IFIT.EQ.20) THEN
         DO IQ2_TEST = 1,NQ2MAX
            DO IZ_TEST  = 1,NZMAX
             BETA_GRID(IZ_TEST)         = TABLE_H1ZEUS97(ILINE, 1) 
             Q2_GRID(IQ2_TEST)          = TABLE_H1ZEUS97(ILINE, 2)
             UDS_GRID(IZ_TEST,IQ2_TEST) = TABLE_H1ZEUS97(ILINE, 3)
             G_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H1ZEUS97(ILINE, 4)
             C_GRID(IZ_TEST,IQ2_TEST)   = TABLE_H1ZEUS97(ILINE, 5)

             UDS_GRID(IZ_TEST,IQ2_TEST)=UDS_GRID(IZ_TEST,IQ2_TEST)/6.
             ILINE=ILINE+1
            ENDDO
         ENDDO
c     table 30
      ELSEIF(IFIT.EQ.30) THEN
         DO IQ2_TEST = 1,NQ2MAX
            DO IZ_TEST  = 1,NZMAX
             BETA_GRID(IZ_TEST)         = TABLE_ZEUS97(ILINE, 1) 
             Q2_GRID(IQ2_TEST)          = TABLE_ZEUS97(ILINE, 2)
             UDS_GRID(IZ_TEST,IQ2_TEST) = TABLE_ZEUS97(ILINE, 3)
             G_GRID(IZ_TEST,IQ2_TEST)   = TABLE_ZEUS97(ILINE, 4)
             C_GRID(IZ_TEST,IQ2_TEST)   = TABLE_ZEUS97(ILINE, 5)

             UDS_GRID(IZ_TEST,IQ2_TEST)=UDS_GRID(IZ_TEST,IQ2_TEST)/6.
             ILINE=ILINE+1
            ENDDO
         ENDDO
      ENDIF



      RETURN
      END
