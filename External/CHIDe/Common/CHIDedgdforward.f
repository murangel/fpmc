       function CHIDedgdforward(xeff,ks,iglu)

	implicit double precision(a-h,k-m,o-z)
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0
        COMMON / INTINIP / IINIP

	data mus2/ 0.5625d0 /


        pi=atan(1d0)*4d0
        IINIP = 0

C     step of differentiation
	kstep=0.0001E0
      
C     ansatz dependent parameters
	if(iglu.eq.1)then  ! GRV2000
           ks0=0.895E0
           kss_hard = sqrt(1. + 0.0018E0*(log(1./xeff))**4)
           kss_soft = 3E0 
           mus1 = 0.01E0
           asoft = 2E0
	     npow = 5
	     nrr = 1
        
	        if(ks.le.ks0) then
		     call GBf(ks0,FB,mus2)
			 call GG(xeff,ks0,FG)
	           CB=FG/FB
		     call GBf(ks,FB,mus2)
			 W=FB*CB
	        else
		     call GG(xeff,ks,W)
			endif
			call GBf(ks,U,mus1)
	  	  
	  elseif(iglu.eq.2)then  ! GRV2002
           ks0=1.45d0
           kss_soft = 1.0E0
           mus1 = 0.01E0
           asoft = 2.66E0
           ksscon1 = 0.4E0
           ksscon2 = 0.245E0
           kss_hard = ksscon1 + ksscon2*log(1./xeff)
           npow = 11
	     nrr = 4

           if(ks.le.ks0) then
              call GBf(ks0,FB,mus2)
              call GG(xeff,ks0,FG)
                CB=FG/FB
              call GBf(ks,FB,mus2)
                W=FB*CB

           else
              call GG(xeff,ks,w)

           endif
           call GBf(ks,U,mus1)
	  	  
	  elseif(iglu.eq.3)then	 ! formula2002

	         ks0=1.4e0
		 kss_soft = 1.0E0
		 mus1 = 0.01E0
		 asoft = 2.6E0
                 ksscon1 = 0.0E0
		 ksscon2 = 0.31E0
		 kss_hard = ksscon1 + ksscon2*log(1./xeff)
		 npow = 11
		 nrr = 4

           if(ks.le.ks0) then
              call GBf(ks0,FB,mus2)
              call GGtype3f(xeff,ks0,FG)
              CB=FG/FB
              
              call GBf(ks,FB,mus2)
              W=FB*CB
           else
              call GGtype3f(xeff,ks,W)
           endif
           call GBf(ks,U,mus1)

	  elseif(iglu.eq.4)then  ! GRV2005

           ks0=0.895e0
           kss_hard = sqrt(0.02E0 + 0.00067E0*(log(1./xeff))**4)
 	   kss_soft = 2.0E0
           asoft = 1.6E0
           mus1 = 0.02E0
	     npow = 5
	     nrr = 3

	        if(ks.le.ks0) then
		     call GBf(ks0,FB,mus2)
                     call GG(xeff,ks0,FG)
                          CB=FG/FB
		     call GBf(ks,FB,mus2)
			  W=FB*CB
	        else
		     call GG(xeff,ks,W)

                  endif
			call GBf(ks,U,mus1)

        else
           print*,'Wrong iglu! iglu = ',iglu
           stop
        end if   

        
        U = U*(1.-xeff)**npow
        dgd_soft=U*asoft*kss_soft**nrr/(kss_soft**nrr+ks**nrr)

        dgd_hard=W*ks**nrr/(ks**nrr+kss_hard**nrr)
        CHIDedgdforward = dgd_soft  + dgd_hard

        return
      end function 

! ===== Dipole-like formfactor ======= ! 
       double precision function ffactorf(k2)
       implicit double precision(a-h,k-z)

C     proton scale, related to the proton radius
       mps=1.0 ! GeV^2
       ffactorf = 1d0/(1d0+k2/mps)**2
        return
       end function


! ===== Born unintegrated DGD ====== !
       subroutine GBf(ks,FB,muss)
	implicit double precision(a-h,k-z)
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0

	FB=4.*alpha_s(ks)/pi*ks**2/(ks+muss)**2*(1d0-ffactorf(3*ks))

	return
      end subroutine


      subroutine GG(xeff,ks,FG)
       implicit double precision(a-h,k-m,o-z)
       dimension QS1(27),QS1mid(26)
! ---- GRV table -- 27 step points
       DATA QS1 / 0.8E0,
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4, 2.2E4, 4.6E4,
     6           1.0E5, 2.2E5, 4.6E5,
     7           1.E6 /  

          DO iq=1,26
             QS1mid(iq) = DSQRT(QS1(iq)*QS1(iq+1))
          END DO
          
          DO iq=1,26
             IF (QS1mid(iq).GT.ks)goto 99
          END DO

 99       a = QS1mid(iq-1)
          b = QS1mid(iq)

	  call GG1(xeff,a,FGa)
	  call GG1(xeff,b,FGb)
       FG = fga + (fgb-fga)/DLOG(b/a)*DLOG(ks/a)
       return
       end


      subroutine GG1(xeff,ks,FG1)
	implicit double precision(a-h,k-m,o-z)
        COMMON/INTINIP/IINIP 
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0

	ksp=ks+kstep

          call GRV98PA(1,xeff,ks, UV, DV, US, DS, SS, GL)
          call GRV98PA(1,xeff,ksp, UV, DV, US, DS, SS, GLp)

	FG1=max((GLp-GL)/kstep*ks,0d0)

	return
      end subroutine


      subroutine GGtype3f(xeff,ks,FG)
       implicit double precision(a-h,k-m,o-z)

       FG = 0.2450D0*(DLOG(ks/0.04d0))**(0.34d0-6d0*DSQRT(xeff))
     .   /xeff**0.40d0
        return
       end


      real*8 function valent(x,Q2)
	implicit double precision(a-h,k-m,o-z)
        COMMON/INTINIP/IINIP 
           call GRV98PA(1,x,Q2, UV, DV, US, DS, SS, GL)
        valent = (4d0/9d0*UV + 1d0/9d0*DV)
       return
      end function


      real*8 function alpha_s(Q2) ! - leading order
	implicit double precision(a-h,k-m,o-z)
        pi=atan(1d0)*4d0
        lqcd2=0.04
        Q2s = max(Q2,0.22d0)
        a_s = 4d0*pi/9d0/log(Q2s/LQCD2)
        a_s = min(0.82d0,a_s)
        alpha_s = a_s
       return 
      end function



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                   *
*     G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S      *
*                                                                   *
*                          1998 UPDATE                              *
*                                                                   *
*                  For a detailed explanation see                   *
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *
*                                                                   *
*   This package contains subroutines returning the light-parton    *
*   distributions in NLO (for the MSbar and DIS schemes) and LO;    *
*   the respective light-parton, charm, and bottom contributions    *
*   to F2(electromagnetic); and the scale dependence of alpha_s.    *
*                                                                   *
*   The parton densities and F2 values are calculated from inter-   *
*   polation grids covering the regions                             *
*         Q^2/GeV^2  between   0.8   and  1.E6 ( 1.E4 for F2 )      *
*            x       between  1.E-9  and   1.                       *
*   Any call outside these regions stops the program execution.     *
*                                                                   *
*   At Q^2 = MZ^2, alpha_s reads  0.114 (0.125) in NLO (LO); the    *
*   heavy quark thresholds, QH^2 = mh^2, in the beta function are   *
*            mc = 1.4 GeV,  mb = 4.5 GeV,  mt = 175 GeV.            *
*   Note that the NLO alpha_s running is different from GRV(94).    *
*                                                                   *
*    Questions, comments etc to:  avogt@physik.uni-wuerzburg.de     *
*                                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
      SUBROUTINE GRV98PA (ISET, X, Q2, UV, DV, US, DS, SS, GL)
*********************************************************************
*                                                                   *
*   THE PARTON ROUTINE.                                             *
*                                     __                            *
*   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         *
*                                                                   *
*   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       *
*            DS = d(bar),  SS = s = s(bar),  GL = gluon.            *
*            Always x times the distribution is returned.           *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIP / IINIP , and the     *
*            integer variable  IINIP  has always to be zero when    *
*            GRV98PA is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   *
*            (1+1809 lines with 6 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1),
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ)
      CHARACTER*80 LINE
      COMMON / INTINIP / IINIP
      SAVE XUVF, XDVF, XDEF, XUDF, XSF, XGF, NA, ARRF
*
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8E0,
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4, 2.2E4, 4.6E4,
     6           1.0E5, 2.2E5, 4.6E5,
     7           1.E6 /
       DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9,
     A           1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8,
     B           1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7,
     C           1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1           1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2           1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3           1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4           1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*
*...CHECK OF X AND Q2 VALUES :
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
C         WRITE (6,91)
C  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C         WRITE (*,*) X
         print*,'PARTON INTERPOLATION: X OUT OF RANGE'
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
C         WRITE(6,92)
C  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
C         WRITE (*,*) Q2
         print*,'PARTON INTERPOLATION: Q2 OUT OF RANGE'
         STOP
      ENDIF
      IF (IINIP .NE. 0) GOTO 16
*
*...INITIALIZATION, IF REQUIRED :
*
*    SELECTION AND READING OF THE GRID :
*    (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
        OPEN (11,FILE='grv98lo.grid',STATUS='old')   ! 7.332E-05
      ELSE IF (ISET .EQ. 2) THEN
        OPEN (11,FILE='grv98nlm.grid',STATUS='old')  ! 1.015E-04
      ELSE IF (ISET .EQ. 3) THEN
        OPEN (11,FILE='grv98nld.grid',STATUS='old')  ! 1.238E-04
      ELSE
C        WRITE(6,93)
C  93    FORMAT (2X,'NO OR INVALID PARTON SET CHOICE')
        print*,'NO OR INVALID PARTON SET CHOICE'
        STOP
      END IF
      IINIP = 1
      READ(11,89,end=100) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1
      DO 15 N = 1, NQ
      READ(11,90,end=100) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M),
     1            PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M)
  90  FORMAT (6(1PE10.3))
  15  CONTINUE
 100  CLOSE(11)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0V = XB(IX)**0.5
        XB0S = XB(IX)**(-0.2)
        XB1 = 1.-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
        XDEF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V)
        XUDF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
  20  CONTINUE
        XUVF(NX,IQ) = 0.E0
        XDVF(NX,IQ) = 0.E0
        XDEF(NX,IQ) = 0.E0
        XUDF(NX,IQ) = 0.E0
        XSF(NX,IQ)  = 0.E0
        XGF(NX,IQ)  = 0.E0
  10  CONTINUE
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
*
*...CONTINUATION, IF INITIALIZATION WAS DONE PREVIOUSLY.
*
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      X1 = 1.- X
      XV = X**0.5
      XS = X**(-0.2)
      UV = FFINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = FFINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = FFINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = FFINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      US = 0.5 * (UD - DE)
      DS = 0.5 * (UD + DE)
      SS = FFINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = FFINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS
*
 60   RETURN
      END subroutine
*
*
*
*
      FUNCTION FFINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FFINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FFINT=FFINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END function



      FUNCTION ALPHAS (Q2, NAORD)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE.                                            *
*                                                                   *
*   INPUT :  Q2    =  scale in GeV**2  (not too low, of course);    *
*            NAORD =  1 (LO),  2 (NLO).                             *
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) for use with the GRV(98) partons.      *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A - Z)
      INTEGER NF, K, I, NAORD
      DIMENSION LAMBDAL (3:6),  LAMBDAN (3:6), Q2THR (3)
*
*...HEAVY QUARK THRESHOLDS AND LAMBDA VALUES :
      DATA Q2THR   /  1.960,  20.25,  30625. /
      DATA LAMBDAL / 0.2041, 0.1750, 0.1320, 0.0665 /
      DATA LAMBDAN / 0.2994, 0.2460, 0.1677, 0.0678 /
*
*...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
      NF = 3
      DO 10 K = 1, 3
      IF (Q2 .GT. Q2THR (K)) THEN
         NF = NF + 1
      ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
         ALP  = 1./(B0 * DLOG (Q2/LAM2))
         GO TO 1
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
C         WRITE (6,91)
C  91     FORMAT ('INVALID CHOICE FOR ORDER IN ALPHA_S')
         print*,'INVALID CHOICE FOR ORDER IN ALPHA_S' 
         STOP
       END IF
*
*...START VALUE FOR NLO ITERATION :
       LQ2 = DLOG (Q2 / LAM2)
       ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)
*
*...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE :
       DO 2 I = 1, 6
       XL  = DLOG (1./(B0*ALP) + B10)
       XLP = DLOG (1./(B0*ALP*1.01) + B10)
       XLM = DLOG (1./(B0*ALP*0.99) + B10)
       Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
       Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP
     1       + 1./ (B0*ALP*0.99) - B10 * XLP) / (0.02D0*ALP)
       ALP = ALP - Y/Y1
  2    CONTINUE
*
*...OUTPUT :
  1    ALPHAS = ALP
       RETURN
       END function
