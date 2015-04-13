      subroutine KMR2_Init(s_, gapsurv_, uscale_, delta_, q2cut_)
      include 'kmr2.inc'
      double precision s_,HWUALF,k, uscale_, gapsurv_, q2cut_
      external HWUALF
      integer alphaSorder,alphaSnfmax,delta_
      integer NN 
      save NN  


C MSTW 2008      
c      double precision distance,tolerance,
c     &     mCharm,mBottom,alphaSQ0,alphaSMZ
c      COMMON/mstwCommon/distance,tolerance,
c     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
c      CALL INITALPHAS(2,1.D0,1.D0,0.45077d0,
c     &     1.4d0,4.75d0,1.D10)
c      prefix = "External/kmrlumi2/Grids/mstw2008nlo" ! prefix for the grid files
c      iset = 0 ! Central PDFs
c      BottomMass = 4.6;
c      CharmMass = 1.42;

      s = s_
      gapsurv = gapsurv_
      uscale = uscale_
      delta = delta_
      q2cut = q2cut_

c      if(abs(sqrt(s)-2000d0).lt.100d0) then
c        gapsurv = 0.1d0
c      elseif(abs(sqrt(s)-14000d0).lt.1000) then
c        gapsurv = 0.03d0
c      else
c        gapsurv = 1d0
c      endif

      pi = atan(1d0)*4d0

      print*, "  I==========================I"
      print*, "  I  Initialising KMR model  I"
      print*, "  I==========================I"
      print*, "  I  sqrt(s)= ", sqrt(s)
      print*, "  I  gap survival= ", gapsurv
      print*, "  I  iset =", iset
      print*, "  I  scale =", uscale
      print*, "  I  delta =", delta
      print*, "  I  q2cut =", q2cut
      print*, "  I==========================I"
      
      end


      function KMR2_IntxPgg(x,kT)
      include 'kmr2.inc'
      integer Nf
      double precision KMR2_IntxPgg,x,kT
      double precision x2,x3,x4,C_A,I1,I2
      C_A = 3d0
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      if(kT .le. CharmMass) then
        Nf = 3;
      elseif(kT .le. BottomMass) then
        Nf = 4;
      else
        Nf = 5;
      endif
      I1 = - x4/4.0 + x3/3.0 - x2 - log(1-x)
      I2 = (0.5*(x-x2) + x3/3.0) !Inter of Pgq
      KMR2_IntxPgg = 2d0*C_A*I1 + Nf*I2
      return
      end
 

      double precision function ExHuME_AlphaS(Q2)
      include 'kmr2.inc'
      double precision Q2
      double precision ASConst, Freeze, LambdaQCD, Q
      ASConst = 12d0*3.14159/22d0/2d0
      LambdaQCD = 80d-3
      Freeze = sqrt(Q2cut) 
      Q = sqrt(Q2)
      if(Q>Freeze) then
        ExHuME_AlphaS = ASConst/log(Q/LambdaQCD)
      else
        ExHuME_AlphaS = ASConst/log(Freeze/LambdaQCD)
      endif
      return
      end
      
      function KMR2_alphas(Q2) 
      include 'kmr2.inc'
      double precision KMR2_alphas, Q2
      double precision ExHuME_alphas
      external ExHuME_alphas
      KMR2_alphas = ExHuME_ALPHAS(Q2)
      return
      end

      function KMR2_xg(x,Q2)
      include 'kmr2.inc'
      double precision KMR2_xg,x,Q2,q
      double precision GetOnePDF
      double precision upv,dnv,usea,dsea,str,chm,bot,glu
      integer mode
      external GetOnePDF
      q = sqrt(Q2)
C MRST 2002      
      mode = 1 !NLO
      call mrst2002(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      KMR2_xg = glu
C MSTW 2008
c      KMR2_xg = GetOnePDF(prefix,iset,x,Q,0)
      return
      end

      function KMR2_dT(lnkT2,par)
      include 'kmr2.inc'
      double precision KMR2_dT, lnkT2, par(10), kT2
      double precision I, mu, kT, AS
      double precision KMR2_alphas, KMR2_IntxPgg
      external KMR2_alphas, KMR2_IntxPgg
      kT2 = exp(lnkT2)
      kT = sqrt(kT2)
      mu = sqrt(par(1))
      if(delta .eq. 1) then
        I = KMR2_IntxPgg(1-kT/(mu+kT),kT) 
      elseif(delta .eq. 2) then
        I = KMR2_IntxPgg(1-kT/mu,kT) 
      else
        I = 0
      endif
      AS = KMR2_alphas(kT2)
      KMR2_dT = AS/(2d0*pi)*I
      return
      end
            
      function KMR2_T(QT2,mu2)
      include 'kmr2.inc'
      double precision KMR2_T, QT2, mu2
      double precision KMR2_integrate2, KMR2_dT, I, par(10)
      external KMR2_integrate2, KMR2_dT
      par(1) = mu2
      I = KMR2_integrate2(KMR2_dT, log(QT2), log(mu2), par, 50) 
      KMR2_T = dexp(-I)
      return
      end
     
      function KMR2_Rg(x,QT2)
      include 'kmr2.inc'
      double precision KMR2_Rg, x, QT2
      double precision q, lambda, x1, x2, g1, g2
      double precision KMR2_xg
      external KMR2_xg
      q = sqrt(QT2)
      if(x.lt.0.002 .and. q.lt.1.5625) then
        KMR2_Rg = 1.0
        return
      endif
      x1 = 0.995 * x;
      x2 = 1.005 * x;
      if(x2.gt.1.0) then
        x1 = 0.99 * x;
        x2 = x;
      endif
      g1 = KMR2_xg(x1,QT2)
      g2 = KMR2_xg(x2,QT2)
      lambda = (g1 - g2) * 100.0 / g1;
      if(lambda.gt.1.0) lambda = 1.0;
      if(lambda.lt.0.0) lambda = 0.0;
      KMR2_Rg = 1.0 + lambda * (0.82 + 0.56 * lambda);
      return
      end      

      function KMR2_f(x,xp,QT2,mu2)
      include 'kmr2.inc'
      double precision KMR2_f,x,xp,QT2,mu2
      double precision sqrtT, grad, ddT,d,mu,q
      double precision xglu
      double precision KMR2_Rg, KMR2_T, KMR2_xg,KMR2_IntxPgg,
     & KMR2_alphas, KMR2_dT
      external KMR2_Rg, KMR2_T,KMR2_xg,KMR2_IntxPgg,
     & KMR2_alphas, KMR2_dT
      d = 1e-3
      mu = sqrt(mu2)
      q = sqrt(QT2)
      sqrtT = sqrt(KMR2_T(QT2,mu2))
c     ddT = KMR2_alphas(QT2)/(4.0*pi)*KMR2_IntxPgg(mu/(q+mu),q)
      ddT = -3d0/2d0/pi*KMR2_alphas(QT2)*log(1d0/3d0*q/(q+mu))
      grad = (KMR2_xg(x,QT2*(1d0+d)) 
     &       - KMR2_xg(x,QT2*(1d0-d)))/(2d0*d)
      xglu = KMR2_xg(x,QT2)
      KMR2_f = KMR2_Rg(x,QT2)*sqrtT*(ddT*xglu + grad)
      return
      end
 
      function KMR2_dL(lnQT2,par)
      include 'kmr2.inc'
      double precision KMR2_dL, x1, x2, Qt2, mu2, par(10)
      double precision x1p, x2p, lnQT2
      double precision KMR2_f
      external KMR2_f
      QT2 = exp(lnQT2)
      x1 = par(1)
      x2 = par(2)
      mu2 = par(3)
      x1p = 1e-5
      x2p = 1e-5
      KMR2_dL = KMR2_f(x1,x1p,QT2,mu2)*KMR2_f(x2,x2p,QT2,mu2)
      KMR2_dL = KMR2_dL/QT2
      return
      end
      
      function KMR2_L(M2,y,mu2)
      include 'kmr2.inc'
      double precision Q2min,Q2max,b
      double precision KMR2_L, x1, x2, mu2, M2,
     & y, sss, par(10)
      double precision I, KMR2_integrate1, KMR2_dL
      external KMR2_integrate1, KMR2_dL
      sss = sqrt(M2/s)
      x1 = dexp(y)*sss
      x2 = dexp(-y)*sss
      Q2min = q2cut
      Q2max = 1d6
      par(1) = x1
      par(2) = x2
      par(3) = mu2
      I = KMR2_integrate1(KMR2_dL,log(Q2min),log(Q2max),par,100)
      KMR2_L = (pi/(8.0) * I)**2
      return
      end 

      function KMR2_GET_LUMI(x1,x2)
      include 'kmr2.inc'
      double precision KMR2_GET_LUMI,x1,x2
      double precision M2,y,mu2
      double precision KMR2_L
      external KMR2_L
      integer N
      save N
      N = N+1
      M2 = s*x1*x2
      y = 0.5*log(x1/x2)
c     mu2 = M2*0.618**2 ! (this is used in ExHuME)
      mu2 = M2*uscale**2 
      KMR2_GET_LUMI = KMR2_L(M2,y,mu2)/16d0 ! 16 = squared b=4.0 from
                                              ! t integration 
                            
      KMR2_GET_LUMI = gapsurv*KMR2_GET_LUMI         
      if(KMR2_GET_LUMI/KMR2_GET_LUMI.ne.1d0) KMR2_GET_LUMI = 0d0
      RETURN
      end


C Two functions due to some strange things that happen when there is
C only one function and there is one integrantion inside another

      
      function KMR2_integrate1(func,a,b,par,N)
      implicit none
      double precision  KMR2_integrate1 
      integer N
      double precision func,a,b,par(10)
      integer i
      double precision x,y,w,d,integ
      if(mod(N,2).eq.1) N = N+1
      d = (b-a)/N

      x = a
      integ = 0
      do i=0,N
        if(i.eq.0 .or. i.eq.N) then
           w = 1./3.
        elseif(mod(i,2).eq.0) then
           w = 2./3.
        else
           w = 4./3.
        endif
        y=func(x,par)
        integ = integ + y*w
        x = x+d
      enddo
      KMR2_integrate1 = integ*d
      return
      end
      

      function KMR2_integrate2(func,a,b,par,N)
      implicit none
      double precision  KMR2_integrate2
      integer N
      double precision func,a,b,par(10)
      integer i
      double precision x,y,w,d,integ
      if(mod(N,2).eq.1) N = N+1
      d = (b-a)/N

      x = a
      integ = 0
      do i=0,N
        if(i.eq.0 .or. i.eq.N) then
           w = 1./3.
        elseif(mod(i,2).eq.0) then
           w = 2./3.
        else
           w = 4./3.
        endif
        y=func(x,par)
        integ = integ + y*w
        x = x+d
      enddo
      KMR2_integrate2 = integ*d
      return
      end
