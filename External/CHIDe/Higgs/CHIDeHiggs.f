C     ====================================================================
C     Analytic part of the integrand
      subroutine CHIDeHiggs(sigma,b1,a3,k1,k3,k,kp)
        implicit none
      integer sudaf,grid,ihist,formfac,event,scale,evolv,interpol,Vgrid
      integer exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      integer nf,iglu
      double precision sigma
      double precision k(2),kp(2),k1(2),k3(2)
      double precision mskk1(2),skpk1(2),skk3(2),mskpk3(2),mk(2),mkp(2)
      double precision skk1(2),skpk3(2),mskk3(2)
      double precision Deno,Denop,Numow1,Numopw1,Numow2,Numopw2
      double precision N1,w1,N2,w2,a,sa,NVertex
      double precision CHIDedot,CHIDedotsum,phi,as,alphas,CHIDedothalf
      double precision pi,s,kmax,ncolor,nb,gq,Gf,mt,mp
      double precision mh
      double precision CHIDedgd2008,dgdforward,prefac
      double precision Sphi,Sphib,Sudaphi
      double precision ephi,phip,phipb
      double precision Ifact,fact,Cf
      double precision Dsuda,Rg
      double precision m,ao,mprime,Qp,Cff,bo,mb,mbprime
      double precision b1,a3
      double precision c,ael
      double precision dgdpol,Dsudapol
      double precision x,xp

      double precision mmu2,S2,regge,a3t,a1t,factq,factps

      common/const/pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      common/higgs/mh
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      common/flavor/nf
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/sudaD/mmu2,s2

 

C     Sum of vectors in phi and dphi
      mskk1(1)=-k(1)-k1(1)
      mskk1(2)=-k(2)-k1(2)
      
      skk1(1)=k(1)+k1(1)
      skk1(2)=k(2)+k1(2)
      skpk1(1)=kp(1)+k1(1)
      skpk1(2)=kp(2)+k1(2)

      skk3(1)=k(1)+k3(1)
      skk3(2)=k(2)+k3(2)
      skpk3(1)=kp(1)+k3(1)
      skpk3(2)=kp(2)+k3(2)
      mskpk3(1)=-kp(1)-k3(1)
      mskpk3(2)=-kp(2)-k3(2)
      mskk3(1)=-k(1)-k3(1)
      mskk3(2)=-k(2)-k3(2)

      mk(1)=-k(1)
      mk(2)=-k(2)
      mkp(1)=-kp(1)
      mkp(2)=-kp(2)


C  --------------------------   Higgs vertex  ----------------------------------   C
C    JR: as(mh**2) is always fixed to 0.1 as a prescription
      w1=sqrt(sqrt(2.)*Gf)*as(mh**2)*mh**2/3./pi
      w2=w1
      a=mh**2/mt**2
      N1=6.*(1.+(1.-4./a)*(atan(sqrt(a/(4.-a))))**2)/a
      N2=-12.*(5.-(1.+4./a)*(atan(sqrt(a/(4.-a))))**2
     &    -4.*atan(sqrt(a/(4.-a)))/sqrt(a/(4.-a)))/a
      

C   ----------------------------------------------------------------------------   C
C   OUR MODEL  k~ki
C   ----------------------------------------------------------------------------   C 
      if(Model.EQ.1.or.Model.EQ.2)then
C     Vertex and Numerator
      Numow1=8.*CHIDedotsum(k,k1,k,k3)
      Numopw1=8.*CHIDedotsum(kp,k1,kp,k3)
      Numow2=8.*CHIDedotsum(k,k1,k,k1)*CHIDedotsum(k,k3,k,k3)
      Numopw2=8.*CHIDedotsum(kp,k1,kp,k1)*CHIDedotsum(kp,k3,kp,k3)

      if(fullvertex.EQ.1)then
       Nvertex=(Numow1*w1*N1/mh**2+Numow2*w2*N2/mh**4)
     &       *(Numopw1*w1*N1/mh**2+Numopw2*w2*N2/mh**4)
      elseif(fullvertex.EQ.2)then
       Nvertex=Numow1*Numopw1*w1**2*N1**2/mh**4
      endif 
      
C     Propagator: denominator!
      Deno=CHIDedot(k,k)*CHIDedotsum(k,k1,k,k1)*CHIDedotsum(k,k3,k,k3)
      Denop=CHIDedot(kp,kp)*CHIDedotsum(kp,k1,kp,k1)
     & *CHIDedotsum(kp,k3,kp,k3)
     
C   -------------------------  Impact Factor  ---------------------------------   C
      if(formfac.EQ.1)then
         phip=ephi(k,mskk1)*ephi(skpk1,mkp)
         phipb=ephi(mk,skk3)*ephi(mskpk3,kp)
       
      elseif(formfac.EQ.2)then

         if(event.EQ.1)then
      phip=CHIDedgd2008(0.41d0*b1,CHIDedot(k,k),
     & CHIDedotsum(k,k1,k,k1),CHIDedot(k1,k1))
     &*CHIDedgd2008(0.41d0*b1,CHIDedot(kp,kp),
     & CHIDedotsum(kp,k1,kp,k1),CHIDedot(k1,k1))
      phipb=CHIDedgd2008(0.41d0*abs(a3),
     & CHIDedot(k,k),CHIDedotsum(k,k3,k,k3)
     &,CHIDedot(k3,k3))*CHIDedgd2008(0.41d0*abs(a3),CHIDedot(kp,kp)
     &,CHIDedotsum(kp,k3,kp,k3),CHIDedot(k3,k3))

        elseif(event.EQ.2)then
C     Exact definition of dgdforward -> Very slow
C       phip=dgdforward(0.41d0*b1,CHIDedothalf(k,k1,k,k1),iglu)
C     &*prefac(k,mskk1,CHIDedot(k1,k1),b1)
C     &*dgdforward(0.41d0*b1,CHIDedothalf(kp,k1,kp,k1),iglu)
C     &*prefac(skpk1,kp,CHIDedot(k1,k1),b1)
C       phipb=dgdforward(0.41d0*abs(a3),CHIDedothalf(k,k3,k,k3),iglu)
C     &*prefac(mk,skk3,CHIDedot(k3,k3),abs(a3)) 
C     &*dgdforward(0.41d0*abs(a3),CHIDedothalf(kp,k3,kp,k3),iglu)
C     &*prefac(mskpk3,kp,CHIDedot(k3,k3),abs(a3))
C    Interpolation of dgdforward
       phip=dgdpol(b1,CHIDedothalf(k,k1,k,k1))
     &      *prefac(k,mskk1,CHIDedot(k1,k1),b1)
     &      *dgdpol(b1,CHIDedothalf(kp,k1,kp,k1))
     &      *prefac(skpk1,mkp,CHIDedot(k1,k1),b1)
      phipb=dgdpol(abs(a3),CHIDedothalf(k,k3,k,k3))
     &      *prefac(mk,skk3,CHIDedot(k3,k3),abs(a3)) 
     &      *dgdpol(abs(a3),CHIDedothalf(kp,k3,kp,k3))
     &      *prefac(mskpk3,kp,CHIDedot(k3,k3),abs(a3))
C      print*,'prefact',prefac(k,mskk1,CHIDedot(k1,k1),b1)

c     phip=1.
c     phipb=1.
        endif

       endif

C   ----------------------  Sudakov Form Factor  -------------------------------   C 
C     Durham Sudakov factor
C     Durham Sudakov factor: Delta=qt/(qt+mu) and mu=0.62mh
C     Related to HiggsDurhamsudatab.d 
      if(sudaf.EQ.1)then
       if(interpol.EQ.2)then
           Sphi=Dsuda((0.62d0*mh)**2/x,CHIDedot(k,k)/xp)
           Sphib=Dsuda((0.62d0*mh)**2/x,CHIDedot(kp,kp)/xp)
       elseif(interpol.EQ.1)then
           Sphi=Dsudapol((0.62d0*mh)**2/x,CHIDedot(k,k)/xp)
           Sphib=Dsudapol((0.62d0*mh)**2/x,CHIDedot(kp,kp)/xp)
       endif  

C     DLA Sudakov form factor
      elseif(sudaf.EQ.2)then
         if(scale.EQ.1)then
       Sphi=Dsuda(mh**2/x,CHIDedot(k,k)/xp)
       Sphib=Dsuda(mh**2/x,CHIDedot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(mh**2/x,CHIDedot(skk1,skk1)/xp))
     &      *sqrt(Dsuda(mh**2/x,CHIDedot(skk3,skk3)/xp))
           Sphib=sqrt(Dsuda(mh**2/x,CHIDedot(skpk1,skpk1)/xp))
     &      *sqrt(Dsuda(mh**2/x,CHIDedot(skpk3,skpk3)/xp))
         endif
                 
C     Naive Sudakov Factor      
      elseif(sudaf.EQ.3)then

       Sphi=Dexp(-3./pi/2.*as(mh**2/x)*
     & ((log(CHIDedot(skk3,skk3)/mh**2/xp)
     &*log(CHIDedot(skk1,skk1)/mh**2/xp))))
       Sphib=Dexp(-3./pi/2.*as(mh**2/x)*
     & ((log(CHIDedot(skpk3,skpk3)/mh**2/xp)
     &*log(CHIDedot(skpk1,skpk1)/mh**2/xp))))

C     Jets Sudakov Factor
      elseif(sudaf.EQ.4)then
         x=2.
         xp=1.
       if(interpol.EQ.2)then
         if(scale.EQ.1)then
           Sphi=Dsuda(mh**2/x,CHIDedot(k,k)/xp)
           Sphib=Dsuda(mh**2/x,CHIDedot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(mh**2/x,CHIDedot(skk1,skk1)/xp))
     &      *sqrt(Dsuda(mh**2/x,CHIDedot(skk3,skk3)/xp))
           Sphib=sqrt(Dsuda(mh**2/x,CHIDedot(skpk1,skpk1)/xp))
     &      *sqrt(Dsuda(mh**2/x,CHIDedot(skpk3,skpk3)/xp))
         endif

       elseif(interpol.EQ.1)then
         if(scale.EQ.1)then
           Sphi=Dsudapol(mh**2/x,CHIDedot(k,k)/xp)
           Sphib=Dsudapol(mh**2/x,CHIDedot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsudapol(mh**2/x,CHIDedot(skk1,skk1)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(skk3,skk3)/xp))
           Sphib=sqrt(Dsudapol(mh**2/x,CHIDedot(skpk1,skpk1)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(skpk3,skpk3)/xp))
         endif
       endif  

C     No Sudakov Factor
      elseif(sudaf.EQ.5)then
       Sphi=1.d0
       Sphib=1.d0 

C    Forshaw Sudakov factor: ~Durham with Delta=qt/mu and mu=mh
C    Related to Higgssudatab.d 
      elseif(sudaf.EQ.6)then
       if(scale.EQ.1)then
        Sphi=Dsudapol(mh**2/x,CHIDedot(k,k)/xp)
        Sphib=Dsudapol(mh**2/x,CHIDedot(kp,kp)/xp)
       elseif(scale.EQ.2)then
        Sphi=sqrt(Dsudapol(mh**2/x,CHIDedot(skk1,skk1)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(skk3,skk3)/xp))
        Sphib=sqrt(Dsudapol(mh**2/x,CHIDedot(skpk1,skpk1)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(skpk3,skpk3)/xp))
       endif  


      endif

C   ----------------------------------------------------------------------------   C
C   DURHAM MODEL  k²>>>ki  
C   ----------------------------------------------------------------------------   C
      elseif(Model.EQ.3)then
C     Vertex and Numerator
      Numow1=8.*CHIDedot(k,k)
      Numopw1=8.*CHIDedot(kp,kp)
      Numow2=8.*CHIDedot(k,k)*CHIDedot(k,k)
      Numopw2=8.*CHIDedot(kp,kp)*CHIDedot(kp,kp)

      if(fullvertex.EQ.1)then
       Nvertex=(Numow1*w1*N1/mh**2+Numow2*w2*N2/mh**4)
     &       *(Numopw1*w1*N1/mh**2+Numopw2*w2*N2/mh**4)
      elseif(fullvertex.EQ.2)then
       Nvertex=Numow1*Numopw1*w1**2*N1**2/mh**4
      endif      


C     Propagator: denominator!
      Deno=CHIDedot(k,k)*CHIDedot(k,k)*CHIDedot(k,k)
      Denop=CHIDedot(kp,kp)*CHIDedot(kp,kp)*CHIDedot(kp,kp)

C   ----------------------- Impact and Sudakov  --------------------------------   C
C    Impact factor
      phip=dgdpol(b1,CHIDedot(k,k))*prefac(k,k,CHIDedot(k1,k1),b1)
     &*dgdpol(b1,CHIDedot(kp,kp))*prefac(kp,kp,CHIDedot(k1,k1),b1)
      phipb=dgdpol(abs(a3),CHIDedot(k,k))*
     & prefac(mk,mk,CHIDedot(k3,k3),abs(a3)) 
     &*dgdpol(abs(a3),CHIDedot(kp,kp))*
     & prefac(mkp,mkp,CHIDedot(k3,k3),abs(a3))

C    Sudakov form factor
       if(moode.EQ.1)then
           Sphi=sqrt(Dsudapol(mh**2/x,CHIDedot(k,k)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(k,k)/xp))
           Sphib=sqrt(Dsudapol(mh**2/x,CHIDedot(kp,kp)/xp))
     &      *sqrt(Dsudapol(mh**2/x,CHIDedot(kp,kp)/xp))
       elseif(moode.EQ.2)then   
C           Rg=1.2d0
           Rg=1.d0
           Sphi=Rg**2*Dsudapol((0.62d0*mh)**2,CHIDedot(k,k))
           Sphib=Rg**2*Dsudapol((0.62d0*mh)**2,CHIDedot(kp,kp))
       endif
      endif   

C   -----------------------------------------------------------------------------  C 

C     Coefficients: 
      gq=sqrt(4.*pi*ael)
      Cf=(2./9.)**2
      fact=Cf*gq**8/4./(2.*pi)**4
      Ifact=pi**2/gq**2 

C     Sudakov and Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif  
      Sudaphi=Sphi*Sphib

C     Coefficient and Phase space
      factps=1./8./(2.*pi)**5/b1/(1.-b1)/(1.+a3)
      factq=81.             !(3x3)² from impact

C     Regge factor
      a1t=0.09d0-0.3d0*CHIDedot(k1,k1)/(1.-b1)
      a3t=0.09d0-0.3d0*CHIDedot(k3,k3)/(1.+a3)

      regge=((1./abs(a3))**a1t*(1./b1)**a3t)**2

C     ----------- Cross section ------------------------------------------------   C
      sigma=factps*factq*nb*fact*NVertex/Deno/Denop
      sigma=s2*sigma

      if(formfac.EQ.1)then
       sigma=sigma*regge
      endif

      sigma=sigma*Sudaphi*phi
  
C      print*,"Vertex and numerator=",NVertex
C      print*,'fact=',fact
C      print*,"Num et deno=",numow1,numopw1,Deno,Denop
C      print*,numow1*numopw1/Deno/Denop
C      print*,"phi=",phi,phip,phipb
C      print*,"Sudaphi=",sudaphi,Sphi,Sphib
C      print*,'In suda: mu=',(0.62d0*mh)**2/x,' l=',CHIDedot(k,k)/xp
C      print*,'sigma=',sigma
C      print*,'******************************'

      end
C   ============================================================================   C
C   SUBROUTINES and FUNCTIONS                                                      C
C   Several functions and definitions                                              C  
C          - Alpha strong                                                          C
C          - LCWF impact factor, elastic and diffractif                            C
C          - UgD forward and prefactor                                             C
C          - Sudakov form factor                                                   C
C          - Interpolation bilinear for UgD (dgdpol) and Sudakov (sudapol)         C
C          - Addition, substraction and multiplication of 4-vectors                C
C   ============================================================================   C

C   ----------------------------------------------------------------------------   C
C   SUBROUTINE ALPHA STRONG
C   Compute as in function of Q²:
C      -> choose a Q²
C      -> compute the corresponding lambda(nf)
C      -> compute as
C   ----------------------------------------------------------------------------   C
      function as(q) 
       implicit none 
       integer i,nf 
      double precision as
      double precision m(2:6)
      double precision smu,smc,sms,smb,smt 
      double precision lambda(2:6),la
      double precision q
      double precision pi,s,kmax,ncolor,gq,nb,Gf,mt,mp

      common/const/pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      common/flavor/nf
            
C     Constituent quarks masses and l(5) en GeV from particle data group 
      m(2)=0.312d0
      m(3)=0.45d0
      m(4)=1.5d0
      m(5)=4.5d0
      m(6)=173.4d0
      smu=m(2)**2
      sms=m(3)**2
      smc=m(4)**2
      smb=m(5)**2
      smt=m(6)**2
      lambda(5)=0.2d0
 
      do i=1,3 
      lambda(5-i)=lambda(6-i)*(m(6-i)
     &/lambda(6-i))**(2./(33.-2.*(5-i)))
      end do
      lambda(6)=lambda(5)**(23./21.)/(m(6)**(2./21.))
      continue

      if(q.GT.smt)then 
        nf=6
      elseif(q.GT.smb)then      
        nf=5
      elseif(q.GT.smc)then 
        nf=4
      elseif(q.GT.sms)then
        nf=3 
      elseif(q.GT.smu)then
        nf=2  
      else
      goto 666
      end if
      continue

      la=lambda(nf)
      as=12.*pi/(33.-2.*nf)/log(q/la**2)

      if(as.LT.0.or.as.GT.0.7d0)as=0.7d0
      
666   end    


C   ----------------------------------------------------------------------------   C
C   LCWF IMPACT FACTOR 
C   Elastic, diffractif, function f(x)~ effective slope B for the proton  
C   ----------------------------------------------------------------------------   C  
C   Elastic Impact factor
C   phi: Chosen from Jr and Hernandez paper
      function ephi(p1,p2)
	implicit none
      integer iglu
      double precision ephi,CHIDedot,CHIDedotsum,f
      double precision p1(2),p2(2),x1,x2,c,ael,gq

      common/Ifactor/c,ael,iglu         
         
      x1=CHIDedotsum(p1,p2,p1,p2)
      x2=CHIDedot(p1,p1)+CHIDedot(p2,p2)-c*CHIDedot(p1,p2)

      ephi=f(x1)-f(x2)
C      if(ephi.LT.10E-8)then
C       ephi=0.
C      endif 

      end
C   ----------------------------------------------------------------------------   C
C   Diffractif Impact factor
C   dphi from C&H p.479
      function dphi(p1,p2,p3,p4)
       implicit none
      double precision p1(2),p2(2),p3(2),p4(2),dphi
      double precision sp1p2(2),sp1p4(2),sp1p3(2)
      double precision f,E3,E2
      double precision x1,x2,x3,x4,x12,x14,x13
      double precision y12,y34,y23,y13,y14,y24
 
       x1=E2(p1)
       x2=E2(p2)
       x3=E2(p3)
       x4=E2(p4)

       sp1p2(1)=p1(1)+p2(1)
       sp1p2(2)=p1(2)+p2(2)
       sp1p4(1)=p1(1)+p4(1)
       sp1p4(2)=p1(2)+p4(2)
       sp1p3(1)=p1(1)+p3(1)
       sp1p3(2)=p1(2)+p3(2)

       x12=E2(sp1p2)
       x14=E2(sp1p4)
       x13=E2(sp1p3)

       y12=E3(p1,p2)
       y34=E3(p3,p4)
       y23=E3(p2,p3)
       y13=E3(p1,p3)
       y14=E3(p1,p4)
       y24=E3(p4,p2)

       dphi=1.d0-f(x1)-f(x2)-f(x3)-f(x4)
     &+2.d0*f(x12)+0.5d0*f(x14)+0.5d0*f(x13)
     &-f(y12)-f(y34)+0.5d0*f(y23)+0.5d0*f(y13)+0.5d0*f(y14)+0.5d0*f(y24)
       dphi=dphi/3.
      end
C  ------------------------------------------------------------------------------  C
C  Function in the Impact factor
      function f(x)
	implicit none
      double precision x,f
      f=(3.53d0+2.79d0*x)/(3.53d0+x)/(1.d0+(x/0.71d0))**2
      if(x.lt.0.)print*,'x<0'
      end

      function E2(p1)
       implicit none
      double precision E2,p1(2),pz(2)
      double precision E3
       
      pz(1)=0.
      pz(2)=0.
      E2=E3(p1,pz)
      end

      function E3(p1,p2)
        implicit none
      integer iglu
      double precision E3,p1(2),p2(2),p3(2)
      double precision c,CHIDedot,ael

      common/Ifactor/c,ael,iglu

      p3(1)=-p1(1)-p2(1)
      p3(2)=-p1(2)-p2(2)

      E3=CHIDedot(p1,p1)+CHIDedot(p2,p2)+CHIDedot(p3,p3)
     &-c*(CHIDedot(p1,p2)+CHIDedot(p2,p3)+CHIDedot(p1,p3))
      end

C   ----------------------------------------------------------------------------   C
C   UgD FORWARD PREFACTOR
C   ----------------------------------------------------------------------------   C
      function prefac(p1,p2,t,x)
         implicit none
      double precision p1(2),p2(2),t,x
      double precision aprim,xo,Bo
      double precision f,prefac,CHIDedot,tfac
      
      Bo=4.d0
      aprim=0.25d0
      xo=3.4d0*10.E-4

      prefac=2.*CHIDedot(p1,p1)*CHIDedot(p2,p2)/
     & (CHIDedot(p1,p1)**2+CHIDedot(p2,p2)**2)
      prefac=prefac*Dexp(-0.5d0*(Bo+2.*aprim*log(xo/x))*t)

C      |t| dependance coming from hep-ph 0802.0177
C       tfac=Dexp(-4.*t)
C       prefac=2.*tfac
C     &*CHIDedot(p1,p1)*CHIDedot(p2,p2)/(CHIDedot(p1,p1)**2+CHIDedot(p2,p2)**2)

      return
      end 

C   ---------------------------------------------------------------------------   C
C   SUDAKOV FORM FACTOR
C   ---------------------------------------------------------------------------   C 
      function Dsuda(mu2,l2)
        implicit none
      integer IER
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol,Vgrid
      integer exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      real*8 Dsuda,mu2,mmu2,l2
      real*8 arg,Iarg
      real*8 ERROR,dcadredo,prec
      double precision as,s2
     
      external dcadredo
      external Iarg

      common/sudaD/mmu2,s2
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps

      prec=0.001d0
      mmu2=mu2

      if(l2.GE.mu2)then
        Dsuda=1.
      elseif(l2.LE.mu2)then
       arg=dcadredo(Iarg,l2,mu2,0.d0,prec,ERROR,IER)
       if(evolv.EQ.1)then 
        Dsuda=Dexp(-arg)
       elseif(evolv.EQ.2)then
        Dsuda=Dexp(-as(mu2)*arg)
       endif
      endif 

      end
C   ---------------------------------------------------------------------------   C
C   Integrand for dcadredo.f
      function Iarg(qt)
        implicit none
      integer nf
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol
      integer Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      real*8 pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      real*8 d,Iarg,Ia,Ib,Ilog,Iconst
      double precision qt,mmu2,s2,nnu2
      real*8 as

      common/const/pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      common/flavor/nf
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      common/sudaD/mmu2,s2

C   Sudakov prescription: 1=Durham, 2=DLA

      if(Model.EQ.3.AND.moode.EQ.2.AND.sudaf.EQ.1)then
       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt))  
       Ia=(-36.0*log(d)-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d+2.0*nf-33.0)/6.0

       Ilog=-6.*log(d)+(2.*nf-33.)/6.
       Iconst=(-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d)/6.0
       
      elseif(sudaf.EQ.1.OR.sudaf.EQ.6)then
       d=sqrt(qt)/sqrt(mmu2)  
       Ia=(-36.0*log(d)-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d+2.0*nf-33.0)/6.0  
 
      elseif(sudaf.EQ.2)then
       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt))  
       Ia=-6.*log(d)
      endif

      if(evolv.EQ.1)then 
       Ib=as(qt)/2./pi/qt
      elseif(evolv.EQ.2)then
       Ib=1./2./pi/qt
      endif
       
       Iarg=Ia*Ib

      end

C   ----------------------------------------------------------------------------   C
C   INTERPOLATION of UgD and SUDAKOV
C   All tabs are made from [interpol.f]
C   ----------------------------------------------------------------------------   C
      function dgdpol(p1,p2)
        implicit none
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol
      integer Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      double precision p1,p2,p1log,p2log,dgdpol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps

      p1log=log(p1)
      p2log=log(p2)
      
C     Dimension and parameter of dgdtab
c       lp1min=log(0.0001d0)                             
c       lp1max=log(0.4d0)
c       step1=(lp1max-lp1min)/1000.
c       lp2min=log(1.E-6)                                
c       lp2max=log(5000d0)                               
c       step2=(lp2max-lp2min)/1000.

       lp1min=-14d0                             
       lp1max=0.d0
       step1=(lp1max-lp1min)/1000.
       lp2min=-9.d0
       lp2max=8.d0
       step2=(lp2max-lp2min)/1000.

      if(p1log.LT.lp1min)p1log=lp1min+step1
      if(p2log.LT.lp2min)p2log=lp2min+step2
      if(p1log.GT.lp1max)p1log=lp1max-step1
      if(p2log.GT.lp2max)p2log=lp2max-step2


C     Bililear interpolation as in Numerical Recipes
      ip1=(p1log-lp1min)/step1
      ip2=(p2log-lp2min)/step2

C     The point to the left is x(1)+ix*stepx
C     The point to the bottom is k(1)+ik*stepk
      lp1left=lp1min+int(ip1)*step1
      lp1right=lp1min+int(ip1+1)*step1
      lp2top=lp2min+int(ip2+1)*step2
      lp2bottom=lp2min+int(ip2)*step2

C      print*,p1log,lp1left,lp1right,lp2top,lp2bottom      

      f1=dgdgrid(int(ip1)+1,int(ip2)+1)
      f2=dgdgrid(int(ip1+1)+1,int(ip2)+1)
      f3=dgdgrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=dgdgrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      dgdpol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
C      print*,f1,f2,f3,f4
   
      return
      end

C   ----------------------------------------------------------------------------   C 
      function Dsudapol(p1,p2)
        implicit none
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol
      integer Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      double precision p1,p2,p1log,p2log,Dsudapol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps

      p1log=log(p1)
      p2log=log(p2)

C     Dimension and parameter of dgdtab
      lp1min=7.82404601d0
      lp1max=11.982929d0
      step1=(lp1max-lp1min)/1000.
      lp2min=-13.81551056d0
      lp2max=9.210340372d0
      step2=(lp2max-lp2min)/1000.
 
      if(p1log.LE.lp1min)p1log=lp1min+step1
      if(p2log.LE.lp2min)p2log=lp2min+step2
      if(p1log.GE.lp1max)p1log=lp1max-step1
      if(p2log.GE.lp2max)p2log=lp2max-step2


C     Bililear interpolation as in Numerical Recipes
      ip1=(p1log-lp1min)/step1
      ip2=(p2log-lp2min)/step2

C     The point to the left is x(1)+ix*stepx
C     The point to the bottom is k(1)+ik*stepk
      lp1left=lp1min+int(ip1)*step1
      lp1right=lp1min+int(ip1+1)*step1
      lp2top=lp2min+int(ip2+1)*step2
      lp2bottom=lp2min+int(ip2)*step2
      
      f1=sudagrid(int(ip1)+1,int(ip2)+1)
      f2=sudagrid(int(ip1+1)+1,int(ip2)+1)
      f3=sudagrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=sudagrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      Dsudapol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
   
      return
      end


