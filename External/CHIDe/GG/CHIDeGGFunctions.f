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
C   SUBROUTINE KINEMATIC
C   ----------------------------------------------------------------------------   C
      subroutine CHIDekinematics(k,kp,k1,k2,k3)
       implicit none 
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1

C   Usefulf sum of vectors
C   k+k1
      skk1(1)=k(1)+k1(1)
      skk1(2)=k(2)+k1(2)
      skpk1(1)=kp(1)+k1(1)
      skpk1(2)=kp(2)+k1(2)

      mskk1(1)=-k(1)-k1(1)
      mskk1(2)=-k(2)-k1(2)
      mskpk1(1)=-kp(1)-k1(1)
      mskpk1(2)=-kp(2)-k1(2)
      
      mk1(1)=-k1(1)
      mk1(2)=-k1(2)

C   k+k3
      skk3(1)=k(1)+k3(1)
      skk3(2)=k(2)+k3(2)
      skpk3(1)=kp(1)+k3(1)
      skpk3(2)=kp(2)+k3(2)

      mskk3(1)=-k(1)-k3(1)
      mskk3(2)=-k(2)-k3(2)
      mskpk3(1)=-kp(1)-k3(1)
      mskpk3(2)=-kp(2)-k3(2)

      mk3(1)=-k3(1)
      mk3(2)=-k3(2)

C   k+k2
      skk2(1)=k(1)+k2(1)
      skk2(2)=k(2)+k2(2)
      skpk2(1)=kp(1)+k2(1)
      skpk2(2)=kp(2)+k2(2)

      mskk2(1)=-k(1)-k2(1)
      mskk2(2)=-k(2)-k2(2)
      mskpk2(1)=-kp(1)-k2(1)
      mskpk2(2)=-kp(2)-k2(2)     

      mk2(1)=-k2(1)
      mk2(2)=-k2(2)

C   k+k2-k3
      skk2k3(1)=k(1)+k2(1)-k3(1)
      skk2k3(2)=k(2)+k2(2)-k3(2)
      skpk2k3(1)=kp(1)+k2(1)-k3(1)
      skpk2k3(2)=kp(2)+k2(2)-k3(2)      

      mskk2k3(1)=-k(1)-k2(1)+k3(1)
      mskk2k3(2)=-k(2)-k2(2)+k3(2)
      mskpk2k3(1)=-kp(1)-k2(1)+k3(1)
      mskpk2k3(2)=-kp(2)-k2(2)+k3(2)      
    
C   k and kp
      mk(1)=-k(1)
      mk(2)=-k(2)
      mkp(1)=-kp(1)
      mkp(2)=-kp(2)
    
      end

C   ----------------------------------------------------------------------------   C
C   SUBROUTINE ALPHA STRONG
C   Compute as in function of Q²:
C      -> choose a Q²
C      -> compute the corresponding lambda(nf)
C      -> compute as
C   ----------------------------------------------------------------------------   C
      function CHIDeas(q) 
       implicit none 
      integer i,nf 
      double precision CHIDeas
      double precision m(2:6)
      double precision smu,smc,sms,smb,smt 
      double precision lambda(2:6),la
      double precision q
      double precision pi,s,couleur,gg,gq,kmax,nb,mp

      common/const/pi,couleur,gg,gq,kmax,nb,mp
      common/flavor/nf
      common/param/s
            
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
      CHIDeas=12.*pi/(33.-2.*nf)/log(q/la**2)

      if(CHIDeas.LT.0.or.CHIDeas.GT.0.7d0)CHIDeas=0.7d0

666   end    

C   ----------------------------------------------------------------------------   C
C   LCWF IMPACT FACTOR 
C   Elastic, diffractif, function f(x)~ effective slope B for the proton  
C   ----------------------------------------------------------------------------   C  
C   Elastic Impact factor
C   phi: Chosen from Jr and Hernandez paper
      function CHIDeephi(p1,p2)
	implicit none
      integer iglu
      double precision CHIDeephi,CHIDedot,CHIDedotsum,CHIDef
      double precision p1(2),p2(2),x1,x2,c,ael

      common/Ifactor/c,ael,iglu         
         
      x1=CHIDedotsum(p1,p2,p1,p2)
      x2=CHIDedot(p1,p1)+CHIDedot(p2,p2)-c*CHIDedot(p1,p2)

      CHIDeephi=CHIDef(x1)-CHIDef(x2)
C      if((CHIDedot(p1,p1)).LT.10E-8.or.(CHIDedot(p2,p2)).LT.10E-8)then
C      CHIDeephi=0.
C      endif 

      end
C   ----------------------------------------------------------------------------   C
C   Diffractif Impact factor
C   dphi from C&H p.479
      function CHIDedphi(p1,p2,p3,p4)
       implicit none
      double precision p1(2),p2(2),p3(2),p4(2),CHIDedphi
      double precision sp1p2(2),sp1p4(2),sp1p3(2)
      double precision CHIDef,CHIDeE3,CHIDeE2
      double precision x1,x2,x3,x4,x12,x14,x13
      double precision y12,y34,y23,y13,y14,y24
 
       x1=CHIDeE2(p1)
       x2=CHIDeE2(p2)
       x3=CHIDeE2(p3)
       x4=CHIDeE2(p4)

       sp1p2(1)=p1(1)+p2(1)
       sp1p2(2)=p1(2)+p2(2)
       sp1p4(1)=p1(1)+p4(1)
       sp1p4(2)=p1(2)+p4(2)
       sp1p3(1)=p1(1)+p3(1)
       sp1p3(2)=p1(2)+p3(2)

       x12=CHIDeE2(sp1p2)
       x14=CHIDeE2(sp1p4)
       x13=CHIDeE2(sp1p3)

       y12=CHIDeE3(p1,p2)
       y34=CHIDeE3(p3,p4)
       y23=CHIDeE3(p2,p3)
       y13=CHIDeE3(p1,p3)
       y14=CHIDeE3(p1,p4)
       y24=CHIDeE3(p4,p2)

       CHIDedphi=1.d0-CHIDef(x1)-CHIDef(x2)-CHIDef(x3)-CHIDef(x4)
     &+2.d0*CHIDef(x12)+0.5d0*CHIDef(x14)+0.5d0*CHIDef(x13)
     &-CHIDef(y12)-CHIDef(y34)+0.5d0*CHIDef(y23)+0.5d0*CHIDef(y13)
     &+0.5d0*CHIDef(y14)+0.5d0*CHIDef(y24)
       CHIDedphi=CHIDedphi/3.
      end
C  ------------------------------------------------------------------------------  C
C  Function in the Impact factor
      function CHIDef(x)
	implicit none
      double precision x,CHIDef
      CHIDef=(3.53d0+2.79d0*x)/(3.53d0+x)/(1.d0+(x/0.71d0))**2
      if(x.lt.0.)print*,'x<0'
      end

      function CHIDeE2(p1)
       implicit none
      double precision CHIDeE2,p1(2),pz(2)
      double precision CHIDeE3
       
      pz(1)=0.
      pz(2)=0.
      CHIDeE2=CHIDeE3(p1,pz)
      end

      function CHIDeE3(p1,p2)
        implicit none
      integer iglu
      double precision CHIDeE3,p1(2),p2(2),p3(2)
      double precision c,CHIDedot,ael

      common/Ifactor/c,ael,iglu

      p3(1)=-p1(1)-p2(1)
      p3(2)=-p1(2)-p2(2)

      CHIDeE3=CHIDedot(p1,p1)+CHIDedot(p2,p2)+CHIDedot(p3,p3)
     &-c*(CHIDedot(p1,p2)+CHIDedot(p2,p3)+CHIDedot(p1,p3))
      end

C   ----------------------------------------------------------------------------   C
C   UgD FORWARD PREFACTOR
C   ----------------------------------------------------------------------------   C
      function CHIDeprefac(p1,p2,t,x)
         implicit none
      double precision p1(2),p2(2),t,x
      double precision aprim,xo,Bo
      double precision CHIDef,CHIDeprefac,CHIDedot,tfac
      
      Bo=4.d0
      aprim=0.25d0
      xo=3.4d0*10.E-4

      CHIDeprefac=2.*CHIDedot(p1,p1)*CHIDedot(p2,p2)
     &        /(CHIDedot(p1,p1)**2+CHIDedot(p2,p2)**2)
      CHIDeprefac=CHIDeprefac*Dexp(-0.5d0*(Bo+2.*aprim*log(xo/x))*t)

C      |t| dependance coming from hep-ph 0802.0177
C       tfac=Dexp(-4.*t)
C       CHIDeprefac=2.*tfac
C     &*CHIDedot(p1,p1)*CHIDedot(p2,p2)/(CHIDedot(p1,p1)**2+CHIDedot(p2,p2)**2)

      return
      end 

C   ---------------------------------------------------------------------------   C
C   SUDAKOV FORM FACTOR
C   ---------------------------------------------------------------------------   C 
      function CHIDeDsuda(mu2,p2)
        implicit none
      integer IER
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      real*8 CHIDeDsuda,mu2,mmu2,p2
      real*8 arg,CHIDeIarg
      real*8 ERROR,dcadredo,prec
      double precision CHIDeas,s2
     
      external dcadredo
      external CHIDeIarg

      common/sudaD/mmu2,s2
      common/switch/event,exp,Sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale

      prec=0.001d0
      mmu2=mu2

      if(p2.GE.mu2)then
        CHIDeDsuda=1.
      elseif(p2.LT.mu2.or.p2.EQ.mu2)then
       arg=dcadredo(CHIDeIarg,p2,mu2,0.d0,prec,ERROR,IER)
      endif 

      if(evolv.EQ.1)then 
       CHIDeDsuda=Dexp(-arg)
      elseif(evolv.EQ.2)then
       CHIDeDsuda=Dexp(-CHIDeas(mu2)*arg)
      endif

      end
C   ---------------------------------------------------------------------------   C
C   Integrand for dcadredo.f
      function CHIDeIarg(qt)
        implicit none
      integer nf,event,exp,Sudaf,grid,formfac,Lcut,splash
      integer scale,Model,evolv,Vgrid,interpol,moode,lhe,nscale
      double precision qt,mmu2,s2,nnu2
      double precision tgg,sgg,ugg,Sudaphi,phi
      real*8 pi,s,ncolor,gg,gq,kmax,nb,mp
      real*8 d,CHIDeIarg,Ia,Ib,Ilog,Iconst
      real*8 CHIDeas

      common/const/pi,ncolor,gg,gq,kmax,nb,mp
      common/flavor/nf
      common/integrand/tgg,sgg,ugg,Sudaphi,phi
      common/switch/event,exp,Sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/sudaD/mmu2,s2
      common/param/s

      if(nscale.EQ.1)then
       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt))
      elseif(nscale.EQ.2)then
       nnu2=sgg/8.
       d=sqrt(qt)/(sqrt(nnu2)+sqrt(qt))
      endif

C   Sudakov prescription: 1=Durham, 2=DDT, 3=DLA
C   Old prescription for delta -> PAPER and sudatab
C       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt))  

      if(sudaf.EQ.1)then
       Ia=(-36.0*log(d)-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d+2.0*nf-33.0)/6.0

       Ilog=-6.*log(d)+(2.*nf-33.)/6.
       Iconst=(-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d)/6.0

      elseif(sudaf.EQ.2)then
       d=qt/mmu2
       Ia=(2.*float(nf)-33.0-36.0*d-36.0*(d+1.0)**2*log(d)
     &+36.0*(d+1.0)**2*log(d+1.0))/6.0

       Ilog=(2.*float(nf)-33.0-36.0*log(d))/6.0
       Iconst=(-36.0*d-36.0*(d**2+2.*d)*log(d)
     &+36.0*(d+1.0)**2*log(d+1.0))/6.0

      elseif(sudaf.EQ.3)then
       Ia=-6.*log(d)
      endif

      if(evolv.EQ.1)then 
       Ib=CHIDeas(qt)/2./pi/qt
      elseif(evolv.EQ.2)then
       Ib=1./2./pi/qt
      endif
       
       CHIDeIarg=Ia*Ib
       
      end

C   ----------------------------------------------------------------------------   C
C   INTERPOLATION of UgD and SUDAKOV
C   ----------------------------------------------------------------------------   C
      function CHIDedgdpol(p1,p2)
        implicit none
      integer pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4    
      double precision p1,p2,p1log,p2log,CHIDedgdpol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid
      common/logs/pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4

      p1log=log(p1)
      p2log=log(p2)

C     Dimension and parameter of dgdtab
      lp1min=-14.d0
      lp1max=0.d0
      step1=(lp1max-lp1min)/1000.
      lp2min=-9.d0
      lp2max=8.d0
      step2=(lp2max-lp2min)/1000.

      if(p1log.LT.lp1min)pout1=pout1+1
      if(p2log.LT.lp2min)pout2=pout2+1
      if(p1log.GT.lp1max)pout3=pout3+1
      if(p2log.GT.lp2max)pout4=pout4+1    

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
      
      f1=dgdgrid(int(ip1)+1,int(ip2)+1)
      f2=dgdgrid(int(ip1+1)+1,int(ip2)+1)
      f3=dgdgrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=dgdgrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      CHIDedgdpol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
   
      return
      end

C   ----------------------------------------------------------------------------   C 
      function CHIDeDsudapol(p1,p2)
        implicit none
      double precision p1,p2,p1log,p2log,CHIDeDsudapol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid

      p1log=log(p1)
      p2log=log(p2)

C     Dimension and parameter of dgdtab
      lp1min=-14.d0
      lp1max=8.d0
      step1=(lp1max-lp1min)/1000.
      lp2min=-14.d0
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
      
      f1=sudagrid(int(ip1)+1,int(ip2)+1)
      f2=sudagrid(int(ip1+1)+1,int(ip2)+1)
      f3=sudagrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=sudagrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      CHIDeDsudapol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
   
      return
      end

