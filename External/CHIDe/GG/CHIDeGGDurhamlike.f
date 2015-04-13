C   ============================================================================   C
C   SUBROUTINE DURHAMLIKE:                                                         C
C   A external subroutine that mimic Durham calculation, it means a different      C 
C   scale in the Sudakov form factor and the approximation that ki << k            C
C       moode = 1: Reference curve with ki << k                                    C
C       moode = 2: Durham like                                                     C
C   ============================================================================   C

C   ============================================================================   C
C   FUNCTION INTEGM1:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 1   ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |   |--                                   --|   |                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      subroutine CHIDedurhamM1(sigmadurhamM1,a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision CHIDedot,sigmadurhamM1
      double precision phi,CHIDeas
      double precision b1,b2,k2(2),a1,a2
      double precision k(2),kp(2),k1(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p
      double precision Ij,C0,C2
      double precision sgg,tgg,ugg
      double precision I1,J1,J2
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi
      double precision phip,phipb,Ifact
      double precision c,ael
      double precision x,xp
      double precision Rg,CHIDeprefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      double precision CHIDedgdpol,CHIDeDsudapol
      double precision mmu2,S2
      double precision fact,factsp,factp,regge,a1t,a3t

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/sudaD/mmu2,S2

C   Common factor
      Ifact=pi**2/gq**2

C   ----------------------  M1: Num and Den  ---------------------------------   C
C     Numerator
      I1=(1./2.)*(1.+(ugg**4/sgg**4)+(tgg**4/sgg**4))*(1./tgg+1./ugg)**2
      J1=(1./2.)*(1.-(ugg**4/sgg**4)-(tgg**4/sgg**4))*(1./tgg+1./ugg)**2
      I2=CHIDedot(k,kp)*CHIDedot(k,kp)
      J2=CHIDedot(k,k)*CHIDedot(kp,kp)-CHIDedot(k,kp)*CHIDedot(k,kp)

C     Angular momentum dependence
C     J=0 then C2=0 or J=2 and then C0=0 
      if(moode.EQ.1)then
       C0=I2+J2  
       C2=(I2-J2)
      elseif(moode.EQ.2)then   
       C0=I2+J2
       C2=0.
      endif

      Ij=C0*(I1+J1)/2.+C2*(I1-J1)/2.

C     Denominator
      I3=CHIDedot(k,k)*CHIDedot(k,k)*CHIDedot(k,k)
      I3p=CHIDedot(kp,kp)*CHIDedot(kp,kp)*CHIDedot(kp,kp)

      if(CHIDedot(k,k).LT.propcut.or.CHIDedot(kp,kp).LT.propcut)goto 666
      if(CHIDedot(k,kp).LT.propcut)goto 666

C   --------------  M1: Impact and Sudakov Factor  --------------------------   C
C    Impact factor
      phip=CHIDedgdpol(b1+b2,CHIDedot(k,k))
     &     *CHIDeprefac(k,k,CHIDedot(k1,k1),b1+b2)
     &     *CHIDedgdpol(b1+b2,CHIDedot(kp,kp))
     &     *CHIDeprefac(kp,kp,CHIDedot(k1,k1),b1+b2)
      phipb=CHIDedgdpol(a1+a2,CHIDedot(k,k))
     &     *CHIDeprefac(mk,mk,CHIDedot(k3,k3),a1+a2) 
     &     *CHIDedgdpol(a1+a2,CHIDedot(kp,kp))
     &     *CHIDeprefac(mkp,mkp,CHIDedot(k3,k3),a1+a2)

C    Sudakov form factor
C    Durham prescription: 0.3844=0.62² cf [arXiv:0802.0177]
      if(moode.EQ.1)then  
       Sphi=CHIDeDsudapol(CHIDedot(k2,k2)/x,CHIDedot(k,k)/xp)
       Sphib=CHIDeDsudapol(CHIDedot(k2,k2)/x,CHIDedot(kp,kp)/xp)

      elseif(moode.EQ.2)then
C       Rg=1.4d0
       Rg=1.d0  
       Sphi=Rg**2*CHIDeDsudapol(0.3844d0*sgg,CHIDedot(k,k))
       Sphib=Rg**2*CHIDeDsudapol(0.3844d0*sgg,CHIDedot(kp,kp))
      endif

C   -------------------------------  Diag1=M1  -----------------------------------  C 
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif  
      Sudaphi=Sphi*Sphib

C     Coefficient, function and phase space
C     First pick up the corresponding value of gg(as) from 
C     subroutine and ael,aDD from sigmael.f
      factsp=1./16./(2.*pi)**8/b1/b2
      factp=81./2.              !3 quarks in p+ and 1/2 because of gluons
      gg=sqrt(4.*pi*CHIDeas(sgg))
      fact=gg**4*gq**8*nb/4./pi**4

C     Regge factor
      if(formfac.EQ.1)then
      a1t=0.09d0-0.3d0*CHIDedot(k1,k1)/(1.-(b1+b2))
      a3t=0.09d0-0.3d0*CHIDedot(k3,k3)/(1.-(a2+a1))
      regge=((1./(a2+a1))**a1t*(1./(b1+b2))**a3t)**2
      else
       regge=1.d0
      endif

C    Color factor
      Cf=(ncolor**2-1.)/ncolor**2

C    Function M1
      sigmadurhamM1=fact*factsp*factp*regge*Cf*Ij*Sudaphi*phi/I3/I3p
      sigmadurhamM1=S2*sigmadurhamM1
C    Debug tools  
C      print*,'In Diag 1 Durham mode:'
C      print*,'phi=',phi,phip,phipb
C      print*,'Fact,Ifact=',Ifact,pi,gq
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo',I3*I3p,Ij
C      print*,'********************************'    

      goto 66
666   sigmadurhamM1=0.
66    end



C   ============================================================================   C
C   FUNCTION INTEGM2:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 2   ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |------                                   ------|                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      subroutine CHIDedurhamM2(sigmadurhamM2,a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision CHIDedot,CHIDedotsum,sigmadurhamM2
      double precision fc,fs
      double precision phi,CHIDeas 
      double precision b1,b2,a1,a2
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p
      double precision I1,Ij
      double precision sgg,tgg,ugg
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi
      double precision phip,phipb,Ifact
      double precision c,ael
      double precision x,xp
      double precision CHIDeprefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      double precision CHIDedgdpol
      double precision mmu2,S2
      double precision fact,factsp,factp,regge,a1t,a3t

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/sudaD/mmu2,S2

C   Common factor
      Ifact=pi**2/gq**2

C   ----------------------  M2: Den and Num  ----------------------------------   C
C     Ij: numerator
      I1=(fc(k,skk2)*fc(kp,skpk2))+fs(k,skk2,kp,skpk2)
      I2=(fc(mk,mskk2)*fc(mkp,mskpk2))+fs(mk,mskk2,mkp,mskpk2)
      Ij=16.*I1*I2/(sqrt(CHIDedot(k2,k2)))**4

C     Denominator!
      I3=CHIDedot(k,k)*CHIDedot(k,k)*CHIDedotsum(k,k2,k,k2)
     &   *CHIDedotsum(k,k2,k,k2)
      I3p=CHIDedot(kp,kp)*CHIDedot(kp,kp)*CHIDedotsum(kp,k2,kp,k2)
     &   *CHIDedotsum(kp,k2,kp,k2)

      if(CHIDedot(k,k).LT.propcut.or.CHIDedot(kp,kp).LT.propcut)goto 666
      if(CHIDedotsum(k,k2,k,k2).LT.propcut)goto 666
      if(CHIDedotsum(kp,k2,kp,k2).LT.propcut)goto 666

C   -----------------  M2: Impact and Sudakov Factor  -------------------------   C
C     Impact factor
       phip=CHIDedgdpol(0.5d0*(b1+b2),CHIDedot(k,k))
     &*CHIDeprefac(k,k,CHIDedot(k1,k1),0.5d0*(b1+b2))
     &*CHIDedgdpol(0.5d0*(b1+b2),CHIDedot(kp,kp))
     &*CHIDeprefac(kp,mkp,CHIDedot(k1,k1),0.5d0*(b1+b2))
      phipb=CHIDedgdpol(0.5d0*(a1+a2),CHIDedotsum(k,k2,k,k2))
     &*CHIDeprefac(mskk2,skk2,CHIDedot(k3,k3),0.5*(a1+a2)) 
     &*CHIDedgdpol(0.5d0*(a1+a2),CHIDedotsum(kp,k2,kp,k2))
     &*CHIDeprefac(mskpk2,skpk2,CHIDedot(k3,k3),0.5*(a1+a2))

C     No Sudakov Factor
       Sphi=1.d0
       Sphib=1.d0 

C   -------------------------------  Dia2=M2  -----------------------------------  C
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif  
      Sudaphi=Sphi*Sphib

C     Coefficient, function and phase space
C     First pick up the corresponding value of gg(as) from 
C     subroutine and ael,aDD from sigmael.f
      factsp=1./16./(2.*pi)**8/b1/b2
      factp=81./2.              !3 quarks in p+ and 1/2 because of gluons
      gg=sqrt(4.*pi*CHIDeas(sgg))
      fact=gg**4*gq**8*nb/4./pi**4

C     Regge factor
      if(formfac.EQ.1)then
      a1t=0.09d0-0.3d0*CHIDedot(k1,k1)/(1.-(b1+b2))
      a3t=0.09d0-0.3d0*CHIDedot(k3,k3)/(1.-(a2+a1))
      regge=((1./(a2+a1))**a1t*(1./(b1+b2))**a3t)**2
      else
       regge=1.d0
      endif

C    Color factor
      Cf=(ncolor**2-1.)/ncolor**2

C    Function INTEG
      sigmadurhamM2=fact*factsp*factp*regge*Cf*Ij*Sudaphi*phi/I3/I3p
      sigmadurhamM2=S2*sigmadurhamM2

C    Debug tools  
C      print*,'In Diag 2:'
C      print*,'phi=',phi,skk2,mk3
C      print*,'Fact,Ifact=',Ifact,pi,gq
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo=',Ij
C      print*,'********************************'
   
      goto 66
666   sigmadurhamM2=0.
66    end

C   ============================================================================   C
C   FUNCTION INTEGM12:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 12  ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |   |--                                   ------|                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      subroutine CHIDedurhamM12(sigmadurhamM12
     &                           ,a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision CHIDedot,CHIDedotsum,sigmadurhamM12
      double precision fc,fs
      double precision b1,b2,a1,a2
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p,I4,I5
      double precision I1,Ij
      double precision sgg,tgg,ugg
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi,CHIDeas
      double precision phi,phip,phipb,Ifact
      double precision c,ael
      double precision x,xp
      double precision Rg,CHIDeprefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      double precision CHIDedgdpol,CHIDeDsudapol
      double precision mmu2,S2
      double precision fact,factsp,factp,regge,a1t,a3t

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/sudaD/mmu2,S2


C   Common factor
      Ifact=pi**2/gq**2

C   ----------------------  M12: Num and Den  --------------------------------   C  
C     Ij: numerator
      I1=(1.+ugg**2/sgg**2+tgg**2/sgg**2)
     &   *(fc(mk2,k)*fc(mk2,k)*fc(kp,skpk2)*fc(mkp,mskpk2)
     &   - fs(mk2,k,mk2,k)*fs(kp,skpk2,mkp,mskpk2))
      I2=(1.-ugg**2/sgg**2-tgg**2/sgg**2)
     &   *(fc(kp,skpk2)*fc(mkp,mskpk2)*fs(mk2,k,mk2,k)
     &   -fc(mk2,k)*fc(mk2,k)*fs(kp,skpk2,mkp,mskpk2))
      I4=(1.+ugg**2/sgg**2-tgg**2/sgg**2)
     &   *(fc(mk2,k)*fc(mkp,mskpk2)*fs(mk2,k,kp,skpk2)
     &   -fc(mk2,k)*fc(kp,skpk2)*fs(mk2,k,mkp,mskpk2))
      I5=(1.-ugg**2/sgg**2+tgg**2/sgg**2)
     &   *(fc(mk2,k)*fc(kp,skpk2)*fs(mk2,k,mkp,mskpk2)
     &   -fc(mk2,k)*fc(mkp,mskpk2)*fs(mk2,k,kp,skpk2))

      Ij=8.d0*(I1+I2+I4+I5)/(sqrt(CHIDedot(k2,k2)))**4

C     Denominator!
      I3=CHIDedot(k,k)*CHIDedot(k,k)*CHIDedot(k,k)
      I3p=CHIDedot(kp,kp)*CHIDedot(kp,kp)*CHIDedotsum(kp,k2,kp,k2)
     &     *CHIDedotsum(kp,k2,kp,k2)

      if(CHIDedot(k,k).LT.propcut.or.CHIDedot(kp,kp).LT.propcut)goto 666
      if(CHIDedotsum(kp,k2,kp,k2).LT.propcut)goto 666

C   -----------------  M12: Impact and Sudakov Factor  -----------------------   C
      phip=CHIDedgdpol((b1+b2),CHIDedot(k,k))
     &      *CHIDeprefac(k,mk,CHIDedot(k1,k1),(b1+b2))
     &      *CHIDedgdpol(0.5d0*(b1+b2),CHIDedot(kp,kp))
     &      *CHIDeprefac(kp,mkp,CHIDedot(k1,k1),0.5d0*(b1+b2))
      phipb=CHIDedgdpol((a1+a2),CHIDedot(k,k))
     &      *CHIDeprefac(mk,k,CHIDedot(k3,k3),(a1+a2)) 
     &      *CHIDedgdpol(0.5d0*(a1+a2),CHIDedotsum(kp,k2,kp,k2))
     &      *CHIDeprefac(mskpk2,skpk2,CHIDedot(k3,k3),0.5d0*(a1+a2))

C    Sudakov form factor
C    Durham prescription: 0.3844=0.62² cf [arXiv:0802.0177]
      if(moode.EQ.1)then  
       Sphi=CHIDeDsudapol(CHIDedot(k2,k2)/x,CHIDedot(k,k)/xp)
       Sphib=1.d0

      elseif(moode.EQ.2)then
       Rg=1.4d0
       Sphi=Rg**2*CHIDeDsudapol(0.3844d0*sgg,CHIDedot(k,k)/xp)
       Sphib=1.d0
      endif

C   -------------------------------  Dia12=M12  ---------------------------------  C
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif 

      Sudaphi=Sphi*Sphib

C     Coefficient, function and phase space
C     First pick up the corresponding value of gg(as) from 
C     subroutine and ael,aDD from sigmael.f
      factsp=1./16./(2.*pi)**8/b1/b2
      factp=81./2.              !3 quarks in p+ and 1/2 because of gluons
      gg=sqrt(4.*pi*CHIDeas(sgg))
      fact=gg**4*gq**8*nb/4./pi**4

C     Regge factor
      if(formfac.EQ.1)then
      a1t=0.09d0-0.3d0*CHIDedot(k1,k1)/(1.-(b1+b2))
      a3t=0.09d0-0.3d0*CHIDedot(k3,k3)/(1.-(a2+a1))
      regge=((1./(a2+a1))**a1t*(1./(b1+b2))**a3t)**2
      else
       regge=1.d0
      endif

C   Color factor
      Cf=-(ncolor**2-1.)/ncolor**2

C   Function INTEG
      sigmadurhamM12=fact*factsp*factp*regge*Cf*Ij*Sudaphi*phi/I3/I3p
      sigmadurhamM12=S2*sigmadurhamM12

C   ------------------------------------------------------------------------------  C 

C    Debug tools  
C      print*,'In Diag 12:'
C      print*,'phi=',phi,phipb
C      print*,'Fact,Ifact=',Ifact,Cf
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo=',Ij,I1,I2,I4,I5,1./I3,1./I3p
C      print*,'M12=',IntegM12
C      print*,'********************************'

      goto 66
666   sigmadurhamM12=0.
66    end
