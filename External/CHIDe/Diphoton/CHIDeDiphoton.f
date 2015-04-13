      subroutine CHIDeDiphoton(sigma,k_, kp_, k1_, k2_,
     & k3_,b1_,b2_,a1_,a2_)
      implicit none
      integer it,ndo  
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,Lcut,nscale
      integer scale,Model,evolv,Vgrid,interpol,moode,lhe,splash
      integer jetveto,update
      integer iglu,nf
      integer ik,ik1,ik2
      double precision k_(2), kp_(2), k1_(2), k2_(2), k3_(2)
      double precision sigma,jacc
      double precision b1_, b2_, abcut_,a1_,a2_, xi1_, xi2_, x1,
     & x2,c1,c2,x1x2,b1_1, b1_2, b2_1, b2_2
      double precision si,swgt,schi,xi
      double precision x(11),CHIDeJetsdsigma,wgt,CHIDedot
      double precision CHIDedotdiff,CHIDedotsum
      double precision sigmaM1,sigmaM2,sigmaM12,sigmaM12p
      double precision sigmadurhamM1,sigmadurhamM2
      double precision sigmadurhamM12,sigmadurhamM12p
      double precision k1(2),ak1,theta1,k2(2),ak2,theta2
      double precision k3(2),ak3,theta3
      double precision k(2),ak,thetak,kp(2),akp,thetap
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision b1,b2,a1,a2,abcut,hbcutd,lbcutd
      double precision sgg,tgg,ugg,Sudaphi,phi
      double precision ael,c,aelm
      double precision kmax,pi,s,ncolor,gg,gq,gelm,nb,mp
      double precision weight,weightsuda,weightphi
      double precision fact,factsp,jac,factp,regge
      double precision CHIDeDiphotonas,Cf
      double precision gap1,gap2,gap3,gap4
      double precision yp1,yp2,hygap,lygap,yj1,yj2,ycut,yclu
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      double precision Ejet1,Ejet2,Etcut,Etran
      double precision a1t,a3t
      double precision mmu2,s2
      double precision xmax,amax,amin
      double precision normk1,normk3,normk,normkp,normk2,xb1b2
      double precision sskk1,sskk2,sskk3,sskk2k3
      double precision Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      double precision avgi,sd
      double precision nj
      double precision mq(1:6),Qup,Qdown

      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist
      common/Energy/Ejet1,Ejet2    
      common/const/pi,ncolor,gg,gq,gelm,kmax,nb,mp
      common/ab/a1,a2,b1,b2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi
      common/Ifactor/c,ael,aelm,iglu
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/switch/event,exp,sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/sudaD/mmu2,s2
      common/BotE/amax,xmax,amin
      common/lhe/k,k1,k2,k3
      common/out/avgi,sd,weight,update
      common/algo/Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/param/s
      common/flavor/mq,Qup,Qdown,nf

      CHIDeJetsdsigma=0.d0
      jac = 0d0

      k(1) = k_(1)
      k(2) = k_(2)
      kp(1) = kp_(1)
      kp(2) = kp_(2)
      k1(1) = k1_(1)
      k1(2) = k1_(2)
      k2(1) = k2_(1)
      k2(2) = k2_(2)
      k3(1) = k3_(1)
      k3(2) = k3_(2)
 
      call CHIDekinematics(k,kp,k1,k2,k3)
      b1 = b1_
      b2 = b2_
      a1 = a1_
      a2 = a2_
            
      if(a1+a2.GT.1.d0.or.b1+b2.GT.1.d0) goto 999

C     Gluon gluon MANDELSTAM variables
      sgg=CHIDedot(k2,k2)*(b1+b2)**2/b1/b2
      tgg=-CHIDedot(k2,k2)*(b1+b2)/b2
      ugg=-CHIDedot(k2,k2)*(b1+b2)/b1

C     Rapidities: antiproton denotes with p2
      yp1=0.5d0*log((1.-b1-b2)**2*s/(CHIDedot(k1,k1)+mp**2))
      yp2=-0.5d0*log((1.-a1-a2)**2*s/(CHIDedot(k3,k3)+mp**2))
      yclu=0.5d0*log((b1+b2)/(a1+a2))
      yj1=0.5d0*log(b1/a1)
      yj2=0.5d0*log(b2/a2)

      if(yp1.LT.0.d0.or.yp2.GT.0.d0) goto 999

C     Transverse energy of the two jets
      Ejet1=sqrt(CHIDedotsum(k1,k2,k1,k2))
      Ejet2=sqrt(CHIDedotsum(k2,k3,k2,k3))

C     Gap
C      gap1=yp1-yj1
C      gap2=yp1-yj2
C      gap3=-yp2+yj1
C      gap4=-yp2+yj2

C     Cuts on transverse energy of the two jets
      if(Ejet1.LT.Etcut.or.Ejet2.LT.Etcut)goto 999
      
C     Cuts on alpha and beta ordered
C      if(a1+a2.GT.hacutd.or.a1+a2.LT.lacutd)goto 999
C      if(b1+b2.GT.hbcutd.or.b1+b2.LT.lbcutd)goto 999 

C     Cuts on rapidity
      if(yj1.GT.ycut.or.yj2.GT.ycut)goto 999 
      if(yj1.LT.ylcut.or.yj2.LT.ylcut)goto 999 


C     Cuts on the gap
c     if(min(gap1,gap2).LT.lygap)goto 999
c     if(min(gap3,gap4).LT.lygap)goto 999

     
C   ============================================================================   C
C                        Cross section = m1-m2+2*m1m2                              C
C   ============================================================================   C
       call CHIDeDiphotonM(sigmaM1,a1,a2,b1,b2,k1,k2,k3,k,kp) 

       CHIDeJetsdsigma=sigmaM1

999    sigma = CHIDeJetsdsigma
       return

      end







