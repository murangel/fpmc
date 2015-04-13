c      subroutine StandaloneCHIDeDiphoton
       program StandaloneCHIDeDiphoton
       implicit none

      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,scale,splash
      integer Model,evolv,Vgrid,Lcut,interpol,moode,lhe,nscale
      integer iglu,j,i
      integer update
      integer verbose
      double precision xl,xu,acc,si,swgt,schi,xi
      double precision alph
      double precision weight,c,ael
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,yclucut,propcut
      double precision Etjetmin,Etmin,CHIDeas
      double precision Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp
      double precision mmu,s2

      double precision avgi, sd, chi2a


      double precision xx(11), ww

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist
      common/out/avgi,sd,weight,update
      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/Ifactor/c,ael,iglu
      common/algo/Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      common/const/pi,ncolor,gg,gq,kmax,nb,mp
      common/sudaD/mmu,s2
      common/param/s
      common/bbb/verbose

      double precision dsigma, CHIDedotdiff
      external CHIDedotdiff
      double precision b1, b2, abcut, a1, a2, xi1, xi2
      double precision k(2), kp(2), k1(2), k2(2), k3(2)
       
       double precision sigma
       external sigma  
       character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       double precision ss,ptmin
       common/aaa/ss,ptmin
       dgdtab4="Data/dgdtab1.d"
       dgdtab4="Data/dgdtab2.d"
       dgdtab4="Data/dgdtab3.d"
       dgdtab4="Data/dgdtab4.d"
       sudatab="Data/ggsudatab.d"
       
        ss=14000.0**2
       ptmin=5.0
       CALL CHIDeDiphotonInit(ss,ptmin,0.03d0,.5d0,1d0,4,
     &                  0.002d0,0.02d0,0.002d0,0.02d0,-2d0,2d0)
       
c      ss=1960.0**2
c      ptmin=5.0
c      CALL CHIDeDiphotonInit(ss,ptmin,0.15d0,.5d0,1d0,4,
c    &                  0.03d0,0.08d0,0.0d0,1.0d0,-1d0,1d0)


      xx(1)=0.5
      xx(2)=0.5
      xx(3)=0.5
      xx(4)=0.5
      xx(5)=0.5
      xx(6)=0.5
      xx(7)=0.5
      xx(8)=0.5
      xx(9)=0.5
      xx(10)=0.5
      xx(11)=0.5
      print*, sigma(xx,ww), ww
      !return 
c      k(1) = 0.5d0
c      k(2) = 1.5d0
c      kp(1) = -1d0
c      kp(2) = 1d0
c      k1(1) = 0.1d0
c      k1(2) = 0.1d0
c      k2(1) = 7d0
c      k2(2) = 0d0
c      k3(1) = -0.5d0
c      k3(2) = -0.3d0
c      b1 = 0.001d0
c      b2 = 0.0004d0
    
c      a1=CHIDedotdiff(k1,k2,k1,k2)/s/b1
c      a2=CHIDedotdiff(k1,k2,k1,k2)/s/b2
     
      
c      call CHIDeDiphoton(dsigma,k,kp,k1,k2,k3,b1,b2,a1,a2)
c      print*, dsigma
c      stop
             
       nprn=0
       ncall=500000 
       itmx=5
       !itmx=10
       verbose=0
       call VEGAS(11,sigma,avgi,sd,chi2a)

       verbose=0
       print*, avgi, " nb +-", sd  
       !ncall=200
       ncall=100000
       itmx=1
       !itmx=10
       call VEGAS1(11,sigma,avgi,sd,chi2a)

       print*,
       print*,
       print*, 1d3*avgi, " fb +-", 1e3*sd  
       end





       
      function sigma(x,wgt)
      implicit none
      integer N,NN,i,ii
      integer verbose
      double precision sigma
      double precision mean
      double precision sum,max
      double precision x(11),wgt
      double precision jac,kmax,pi
      double precision thetak, thetap, theta1, theta2, theta3
      double precision ak, akp, ak1, ak2, ak3
      double precision b1,b2,a1,a2
      double precision k(2), kp(2), k1(2), k2(2), k3(2)
      double precision AA, BB, CC, DD, bmin,bmax
      
      double precision ss,ptmin
      common/aaa/ss,ptmin
      common/bbb/verbose

      
      double precision  CHIDedotdiff, CHIDedotsum
      external  CHIDedotdiff, CHIDedotsum

      kmax = 10d0
      pi=3.14159d0
      jac = 1d0
      bmin=0.0001d0
      bmax=1d0
      
      ak=kmax*x(1)
      thetak=2.*pi*x(2)
      jac=jac*2.*pi*ak*kmax
      k(1)=ak*cos(thetak)
      k(2)=ak*sin(thetak)
      
      akp=kmax*x(3)
      thetap=2.*pi*x(4)
      jac=jac*2.*pi*akp*kmax
      kp(1)=akp*cos(thetap)
      kp(2)=akp*sin(thetap)
 
      ak1=kmax*x(5)
      theta1=2.*pi*x(6)
      jac=jac*2.*pi*ak1*kmax
      k1(1)=ak1*cos(theta1)
      k1(2)=ak1*sin(theta1)

      ak2=5.+10.*x(7)*kmax
      theta2=2.*pi*0d0
      jac=jac*2.*pi*ak2*kmax*10.
      k2(1)=ak2*cos(theta2)
      k2(2)=ak2*sin(theta2)

      ak3=kmax*x(8)
      theta3=2.*pi*x(9)
      jac=jac*2.*pi*ak3*kmax
      k3(1)=ak3*cos(theta3)
      k3(2)=ak3*sin(theta3)
      
c     b1 = x(10)    
c     b2 = x(11)
      b1 = (bMIN/bmAX)**x(10)*bMAX
      b2 = (bMIN/bMAX)**x(11)*bMAX

      jac = jac*DLOG(bMAX/(bMIN))*B1
      jac = jac*DLOG(bMAX/(bMIN))*B2

      a1=CHIDedotsum(k1,k2,k1,k2)/b1/ss
      a2=CHIDedotsum(k2,k3,k2,k3)/b2/ss

      call CHIDeDiphoton(sigma,k,kp,k1,k2,k3,b1,b2,a1,a2)
      
      sigma = sigma*jac
      if(abs(sigma).LT.1.d-25)sigma=0.d0
      if(verbose.eq.1.and.abs(sigma).GT.1d-30) print*,
     & "AAAA", a1, a2, b1, b2, ak1, ak2, ak3, wgt*sigma

      return
      end
