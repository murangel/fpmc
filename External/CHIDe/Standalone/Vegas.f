
       program StandaloneCHIDeGG
       implicit none

      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,scale,splash
      integer Model,evolv,Vgrid,Lcut,interpol,moode,lhe,nscale
      integer iglu,j,i
      integer update
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

      double precision CHIDeJetsdsigma

      external CHIDeJetsdsigma

     
      double precision CSGG, CSH
       double precision sigmaGG, sigmaHiggs
       external sigmaGG, sigmaHiggs
       character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       dgdtab4="Data/Higgsdgdtab1.d"
       dgdtab4="Data/Higgsdgdtab2.d"
       dgdtab4="Data/Higgsdgdtab3.d"
       dgdtab4="Data/Higgsdgdtab4.d"
       sudatab="Data/Higgssudatab.d"
       
       CALL CHIDeHiggsInit
     & (120.0d0,173.3d0,14000.0d0**2,4,1d0,1d0,0.075d0)

       nprn=0
       ncall=100000 
       itmx=10
       call VEGAS(10,sigmaHiggs,avgi,sd,chi2a)
       print*, avgi, "+-", sd  
       ncall=200000
       itmx=10
       call VEGAS1(10,sigmaHiggs,avgi,sd,chi2a)
       print*, avgi, "+-", sd  
       CSH = avgi
       

       dgdtab4="Data/ggdgdtab1.d"
       dgdtab4="Data/ggdgdtab2.d"
       dgdtab4="Data/ggdgdtab3.d"
       dgdtab4="Data/ggdgdtab4.d"
       sudatab="Data/ggsudatab.d"
       
      
       CALL CHIDeGGInit(10d0,1960.0d0**2,0.15d0)
       nprn=0
       ncall=100000
       itmx=10
       call VEGAS(11,sigmaGG,avgi,sd,chi2a)

       print*, avgi, "+-", sd  
       ncall=200000
       itmx=10
       call VEGAS1(11,sigmaGG,avgi,sd,chi2a)

       print*, avgi, "+-", sd  
       CSGG = avgi




       print*, "GG:",CSGG
       print*, "H: ",CSH
       
       end





       
      function sigmaGG(x,wgt)
      implicit none
      integer N,NN,i,ii
      double precision sigmaGG
      double precision mean
      double precision sum,max
      double precision x(11),wgt
      double precision jac,kmax,pi
      double precision thetak, thetap, theta1, theta2, theta3
      double precision ak, akp, ak1, ak2, ak3
      double precision b1,b2,a1,a2
      double precision k(2), kp(2), k1(2), k2(2), k3(2)
      double precision AA, BB, CC, DD
      
      
      double precision  CHIDedotdiff
      external  CHIDedotdiff

      AA = 10d0
      BB = 2d0
      CC = 14d0
      DD = 8d0
      pi=3.14159d0
      jac = 1d0
      
      ak=AA*x(1)
      thetak=2.*pi*x(2)
      jac=jac*2.*pi*ak*AA
      k(1)=ak*cos(thetak)
      k(2)=ak*sin(thetak)
      
      akp=AA*x(3)
      thetap=2.*pi*x(4)
      jac=jac*2.*pi*akp*AA
      kp(1)=akp*cos(thetap)
      kp(2)=akp*sin(thetap)
 
      ak1=BB*x(5)
      theta1=2.*pi*x(6)
      jac=jac*2.*pi*ak1*BB
      k1(1)=ak1*cos(theta1)
      k1(2)=ak1*sin(theta1)

      ak2=CC*x(7)+DD
      theta2=2.*pi*0d0
      jac=jac*2.*pi*ak2*CC
      k2(1)=ak2*cos(theta2)
      k2(2)=ak2*sin(theta2)

      ak3=BB*x(8)
      theta3=2.*pi*x(9)
      jac=jac*2.*pi*ak3*BB
      k3(1)=ak3*cos(theta3)
      k3(2)=ak3*sin(theta3)
      
      b1 = x(10)    
      b2 = x(11)
      
      a1=CHIDedotdiff(k1,k2,k1,k2)/1960d0**2/b1
      a2=CHIDedotdiff(k3,k2,k3,k2)/1960d0**2/b2

      call CHIDeGG(sigmaGG,k,kp,k1,k2,k3,b1,b2,a1,a2)
      
      sigmaGG = sigmaGG*jac
      
      return
      end



          function sigmaHiggs(x,wgt)
      implicit none
      integer n,nn,i,ii
      double precision sigmahiggs
      double precision mean
      double precision sum,max
      double precision x(11),wgt
      double precision jac,kmax,pi
      double precision thetak, thetap, theta1, theta2, theta3
      double precision ak, akp, ak1, ak2, ak3
      double precision b1,a3
      double precision k(2), kp(2), k1(2), k2(2), k3(2)
      double precision aa, bb, cc, dd
      
      
      double precision  chidedotdiff
      external  chidedotdiff

      aa = 20d0
      bb = 10d0
      pi=3.14159d0
      jac = 1d0
      
      ak=aa*x(1)
      thetak=2.*pi*x(2)
      jac=jac*2.*pi*ak*aa
      k(1)=ak*cos(thetak)
      k(2)=ak*sin(thetak)
      
      akp=aa*x(3)
      thetap=2.*pi*x(4)
      jac=jac*2.*pi*akp*aa
      kp(1)=akp*cos(thetap)
      kp(2)=akp*sin(thetap)
 
      ak1=bb*x(5)
      theta1=2.*pi*x(6)
      jac=jac*2.*pi*ak1*bb
      k1(1)=ak1*cos(theta1)
      k1(2)=ak1*sin(theta1)

      ak3=bb*x(8)
      theta3=2.*pi*x(9)
      jac=jac*2.*pi*ak3*bb
      k3(1)=ak3*cos(theta3)
      k3(2)=ak3*sin(theta3)
      
      b1 = x(9)    
      a3 = x(10)
      
c      a3=0.01
c      b1=0.01
c      k1(1)=1d0
c      k1(2)=1d0
c      k3(1)=1d0
c      k3(2)=1d0
c      k(1)=2d0
c      k(2)=2d0
c      kp(1)=2d0
c      kp(2)=2d0
c                                     
c                   call chidehiggs(sigma,a3,b1,
c     &                             k1,k3,k,kp)
c                 print*, "sigma =", sigma
c                   stop
c  
      call chidehiggs(sigmahiggs,a3,b1,k1,k3,k,kp)
c     print*, "sigma =", sigma

c      stop
          
      sigmahiggs = sigmahiggs*jac
      
      return
      end
