C   ============================================================================   C  
C   SUBROUTINE INIT:                                                               C
C   Initialization of parameters and constants                                     C
C   Content: value for gap survival, possibily to change the scales in the         C
C            the Sudakov form factor by changing the value of x, values for cuts,  C
C            parameters for VEGAS,...                                              C
C   =============================================================================  C
      subroutine CHIDeDiphotoninit(s_,Etmin_,s2_,x_,xp_,iglu_,
     &                       xi1min,xi1max,xi2min,xi2max,
     &                       ymin,ymax)
      implicit none
      integer iglu_
      double precision Etmin_,s_,s2_,x_,xp_
      double precision xi1min,xi1max,xi2min,xi2max
      double precision ymin,ymax
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer ncall,itmx,nprn,ndev
      integer event,exp,sudaf,grid,formfac,scale,Vgrid,Lcut
      integer Model,evolv,interpol,moode,lhe,splash,nscale
      integer i,j,iglu,nf
      double precision xl,xu,acc 
      double precision pi,s,ncolor,gg,gelm,gq,kmax,nb,mp
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      double precision facwgt,ael,aelm
      double precision etcut,hygap,lygap,hbcutd,lbcutd,ycut,yclucut
      double precision ylcut,hacutd,lacutd,Etmin,propcut
      double precision c
      double precision s2
      double precision xmax,amax,amin
      double precision x,xp,mmu
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4
      double precision Etjetmin,gapsurv
      double precision mq(1:6),Qup,Qdown

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist 
      common/histo/nbhis,y,xhis,count,facwgt,ave
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &ylcut,hacutd,lacutd,Etmin,yclucut,propcut
      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,aelm,iglu
      common/BotE/amax,xmax,amin
      common/uncert/x,xp
      common/cinterpol/dgdgrid,sudagrid
      common/log/pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4
      common/const/pi,ncolor,gg,gq,gelm,kmax,nb,mp
      common/flavor/mq,Qup,Qdown,nf
      
      common/sudaD/mmu,s2
      double precision Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      common/param/s
      common/algo/Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj

      character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
      common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab

      
      s=s_
      exp=2
      Vgrid=1
      model=2
      splash=1
      Etjetmin=Etmin_
      Etcut=Etjetmin
      Etjetveto=10.
      
      s2=s2_

      propcut=-10.
      formfac=2
      event=3
      iglu=iglu_
      sudaf=1
      evolv=1
      interpol=1
      scale=2
      nscale=1
      
      
      pi=atan(1.d0)*4

      nb=0.3894*1.E6
      mp=0.938d0
      ncolor=3.
      pi=atan(1.d0)*4



      ael=1.d0
      gq=sqrt(4.*pi*ael)
      aelm=1./137.036d0
      gelm=sqrt(4.*pi*aelm)
      

C  Quarks masse and charge
      mq(2)=0.312d0
      mq(3)=0.45d0
      mq(4)=1.5d0
      mq(5)=4.5d0
      mq(6)=173.4d0
      mq(1)=mq(2)
      Qup=2./3.d0
      Qdown=-1./3.d0

      kmax=50.

      x=x_
      xp=xp_

      ylcut=ymin
      ycut=ymax                           ! Cut on rapidity: ycut on jets = 2.5
c     lygap=gapmin                        ! Size of the gap: jets
c     hygap=gapmax                        !    3.6 < eta(gap) < 5.9
      hacutd=xi2max                       ! Higher cut on transfered energy from the p+
      lacutd=xi2min                       ! Lower cut on transfered energy from the p+
      hbcutd=xi1max                       ! Higher cut on transfered energy from the p-
      lbcutd=xi1min                       ! Lower cut on transfered energy from the p-



C     Chose of dgdtab as a function of iglu
       if(iglu.EQ.1)then
        open(unit = 99, file = dgdtab1, status = "old")
       elseif(iglu.EQ.2)then
        open(unit = 99, file = dgdtab2, status = "old")
       elseif(iglu.EQ.3)then
        open(unit = 99, file = dgdtab3, status = "old") 
       elseif(iglu.EQ.4)then
        open(unit = 99, file = dgdtab4, status = "old")
       endif
C     Open sudatab 
      open(unit = 98, file = sudatab, status = "old")
      
       
       
      print*, ""
      print*, ""
      print*, "  __________________________________ "
      print*, " |                                  |"
      print*, " |   Initialisation of CHIDe model  |"
      print*, " |__________________________________|"
      print*, " |                                  |"
      print*, " |    pp -> p + gammagamma + p      |"
      print*, " |__________________________________|" 
      print*, " |                                  |"
      print*, " |  Using parameters:               |" 
      write(*,3000) 'sqrt(s)', sqrt(s)
      write(*,3000) 'Etjetmin', Etjetmin
      write(*,4000) 'xi_1 min', lbcutd
      write(*,4000) 'xi_1 max', hbcutd
      write(*,4000) 'xi_2 min', lacutd
      write(*,4000) 'xi_2 max', hacutd
      write(*,5000) 'y min', ylcut
      write(*,5000) 'y max', ycut
      write(*,1000) 'gluon density', iglu
      write(*,2000) 'gap survival', s2
      write(*,2000) 'upper sudakov limit', x
      write(*,2000) 'lower sudakov limit', xp
      print*, " |__________________________________|" 
      print*, ""
      print*, ""
1000  format ("  |   ",A20," = ", I2, "      |")
2000  format ("  |   ",A20," = ", F6.3, "  |")
3000  format ("  |   ",A20," = ", F6.0, "  |")
4000  format ("  |   ",A20," = ", F6.4, "  |")
5000  format ("  |   ",A20," = ", F6.2, "  |")



C     Read the tab and put it in the common
       do i = 0, 1000
        do j = 0, 1000
         read(99,*) dgdgrid(i+1,j+1)
         read(98,*) sudagrid(i+1,j+1)
        end do
       end do
      continue



      end
