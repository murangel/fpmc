C   ============================================================================   C  
C   SUBROUTINE INIT:                                                               C
C   Initialization of parameters and constants                                     C
C   Content: value for gap survival, possibily to change the scales in the         C
C            the Sudakov form factor by changing the value of x, values for cuts,  C
C            parameters for VEGAS,...                                              C
C   =============================================================================  C
      subroutine CHIDeHiggsInit(mH_,mt_,s_,iglu_,x_,xp_,s2_)
      implicit none
      double precision mH_,mt_,s_,x_,xp_,s2_
      integer iglu_ 
      integer ncall1,ncall2,itmx1,itmx2
      integer ncall,itmx,nprn,ndev
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol
      integer Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      integer i,j,iglu
      double precision mh
      double precision xl,xu,acc
      double precision kmax,pi,s,ncolor,nb,Gf,mt,mp
      double precision c,ael,gq
      double precision acut,bcut
      double precision hacutd,lacutd,hbcutd,lbcutd,ycut
      double precision hygap,lygap
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      double precision facwgt
      double precision mmu2,s2
      double precision x,xp
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
    

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/paramint/ncall1,ncall2,itmx1,itmx2   
      common/const/pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      common/Ifactor/c,ael,iglu
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps 
      common/histo/nbhis,y,xhis,count,facwgt,ave
      common/cuts/acut,bcut,hacutd,lacutd,hbcutd,lbcutd,ycut,
     &hygap,lygap
      common/uncert/x,xp
      common/sudaD/mmu2,s2
      common/cinterpol/dgdgrid,sudagrid
      common/higgs/mh
      character*500  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
      common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab

      mh=mH_

      exp=2
      Lenergy=1
      Lcut=3
      Model=2

      fullvertex=2
      formfac=2
      event=2
      iglu=iglu_       
      sudaf=6
      evolv=1
      interpol=1
      scale=2
      gaps=1

C     Constants
      pi=atan(1.d0)*4
      s=s_

      ncolor=3.
      nb=389379.D0
      Gf=1.16639E-5
      mt=mt_
      mp=0.938d0

      s2 = s2_ 

C   Scale in Sudakov: uncertainties evaluation
      x=x_
      xp=xp_ 

C   Cuts for LCWF
C   Typical cuts alpha,beta < 0.1 lost by the p+: 3 quarks so 0.3
      bcut=0.3d0
      acut=0.3d0

C     Cuts (no cuts)
      hacutd=1.d0
      lacutd=0.d0
      hbcutd=1.d0
      lbcutd=0.d0
      ycut=20.d0          

C     Parameters from CH paper
      c=-0.41d0
      ael=0.86d0
      

C   All tabs are made from [interpol.f]
      if(iglu.EQ.1)then
         open(unit = 99, file = dgdtab1, status = "old")
      elseif(iglu.EQ.2)then
         open(unit = 99, file = dgdtab2, status = "old")
      elseif(iglu.EQ.3)then
         open(unit = 99, file = dgdtab3, status = "old")
      elseif(iglu.EQ.4)then
         open(unit = 99, file = dgdtab4, status = "old")
      endif
      
      open(unit = 98, file = sudatab, status = "old")
      
      print*, ""
      print*, ""
      print*, "  __________________________________ "
      print*, " |                                  |"
      print*, " |   Initialisation of CHIDe model  |"
      print*, " |__________________________________|"
      print*, " |                                  |"
      print*, " |         pp -> p + H + p          |"
      print*, " |__________________________________|" 
      print*, " |                                  |"
      print*, " |  Using parameters:               |" 
      write(*,4000) 'Higgs mass', mh
      write(*,4000) 'Top mass', mt
      write(*,3000) 'sqrt(s)', sqrt(s)
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
4000  format ("  |   ",A20," = ", F6.2, "  |")


C     Read the tab and put it in the common
      do i = 0, 1000
       do j = 0, 1000
        read(99,*,end=999) dgdgrid(i+1,j+1)
        read(98,*,end=999) sudagrid(i+1,j+1)
       end do
      end do

 999  continue  

      end
