c////////////////////////////////////////////////
c// @purpose: Simple ntuple
c//
c// @date:    20/02/2011 
c// @author:  O. Kepka (kepkao@fzu.cz)
c////////////////////////////////////////////////


c______________________________________________________________________________ 
c___ create ntuple

      SUBROUTINE NTINIT
c______________________________________________________________________________ 

      implicit none
      INCLUDE '../Examples/ffcard.inc'
      REAL hmemor
      COMMON/pawc/hmemor(1000000)
 
***** common block - generator level particles 
       integer ngenmax
       parameter(ngenmax=1000)
       integer ngen
       real px(ngenmax),py(ngenmax),pz(ngenmax)
       real e(ngenmax),m(ngenmax)
       integer id(ngenmax)
       integer ii(ngenmax)
       integer ist(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,m,id,ii,ist

***** common block - bjorken-x of the parton from pomeron         
c       common /remnant/xg1b,xg2b


c***** jets 
      integer njetmax
      parameter(njetmax=30)
      integer njets
      real
     &     pxjet(njetmax),pyjet(njetmax),
     &     pzjet(njetmax),ejet(njetmax)
      common/jets/njets,
     &     pxjet,pyjet,
     &     pzjet,ejet


      INTEGER MODE, ICYCLE, ISTAT


*
      print '(A)',' hropen ...'


      call hbnt(777,'ntuple',' ')
      
      print *,'Creating the ntuple 777'


c       call hbname(777,'remnant',xg1b,'xg1,'//
c     & 'xg2')

       
     
       call hbname(777,'gener',ngen,'ngen[0,1000]:I,'//
     &     'px(ngen),py(ngen),pz(ngen),e(ngen),rm(ngen),id(ngen),'//
     &     'ii(ngen), ist(ngen)')

       if(UOUTPUT.ge.2)
     &  call hbname(777,'jets',njets,'njets[0,30]:I,'//            
     & 'pxjet(njets),pyjet(njets),'//              
     & 'pzjet(njets),pejet(njets)')



      END 
      
c______________________________________________________________________________ 
      
      SUBROUTINE NTEND
c______________________________________________________________________________ 

*---Close the ntuple:
      INTEGER idbg,igener,irad,ifrad
      COMMON /CONST1/ idbg,igener,irad,ifrad

      call hrout(777,icycle,' ')

      END

