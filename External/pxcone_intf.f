* OK: 1/12/2006 interface to the pxcone jet algo. 
* 1) convert hep to structure of pxcone(pythia typ)
* 2) do the clustering
* 3) fill the appropriate structure in ntuples      

      SUBROUTINE CONE(IERR)
      !IERR = 0 OK
      !IERR!= 0 error

      include 'herwig6500.inc'
      integer nfoundjet

      
C ... PXCONE variables, two pxcone algorithms 
      INTEGER  ITKDM,MXTRK
      ! ... ctyrvektor, max. pocet castic v jednom eventu
      PARAMETER  (ITKDM=4,MXTRK=4000)
      INTEGER  MXJET, MXTRAK, MXPROT
      PARAMETER  (MXJET=20,MXTRAK=4000,MXPROT=500)
      INTEGER  IPASS (MXTRAK),IJMUL (MXJET)
      INTEGER  NTRAK,MODE, NJET, IERR
      DOUBLE PRECISION  PTRAK(ITKDM,MXTRK),PJET(5,MXJET)
      DOUBLE PRECISION  CONER, EPSLON, OVLIM

      
      
c ... particles on generator level     
       integer ngenmax
       parameter(ngenmax=1000)
       integer ngen
       real px(ngenmax),py(ngenmax),pz(ngenmax)
       real e(ngenmax),rm(ngenmax)
       integer id(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,rm,id

c ... jet structure (as in simul)
      integer nrecmax
      parameter(nrecmax=100)
      integer nrec
      real typrec(nrecmax),pxrec(nrecmax),pyrec(nrecmax)
      real pzrec(nrecmax)
      real erec(nrecmax),qrec(nrecmax),eemrec(nrecmax),
     &ehadrec(nrecmax)
      real etrarec(nrecmax),widrec(nrecmax),var6(nrecmax)
      real ctagrec(nrecmax),btag1rec(nrecmax),btag2rec(nrecmax)
      integer ntrarec(nrecmax)
      real sum15ec(nrecmax),sum15hc(nrecmax), sum40(nrecmax)
      real sum40ec(nrecmax)
      
      common /recons1/ nrec,
     & typrec,pxrec,pyrec,pzrec,erec,qrec,eemrec,ehadrec,
     & etrarec,widrec,var6,
     & ctagrec,btag1rec,btag2rec,
     & ntrarec,
     & sum15ec,sum15hc, sum40, sum40ec

      data nevhep /0/

c ... loop variables
      INTEGER I, J, N, IPART
      DOUBLE PRECISION RAP


C...HEPEVT commonblock.
c      PARAMETER (NMXHEP=4000)
c      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
c      DOUBLE PRECISION PHEP,VHEP
c      SAVE /HEPEVT/
      
c------------------------------------------------------------------- 

      
C ... PXCONE parameters ...
      MODE    = 2        !Snow mass scheme
      CONER   = 0.7
      EPSLON  = 8
      OVLIM   = 0.5d0

c------------------------------------------------------------------- 
c particles on the generator level
      call vzero(px,ngenmax)
      call vzero(py,ngenmax)
      call vzero(pz,ngenmax)
      call vzero(e,ngenmax)
      call vzero(rm,ngenmax)
      call vzero(id,ngenmax)
        
        
      N=NHEP
      IPART=0
      do 1515 I=1,N
          IF(ISTHEP(I).EQ.1) THEN
             IPART=IPART+1
             px(IPART)=sngl(PHEP(1,I))
             py(IPART)=sngl(PHEP(2,I))
             pz(IPART)=sngl(PHEP(3,I))
             e(IPART) =sngl(PHEP(4,I))
c         rm(i)=sngl(PTRAK(5,i))
c          print '(A,4F8.2)','gen:',px(i),py(i),pz(i),e(i)
c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
          ENDIF
1515  continue
      ngen=IPART
          
c-------------------------------------------------------------------          
      NTRAK=0
      N=NHEP
c     print '(A,I6)', 'NHEP', NHEP
      DO 180 I=1,N
          !if particle in the final state
          IF(ISTHEP(I).EQ.1) THEN
c             RAP = DABS( 0.5*DLOG( (PHEP(4,I)+PHEP(3,I))/
c     .                (PHEP(4,I)-PHEP(3,I)) ) )
c             print*,rap
c             IF(RAP.LT.4) THEN
                NTRAK=NTRAK+1
                DO 160 J=1,5
                  PTRAK(J,NTRAK)=PHEP(J,I)
 160            CONTINUE
c             ENDIF     
               
c          print '(A,4F8.2)','phep:',PHEP(1,i),PHEP(2,i),
c     .                              PHEP(3,i),PHEP(4,i)
               
          ENDIF     
  180 CONTINUE       


c-------------------------------------------------------------------           
      !clustering
      
      NJET=0                          !midpoint alg.
      CALL PXCONE(MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET,
     +            NJET,PJET,IPASS,IJMUL,IERR)

      IF(IERR.NE.0) THEN
         print *, 'PXCONE did not converge'
         RETURN
      ENDIF   


      print *,' *** cone_interface nevhep :', nevhep

c...jets common blocks
      call vzero(typrec,nrecmax)
      call vzero(pxrec,nrecmax)
      call vzero(pyrec,nrecmax)
      call vzero(pzrec,nrecmax)
      call vzero(erec,nrecmax)
      call vzero(qrec,nrecmax)
      call vzero(eemrec,nrecmax)
      call vzero(ehadrec,nrecmax)
      call vzero(etrarec,nrecmax)
      call vzero(widrec,nrecmax)
      call vzero(ctagrec,nrecmax)
      call vzero(btag1rec,nrecmax)
      call vzero(btag2rec,nrecmax)
      call vzero(ntrarec,nrecmax)
      call vzero(sum15ec,nrecmax)
      call vzero(sum15hc,nrecmax)
      call vzero(sum40,nrecmax)
      call vzero(sum40ec,nrecmax)

      !assumtion all jet hadronic
      print '(A,I6)', 'njet', njet
      do 1550 i=1, njet
         typrec(i)= 4
         eemrec(i)= 0
         pxrec(i) = PJET(1,i)
         pyrec(i) = PJET(2,i)
         pzrec(i) = PJET(3,i)
         erec(i)  = PJET(4,i)
c         print '(A, 4F8.2)', 'pjet x y z e :', PJET(1,i),PJET(2,i),
c     .                                         PJET(3,i),PJET(4,i)
c         print *,''
 1550 continue         
      nrec = njet
      
      
       
       RETURN
       END

      

