* OK: 1/12/2006 interface to the pxcone jet algo. 
* 1) convert hep to structure of pxcone(pythia typ)
* 2) do the clustering
* 3) fill the appropriate structure in ntuples      
      SUBROUTINE MYRECO(IERR)
      !IERR = 0 OK
      !IERR!= 0 error

      include 'HERWIG65.INC'
      integer nfoundjet

      
C ... PXCONE variables, two pxcone algorithms for jet reconstructions
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
       real e(ngenmax),m(ngenmax)
       integer id(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,m,id

c ... jets 
      integer njetmax
      parameter(njetmax=30)
      integer njets
      real
     &     pxjet(njetmax),pyjet(njetmax),
     &     pzjet(njetmax),ejet(njetmax)
      common/jets/njets,
     &     pxjet,pyjet,
     &     pzjet,ejet



      data nevhep /0/

c ... variables for selecting leptons
      real*8 pt, pttresh
      integer absid
      integer nleptons
      INTEGER MXLEP
      PARAMETER (MXLEP=20)
      REAL*8 PLEP(5,MXLEP)
c ... local variables for misset
      real*8 misspx, misspy


c ... loop variables
      INTEGER I, J, N, IPART
      DOUBLE PRECISION RAP

      integer idiv

      
c------------------------------------------------------------------- 

C ... PXCONE parameters ...
      MODE    = 2        !Snow mass scheme
      CONER   = 0.7
      EPSLON  = 4
c     EPSLON  = 8
      OVLIM   = 0.5d0

c------------------------------------------------------------------- 
c particles on the generator level
      call vzero(px,ngenmax)
      call vzero(py,ngenmax)
      call vzero(pz,ngenmax)
      call vzero(e,ngenmax)
      call vzero(m,ngenmax)
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
             m(IPART)=sngl(PHEP(5,I))
             id(IPART)=IDHEP(I)
c         rm(i)=sngl(PTRAK(5,i))
c          print '(A,4F8.2)','gen:',px(i),py(i),pz(i),e(i)
c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
          ENDIF
1515  continue
      ngen=IPART


      idiv = MOD(nevhep, 1000)
      if(nevhep.lt.100) then 
        print *,' *** cone_interface nevhep :', nevhep
      else  if ( idiv.EQ.0 ) then
            print *,' ***  cone_interface nevhep :', nevhep
      ENDIF
      
          
c-------------------------------------------------------------------          
      NTRAK=0
      N=NHEP
c     print '(A,I6)', 'NHEP', NHEP
      DO 180 I=1,N
          !if particle in the final state
          IF(ISTHEP(I).EQ.1) THEN
             RAP = DABS( 0.5*DLOG( (PHEP(4,I)+PHEP(3,I))/
     .                (PHEP(4,I)-PHEP(3,I)) ) )
c             print*,rap
             IF(RAP.LT.4) THEN
                NTRAK=NTRAK+1
                DO 160 J=1,5
                  PTRAK(J,NTRAK)=PHEP(J,I)
 160            CONTINUE
             ENDIF     
               
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





c       fill jets
      do 1550 i=1, njet
         pxjet(i) = PJET(1,i)
         pyjet(i) = PJET(2,i)
         pzjet(i) = PJET(3,i)
         ejet(i)  = PJET(4,i)
c         print '(A, 4F8.2)', 'pjet x y z e :', PJET(1,i),PJET(2,i),
c     .                                         PJET(3,i),PJET(4,i)
c         print *,''
 1550 continue         
       njets = njet
  
       
       RETURN
       END

      


**********************************************
* R.S. 05.06.2012                            *
* MYRECO copy with additional parameter that *
* decides whether to run jet cone algorithm. *
**********************************************
      
      SUBROUTINE FPMC_RECO(output,IERR)
      !IERR = 0 OK
      !IERR!= 0 error

      include 'HERWIG65.INC'
      integer output
      integer nfoundjet

C ... PXCONE variables, two pxcone algorithms for jet reconstructions
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
      real e(ngenmax),m(ngenmax)
      integer id(ngenmax)
      integer ii(ngenmax)
      integer ist(ngenmax)
       
      common /gener/ngen,
     &     px,py,pz,e,m,id,ii,ist

c ... jets 
      integer njetmax
      parameter(njetmax=30)
      integer njets
      real
     &     pxjet(njetmax),pyjet(njetmax),
     &     pzjet(njetmax),ejet(njetmax)
      common/jets/njets,
     &     pxjet,pyjet,
     &     pzjet,ejet

      data nevhep /0/

c ... variables for selecting leptons
      real*8 pt, pttresh
      integer absid
      integer nleptons
      INTEGER MXLEP
      PARAMETER (MXLEP=20)
      REAL*8 PLEP(5,MXLEP)
c ... local variables for misset
      real*8 misspx, misspy

c ... loop variables
      INTEGER I, J, N, IPART
      DOUBLE PRECISION RAP

      integer idiv

c------------------------------------------------------------------- 

C ... PXCONE parameters ...
      MODE    = 2        !Snow mass scheme
      CONER   = 0.7
      EPSLON  = 4
c     EPSLON  = 8
      OVLIM   = 0.5d0

c------------------------------------------------------------------- 
c particles on the generator level
      call vzero(px,ngenmax)
      call vzero(py,ngenmax)
      call vzero(pz,ngenmax)
      call vzero(e,ngenmax)
      call vzero(m,ngenmax)
      call vzero(id,ngenmax)
        
      N=NHEP
      IPART=0
      do 1515 I=1,N
          IF((ISTHEP(I).EQ.1.OR.IDHEP(I).EQ.5.OR.IDHEP(I).EQ.-5)) THEN
             IPART=IPART+1
             px(IPART)=sngl(PHEP(1,I))
             py(IPART)=sngl(PHEP(2,I))
             pz(IPART)=sngl(PHEP(3,I))
             e(IPART) =sngl(PHEP(4,I))
             m(IPART)=sngl(PHEP(5,I))
             id(IPART)=IDHEP(I)
             ii(IPART) = I
             ist(IPART) = ISTHEP(I)
c         rm(i)=sngl(PTRAK(5,i))
c          print '(A,4F8.2)','gen:',px(i),py(i),pz(i),e(i)
c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
          ENDIF
1515  continue
      ngen=IPART


c      idiv = MOD(nevhep, 1000)
c      if(nevhep.lt.100) then 
c        print *,' *** cone_interface nevhep :', nevhep
c      else  if ( idiv.EQ.0 ) then
c            print *,' ***  cone_interface nevhep :', nevhep
c      ENDIF
      
          
c-------------------------------------------------------------------          
      NTRAK=0
      N=NHEP
c     print '(A,I6)', 'NHEP', NHEP
      DO 180 I=1,N
          !if particle in the final state
          IF(ISTHEP(I).EQ.1 .and. ABS(IDPDG(I)).GT.100 ) THEN
            RAP = DABS( 0.5*DLOG( (PHEP(4,I)+PHEP(3,I))/
     .                (PHEP(4,I)-PHEP(3,I)) ) )
            IF(RAP.LT.5) THEN
                NTRAK=NTRAK+1
                DO 160 J=1,5
                  PTRAK(J,NTRAK)=PHEP(J,I)
 160            CONTINUE
            ENDIF     
               
c          print '(A,4F8.2)','phep:',PHEP(1,i),PHEP(2,i),
c     .                              PHEP(3,i),PHEP(4,i)
               
          ENDIF     
  180 CONTINUE       


c-------------------------------------------------------------------           
      !clustering
      if(output.ge.2) then
      
      NJET=0                          !midpoint alg.
      CALL PXCONE(MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET,
     +            NJET,PJET,IPASS,IJMUL,IERR)

      IF(IERR.NE.0) THEN
         print *, 'PXCONE did not converge'
         RETURN
      ENDIF   

c       fill jets
      do 1550 i=1, njet
         pxjet(i) = PJET(1,i)
         pyjet(i) = PJET(2,i)
         pzjet(i) = PJET(3,i)
         ejet(i)  = PJET(4,i)
 1550 continue         
         njets = njet
      
      endif
 
       
       RETURN
       END

      

