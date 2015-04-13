c     rutine to initialize ntuples for detector simulation

      SUBROUTINE NTINIT


      REAL hmemor
      COMMON/pawc/hmemor(8000000)
 
c______________________________________________________________________________       
c     common block for ntuples
c     exactly the same structure as simul_interf_cms.f


***** gener common block
       integer ngenmax
       parameter(ngenmax=1000)
       integer ngen
       real px(ngenmax),py(ngenmax),pz(ngenmax)
       real e(ngenmax),rm(ngenmax)
       integer id(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,rm,id

***** remnant common block        
       common /remnant/xg1b,xg2b


***** renostruction common block      
C CHR 25/05/99 modif pour ntuple reconstruit
      integer nrecmax
      parameter(nrecmax=100)
      integer nrec
      real typrec(nrecmax),pxrec(nrecmax),pyrec(nrecmax)
      real pzrec(nrecmax)
       real erec(nrecmax),qrec(nrecmax),eemrec(nrecmax),
     & ehadrec(nrecmax)
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
     

***** global common block      
      real circul,transen,transmass,ptmisstx,ptmissty
      real ptmissx,ptmissy,pznu1,pznu2
       
      common /global/ circul,transen,transmass,ptmisstx,ptmissty,
     . ptmissx,ptmissy,pznu1,pznu2


******** track common block
*       integer ntramax
*       parameter(ntramax=150)
*       integer ntra
*       real pxtra(ntramax),pytra(ntramax),pztra(ntramax)
*       real qtra(ntramax)
       
*       common /track/ ntra,pxtra,pytra,pztra,qtra


***** cluster common block      
c       integer nclmax
c       parameter(nclmax=500)
c       integer nclu
c       real pxrclu(nclmax),pyrclu(nclmax),pzrclu(nclmax)
c       real eclu(nclmax)
c       real eemclu(nclmax)
c       real ehadclu(nclmax),widtclu(nclmax),rmisset,phimet
c       integer multclu(nclmax)
       
c       common /cluster/ nclu,pxrclu,pyrclu,pzrclu,eclu,eemclu,
c     &      ehadclu,widtclu,multclu,rmisset,phimet
     
     
***** trigger common block
c       integer ntrimax
c       parameter(ntrimax=50)
c       integer ntri
c       real typtri(ntrimax),etatri(ntrimax),phitri(ntrimax)
c       real ettri(ntrimax)
c       real clutri(ntrimax)
c       real pttri(ntrimax),tratri(ntrimax)
       
c       common /trigger/ ntri,typtri,etatri,phitri,ettri,clutri,pttri,
c     &     tratri
     
***** cell common block      
c       integer ncellmax
c       parameter(ncellmax=5000)
c       integer ncell
c       real pxcell(ncellmax),pycell(ncellmax)
c       real etsheta(ncellmax)
c       real ecell(ncellmax),etcell(ncellmax)
c       real etacell(ncellmax)
c       real phicell(ncellmax),rmrkcell(ncellmax)
       
c       common /cells/ ncell, pxcell, pycell, etsheta, ecell, etcell,
c     +     etacell,phicell,rmrkcell    


c OK:10/06: following is not used also. 
***** jets : first jets from LUCELL, then LUCLUS
      integer njetmax
      parameter(njetmax=30)
      integer njetscel
      real
     &     pxjet(njetmax),pyjet(njetmax),
     &     pzjet(njetmax),pejet(njetmax),
     &     pmjet(njetmax)
      common/celljets/
     &     njetscel,
     &     pxjet,pyjet,
     &     pzjet,pejet,
     &     pmjet

      
      integer nclumax
      parameter(nclumax=30)
      integer njetsclu
      real
     &     pxclu(nclumax),pyclu(nclumax),
     &     pzclu(nclumax),peclu(nclumax),
     &     pmclu(nclumax)
      common/clujets/
     &     njetsclu,
     &     pxclu,pyclu,
     &     pzclu,peclu,
     &     pmclu

      INTEGER i
*...Max. number of preclusters/ clusters/jets
      INTEGER    NCLUS
      PARAMETER (NCLUS = 50)
*KEEP,PTMIS.
*...Et missing true (neutrinos), Et miss. measured, reconstructed Z-compon. of
*...neutrino from W -> l nu decay (if any).
      REAL            PTMT,    PTM,    PZNU
      COMMON /PTMIS/  PTMT(2), PTM(2), PZNU(2)
*KEEP,TRANCH.
*...Circularity, transverse energy and mass values calculated in sub. TRAPAR
      REAL            CIRC, TREN, TRMA
      COMMON /TRANCH/ CIRC, TREN, TRMA
*KEEP,LEPTS.
*...Desired leptons found by user (PROLEP, then after CISOL check and,
*...possibly, after GLOBJF),
*...PL - 4-VECTORS,  KL - code :  +-1(2) for e(mu)-+ ;
*...INDLEP - number of each lepton in /PLIST/, LHIT - in /HITPAR/,
*...MRL - origin (code of parent, =999 for SUSY particles),
*...LISOL - isolation mark (0 - isolated, 1 - non-isolated in tracker,
*... 2 - non-isolated in calo, 3 - non-isolated in both tracker and calo,
*... 4 - isolated in jet,),   NLE - number of leptons.
      INTEGER    NLMAX, KL, INDLEP, LHIT, MRL, LISOL, NLE
      PARAMETER (NLMAX = 100)
      REAL            PL
      COMMON /LEPTS/  PL(4,NLMAX),  KL(NLMAX), INDLEP(NLMAX),
     &                LHIT(NLMAX), MRL(NLMAX),  LISOL(NLMAX), NLE
*KEEP,GAMMAS.
*...Gammas found (PROGAM, then CISOL and possibly GLOBJF),
*...PG     - 4-VECTORS,
*...INDGAM - number of each lepton in /PLIST/,
*...IGHIT  - number in /HITPAR/,
*...LISGAM - isolation mark : 0 - isolated,
*...         2 - non-isolated in calo, 4 - isolated in jet;
*...NGAM   - number of gammas.
      INTEGER NGMAX, INDGAM, IGHIT, LISGAM, NGAM
      REAL    PG
      PARAMETER (NGMAX = 100)
      COMMON /GAMMAS/  PG(4,NGMAX), INDGAM(NGMAX),
     &                IGHIT(NGMAX), LISGAM(NGMAX), NGAM
*KEEP,JETOUT.
*...Jets found by subroutine GLOBJF
      INTEGER                       NJG
      REAL            PJG
      COMMON /JETOUT/ PJG(4,NCLUS), NJG
*KEEP,JEMARK.
*...Jet marks filled by sub. QJMAT : 0 - quark/gluon not assigned,
*... > 0 - PYTHIA`s code of the quark/gluon (without sign).
*...CALL TAUTAG can add 1000 to this mark of jet if recognized as 1-prong tau.
*...BTAG contains a likelihood of jet b-tagging (if IBTAG =/= 0)
      INTEGER         MRKJET
      REAL                           BTAG
      COMMON /JEMARK/ MRKJET(NCLUS), BTAG(NCLUS)
*KEEP,DNAMES.
*...Names and LUNs of the I/O files.
      INTEGER          LDAT1,  LDAT2,  LDAT3,  LDAT4,  LDAT5,  LDAT6,
     &                LRNDM1, LRNDM2, IDNTIN, IDNTOUT
      CHARACTER*80    DATNA1, DATNA2, DATNA3, DATNA4, DATNA5, DATNA6,
     &                RNDMN1, RNDMN2, HISTNA, HISTNB
      COMMON /DNAMES/  LDAT1,  LDAT2,  LDAT3,  LDAT4,  LDAT5,  LDAT6,
     &                LRNDM1, LRNDM2, IDNTIN, IDNTOUT,
     &                DATNA1, DATNA2, DATNA3, DATNA4, DATNA5, DATNA6,
     &                RNDMN1, RNDMN2, HISTNA, HISTNB
*KEND.
*............

      INTEGER          NLE_CWN, KL_CWN, MRL_CWN, LISOL_CWN
      REAL             PL_CWN
      COMMON /LEP_CWN/ NLE_CWN, PL_CWN(4,NLMAX), KL_CWN(NLMAX),
     &                 MRL_CWN(NLMAX), LISOL_CWN(NLMAX)

      INTEGER          NGAM_CWN, LISGAM_CWN
      REAL             PG_CWN
      COMMON /GAM_CWN/ NGAM_CWN, PG_CWN(4,NGMAX), LISGAM_CWN(NGMAX)

      INTEGER          NJG_CWN, MRKJET_CWN
      REAL             PJG_CWN,                   BTAG_CWN
      COMMON /JET_CWN/ NJG_CWN, PJG_CWN(4,NCLUS), MRKJET_CWN(NCLUS),
     &                                            BTAG_CWN(NCLUS)


      INTEGER MODE, ICYCLE, ISTAT



c______________________________________________________________________________ 
c create ntuple
*
      print '(A)',' hropen ...'

*      call hropen(31,'pythia','higgs.ntp','n',1024,istat)

c      IF (ISTAT.NE. 0) THEN
c       PRINT * ,'*** ERROR ',ISTAT,' DURING OPEN OF FILE 31'
c       CLOSE(31)
c       RETURN
c      ENDIF

c      print '(A)',' hbnt ...'

c      call hmdir('//pythia','s')
c      call hcdir('//pythia',' ')
      call hbnt(777,'ntuple',' ')
      
      print *,'on vient de creer le ntuple 777'

      

c      call hbname(777,'LUCELL',njetscel,'njetscel[0,30]:U,'//
c     &'pxjet(njetscel),pyjet(njetscel),'//
c     &'pzjet(njetscel),pejet(njetscel),pmjet(njetscel)')

c      call hbname(777,'LUCLUS',njetsclu,'njetsclu[0,30]:U,'//
c     &'pxclu(njetsclu),pyclu(njetsclu),'//
c     &'pzclu(njetsclu),peclu(njetsclu),pmclu(njetsclu)')

      call hbname(777,'recons1',nrec,'nrec[0,100]:I,'//
     & 'typrec(nrec),'//
     & 'pxrec(nrec),pyrec(nrec),pzrec(nrec),'//
     & 'erec(nrec),'//
     & 'qrec(nrec),eemrec(nrec),ehadrec(nrec),'//
     & 'etrarec(nrec),widrec(nrec),var6(nrec),'//
     & 'ctagrec(nrec),btag1rec(nrec),btag2rec(nrec),'//
     & 'ntrarec(nrec),'//
     & 'sum15ec(nrec),sum15hc(nrec), sum40(nrec), sum40ec(nrec)')
       
c       call hbname(777,'track',ntra,'ntra[0,100]:I,'//
c     & 'pxtra(ntra),pytra(ntra),pztra(ntra),qtra(ntra)')

       
        call hbname(777,'global',circul,'circul,transen,transmass,'//
     & 'ptmisstx,ptmissty,'//
     & 'ptmissx,ptmissy,pznu1,pznu2')


* do not put cells in root tuple       
c        call hbname(777,'cells',ncell,'ncell[0,5000]:I,'// 
c     & 'pxcell(ncell),'//
c     & 'pycell(ncell), etsheta(ncell), ecell(ncell), etcell(ncell),'//
c     & 'etacell(ncell),phicell(ncell),rmrkcell(ncell)')
   
     


       call hbname(777,'remnant',xg1b,'xg1,'//
     & 'xg2')

       
c       call hbname(777,'cluster',nclu,'nclu[0,500]:I,'//
c     &  'pxrclu(nclu),pyrclu(nclu),pzrclu(nclu),eclu(nclu),'//
c     &  'eemclu(nclu),ehadclu(nclu),'//
c     &  'widtclu(nclu),multclu(nclu),rmisset,phimet')
     
c       call hbname(777,'trigger',ntri,'ntri[0,50]:I,'//
c     & 'typtri(ntri),etatri(ntri),phitri(ntri),ettri(ntri),'//
c     &  'clutri(ntri),pttri(ntri),tratri(ntri)')
     
       call hbname(777,'gener',ngen,'ngen[0,1000]:I,'//
     &     'px(ngen),py(ngen),pz(ngen),e(ngen),rm(ngen),id(ngen)')

c        CALL HBNAME(777,'TRANCH',CIRC,'CIRC,'//
c     &                   'TREN,TRMA')
c        CALL HBNAME(777,'PTMIS', PTMT,'PTMT(2),'//
c     &                   'PTM(2),PZNU(2)')
c        CALL HBNAME(777,'LEPTONS',NLE_CWN,'NLE_CWN[0,100],'//
c     &                   'PL_CWN(4,NLE_CWN),KL_CWN(NLE_CWN)[-2,2],'//
c     &                   'MRL_CWN(NLE_CWN),LISOL_CWN(NLE_CWN)[0,4]')
c        CALL HBNAME(777,'GAMMAS',NGAM_CWN,'NGAM_CWN[0,100],'//
c     &                   'PG_CWN(4,NGAM_CWN),'//
c     &                   'LISGAM_CWN(NGAM_CWN)[0,4]')
c        CALL HBNAME(777,'JETS',NJG_CWN,'NJG_CWN[0,50],'//
c     &                   'PJG_CWN(4,NJG_CWN),'//
c     &                   'MRKJET_NJG(NJG_CWN)[0,2000],'//
c     &                   'BTAG_NJG(NJG_CWN)')

      return
      end
      
      
      
      SUBROUTINE NTEND

*KEND.
* close the ntuple:
    
*KEEP,CONST.
      INTEGER idbg,igener,irad,ifrad
      COMMON /CONST1/ idbg,igener,irad,ifrad



       call hrout(777,icycle,' ')


      END

