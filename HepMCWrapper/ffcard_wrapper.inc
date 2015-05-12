C-----------------------------------------------------------------------
C * 13/03/2008 Oldrich Kepka 
C * Common blocks and declarations needed for the use of Format 
C * Free data cards. They are used in the module.f, module_cms.f, etc.
C----------------------------------------------------------------------- 
      integer 
     $ UOUTPUT,
     $ UMAXEV, UIPROC, UNFLUX, UNRN1, UNRN2, UIFIT, UISOFTM,
     $ UAAANOM, UZION, UAION, UAAEXOTIC,
     $ UCHIDeIGLU, UKMR2DELTA
     
      real URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
     $     UYJMIN,UYJMAX,UPTMIN,UPTMAX,UEMMIN, UEMMAX, UDKAPPA,
     $   UDLAMBDA,
     $   UANOMCUTOFF, UQ2WWMN, UQ2WWMX,UYWWMIN, UYWWMAX, UA0W, UACW,
     $   UA0Z, UACZ, UA1A, UA2A,
     $   UAAM, UAAQ, UAAN,
     $   UCHIDeX, UCHIDeXP, UCHIDeS2, ! R.S.
     $   UXI1Min, UXI1Max, UXI2Min, UXI2Max,
     $   UCHIDeGapMin, UCHIDeGapMax,
     $   UKMR2Q2CUT, UKMR2SURV, UKMR2SCALE,UBMIN
     
      COMMON /MYFFREAD1/ URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
     $     UYJMIN, UYJMAX,UPTMIN, UPTMAX, UEMMIN, UEMMAX, UDKAPPA,
     $     UACW, UA0W, UA0Z, UACZ, UA1A, UA2A,
     $     UAAM, UAAQ, UAAN,
     $     UCHIDeX, UCHIDeXP, UCHIDeS2, ! R.S.
     $   UXI1Min, UXI1Max, UXI2Min, UXI2Max,
     $   UCHIDeGapMin, UCHIDeGapMax,
     $    UKMR2Q2CUT, UKMR2SURV, UKMR2SCALE
 

     
      COMMON /MYFFREAD2/ UDLAMBDA, UANOMCUTOFF, UYWWMIN, UYWWMAX,
     $ UQ2WWMN, UQ2WWMX

      COMMON /MYFFREAD3/ UOUTPUT,UMAXEV,UIPROC,UNFLUX,UNRN1,UNRN2,UIFIT, 
     $     UISOFTM, UZION, UAION, UBMIN, UAAANOM, UAAEXOTIC,
     $     UCHIDeIGLU, UKMR2DELTA ! R.S.

      character*1 UHADR
      COMMON /CC0/ UHADR

      character*3 UTYPEPR, UTYPINT
      COMMON /CC1/ UTYPEPR
      COMMON /CC2/ UTYPINT

      character*4 UPART1, UPART2
      COMMON /CC3/ UPART1
      COMMON /CC4/ UPART2

      INTEGER UMODPDF1, UMODPDF2 
      COMMON /CC5/ UMODPDF1, UMODPDF2

      character*32 UNTNAME, UCHIDePATH
      COMMON /CYFFLONG1/ UNTNAME, UCHIDePATH




