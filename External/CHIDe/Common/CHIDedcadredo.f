C-----------------------------------------------------------------------
C                                                                       
      real*8 FUNCTION DCADREDO(F,A,B,AERR,RERR,ERROR,IER)                   
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            IER                                            
      real*8               F,A,B,AERR,RERR,ERROR                          
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            IBEGS(30),MAXTS,MAXTBL,MXSTGE,IBEG,II,NNLEFT   
      INTEGER            I,N2,III,ISTEP2,IEND,ISTEP,L,LM1,IT,ISTAGE,N   
      real*8               T(10,10),R(10),AIT(10),DIF(10),RN(4),TS(2049)  
      real*8               BEGIN(30),FINIS(30),EST(30)                    
      real*8               H2TOL,AITTOL,LENGTH,JUMPTL,ZERO,P1,HALF,ONE    
      real*8               TWO,FOUR,FOURP5,TEN,HUN,CADRE,AITLOW           
      real*8               STEPMN,STEPNM,STAGE,CUREST,FNSIZE,HRERR        
      real*8               PREVER,BEG,FBEG,EDN,FEND,STEP,ASTEP,TABS,HOVN  
      real*8               FN,SUM,SUMABS,VINT,TABTLM,ERGL,ERGOAL          
      real*8               ERRA,ERRR,FEXTRP,ERRER,DIFF,SING,FEXTM1        
      real*8               H2NEXT,SINGNX,SLOPE,FBEG2,ERRET,H2TFEX,FI      
      LOGICAL            H2CONV,AITKEN,RIGHT,REGLAR,REGLSV(30)          
      DATA               AITLOW,H2TOL,AITTOL,JUMPTL,MAXTS,MAXTBL,MXSTGE 
     1                   /1.1,.15,.1,.01,2049,10,30/                    
      DATA               RN(1),RN(2),RN(3),RN(4)/                       
     1                   .7142005,.3466282,.843751,.1263305/            
      DATA               ZERO,P1,HALF,ONE,TWO,FOUR,FOURP5,TEN,HUN       
     1                   /0.0,0.1,0.5,1.0,2.0,4.0,                      
     2                   4.5,10.0,100.0/                                
C     FIRST EXECUTABLE STATEMENT 

      IER = 0                                                           
      CADRE = ZERO                                                      
      ERROR = ZERO                                                      
      CUREST = ZERO                                                     
      VINT = ZERO                                                       
      LENGTH = ABS(B-A)                                                 
      IF (LENGTH .EQ. ZERO) GO TO 215                                   
      IF (RERR .GT. P1 .OR. RERR .LT. ZERO) GO TO 210                   
      HRERR = RERR+HUN                                                  
      IF (AERR .EQ. ZERO .AND. HRERR .LE. HUN) GO TO 210                
      ERRR = RERR                                                       
      ERRA = ABS(AERR)                                                  
      STEPMN = LENGTH/(TWO**MXSTGE)                                     
      STEPNM = DMAX1(LENGTH,ABS(A),ABS(B))*TEN                          
      STAGE = HALF                                                      
      ISTAGE = 1                                                        
      FNSIZE = ZERO                                                     
      PREVER = ZERO                                                     
      REGLAR = .FALSE.                                                  
C                                  THE GIVEN INTERVAL OF INTEGRATION    
C                                    IS THE FIRST INTERVAL CONSIDERED.  
      BEG = A                                                           
      FBEG = F(BEG)*HALF                                                
      TS(1) = FBEG                                                      
      IBEG = 1                                                          
      EDN = B                                                           
      FEND = F(EDN)*HALF                                                
      TS(2) = FEND                                                      
      IEND = 2                                                          
    5 RIGHT = .FALSE.                                                   
C                                  INVESTIGATION OF A PARTICULAR        
C                                    SUBINTERVAL BEGINS AT THIS POINT. 
   10 STEP = EDN - BEG                                                  
      ASTEP =  ABS(STEP)                                                
      IF (ASTEP .LT. STEPMN) GO TO 205                                  
      HRERR = STEPNM+ASTEP                                              
      IF (HRERR .EQ. STEPNM) GO TO 205                                  
      T(1,1) = FBEG + FEND                                              
      TABS = ABS(FBEG) + ABS(FEND)                                      
      L = 1                                                             
      N = 1                                                             
      H2CONV = .FALSE.                                                  
      AITKEN = .FALSE.                                                  
   15 LM1 = L                                                           
      L = L + 1   
C                                  CALCULATE THE NEXT TRAPEZOID SUM,    
C                                    T(L,1), WHICH IS BASED ON *N2* + 1 
C                                    EQUISPACED POINTS. HERE,           
C                                    N2 = N*2 = 2**(L-1).               
      N2 = N+N                                                          
      FN = N2                                                           
      ISTEP = (IEND - IBEG)/N                                           
      IF (ISTEP .GT. 1) GO TO 25                                        
      II = IEND                                                         
      IEND = IEND + N                                                   
      IF (IEND .GT. MAXTS) GO TO 200                                    
      HOVN = STEP/FN                                                    
      III = IEND                                                        
      FI = ONE                                                          
      DO 20 I=1,N2,2                                                    
         TS(III) = TS(II)                                               
         TS(III-1) = F(EDN - FI * HOVN)                                 
         FI = FI+TWO                                                    
         III = III-2                                                    
         II = II-1                                                      
   20 CONTINUE                                                          
      ISTEP = 2                                                         
   25 ISTEP2 = IBEG + ISTEP/2                                           
      SUM = ZERO                                                        
      SUMABS = ZERO                                                     
      DO 30 I=ISTEP2,IEND,ISTEP                                         
         SUM = SUM + TS(I)                                              
         SUMABS = SUMABS + ABS(TS(I))                                   
   30 CONTINUE                                                          
      T(L,1) = T(L-1,1)*HALF+SUM/FN                                     
      TABS = TABS*HALF+SUMABS/FN                                        
      N = N2                                                            
C                                  GET PRELIMINARY VALUE FOR *VINT*     
C                                    FROM LAST TRAPEZOID SUM AND UPDATE 
C                                    THE ERROR REQUIREMENT *ERGOAL*     
C                                    FOR THIS SUBINTERVAL. 

          
      IT = 1                                                            
      VINT = STEP*T(L,1)                                                
      TABTLM = TABS*TEN                                                 
      FNSIZE = DMAX1(FNSIZE,ABS(T(L,1)))                                
      ERGL = ASTEP*FNSIZE*TEN                                           
      ERGOAL = STAGE*DMAX1(ERRA,ERRR*ABS(CUREST+VINT))                  
C                                  COMPLETE ROW L AND COLUMN L OF *T*   
C                                    ARRAY.                             
      FEXTRP = ONE                                                      
      DO 35 I=1,LM1                                                     
         FEXTRP = FEXTRP*FOUR                                           
         T(I,L) = T(L,I) - T(L-1,I)                                     
         T(L,I+1) = T(L,I) + T(I,L)/(FEXTRP-ONE)                        
   35 CONTINUE                                                          
      ERRER = ASTEP*ABS(T(1,L))                                         
C                                  PRELIMINARY DECISION PROCEDURE       
C                                    IF L = 2 AND T(2,1) = T(1,1),      
C                                    GO TO 135 TO FOLLOW UP THE         
C                                    IMPRESSION THAT INTERGRAND IS      
C                                    STRAIGHT LINE.                     
      IF (L .GT. 2) GO TO 40                                            
      HRERR = TABS+P1*ABS(T(1,2))                                       
      IF (HRERR .EQ. TABS) GO TO 135                                    
      GO TO 15                                                          
C                                  CACULATE NEXT RATIOS FOR             
C                                    COLUMNS 1,...,L-2 OF T-TABLE       
C                                    RATIO IS SET TO ZERO IF DIFFERENCE 
C                                    IN LAST TWO ENTRIES OF COLUMN IS   
C                                    ABOUT ZERO                         
   40 DO 45 I=2,LM1                                                     
         DIFF = ZERO                                                    
         HRERR = TABTLM+ABS(T(I-1,L))                                   
         IF (HRERR .NE. TABTLM) DIFF = T(I-1,LM1)/T(I-1,L)              
         T(I-1,LM1) = DIFF                                              
   45 CONTINUE                                                          
      IF (ABS(FOUR-T(1,LM1)) .LE. H2TOL) GO TO 60                       
      IF (T(1,LM1) .EQ. ZERO) GO TO 55                                  
      IF (ABS(TWO-ABS(T(1,LM1))) .LT. JUMPTL) GO TO 130                 
      IF (L .EQ. 3) GO TO 15                                            
      H2CONV = .FALSE.                                                  
      IF (ABS((T(1,LM1)-T(1,L-2))/T(1,LM1)) .LE. AITTOL) GO TO 75       
   50 IF (REGLAR) GO TO 55                                              
      IF (L .EQ. 4) GO TO 15                                            
      HRERR = ERGL+ERRER                                                
   55 IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175            
      GO TO 145                                                         
C                                  CAUTIOUS ROMBERG EXTRAPOLATION       
   60 IF (H2CONV) GO TO 65                                              
      AITKEN = .FALSE.                                                  
      H2CONV = .TRUE.                                                   
   65 FEXTRP = FOUR                                                     
   70 IT = IT + 1                                                       
      VINT = STEP*T(L,IT)                                               
      ERRER = ABS(STEP/(FEXTRP-ONE)*T(IT-1,L))                          
      IF (ERRER .LE. ERGOAL) GO TO 160                                  
      HRERR = ERGL+ERRER                                                
      IF (HRERR .EQ. ERGL) GO TO 160                                    
      IF (IT .EQ. LM1) GO TO 125                                        
      IF (T(IT,LM1) .EQ. ZERO) GO TO 70                                 
      IF (T(IT,LM1) .LE. FEXTRP) GO TO 125                              
      IF (ABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP .LT. AITTOL)                
     1       FEXTRP = FEXTRP*FOUR                                       
      GO TO 70                                                          
C                                  INTEGRAND MAY HAVE X**ALPHA TYPE     
C                                    SINGULARITY                        
C                                    RESULTING IN A RATIO OF *SING*  =  
C                                    2**(ALPHA + 1)                     
   75 IF (T(1,LM1) .LT. AITLOW) GO TO 175                               
      IF (AITKEN) GO TO 80                                              
      H2CONV = .FALSE.                                                  
      AITKEN = .TRUE.                                                   
   80 FEXTRP = T(L-2,LM1)                                               
      IF (FEXTRP .GT. FOURP5) GO TO 65                                  
      IF (FEXTRP .LT. AITLOW) GO TO 175                                 
      IF (ABS(FEXTRP-T(L-3,LM1))/T(1,LM1) .GT. H2TOL) GO TO 175         
      SING = FEXTRP                                                     
      FEXTM1 = ONE/(FEXTRP - ONE)                                       
      AIT(1) = ZERO                                                     
      DO 85 I=2,L                                                       
         AIT(I) = T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1                     
         R(I) = T(1,I-1)                                                
         DIF(I) = AIT(I) - AIT(I-1)                                     
   85 CONTINUE                                                          
      IT = 2                                                            
   90 VINT = STEP*AIT(L)                                                
      ERRER = ERRER*FEXTM1                                              
      HRERR = ERGL+ERRER                                                
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 95             
      IER = MAX0(IER,65)                                                
      GO TO 160                                                         
   95 IT = IT + 1                                                       
      IF (IT .EQ. LM1) GO TO 125                                        
      IF (IT .GT. 3) GO TO 100                                          
      H2NEXT = FOUR                                                     
      SINGNX = SING+SING                                                
  100 IF (H2NEXT .LT. SINGNX) GO TO 105                                 
      FEXTRP = SINGNX                                                   
      SINGNX = SINGNX+SINGNX                                            
      GO TO 110                                                         
  105 FEXTRP = H2NEXT                                                   
      H2NEXT = FOUR*H2NEXT                                              
  110 DO 115 I=IT,LM1                                                   
         R(I+1) = ZERO                                                  
         HRERR = TABTLM+ABS(DIF(I+1))                                   
         IF (HRERR .NE. TABTLM) R(I+1) = DIF(I)/DIF(I+1)                
  115 CONTINUE                                                          
      H2TFEX = -H2TOL*FEXTRP                                            
      IF (R(L) - FEXTRP .LT. H2TFEX) GO TO 125                          
      IF (R(L-1)-FEXTRP .LT. H2TFEX) GO TO 125                          
      ERRER = ASTEP*ABS(DIF(L))                                         
      FEXTM1 = ONE/(FEXTRP - ONE)                                       
      DO 120 I=IT,L                                                     
         AIT(I) = AIT(I) + DIF(I)*FEXTM1                                
         DIF(I) = AIT(I) - AIT(I-1)                                     
  120 CONTINUE                                                          
      GO TO 90                                                          
C                                  CURRENT TRAPEZOID SUM AND RESULTING  
C                                    EXTRAPOLATED VALUES DID NOT GIVE   
C                                    A SMALL ENOUGH *ERRER*.            
C                                    NOTE -- HAVING PREVER .LT. ERRER   
C                                    IS AN ALMOST CERTAIN SIGN OF       
C                                    BEGINNING TROUBLE WITH IN THE FUNC-
C                                    TION VALUES. HENCE, A WATCH FOR,   
C                                    AND CONTROL OF, NOISE SHOULD       
C                                    BEGIN HERE.                        
  125 FEXTRP = DMAX1(PREVER/ERRER,AITLOW)                               
      PREVER = ERRER                                                    
      IF (L .LT. 5) GO TO 15                                            
      IF (L-IT .GT. 2 .AND. ISTAGE .LT. MXSTGE) GO TO 170               
      ERRET = ERRER/(FEXTRP**(MAXTBL-L))                                
      HRERR = ERGL+ERRET                                                
      IF (ERRET .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170            
      GO TO 15                                                          
C                                  INTEGRAND HAS JUMP (SEE NOTES)       
  130 HRERR = ERGL+ERRER                                                
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170            
C                                  NOTE THAT  2*FN = 2**L               
      DIFF = ABS(T(1,L))*(FN+FN)                                        
      GO TO 160                                                         
C                                  INTEGRAND IS STRAIGHT LINE           
C                                    TEST THIS ASSUMPTION BY COMPARING  
C                                    THE VALUE OF THE INTEGRAND AT      
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH 
C                                    THE VALUE OF THE STRAIGHT LINE     
C                                    INTERPOLATING THE INTEGRAND AT THE 
C                                    TWO END POINTS OF THE SUB-INTERVAL.
C                                    IF TEST IS PASSED, ACCEPT *VINT*   
  135 SLOPE = (FEND-FBEG)*TWO                                           
      FBEG2 = FBEG+FBEG                                                 
      DO 140 I=1,4                                                      
         DIFF = ABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)              
         HRERR = TABTLM+DIFF                                            
         IF(HRERR .NE. TABTLM) GO TO 155                                
  140 CONTINUE                                                          
      GO TO 160                                                         
C                                  NOISE MAY BE DOMINANT FEATURE        
C                                    ESTIMATE NOISE LEVEL BY COMPARING  
C                                    THE VALUE OF THE INTEGRAND AT      
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH 
C                                    THE VALUE OF THE STRAIGHT LINE     
C                                    INTERPOLATING THE INTEGRAND AT THE 
C                                    TWO ENDPOINTS. IF SMALL ENOUGH,    
C                                    ACCEPT *VINT*                      
  145 SLOPE = (FEND-FBEG)*TWO                                           
      FBEG2 = FBEG+FBEG                                                 
      I = 1                                                             
  150 DIFF = ABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)                 
  155 ERRER = DMAX1(ERRER,ASTEP*DIFF)                                   
      HRERR = ERGL+ERRER                                                
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175            
      I = I+1                                                           
      IF (I .LE. 4) GO TO 150                                           
      IER = 66
                                                          
C                                  INTERGRATION OVER CURRENT SUB-       
C                                    INTERVAL SUCCESSFUL                
C                                    ADD *VINT* TO *DCADRE* AND *ERRER* 
C                                    TO *ERROR*, THEN SET UP NEXT SUB-  
C                                    INTERVAL, IF ANY.                  
  160 CADRE = CADRE + VINT                                              
      ERROR = ERROR + ERRER                                             
      IF (RIGHT) GO TO 165                                              
      ISTAGE = ISTAGE - 1                                               
      IF (ISTAGE .EQ. 0) GO TO 220                                      
      REGLAR = REGLSV(ISTAGE)                                           
      BEG = BEGIN(ISTAGE)                                               
      EDN = FINIS(ISTAGE)                                               
      CUREST = CUREST - EST(ISTAGE+1) + VINT                            
      IEND = IBEG - 1                                                   
      FEND = TS(IEND)                                                   
      IBEG = IBEGS(ISTAGE)                                              
      GO TO 180                                                         
  165 CUREST = CUREST + VINT                                            
      STAGE = STAGE+STAGE                                               
      IEND = IBEG                                                       
      IBEG = IBEGS(ISTAGE)                                              
      EDN = BEG                                                         
      BEG = BEGIN(ISTAGE)                                               
      FEND = FBEG                                                       
      FBEG = TS(IBEG)                                                   
      GO TO 5                                                           
C                                  INTEGRATION OVER CURRENT SUBINTERVAL 
C                                    IS UNSUCCESSFUL. MARK SUBINTERVAL  
C                                    FOR FURTHER SUBDIVISION. SET UP    
C                                    NEXT SUBINTERVAL.                  
  170 REGLAR = .TRUE.                                                   
  175 IF (ISTAGE .EQ. MXSTGE) GO TO 205                                 
      IF (RIGHT) GO TO 185                                              
      REGLSV(ISTAGE+1) = REGLAR                                         
      BEGIN(ISTAGE) = BEG                                               
      IBEGS(ISTAGE) = IBEG                                              
      STAGE = STAGE*HALF                                                
  180 RIGHT = .TRUE.                                                    
      BEG = (BEG+EDN)*HALF                       
      IBEG = (IBEG+IEND)/2                                              
      TS(IBEG) = TS(IBEG)*HALF  
      FBEG = TS(IBEG)                                                   
      GO TO 10                                                            
  185 NNLEFT = IBEG - IBEGS(ISTAGE)                                     
      IF (IEND+NNLEFT .GE. MAXTS) GO TO 200                             
      III = IBEGS(ISTAGE)                                               
      II = IEND                                                         
      DO 190 I=III,IBEG                                                 
         II = II + 1                                                    
         TS(II) = TS(I)                                                 
  190 CONTINUE                                                          
      DO 195 I=IBEG,II                                                  
         TS(III) = TS(I)                                                
         III = III + 1                                                  
  195 CONTINUE                                                          
      IEND = IEND + 1                                                   
      IBEG = IEND - NNLEFT                                              
      FEND = FBEG                                                       
      FBEG = TS(IBEG)                                                   
      FINIS(ISTAGE) = EDN                                               
      EDN = BEG                                                         
      BEG = BEGIN(ISTAGE)                                               
      BEGIN(ISTAGE) = EDN                                               
      REGLSV(ISTAGE) = REGLAR                                           
      ISTAGE = ISTAGE + 1                                               
      REGLAR = REGLSV(ISTAGE)                                           
      EST(ISTAGE) = VINT                                                
      CUREST = CUREST + EST(ISTAGE)                                     
      GO TO 5                                                           
C                                  FAILURE TO HANDLE GIVEN INTEGRA-     
C                                    TION PROBLEM                       
  200 IER = 131                                                         
      GO TO 215                                                         
  205 IER = 132                                                         
      GO TO 215                                                         
  210 IER = 133                                                         
  215 CADRE = CUREST + VINT                                             
  220 DCADREDO= CADRE
                                                     
 9000 CONTINUE                                                          
C      IF (IER .NE. 0) CALL UERTST (IER,6HDCADRE)                       
 9005 RETURN  

                                                     
      END                                                               
