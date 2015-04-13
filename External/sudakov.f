CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C This s numerical approximation to the sudakov form factor
C as obtained by A.Bzdak
C Returns the probability not to radiate vs. the hard 
C mass MJJ and transverse energy ET
C
C MB, 19/03/05
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION SUDAKOV_NEW(MJJ)
      DOUBLE PRECISION MJJ
      
      SUDAKOV_NEW = 
     & (1.7689             * UNITSTEP(MJJ,  7d0,  15d0)) / MJJ**2 + 
     & (5.051642565797349  * UNITSTEP(MJJ, 15d0,  30d0)) / MJJ**2.4 + 
     & (13.917570348992509 * UNITSTEP(MJJ, 30d0,  60d0)) / MJJ**2.7 + 
     & (45.576             * UNITSTEP(MJJ, 60d0, 120d0)) / MJJ**3 + 
     & (230.7794760047412  * UNITSTEP(MJJ,120d0, 240d0)) / MJJ**3.35 + 
     & (900.9984435801367  * UNITSTEP(MJJ,240d0, 480d0)) / MJJ**3.6 + 
     & (3683.4738353490006 * UNITSTEP(MJJ,480d0, 960d0)) / MJJ**3.83 + 
     & (14388.742638648244 * UNITSTEP(MJJ,960d0,1909d0)) / MJJ**4.03      

      RETURN
      END

C------------------------------------------------------------

      DOUBLE PRECISION FUNCTION UNITSTEP(MJJ,MMIN,MMAX)
      DOUBLE PRECISION MJJ,MMIN,MMAX
      
      UNITSTEP = 0d0
      IF(MJJ.GE.MMIN.AND.MJJ.LE.MMAX) UNITSTEP = 1d0

      RETURN
      END

C------------------------------------------------------------

      DOUBLE PRECISION FUNCTION SUDAKOV(MJJ,ET)
      DOUBLE PRECISION MJJ,ET,C1,C2,C3,C4
      PARAMETER(C1 = 0.22d0)
      PARAMETER(C2 = 1.657d0)
      PARAMETER(C3 = 24.062d0)
      PARAMETER(C4 = 6.658d-4)
      
      T1 = MJJ**2/0.16d0
      T2 = 4d0/3d0*DLOG(DLOG(ET**2/0.04d0))
      T3 = C4*(DLOG(MJJ/2d0)-C3)**2
      T4 = (DLOG(ET**2/0.04d0))**2
      T5 = (ET**2)**(4d0/3d0)

      SUDAKOV = C1 * T1**(C2-T2-T3) * T4 * T5

      RETURN
      END
