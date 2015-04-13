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
      FUNCTION SUDAKOVSP(MJJ,ET)
      REAL*4 MJJ,ET,C1,C2,C3,C4
      REAL*4 T1,T2,T3,T4,T5
      PARAMETER(C1 = 0.22d0)
      PARAMETER(C2 = 1.657d0)
      PARAMETER(C3 = 24.062d0)
      PARAMETER(C4 = 6.658d-4)
      
      T1 = MJJ**2/0.16d0
      T2 = 4d0/3d0*DLOG(DLOG(ET**2/0.04d0))
      T3 = C4*(DLOG(MJJ/2d0)-C3)**2
      T4 = (DLOG(ET**2/0.04d0))**2
      T5 = (ET**2)**(4d0/3d0)

      SUDAKOVSP = C1 * T1**(C2-T2-T3) * T4 * T5

      RETURN
      END
