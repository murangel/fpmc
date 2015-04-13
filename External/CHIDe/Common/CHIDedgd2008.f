       double precision function CHIDedgd2008(xeff,k12,k22,delta2)

	implicit double precision(a-h,k-m,o-z)
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0
	data mus2/0.5625d0/


        pi=atan(1d0)*4d0

c 	  ks = (k12+k22)/2d0
 	  ks = DSQRT(k12*k22)
 
           ks0=1.4e0
           kss_soft = 1.0E0
           mus1 = 0.01E0
           asoft = 2.6E0
           ksscon1 = 0.05E0
           ksscon2 = 0.31E0
           kss_hard = ksscon1 + ksscon2*log(1./xeff)
           if(ks.le.ks0) then
              call GB_08(ks0,ks0,FB,mus2)
              call GGtype3_08(xeff,ks0,FG)
              CB=FG/FB
              call GB_08(k12,k22,FB,mus2)
              W=FB*CB
           else
              call GGtype3_08(xeff,ks,W)
           endif
           call GB_08(k12,k22,U,mus1)
           U = U*(1-xeff)**5

      
        W=W*ffactor_08(delta2,delta2)
        U=U*ffactor_08(delta2,delta2)
        
        npower=4
	  dgd_soft=U*asoft*kss_soft**npower/(kss_soft**npower+ks**npower)
        dgd_hard=W*ks**npower/(ks**npower+kss_hard**npower)
        
        CHIDedgd2008 = dgd_soft  + dgd_hard

      end function 

! ===== Dipole-like formfactor ======= ! 
       double precision function ffactor_08(k12,k22)
       implicit double precision(a-h,k-z)

C     proton scale, related to the proton radius
       mps=1.0 ! GeV^2
       ffactor_08 = 1d0/(1d0+k12/mps)/(1d0+k22/mps)
       
       end function


! ===== Born unintegrated DGD ====== !
       subroutine GB_08(k12,k22,FB,muss)
	implicit double precision(a-h,k-z)
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0
C        common/alphaborndgd/alphas

        alphas=0.88d0

	  FB=4d0*alphas/pi * k12/(k12+muss) * k22/(k22+muss) 
     .      *(1d0-ffactor_08(3*k12,3*k22))
c     .      *(1d0-ffactor_08(k12,k22))
	return
      end subroutine

      subroutine GGtype3_08(xeff,ks,FG)
       implicit double precision(a-h,k-m,o-z)

       FG=0.2450D0*(DLOG((ks+0.04d0)/0.04d0))**(0.34d0-6d0*DSQRT(xeff))
     .   /xeff**0.40d0

       end

