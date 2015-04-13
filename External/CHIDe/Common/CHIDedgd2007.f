       double precision function CHIDedgd07(xeff,ks,delta2,ksdel,iglu)

	implicit double precision(a-h,k-m,o-z)
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0
	data mus2/0.5625d0/


        pi=atan(1d0)*4d0
        IINIP = 0

           ks0=1.4e0
           kss_soft = 1.0E0
           mus1 = 0.01E0
           asoft = 2.6E0
           ksscon1 = 0.0E0
           ksscon2 = 0.31E0
           kss_hard = ksscon1 + ksscon2*log(1./xeff)
           if(ks.le.ks0) then
              call GB(ks0,FB,mus2)
              call GGtype3(xeff,ks0,FG)
              CB=FG/FB
              call GB(ks,FB,mus2)
              W=FB*CB
           else
              call GGtype3(xeff,ks,W)
           endif
           call GB(ks,U,mus1)
           U = U*(1.-xeff)**5

      
        W=W*ffactor07(delta2)
        U=U*ffactor07(delta2)
        
        dgd_soft=U*asoft*kss_soft**4/(kss_soft**4+ks**4)
        dgd_hard=W*ks**4/(ks**4+kss_hard**4)
        
        CHIDedgd07 = dgd_soft  + dgd_hard
C	dgd07=dgd07*fJR
C        print*,'fJr=',fJR

      end function 

! ===== Dipole-like formfactor ======= ! 
       double precision function ffactor07(k2)
       implicit double precision(a-h,k-z)

C     proton scale, related to the proton radius
       mps=1.0 ! GeV^2
       ffactor07 = 1d0/(1d0+k2/mps)**2
       
       end function

! ===== Born unintegrated DGD ====== !
       subroutine GB(ks,FB,muss)
	implicit double precision(a-h,k-z)
        common/del/del2
	common/cons/pi,LQCDs,mus,mps,kstep,kss,ks0

        FB=4d0*0.7d0/pi*ks**2/((ks+muss)**2+muss*del2/2.0)
     .      *(1d0-ffactor07(3*ks))
	return
      end subroutine

      subroutine GGtype3(xeff,ks,FG)
       implicit double precision(a-h,k-m,o-z)

       FG=0.2450D0*(DLOG((ks+0.04d0)/0.04d0))**(0.34d0-6d0*DSQRT(xeff))
     .   /xeff**0.40d0

       end


