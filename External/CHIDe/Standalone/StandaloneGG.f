       program StandaloneCHIDeHiggs
       implicit none
       
       character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       dgdtab4="Data/ggdgdtab1.d"
       dgdtab4="Data/ggdgdtab2.d"
       dgdtab4="Data/ggdgdtab3.d"
       dgdtab4="Data/ggdgdtab4.d"
       sudatab="Data/ggsudatab.d"
       
       CALL CHIDeGGInit(14000.0d0**2,10d0,0.15d0,2d0,1d0,4,
     &                  0.03d0,0.08d0,0d0,1d0,-2.5d0,2.5d0)
       call integrate(10000000) 
       end

      subroutine integrate(N)
      implicit none
      integer N,NN,i,ii
      double precision sigma
      double precision mean
      double precision sum,max
      double precision x(12)
      double precision jac,kmax,pi,kmax1,kmax2
      double precision kmin,kmin1,kmin2
      double precision bmin,bmax,b,c
      double precision thetak, thetap, theta1, theta2, theta3
      double precision ak, akp, ak1, ak2, ak3
      double precision b1,b2,a1,a2
      double precision k(2), kp(2), k1(2), k2(2), k3(2)
  
      double precision  CHIDedotdiff
      external  CHIDedotdiff

      sum = 0
      mean = 0
      ii=0
      NN=N/100
      print*, "integrate, N=",N


      do i = 1, N

         x(1)=rand()
         x(2)=rand()
         x(3)=rand()
         x(4)=rand()
         x(5)=rand()
         x(6)=rand()
         x(7)=rand()
         x(8)=rand()
         x(9)=rand()
         x(10)=rand()
         x(11)=rand()
         x(12)=rand()
       

         kmin = 0d0
         kmin1 = 0.010d0
         kmin2 = 0d0
         kmax = 10d0
         kmax1 = 2d0
         kmax2 = 14d0
         bmin = 0.00001d0
         bmax = 0.1d0
         pi=3.14159d0
         jac = 1d0
         ak=kmax*x(1)
         thetak=2.*pi*x(2)
         jac=jac*2.*pi*ak*kmax
         k(1)=ak*cos(thetak)
         k(2)=ak*sin(thetak)
      
         akp=kmax*x(3)
         thetap=2.*pi*x(4)
         jac=jac*2.*pi*akp*kmax
         kp(1)=akp*cos(thetap)
         kp(2)=akp*sin(thetap)

         ak2=8d0+kmax2*x(7)
         theta1=2.*pi*x(8)
         jac=jac*2.*pi*ak2*kmax2
         k2(1)=ak2*cos(theta2)
         k2(2)=ak2*sin(theta2)


         b = 4.d0
         C = B/(DEXP(-B*kMIN1)-DEXP(-B*kMAX1))
         ak1 = -1.0/B*DLOG(DEXP(-B*kMIN1)
     &            - B/C*x(5))
         ak3 = -1.0/B*DLOG(DEXP(-B*kMIN1)
     &            - B/C*x(9))
         jac = pi*jac/C/DEXP(-B*ak1)
         jac = pi*jac/C/DEXP(-B*ak3)
         ak1 = sqrt(ak1)
         ak3 = sqrt(ak3)

c        ak1=kmax1*x(5)
         theta1=2.*pi*x(6)
c        jac=jac**pi*ak1*kmax1
         k1(1)=ak1*cos(theta1)
         k1(2)=ak1*sin(theta1)


c        ak3=kmax1*x(9)
         theta3=2.*pi*x(10)
c        jac=jac*2.*pi*ak3*kmax1
         k3(1)=ak3*cos(theta3)
         k3(2)=ak3*sin(theta3)
      
         b1=(bmin/bmax)**x(11)*bmax
         b2=(bmin/bmax)**x(12)*bmax
         jac = jac*b1*b2*DLOG(bmax/bmin)**2
      
         a1=CHIDedotdiff(k1,k2,k1,k2)/1960d0**2/b1
         a2=CHIDedotdiff(k3,k2,k3,k2)/1960d0**2/b2

         call CHIDeGG(sigma,k,kp,k1,k2,k3,b1,b2,a1,a2)
         sum = sum + sigma*jac
         mean = sum/i
      if(MOD(i,NN).eq.0) print*, i," /",N,"\t cs : ",sum/i
      enddo
      sum = sum/N
      print*, "Cross section = ", sum

      end
