       program StandaloneCHIDeDiphoton
       implicit none
      
       double precision ss,ptmin
       common/aaa/ss,ptmin
       
       character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       dgdtab4="Data/dgdtab1.d"
       dgdtab4="Data/dgdtab2.d"
       dgdtab4="Data/dgdtab3.d"
       dgdtab4="Data/dgdtab4.d"
       sudatab="Data/ggsudatab.d"

       ss=14000.0**2
       ptmin=5.0
       
       CALL CHIDeDiphotonInit(ss, ptmin, 0.03d0, 0.5d0, 1d0,4,
     &                  0.002d0, 0.02d0, 002d0, 0.02d0, -2d0, 2d0)

       call integrate(100000000) 

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
      
      double precision ss,ptmin
      common/aaa/ss,ptmin

      double precision  CHIDedotdiff, CHIDedotsum
      external  CHIDedotdiff, CHIDedotsum

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
       

         kmax = 10d0
         pi=3.14159d0
         jac = 1d0
         bmin = 0.0001d0
         bmax = 1d0
         
         ak=0.0000000001+kmax*x(1)
         thetak=2.*pi*x(2)
         jac=jac*2.*pi*ak*kmax
         k(1)=ak*cos(thetak)
         k(2)=ak*sin(thetak)
      
         akp=0.000000001 + kmax*x(3)
         thetap=2.*pi*x(4)
         jac=jac*2.*pi*akp*kmax
         kp(1)=akp*cos(thetap)
         kp(2)=akp*sin(thetap)

 
         ak1=kmax*x(5)
         theta1=2.*pi*x(6)
         jac=jac*2.*pi*ak1*kmax
         k1(1)=ak1*cos(theta1)
         k1(2)=ak1*sin(theta1)

         ak2=5.+10.*x(7)*kmax
         theta2=2.*pi*0d0
         jac=jac*2.*pi*ak2*kmax*10.
         k2(1)=ak2*cos(theta2)
         k2(2)=ak2*sin(theta2)

         ak3=kmax*x(8)
         theta3=2.*pi*x(9)
         jac=jac*2.*pi*ak3*kmax
         k3(1)=ak3*cos(theta3)
         k3(2)=ak3*sin(theta3)

c         C = B/(DEXP(-B*kMIN1)-DEXP(-B*kMAX1))
c         ak1 = -1.0/B*DLOG(DEXP(-B*kMIN1)
c     &            - B/C*x(5))
c         ak3 = -1.0/B*DLOG(DEXP(-B*kMIN1)
c     &            - B/C*x(9))
c         jac = pi*jac/C/DEXP(-B*ak1)
c         jac = pi*jac/C/DEXP(-B*ak3)
c         ak1 = sqrt(ak1)
c         ak3 = sqrt(ak3)

c         theta1=2.*pi*x(6)
c         k1(1)=ak1*cos(theta1)
c         k1(2)=ak1*sin(theta1)
c
c
c         theta3=2.*pi*x(10)
c        k3(1)=ak3*cos(theta3)
c         k3(2)=ak3*sin(theta3)
      
         b1=(bmin/bmax)**x(11)*bmax
         b2=(bmin/bmax)**x(12)*bmax
 
         jac = jac*DLOG(bMAX/(bMIN))*B1
         jac = jac*DLOG(bMAX/(bMIN))*B2

         a1=CHIDedotsum(k1,k2,k1,k2)/b1/ss
         a2=CHIDedotsum(k2,k3,k2,k3)/b2/ss

         call CHIDeDiphoton(sigma,k,kp,k1,k2,k3,b1,b2,a1,a2)
        
         if(sigma-sigma.NE.0d0) then
          print*, "x = ", x
          print*, "sigma = ", sigma
          print*, "k=", k
          print*, "kp  = ", kp
          print*, "k1  = ", k1
          print*, "k2  = ", k2
          print*, "k3  = ", k3
          print*, "b1  = ", b1
          print*, "b2  = ", b2
          print*, "a1  = ", a1
          print*, "a2  = ", a1
         endif
c        print*, sigma, jac
         sum = sum + sigma*jac
         mean = sum/i
      if(MOD(i,NN).eq.0) print*, i," /",N,"   cs : ",sum/i
      enddo
      sum = sum/N
      print*, "Cross section = ", sum

      end
