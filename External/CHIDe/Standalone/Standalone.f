       program StandaloneCHIDeHiggs
       implicit none
       
       character*50  dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
       common/CHIDePATH/ dgdtab1, dgdtab2, dgdtab3, dgdtab4, sudatab
      
       dgdtab4="Data/Higgsdgdtab1.d"
       dgdtab4="Data/Higgsdgdtab2.d"
       dgdtab4="Data/Higgsdgdtab3.d"
       dgdtab4="Data/Higgsdgdtab4.d"
       sudatab="Data/Higgssudatab.d"
       
       CALL CHIDeHiggsInit
     & (120.0d0,173.3d0,14000.0d0**2,4,1d0,1d0,0.075d0)
       call integrateHiggs(10000) 


       dgdtab4="Data/ggdgdtab1.d"
       dgdtab4="Data/ggdgdtab2.d"
       dgdtab4="Data/ggdgdtab3.d"
       dgdtab4="Data/ggdgdtab4.d"
       sudatab="Data/ggsudatab.d"
       
       CALL CHIDeGGInit(10.0d0,1960.0d0**2,0.15d0)
       call integrateGG(10000) 


       end

      subroutine integrateGG(N)
      a = 2
      b = 3
      end

       
      subroutine integrateHiggs(N)
      implicit none
      integer N,NN,i,ii
      double precision sigma
      double precision mean
      double precision sum,max
      double precision x(11)
      double precision jac,kmax,pi
      double precision thetak, thetap, theta1, theta3
      double precision ak, akp, ak1, ak3
      double precision a3,b1
      double precision k(2), kp(2), k1(2), k3(2)
  
      double precision dsigma
      external dsigma

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
        
         kmax = 10.d0
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
 
         ak1=kmax*x(5)
         theta1=2.*pi*x(6)
         jac=jac*2.*pi*ak1*kmax
         k1(1)=ak1*cos(theta1)
         k1(2)=ak1*sin(theta1)

         ak3=kmax*x(7)
         theta3=2.*pi*x(8)
         jac=jac*2.*pi*ak3*kmax
         k3(1)=ak3*cos(theta3)
         k3(2)=ak3*sin(theta3)
      
         a3=0.02*x(9)    
         b1=0.02*x(10)

         jac=jac*0.02**2
         
         sigma=0d0

         call CHIDeHiggs(sigma,a3,b1,k1,k3,k,kp)
         sum = sum + sigma*jac
         mean = sum/i
      if(MOD(i,NN).eq.0) print*, i," /",N,"\t cs : ",sum/i
      enddo
      sum = sum/N
      print*, "Cross section = ", sum

      end
