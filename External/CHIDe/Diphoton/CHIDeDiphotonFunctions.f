
C   ----------------------------------------------------------------------------   C
C   SUBROUTINE ALPHA STRONG
C   Compute as in function of Q²:
C      -> choose a Q²
C      -> compute the corresponding lambda(nf)
C      -> compute as
C   ----------------------------------------------------------------------------   C
      function CHIDeDiphotonas(q) 
       implicit none 
      integer i,nf 
      double precision CHIDeDiphotonas
      double precision mq(1:6)
      double precision smu,smc,sms,smb,smt 
      double precision lambda(1:6),la
      double precision q
      double precision pi,s,ncolor,gg,gq,gelm,kmax,nb,mp
      double precision Qup,Qdown

      common/const/pi,ncolor,gg,gq,gelm,kmax,nb,mp
      common/flavor/mq,Qup,Qdown,nf
      common/param/s
            
C     Constituent quarks masses and l(5) en GeV from particle data group 
      smu=mq(2)**2
      sms=mq(3)**2
      smc=mq(4)**2
      smb=mq(5)**2
      smt=mq(6)**2
      lambda(5)=0.2d0
 
      do i=1,3 
      lambda(5-i)=lambda(6-i)*(mq(6-i)
     &/lambda(6-i))**(2./(33.-2.*(5-i)))
      end do
      lambda(6)=lambda(5)**(23./21.)/(mq(6)**(2./21.))
      continue

      if(q.GT.smt)then 
        nf=6
      elseif(q.GT.smb)then      
        nf=5
      elseif(q.GT.smc)then 
        nf=4
      elseif(q.GT.sms)then
        nf=3 
      elseif(q.GT.smu)then
        nf=2  
      else
      goto 666
      end if
      continue

      la=lambda(nf)
      CHIDeDiphotonas=12.*pi/(33.-2.*nf)/log(q/la**2)

      if(CHIDeDiphotonas.LT.0.or.CHIDeDiphotonas.GT.0.7d0)
     &   CHIDeDiphotonas=0.7d0

666   end    

