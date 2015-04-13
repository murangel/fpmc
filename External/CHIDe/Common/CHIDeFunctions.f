C   ----------------------------------------------------------------------------   C  
C     OTHER FUNCTION for 4-VECTORS
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product
      function CHIDedot(p1,p2)
	implicit none
      double precision p1(2),p2(2),CHIDedot
      CHIDedot=p1(1)*p2(1)+p1(2)*p2(2)
      end
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product of two sum of two vectors
      function CHIDedotsum(p1,p2,p3,p4)
	implicit none
      double precision p1(2),p2(2),p3(2),p4(2),CHIDedotsum
      CHIDedotsum=(p1(1)+p2(1))*(p3(1)+p4(1))
     &           +(p1(2)+p2(2))*(p3(2)+p4(2))
      end
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product of two difference of two vectors
      function CHIDedotdiff(p1,p2,p3,p4)
	implicit none
      double precision p1(2),p2(2),p3(2),p4(2),CHIDedotdiff
      CHIDedotdiff=(p1(1)-p2(1))*(p3(1)-p4(1))
     &            +(p1(2)-p2(2))*(p3(2)-p4(2))
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the scalar product of two difference in unintegrated 
C     distribution from Igor
      function CHIDedothalf(p1,p2,p3,p4)
        implicit none
      double precision p1(2),p2(2),p3(2),p4(2),CHIDedothalf
      CHIDedothalf=(p1(1)+0.5d0*p2(1))*(p3(1)+0.5d0*p4(1))
     &+(p1(2)+0.5d0*p2(2))*(p3(2)+0.5d0*p4(2))
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the ~cos of the angle between two 4-vectors
C     [cf Dijet2,b3 and g3]
      function fc(p1,p2)
        implicit none
      double precision p1(2),p2(2),fc
      double precision CHIDedot
      fc=CHIDedot(p1,p2)
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the ~sin of the angle between two 4-vectors
C     [cf Dijet2,b3 and g3]
      function fs(p1,p2,p3,p4)
        implicit none
      double precision p1(2),p2(2),p3(2),p4(2),fs
      double precision CHIDedot
      fs=CHIDedot(p1,p3)*CHIDedot(p2,p4)
     &   -(CHIDedot(p1,p4)*CHIDedot(p2,p3))
      end
C   ============================================================================   C

