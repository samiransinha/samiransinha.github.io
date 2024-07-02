 



C        Test for spline basis

         subroutine splinebasis(d, n, m, m1, k, x, innerknots, 
     *   boundaryknots, basis)
C        This subroutine generates Bspline basis functions.      
C        x(n) is a n by 1 input vector for which B-spline
C        basis function will be evaluated. 
C        innerknots(m1) set of m1 innerknot points.
C        newknots is the entire set of knots, of length m=m1+2(d+1)
C        where d is the degree of the splines.
C        k=number of spline basis=m1+d+1  
         IMPLICIT NONE
         integer*4 d, k, m,  m1, n
         double precision x(n), innerknots(m1), boundaryknots(2)
         double precision newknots(m), basis(n, k), result
         external b
         integer*4 i1, i, j
         do i1=1, (d+1)
          newknots(i1)=boundaryknots(1)
         end do
         do i1=(d+2), (m1+d+1)
          newknots(i1)=innerknots(i1-d-1)
         end do
         do i1=(m1+d+2), m
          newknots(i1)=boundaryknots(2)
         end do
         do i=1, n
          if(x(i).eq.boundaryknots(2)) then 
           basis(i, k)=1.d0
           do j=1, (k-1)
            basis(i, j)=0.d0
           end do
          else 
           do j=1, k
            call b(m, j, (d+1), x(i), newknots, result, b)
            basis(i, j)=result  
           end do
          endif
         end do
         return 
         end
C ---------------- 
         subroutine b(i1, i2, i3, y, newknots, result, dumsub)
C        This subroutine calculates the i2 th basis of the spline of degree (i3-1)
         IMPLICIT NONE
         integer*4 i1, i2, i3
         double precision y, newknots(i1), temp1, temp2, result, 
     *   result1, result2
         external dumsub
         if(i3.eq.1) then 
          if((y.ge.newknots(i2)).and.(y.lt.newknots(i2+1))) then 
           result=1.d0
          else 
           result=0.d0
          endif
         else
          call dumsub(i1, i2, (i3-1), y, newknots, result1, dumsub)
          temp1=(y-newknots(i2))*result1/(newknots(i2+i3-1)-
     *    newknots(i2))
          if(temp1.ne.temp1) temp1=0.d0
          call dumsub(i1, (i2+1), (i3-1), y, newknots, result2, dumsub)
          temp2=(newknots(i2+i3)-y)*result2/(
     *    newknots(i2+i3)-newknots(i2+1))
          if(temp2.ne.temp2) temp2=0.d0
          result=temp1+temp2
         endif
         return 
         end 
