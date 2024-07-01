
C ------------------------------------------------------------------
       SUBROUTINE gaussj(a,n,np)!
       INTEGER*4 m,mp,n,np,NMAX
       DOUBLE PRECISION a(np,np)! 
       PARAMETER (NMAX=50) 
C      Linear equation solution by Gauss-Jordan elimination, 
c      equation (2.1.1) above. 
c      a(1:n,1:n) is an input matrix stored in an array of 
c      physical dimensions np by np.
c      b(1:n,1:m) is an input matrix containing the 
c      m right-hand side vectors, 
c      stored in an array of physical dimensions np by mp. 
c      On output, a(1:n,1:n) is 
c      replaced by its matrix inverse, and b(1:n,1:m) is 
c      replaced by the corresponding 
c      set of solution vectors. Parameter: NMAX is the largest 
c      anticipated value of n. 
        INTEGER*4 i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     +  ipiv(NMAX)   ! The integer arrays ipiv, indxr, and indxc are used 
         double precision big,dum,pivinv    !for bookkeeping on the pivoting.
        do 11 j=1,n 
          ipiv(j)=0
 11     continue 
         do 22 i=1,n !This is the main loop over the columns to be reduced. 
          big=0.0d0
          do 13 j=1,n! This is the outer loop of the search for aivot element. 
            if(ipiv(j).ne.1) then 
              do 12 k=1,n 
               if (ipiv(k).eq.0) then 
                 if (abs(a(j,k)).ge.big)then 
                   big=abs(a(j,k)) 
                   irow=j 
                   icol=k 
                 endif 
               endif 
 12           continue
            endif 
 13       continue
          ipiv(icol)=ipiv(icol)+1
c       We now have the pivot element, so we interchange rows, 
c       if needed, to put the pivot element on the diagonal. The columns 
c       are not physically interchanged, only relabeled:
c        indxc(i), the column of the ith pivot element, is the ith column that 
c       is reduced, while indxr(i) is the row in which that pivot element 
c       was originally located. 
c      If indxr(i)  = indxc(i) there is an implied column
c       interchange. With this form of bookkeeping, the 
c       solution b s will end up in the correct order, and 
c       the inverse matrix will be scrambled by columns.
         if (irow.ne.icol) then 
           do 14 l=1,n 
               dum=a(irow,l) 
               a(irow,l)=a(icol,l) 
               a(icol,l)=dum 
 14        continue 
c           do 15 l=1,m 
c               dum=b(irow,l) 
c               b(irow,l)=b(icol,l) 
c               b(icol,l)=dum 
c 15        continue 
         endif 
         indxr(i)=irow   ! We are now ready to divide the pivot 
         indxc(i)=icol   !row by the pivot element, located at irow and icol. 
         if (a(icol,icol).eq.0.) return! pause  'singular matrix in gaussj'
         pivinv=1.0d0/a(icol,icol) 
         a(icol,icol)=1.0d0
         do 16 l=1,n 
             a(icol,l)=a(icol,l)*pivinv 
 16      continue
c         do 17 l=1,m 
c             b(icol,l)=b(icol,l)*pivinv
c 17      continue
         do 21 ll=1,n               ! Next, we reduce the rows... 
              if(ll.ne.icol)then    !...except for the pivot one, of course. 
                  dum=a(ll,icol) 
                  a(ll,icol)=0.0d0
                  do 18 l=1,n 
                      a(ll,l)=a(ll,l)-a(icol,l)*dum 
 18               continue 
c                  do 19 l=1,m 
c                      b(ll,l)=b(ll,l)-b(icol,l)*dum 
c 19               continue
              endif 
 21       continue
 22     continue ! This is the end of the main 
c                ! loop over columns of the reduction.
        do 24 l=n,1,-1        !It only remains to unscramble the 
c                             !solution in view of the column interchanges. 
           if(indxr(l).ne.indxc(l))then      !We do this by 
c                                            !interchanging pairs of columns 
c        in the reverse order that the permutation was built up. 
               do 23 k=1,n 
                  dum=a(k,indxr(l)) 
                  a(k,indxr(l))=a(k,indxc(l)) 
                  a(k,indxc(l))=dum 
 23            continue
           endif 
 24     continue 
        return                      ! And we are done. 
        END
C ------------    
       double precision function indicator(i1, i2)
       implicit none
       integer*4 i1, i2
       if(i1.eq.i2) then 
        indicator=1.d0
       else 
        indicator=0.d0
       endif
       return 
       end
C ============================
C       The following  method is for  many subtypes along with standard error calculation
C       This method can handle a scalar covariate x.
         subroutine eewtsd(nf, npsi, nalpha, ndim,
     *   x, delta, initialpsi, psi, alpha, z, zn, p1, ndimz, ndimzn,
     *   numlevels,   der_mat, t1, n, totsubtype, nrep, newsize1,
     *   newmat, ss_score, store_eps, store_iteration, store_u)
        
         implicit none
        
         integer*4  nalpha, ndim, ndimz, ndimzn,  npsi, nf, 
     *   n, totsubtype, irep, nrep, store_iteration
         integer*4 numlevels(npsi), t1(nf), n1, iu 
         double precision  x(n), delta(n), theta(ndimz),   
     *   alpha(nalpha),  z(totsubtype, ndimz), 
     *   zn(totsubtype, nalpha), psi(npsi)
         double precision  p1(nf, totsubtype), initialpsi(npsi)
         double precision  beta(totsubtype), xi(totsubtype)
         double precision w_u(nf, totsubtype),
     *   w(nf, totsubtype), wx(nf, totsubtype), new_u(nf, ndim),
     *   u_alpha(nf, totsubtype), sum_1, 
     *   u_theta(nf, totsubtype), !anewvar(n, n), 
     *   s_new_u(ndim) 
         double precision s_1(nf, totsubtype), s_u_1(nf, totsubtype), 
     *   s_2(nf, totsubtype), indicator 
         double precision s_0(nf), s_w(nf), s_u_0(nf), s_w_u(nf)
         integer*4 i, i1, i2, j, j1, j2,  k, k1, k2, jnew, 
     *   new_index,  newsize1(n), newmat(n, nf)
         double precision der_mat(ndim, ndim),
     *   s_1nox(nf, totsubtype), s_u_1withx(nf, totsubtype), 
     *   s2_new_u(ndim, ndim), eps, old_psi(npsi), old_alpha(nalpha), 
     *   oldpara(ndim), inc(ndim), tempo_w
C         double precision newmat(n, nf)
         double precision ss_score(ndim, ndim),
     *   correction_theta(totsubtype), correction_alpha(totsubtype), 
     *   omega, omega_u, store_eps,  store_u(ndim)

C       totsubtype=total number of disease subtypes
C       ndischar=number of disease characteristics
C       numlevels=a vector of length ndischar which gives the number
C       of levels of each disease characteristics
C       ndimz= the number columns of Z, that depends on the type of model 
C       we use for this purpose. 
C       ndimzn= the number columns of Zn, that depends on the type of model 
C       we use for this purpose. 
C       score4theta and score4alpha =user supplied score vector for theta and alpha   
C       Defining  theta and alpha
C ------------------------------------------------------
C       npsi=length of psi vector=ndischar+1
C       nalpha=
C       nf=number of failures
C ------------------
        do i=1, npsi
         psi(i)=initialpsi(i)
        end do
        do i=1, nalpha
         alpha(i)=0.01d0
        end do
C        print*, 'check 0'
        eps=2.d0
        irep=1
        do while((eps.gt.0.00001d0).and.(irep.lt.nrep))
         irep=irep+1 
C ------------------
        do j=1, totsubtype
         beta(j)=0.d0
         xi(j)=0.d0
         do k=1, ndimz
          beta(j)=beta(j)+z(j, k)*psi(k)
         end do
         do k1=1, ndimzn
          xi(j)=xi(j)+zn(j, k1)*alpha(k1)
         end do
        end do
C        print*, '*************************'
C        print*, 'alpha= ', alpha
C        print*, 'psi= ', psi
C ---
         do k1=1, ndim 
          do k2=1, ndim 
           der_mat(k1, k2)=0.d0
          end do
         end do
C        print*, 'hello'
C ===== The end of initialization 
         do i1=1, nf
C          print*, 'i1= ', i1
          s_w(i1)=0.d0
          s_0(i1)=0.d0
          s_u_0(i1)=0.d0
          n1=t1(i1)
          
          do j=1, totsubtype
           w(i1, j)=exp(xi(j)+x(n1)*beta(j))*p1(i1, j)
           wx(i1, j)=w(i1, j)*x(n1)
           s_w(i1)=s_w(i1)+w(i1, j)
C ------- initialization of S_2, S_1 and S_0
           s_2(i1, j)=0.d0
           s_1(i1, j)=0.d0
           s_u_1withx(i1, j)=0.d0
C-------- initialization of S_U_1 and S_U_0
           s_1nox(i1, j)=0.d0          
           s_u_1(i1, j)=0.d0
           do i2=n1, n
            tempo_w=exp(xi(j)+x(i2)*beta(j))
            s_1(i1, j)=s_1(i1, j)+tempo_w*p1(i1, j)*x(i2)
            s_2(i1, j)=s_2(i1, j)+tempo_w*p1(i1, j)*x(i2)**2
            s_1nox(i1, j)=s_1nox(i1, j)+tempo_w*p1(i1, j)
            s_0(i1)=s_0(i1)+tempo_w*p1(i1, j) 
            s_u_0(i1)=s_u_0(i1)+tempo_w
            s_u_1(i1, j)=s_u_1(i1, j)+tempo_w
            s_u_1withx(i1, j)=s_u_1withx(i1, j)+tempo_w*x(i2) 
           end do 
          end do
C--- the score functions for theta
          do j=1, totsubtype
           u_theta(i1, j)=0.d0 
           u_theta(i1, j)=wx(i1, j)/s_w(i1)-(s_1(i1, j)/s_0(i1))
C --- the score functions for alpha
           u_alpha(i1, j)=0.d0 
           u_alpha(i1, j)=w(i1, j)/s_w(i1)-(s_u_1(i1, j)/s_u_0(i1))
          end do
C          print*, 'check post -I'
C ---- calculation of the hessian matrix
          do k1=1, npsi 
           do k2=1, npsi   
            do j1=1, totsubtype
             do j2=1, totsubtype  
              der_mat(k1, k2)=der_mat(k1, k2)+
     *        (w(i1, j1)*x(n1)**2*
     *        indicator(j1, j2)/s_w(i1)-
     *        wx(i1, j1)*wx(i1, j2)/s_w(i1)**2-
     *        s_2(i1, j1)*indicator(j1, j2)/s_0(i1)+ 
     *        s_1(i1, j1)*s_1(i1, j2)/s_0(i1)**2 )*
     *        z(j1, k1)*z(j2, k2)
             end do
            end do
           end do
          end do  
C ----          
          do k1=1, npsi
           do k2=(npsi+1), ndim 
            do j1=1, totsubtype
             do j2=1, totsubtype  
              der_mat(k1, k2)=der_mat(k1, k2)+(w(i1, j1)*x(n1)*
     *        indicator(j1, j2)/s_w(i1)-wx(i1, j1)*w(i1, j2)/s_w(i1)**2-
     *        s_1(i1, j1)*indicator(j1, j2)/s_0(i1)+ 
     *        s_1(i1, j1)*s_1nox(i1, j2)/s_0(i1)**2)*
     *        z(j1, k1)*zn(j2, (k2-npsi))
             end do
            end do
           end do
          end do  
C ----
          do k1=(npsi+1), ndim 
           do k2=1, npsi 
            do j1=1, totsubtype
             do j2=1, totsubtype 
              der_mat(k1, k2)=der_mat(k1, k2)+(w(i1, j1)*x(n1)*
     *        indicator(j1, j2)/s_w(i1)-w(i1, j1)*wx(i1, j2)/s_w(i1)**2-
     *        s_u_1withx(i1, j1)*indicator(j1, j2)/s_u_0(i1)+ 
     *        s_u_1(i1, j1)*s_u_1withx(i1, j2)/s_u_0(i1)**2)*
     *        z(j2, k2)*zn(j1, (k1-npsi))
             end do
            end do
           end do
          end do
           
C -------  
          do k1=(npsi+1), ndim 
           do k2=(npsi+1), ndim  
            do j1=1, totsubtype
             do j2=1, totsubtype 
              der_mat(k1, k2)=der_mat(k1, k2)+(w(i1, j1)*
     *        indicator(j1, j2)/s_w(i1)-w(i1, j1)*w(i1, j2)/s_w(i1)**2-
     *        s_u_1(i1, j1)*indicator(j1, j2)/s_u_0(i1)+ 
     *        s_u_1(i1, j1)*s_u_1(i1, j2)/s_u_0(i1)**2)*
     *        zn(j2, (k2-npsi))*zn(j1, (k1-npsi))
             end do
            end do
           end do
          end do
C ---------------
C          print*, 'check-post-III'
C          print*, 'der_mat= '
C          print10, ((der_mat(j1, j2), j2=1, ndim), j1=1, ndim)
 10       format(12f8.4)
         end do
C -----
         do i1=1, nf
          do k=1, npsi
           new_u(i1, k)=0.d0
           do j=1, totsubtype
            new_u(i1, k)=new_u(i1, k)+u_theta(i1, j)*z(j, k)
           end do
          end do
         end do
         do i=1, nf
          do k=(npsi+1), ndim
           new_u(i, k)=0.d0
           do j=1, totsubtype
            new_u(i, k)=new_u(i, k)+u_alpha(i, j)*zn(j, (k-npsi))
           end do
          end do
         end do
         do k=1, ndim
          s_new_u(k)=0.0d0
          do i1=1, nf
           s_new_u(k)=s_new_u(k)+new_u(i1, k)       
          end do
         end do

C ----- 
        call gaussj(der_mat, ndim, ndim)
        do k=1, ndim 
         inc(k)=0.d0
         do k1=1, ndim
          inc(k)=inc(k)-der_mat(k, k1)*s_new_u(k1)
         end do
        end do
        eps=0.d0
        do k=1, npsi
         oldpara(k)=psi(k)
         psi(k)=psi(k)+inc(k)
         eps=eps+abs(inc(k))/max(abs(oldpara(k)), 0.001d0)
         if(abs(psi(k)).gt. 5.d0) psi(k)=0.5d0*psi(k)/abs(psi(k))
        end do
        do k=(npsi+1), ndim
         oldpara(k)=alpha(k-npsi)
         alpha(k-npsi)=alpha(k-npsi)+inc(k)
         eps=eps+abs(inc(k))/max(abs(oldpara(k)), 0.001d0)
	 if(abs(alpha(k-npsi)).gt. 5.d0) alpha(k-npsi)=
     *   0.50d0*alpha(k-npsi)/abs(alpha(k-npsi))
        end do

C        print*, 'eps= ', eps
C        print*, 'alpha= ', alpha
C        print*, 'theta= ', theta
C        print*, 'irep= ', irep
C        print*, 'u(k)= ', s_new_u
       end do
       store_eps=eps
       store_iteration=irep
       do k1=1, ndim 
        store_u(k1)=s_new_u(k1)
       end do

C ==== This is the end of estimation 
C ==== The is the begining of the calculation of the std error.         
       new_index=0
       do k1=1, ndim 
        do k2=1, ndim
         ss_score(k1, k2)=0.d0 
        end do
       end do
C	   print*, 'the begining of the final check point'
       do j=1, n
C        print*, 'j= ', j
        if(delta(j).eq.1.d0) then 
         new_index=new_index+1
         do j1=1, totsubtype
          correction_theta(j1)= u_theta(new_index, j1)
          correction_alpha(j1)= u_alpha(new_index, j1)
         end do
        else 
         do j1=1, totsubtype
          correction_theta(j1)=0.d0
          correction_alpha(j1)=0.d0
         end do
        endif
C -------
        if(newsize1(j).gt.0) then 
         do i1=1, newsize1(j)
          omega=0.d0
          omega_u=0.d0
          do j1=1, totsubtype
           omega=omega+exp(xi(j1)+x(j)*beta(j1))*
     *     p1(newmat(j, i1), j1) 
           omega_u=omega_u+exp(xi(j1)+x(j)*beta(j1))   
          end do
          do j1=1, totsubtype
           correction_theta(j1)=correction_theta(j1)-(omega/s_0(i1))*
     *     ((exp(xi(j1)+x(j)*beta(j1))*
     *     p1(newmat(j, i1), j1)*x(j)/omega)-
     *      (s_1(i1, j1)/s_0(i1)))
           
           correction_alpha(j1)=correction_alpha(j1)-
     *     (omega_u/s_u_0(i1))*((exp(xi(j1)+x(j)*beta(j1))/
     *      omega_u)- 
     *     (s_u_1(i1, j1)/s_u_0(i1)))
          end do
         end do
        endif
C		print*, 'check post'
        do k1=1, npsi
         do k2=1, npsi
          do j1=1, totsubtype
           do j2=1, totsubtype
            ss_score(k1, k2)=ss_score(k1, k2)+correction_theta(j1)*
     *      correction_theta(j2)*z(j1, k1)*z(j2, k2)
           end do
          end do 
         end do
        end do
        do k1=1, npsi
         do k2=(npsi+1), ndim
          do j1=1, totsubtype
           do j2=1, totsubtype
            ss_score(k1, k2)=ss_score(k1, k2)+correction_theta(j1)*
     *      correction_alpha(j2)*z(j1, k1)*zn(j2, (k2-npsi))
           end do
          end do 
         end do
        end do
        do k1=(npsi+1), ndim
         do k2=(npsi+1), ndim
          do j1=1, totsubtype
           do j2=1, totsubtype
            ss_score(k1, k2)=ss_score(k1, k2)+correction_alpha(j1)*
     *      correction_alpha(j2)*zn(j1, (k1-npsi))*zn(j2, (k2-npsi))
           end do
          end do 
         end do
        end do
C -----------		
       end do
C	   print*, 'last check point'
C =============       
        do k1=(npsi+1), ndim
         do k2=1, npsi
          ss_score(k1, k2)=ss_score(k2, k1)
         end do
        end do
       return 
       end 

