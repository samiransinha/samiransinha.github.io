        subroutine calhapprob(vecdip, ncount, nr, ndip, prob, pi, mpi,
     *   nref)
        IMPLICIT NONE
        integer*4 nr, ndip, i1, i2, mpi, i, nref
        integer*4 ncount(nr), vecdip(ndip, 2)
        double precision prob(ndip), pi(mpi), s
        integer*4 j
        do i=1, nr
         if (i .eq. 1) then 
          i1=1
         else
          i1=ncount(i-1)+1
         endif
         i2=ncount(i)
         s=0.d0
          do j=i1, i2
           if(vecdip(j, 1) .eq. vecdip(j, 2)) then 
            prob(j)=pi(vecdip(j, 2))*pi(vecdip(j, 2))
           else
            prob(j)=2*pi(vecdip(j, 1))*pi(vecdip(j, 2))
           endif
           s=s+prob(j)
          end do
          do j=i1, i2
           prob(j)=prob(j)/s
           if(i.lt.10) then 
           print*,' prob= ','i= ', i, prob(j)
           endif
          end do
        end do
        return 
        end
C *************************************************************
        integer function identical(a, b)
        integer*4 a, b
        if(a.eq.b) then 
         identical=1
        else
         identical=0
        endif
        return
        end
C *************************************************************
C  Determination of pi for the proposed method  using additive model. 
        subroutine determinepi(nr, nstr, ncnt, nl, nonzeropi, pcrs,  
     *  ndip, vecdip, beta, pi,  new_pi, mpi, nbet, nref, m)    
        IMPLICIT NONE  
        integer*4 nr, nl, ndip, nstr, identical, mpi, temp2, nbet, m
        integer*4 nonzeropi(nl), ncnt(nr), vecdip(ndip, 2), nref
        double precision epsilon, pi(mpi), old_pi(mpi), 
     *  new_pi(mpi), crs(ndip), beta(nbet)
        double precision store_s1(nl), store_qnty(nl), a, al, 
     *  lambda_u, lambda_l, lb, ub, pcrs(ndip), temp
        integer*4 iteration, kount, u, k1, k2, j, i, k, j0
        double precision a1, a1p, a2, a2p, b1, b1p, b2, b2p, sum_pi
        double precision qnty4, qnty5, s1, tempo_1, tempo_2
C        print*, 'step------0'
        do u=1, mpi
         new_pi(u)=0.d0
        end do
        epsilon=0.20d0
        iteration=1
        do while((epsilon.gt .0.0001d0).and.(iteration.lt.200))
         iteration=iteration+1
         do u=1, nl
          s1=0.d0
          do j=1, nr
           if( j.eq.1)then 
            k1=1 
           else 
            k1=ncnt(j-1)+1
           endif  
           k2=ncnt(j)
           do k=k1, k2
            s1=s1+pcrs(k)*(identical(vecdip(k, 1), nonzeropi(u))+
     *      identical(vecdip(k, 2), nonzeropi(u)))
           end do         
          end do
          qnty4=0.0d0
          qnty5=0.0d0
          do i=1, nstr
           j=(m+1)*(i-1)+1
           if(j.eq. 1) then
            k1=1 
           else 
            k1=ncnt(j-1)+1
           endif
           k2=ncnt(j)
           a1=0.0d0
           a1p=0.0d0 
           a2=0.0d0
           a2p=0.0d0
           do k= k1, k2
            temp=0.d0
            if(vecdip(k, 1).ne. nref) temp=beta(vecdip(k, 1)) 
            if(vecdip(k, 2).ne. nref) temp=temp+beta(vecdip(k, 2)) 
            a1=a1+exp(temp)*(-identical(vecdip(k, 1), 
     *      vecdip(k, 2))+2)*pi(vecdip(k, 1))*pi(vecdip(k, 2))
            a1p=a1p+exp(temp)*(-identical(vecdip(k, 1), 
     *      vecdip(k, 2))+2)*
     *      (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *      pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
            a2=a2+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *      pi(vecdip(k, 1))*pi(vecdip(k, 2))
            a2p=a2p+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *      (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *      pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
           end do
           qnty4=qnty4+a2p/a2
           tempo_1=(a2*a1p-a1*a2p)/(a2*a2)
           tempo_2=(a1/a2)
           do j0=1, m
            j=(m+1)*(i-1)+(1+j0)
            k1=ncnt(j-1)+1
            k2=ncnt(j)
            b1=0.0d0
            b1p=0.0d0
            b2=0.0d0
            b2p=0.0d0
            do k= k1, k2
             temp=0.d0
             if(vecdip(k,1).ne. nref) temp=beta(vecdip(k, 1)) 
             if(vecdip(k,2).ne. nref) temp=temp+beta(vecdip(k, 2))
             b1=b1+exp(temp)*(-identical(vecdip(k, 1), 
     *       vecdip(k, 2))+2)*pi(vecdip(k, 1))*pi(vecdip(k, 2))
             b1p=b1p+exp(temp)*(-identical(vecdip(k, 1), 
     *       vecdip(k, 2))+2)*
     *       (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *       pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
             b2=b2+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *       pi(vecdip(k, 1))*pi(vecdip(k, 2))
             b2p=b2p+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *       (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *       pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
            end do
            qnty4=qnty4+(b2p/b2)
            tempo_1=tempo_1+(b2*b1p-b1*b2p)/(b2*b2)
            tempo_2=tempo_2+(b1/b2)
           end do ! for j0
           qnty5=qnty5+tempo_1/tempo_2
          end do ! for i
          store_s1(u)=s1
          store_qnty(u)= qnty4+qnty5
         end do
         do u=1, nl
          old_pi(nonzeropi(u))=pi(nonzeropi(u))
         end do
          a=0.100d0
          lambda_u=5.d0
          lambda_l=-5.d0
          lb=1.d0
          ub=1.d0
          do u=1, nl
           lb=lb-store_s1(u)/(store_qnty(u)+lambda_l)
           ub=ub-store_s1(u)/(store_qnty(u)+lambda_u)
          end do
          kount=0
          do while((abs(a).gt. 0.000001d0).and.(kount.lt.100))
           kount=kount+1
           al=(lambda_l+lambda_u)/2.d0
           a=1.d0
           do u=1, nl
            a=a-store_s1(u)/(store_qnty(u)+al)
           end do
           if(a.gt.0.d0) then 
            lambda_u=al 
           else 
            lambda_l=al
           endif
          end do
C         print*, 'a= ', a, ' kount= ', kount
          sum_pi=0.d0
          epsilon=0.d0
          do u=1, nl
           pi(nonzeropi(u))=store_s1(u)/(store_qnty(u)+al)
           sum_pi=sum_pi+pi(nonzeropi(u))
          end do
          do u=1, nl
           pi(nonzeropi(u))=pi(nonzeropi(u))/sum_pi
           epsilon=epsilon+abs(old_pi(nonzeropi(u))-
     *     pi(nonzeropi(u)))/old_pi(nonzeropi(u))
          end do
        end do
         do u=1, mpi
         new_pi(u)=pi(u)
         end do
         return
         end
C -----------------------------------------------------
C  Determination of pi for the proposed method  using dominant  model. 
        subroutine determinepi2(nr, nstr, ncnt, nl, nonzeropi, pcrs,  
     *  ndip, vecdip, beta, pi,  new_pi, mpi)  
        integer*4 nr, nl, ndip, nstr, identical, mpi, temp2
        integer*4 nonzeropi(nl), ncnt(nr), vecdip(ndip, 2)
        double precision epsilon, pi(mpi), old_pi(mpi), 
     *  new_pi(mpi), crs(ndip), beta
        double precision store_s1(nl), store_qnty(nl), a, al, 
     *  lambda_u, lambda_l, lb, ub, pcrs(ndip), temp
        integer*4 iteration, kount, u, k1, k2, j, i, k
        double precision a1, a1p, a2, a2p, b1, b1p, b2, b2p, sum_pi
        double precision qnty4, qnty5, s1
C        print*, 'step------0'
        do u=1, mpi
        new_pi(u)=0.d0
        end do
        epsilon=0.20d0
        iteration=1
        do while((epsilon.gt .0.0001d0).and.(iteration.lt.200))
        iteration=iteration+1
C       store_s1=NULL; store_qnty=NULL;
C       newpi=rep(0, 32)
C       nonzeropi=(1:32)[pi!=0]
C       nl=length(nonzeropi)
        do u=1, nl
C        print*, 'step---------I'
        s1=0.d0
          do j=1, nr
            if( j.eq.1)then 
               k1=1 
            else 
               k1=ncnt(j-1)+1
            endif  
            k2=ncnt(j)
           do k=k1, k2
            s1=s1+pcrs(k)*(identical(vecdip(k, 1), nonzeropi(u))+
     *      identical(vecdip(k, 2), nonzeropi(u)))
           end do         
          end do
          qnty4=0.0d0
          qnty5=0.0d0
          do i=1, nstr
C           print*, 'i ===', i
C           print*, 'check-----I'
           j=2*(i-1)+1
C             print*, 'check-----II'
           if(j.eq. 1) then
            k1=1 
           else 
            k1=ncnt(j-1)+1
           endif
           k2=ncnt(j)
           a1=0.0d0
           a1p=0.0d0 
           a2=0.0d0
           a2p=0.0d0
           do k= k1, k2
C           if((vecdip(k, 1).eq. 3).or.(vecdip(k, 2).eq. 3)) then ! dominnat model ! for example 1
C            if((vecdip(k, 1).eq. 3).and .(vecdip(k, 2).eq. 3)) then ! for recissive model ! for example 1
           if((vecdip(k, 1).eq. 1).and .(vecdip(k, 2).eq. 1)) then 
            temp=beta
           else
            temp=0.d0
           endif
C          print*, 'check--------III'
           a1=a1+exp(temp)*(-identical(vecdip(k, 1), 
     *     vecdip(k, 2))+2)*pi(vecdip(k, 1))*pi(vecdip(k, 2))
           a1p=a1p+exp(temp)*(-identical(vecdip(k, 1), 
     *     vecdip(k, 2))+2)*
     *     (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *     pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
           a2=a2+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *     pi(vecdip(k, 1))*pi(vecdip(k, 2))
           a2p=a2p+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *     (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *     pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
          end do
          j=2*(i-1)+2
          k1=ncnt(j-1)+1
          k2=ncnt(j)
          b1=0.0d0
          b1p=0.0d0
          b2=0.0d0
          b2p=0.0d0
         
           do k= k1, k2
C           print*, 'oooooooooooooooo=', vecdip(k, 1), vecdip(k, 2)
C          if((vecdip(k,1).eq. 3).or. (vecdip(k, 2).eq. 3)) then !
C          if((vecdip(k, 1).eq. 3).and .(vecdip(k, 2).eq. 3)) then ! for recissive model 
           if((vecdip(k, 1).eq. 1).and.(vecdip(k, 2).eq. 1)) then 
            temp=beta
           else
            temp=0.d0
           endif
            b1=b1+exp(temp)*(-identical(vecdip(k, 1), 
     *     vecdip(k, 2))+2)*pi(vecdip(k, 1))*pi(vecdip(k, 2))
           b1p=b1p+exp(temp)*(-identical(vecdip(k, 1), 
     *     vecdip(k, 2))+2)*
     *     (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *     pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
           b2=b2+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *     pi(vecdip(k, 1))*pi(vecdip(k, 2))
           b2p=b2p+(-identical(vecdip(k, 1), vecdip(k, 2))+2)*
     *     (pi(vecdip(k, 1))*identical(vecdip(k, 2), nonzeropi(u))+
     *     pi(vecdip(k, 2))*identical(vecdip(k, 1), nonzeropi(u)))
           end do
           qnty4=qnty4+(a2p/a2)+(b2p/b2)
           qnty5=qnty5+ (((a2*a1p-a1*a2p)/(a2*a2))+((b2*b1p-b1*b2p)
     *     /(b2*b2)))/((a1/a2)+(b1/b2))
          end do
           store_s1(u)=s1
           store_qnty(u)= qnty4+qnty5
        end do
C        print*, 'step-------II'
          do u=1, nl
          old_pi(nonzeropi(u))=pi(nonzeropi(u))
          end do
          a=0.100d0
          lambda_u=5.d0
          lambda_l=-5.d0
          lb=1.d0
          ub=1.d0
          do u=1, nl
          lb=lb-store_s1(u)/(store_qnty(u)+lambda_l)
          ub=ub-store_s1(u)/(store_qnty(u)+lambda_u)
          end do
C          print*, lb, ub, '(((((((((('
          kount=0
          do while((abs(a).gt. 0.000001d0).and.(kount.lt.100))
            kount=kount+1
            al=(lambda_l+lambda_u)/2.d0
            a=1.d0
            do u=1, nl
             a=a-store_s1(u)/(store_qnty(u)+al)
            end do
            if(a.gt.0.d0) then 
              lambda_u=al 
            else 
              lambda_l=al
            endif
            
          end do
C          print*, 'a= ', a, ' kount= ', kount
          sum_pi=0.d0
          epsilon=0.d0
          do u=1, nl
           pi(nonzeropi(u))=store_s1(u)/(store_qnty(u)+al)
           sum_pi=sum_pi+pi(nonzeropi(u))
          end do
          do u=1, nl
           pi(nonzeropi(u))=pi(nonzeropi(u))/sum_pi
           epsilon=epsilon+abs(old_pi(nonzeropi(u))-
     *     pi(nonzeropi(u)))/old_pi(nonzeropi(u))
          end do
C          print*, 'step----------III'
         end do
         do u=1, mpi
         new_pi(u)=pi(u)
         end do
         return
         end
C **************************************************************
C Calculation of c_rs of the ECM algorithm with additive model
         subroutine estepofem(nr, nstr, ndip, ncnt, vecdip, crs, pi, 
     *   beta, mpi, nbet, nref, m)
         IMPLICIT NONE
         integer*4 nr, nstr, ndip, identical, mpi, nbet, nref, m
         integer*4 ncnt(nr), vecdip(ndip, 2), j0
         double precision pi(mpi), crs(ndip), s, beta(nbet), temp
         integer*4 i, j, k1, k2, j1
         do i=1, ndip
         crs(i)=0.d0
         end do
         do i=1, nstr
          j1=(m+1)*(i-1)+1
          if( j1.eq.1)then 
           k1=1 
          else 
           k1=ncnt(j1-1)+1
          endif
          k2=ncnt(j1)
          s=0.d0
          do j=k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne. nref) temp=beta(vecdip(j, 1))
           if(vecdip(j, 2).ne. nref) temp=temp+beta(vecdip(j, 2))
           s=s+exp(temp)*
     *     (-identical(vecdip(j, 1), 
     *     vecdip(j, 2))+2)*pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          do j=k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne. nref) temp=beta(vecdip(j, 1))
           if(vecdip(j, 2).ne. nref) temp=temp+beta(vecdip(j, 2))
C ---------------------
           crs(j)=exp(temp)*(-identical(vecdip(j, 1), 
     *     vecdip(j, 2))+2)*pi(vecdip(j, 1))*pi(vecdip(j, 2))/s
          end do
C------- Folowwing part is for the controls 
          do j0=1, m
           j1=(m+1)*(i-1)+(j0+1)
           k1=ncnt(j1-1)+1
           k2=ncnt(j1)
           s=0.d0
           do j=k1, k2
            s=s+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *      pi(vecdip(j, 1))*pi(vecdip(j, 2))
           end do
           do j=k1, k2
            crs(j)=(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *      pi(vecdip(j, 1))*pi(vecdip(j, 2))/s
           end do
          end do
         end do
         return 
         end
C ************************************************************
C  Function of beta for the proposed method and for the additive model
         subroutine funcbeta(nr, nstr, ndip, crs, vecdip, ncnt, beta, 
     *   pi, loglk, mpi, nbet, nref, m)
         IMPLICIT NONE 
         integer*4 nr, nstr, ndip, identical, mpi, nbet, nref, m
         integer*4 ncnt(nr), vecdip(ndip, 2), j0
         double precision beta(nbet), pi(mpi), crs(ndip), temp, tempo
         double precision qnty1, qnty2, qnty3, s1, s2, s3, s4, loglk
         integer*4 i, j1, j, k1, k2
C         print*, 'beta= ', beta
         qnty1=0.d0
         qnty2=0.d0
         qnty3=0.d0
C         print*, 'check---------O'
         do i=1, nstr
          j1=(m+1)*(i-1)+1
          if(j1.eq.1) then 
           k1=1 
          else 
           k1=ncnt(j1-1)+1
          endif
          k2=ncnt(j1)
          s1=0.d0
          s3=0.d0
          
          do j=k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne. nref) temp=beta(vecdip(j, 1)) 
           if(vecdip(j, 2).ne. nref) temp=temp+beta(vecdip(j, 2))
           s1=s1+exp(temp)*
     *     (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
           s3=s3+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          
          qnty1=qnty1+log(s1/s3)
          
C -----------------------
          do j=k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne. nref) temp=beta(vecdip(j, 1)) 
           if(vecdip(j, 2).ne. nref) temp=temp+beta(vecdip(j, 2))
           qnty3=qnty3+crs(j)*(temp+log(
     *     (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2)))-log(s1))
          end do
C -----------------------
          tempo=s1/s3
          do j0=1, m
           j1=(m+1)*(i-1)+(j0+1)
           k1=ncnt(j1-1)+1
           k2=ncnt(j1)
           s2=0.d0
           s4=0.d0
           do j = k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne. nref) temp=beta(vecdip(j, 1)) 
           if(vecdip(j, 2).ne. nref) temp=temp+beta(vecdip(j, 2))
            s2=s2+exp(temp)*
     *      (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *      pi(vecdip(j, 1))*pi(vecdip(j, 2))
            s4=s4+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *      pi(vecdip(j, 1))*pi(vecdip(j, 2))
           end do
           tempo=tempo+s2/s4           
          end do
          qnty2=qnty2-log(tempo)
         end do
         loglk=-(qnty1+qnty2+qnty3) ! this gives negative of log likelihood 
C         if(qnty1.ne. qnty1)print*, qnty1, qnty2, qnty3
         return
         end
C -----------------------------------------------------------
C Calculation of c_rs of the ECM algorithm using dominant model
         subroutine estepofem2(nr, nstr, ndip, ncnt, vecdip, crs, pi, 
     *   beta, mpi)
         integer*4 nr, nstr, ndip, identical, mpi
         integer*4 ncnt(nr), vecdip(ndip, 2)
         double precision pi(mpi), crs(ndip), s, beta, temp
         integer*4 i, j, k1, k2, j1
         do i=1, ndip
         crs(i)=0.d0
         end do
         do i=1, nstr
          j1=2*(i-1)+1
          if( j1.eq.1)then 
           k1=1 
          else 
           k1=ncnt(j1-1)+1
          endif
          k2=ncnt(j1)
          s=0.d0
          do j=k1, k2
C           if((vecdip(j,1).eq. 3).or. (vecdip(j, 2).eq. 3)) then !for dominant model! for example 1 
C           if((vecdip(j,1).eq. 3).and. (vecdip(j, 2).eq. 3)) then ! for recessive model ! for example 2
            if((vecdip(j,1).eq. 1).and. (vecdip(j, 2).eq. 1)) then 

            temp=beta
           else
            temp=0.d0
           endif
           s=s+exp(temp)*
     *     (-identical(vecdip(j, 1), 
     *     vecdip(j, 2))+2)*pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          do j=k1, k2
C           if((vecdip(j, 1).eq. 3).or.(vecdip(j, 2).eq. 3)) then !for dominant model ! for example 1 
C           if((vecdip(j,1).eq. 3).and. (vecdip(j, 2).eq. 3)) then ! for recessive model for example 2
             if((vecdip(j, 1).eq. 1).and.(vecdip(j, 2).eq. 1)) then 
            temp=beta
           else
            temp=0.d0
           endif
           crs(j)=exp(temp)*(-identical(vecdip(j, 1), 
     *     vecdip(j, 2))+2)*pi(vecdip(j, 1))*pi(vecdip(j, 2))/s
          end do
          j1=2*(i-1)+2
          k1=ncnt(j1-1)+1
          k2=ncnt(j1)
          s=0.d0
          do j=k1, k2
           s=s+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          do j=k1, k2
           crs(j)=(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))/s
          end do
         end do
         return 
         end
C ************************************************************
C        Function of beta for Chen and Chatterjee
         subroutine funcbetacc(nr, nstr, ndip, vecdip, ncnt, 
     *   beta, pi, loglk, nref, nbet, mupi, m)
         IMPLICIT NONE 
         integer*4 nr, nstr, ndip, identical, nref, nbet, mupi, m
         integer*4 ncnt(nr), vecdip(ndip, 2)
         double precision beta(nbet), pi(mupi)
         double precision qnty1, qnty2, s1, s2, s3, s4, loglk,
     *   temp, tempo
         integer*4 i, j1, j, k1, k2, j0
         qnty1=0.d0
         qnty2=0.d0
C         print*, 'hii'
         do i=1, nstr
         j1=(m+1)*(i-1)+1
           if(j1.eq.1) then 
             k1=1 
           else 
             k1=ncnt(j1-1)+1
           endif
           k2=ncnt(j1)
           s1=0.d0
           s3=0.d0
          do j=k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne.nref) temp=beta(vecdip(j, 1))
           if(vecdip(j, 2).ne.nref) temp=temp+beta(vecdip(j, 2))
          s1=s1+exp(temp)*
     *    (-identical(vecdip(j, 1), vecdip(j, 2))+2.d0)*
     *    pi(vecdip(j, 1))*pi(vecdip(j, 2))
          s3=s3+(-identical(vecdip(j, 1), vecdip(j, 2))+2.d0)*
     *    pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
C          print*, s1, s3
          if( (s1.ne.0.d0).and.(s3.ne.0.d0) ) then
             qnty1=qnty1+log(s1/s3)
C           print*, 'qnty1= ', qnty1
C -----------------------
          tempo=s1/s3
         do j0=1, m          
          j1=(m+1)*(i-1)+(j0+1)
          k1=ncnt(j1-1)+1
          k2=ncnt(j1)
          s2=0.d0
          s4=0.d0
          do j = k1, k2
           temp=0.d0
           if(vecdip(j, 1).ne.nref) temp=beta(vecdip(j, 1))
           if(vecdip(j, 2).ne.nref) temp=temp+beta(vecdip(j, 2))
           s2=s2+exp(temp)*
     *     (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
           s4=s4+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          tempo=tempo+s2/s4
         end do
          qnty2=qnty2-log(tempo)
         endif
         end do
         loglk=-(qnty1+qnty2)
C         print*,'*********', qnty1, qnty2
         return
         end
C--------------------------------------------------------------
C  Function of beta for the proposed method and for the dominant model
         subroutine funcbeta2(nr, nstr, ndip, crs, vecdip, ncnt, beta, 
     *   pi, loglk, mpi) 
         integer*4 nr, nstr, ndip, identical, mpi
         integer*4 ncnt(nr), vecdip(ndip, 2)
         double precision beta, pi(mpi), crs(ndip), temp
         double precision qnty1, qnty2, qnty3, s1, s2, s3, s4, loglk
         integer*4 i, j1, j, k1, k2
         qnty1=0.d0
         qnty2=0.d0
         qnty3=0.d0
C         print*, 'check---------O'
         do i=1, nstr
         
         j1=2*(i-1)+1
           if(j1.eq.1) then 
             k1=1 
           else 
             k1=ncnt(j1-1)+1
           endif
           k2=ncnt(j1)
           s1=0.d0
           s3=0.d0
          do j=k1, k2
C          if((vecdip(j, 1).eq. 3).or.(vecdip(j, 2).eq. 3)) then !for the dominant model ! for example 1
C          if((vecdip(j, 1).eq. 3). and.(vecdip(j, 2).eq. 3)) then ! for the recissive model ! for example 1
           if((vecdip(j, 1).eq. 1).and.(vecdip(j, 2).eq. 1)) then
           temp=beta
          else
           temp=0.d0
          endif
C          print*, 'check----------I', ' i= ', i
          s1=s1+exp(temp)*
     *    (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *    pi(vecdip(j, 1))*pi(vecdip(j, 2))
          s3=s3+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *    pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
C          print*, s1, s3
          qnty1=qnty1+log(s1/s3)
C -----------------------
          do j=k1, k2
C          if((vecdip(j, 1).eq. 3).or.(vecdip(j, 2).eq. 3)) then !for the dominant model
C          if((vecdip(j, 1).eq. 3). and.(vecdip(j, 2).eq. 3)) then ! for the recissive model
           if((vecdip(j, 1).eq. 1).and.(vecdip(j, 2).eq. 1)) then
           temp=beta
          else
           temp=0.d0
          endif
C          print*,'check---------II' 
          qnty3=qnty3+crs(j)*(temp+log(
     *    (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *    pi(vecdip(j, 1))*pi(vecdip(j, 2)))-log(s1))
          end do
C -----------------------
          j1=2*(i-1)+2
          k1=ncnt(j1-1)+1
          k2=ncnt(j1)
          s2=0.d0
          s4=0.d0
          do j = k1, k2
C          if((vecdip(j, 1).eq. 3).or.(vecdip(j, 2).eq.3)) then ! for the dominant model
C          if((vecdip(j, 1).eq. 3). and.(vecdip(j, 2).eq. 3)) then ! for the recissive model
          if((vecdip(j, 1).eq. 1).and.(vecdip(j, 2).eq. 1)) then !  for the second example
           temp=beta
          else
           temp=0.d0
          endif
C          print*,'check---------III' 
           s2=s2+exp(temp)*
     *     (-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
           s4=s4+(-identical(vecdip(j, 1), vecdip(j, 2))+2)*
     *     pi(vecdip(j, 1))*pi(vecdip(j, 2))
          end do
          qnty2=qnty2-log((s1/s3)+(s2/s4))
C          print*, 'check---------IV'
         end do
C          print*, 'qntys:= ', qnty1, qnty2, qnty3
         loglk=-(qnty1+qnty2+qnty3)
         print*, '*****= ', qnty1+qnty2
C          print*, 'loglk= ', loglk
C         print*,'*********', qnty1, qnty2, qnty3
         return
         end
C -------------------------------------------------------------
C       This subroutine calculates hap -frequencies by the EM algo. 
C       using only the control sample.
        subroutine expem(nr, ns, ncount, vecdip, ndip, pi, mpi, nref, m)
        integer*4  nr, ns, ndip, i1, i2, i0, j, mpi, m
        integer*4 vecdip(ndip, 2), ncount(nr), nref
        double precision  pi(mpi), oldpi(mpi), l, s(ns), temp_l, 
     *  prold(ndip), pr, eps !, prob(16)
        do i0=1, mpi
         oldpi(i0)=1.d0/mpi
        end do
        eps=2.d0
        k=0
        do while((eps. gt. 0.000001d0).and.(k.lt.500))
         k=k+1
         do i11=1,  mpi
          pi(i11)=0.d0
         end do
         do i0=1, ns
          do j0=1, m
           i=(1+m)*(i0-1)+(j0+1)
           i1=ncount(i-1)+1
           i2=ncount(i)
           s(i0)=0.d0
           do j=i1, i2
            if(vecdip(j, 1).eq.vecdip(j, 2)) then
             prold(j)=oldpi(vecdip(j, 1))*oldpi(vecdip(j, 1))
             s(i0)=s(i0)+prold(j)
            else
             prold(j)=2.d0*oldpi(vecdip(j, 1))*oldpi(vecdip(j, 2))
             s(i0)=s(i0)+prold(j)
            endif
           end do
           do j=i1, i2
            prold(j)=prold(j)/s(i0)
            pi(vecdip(j, 1))=pi(vecdip(j, 1))+prold(j)
            pi(vecdip(j, 2))=pi(vecdip(j, 2))+prold(j) 
           end do
          end do
         end do
C         print*, '************'
          eps=0.0d0
         do ij=1, mpi
          pi(ij)=pi(ij)/(2*ns)
          if(pi(ij). lt. 0.000001d0) pi(ij)=0.d0
          if( pi(ij).ne.0.d0) eps=eps+abs(pi(ij)-oldpi(ij))/oldpi(ij)
          oldpi(ij)=pi(ij)
         end do
        end do ! for while 
        return
        end
C ************************************************************
C       This subroutine calculates hap -frequencies in the general popl. 
C       by the EM algo using cases and controls.
        subroutine expemcc(nr, ns, ncount, vecdip, ndip, pi, mpi, nref)
        integer*4  nr, ns, ndip, i1, i2, i0, j, mpi, nref
        integer*4 vecdip(ndip, 2), ncount(nr)
        double precision  pi(mpi), oldpi(mpi), l, s(nr), temp_l, 
     *  prold(ndip), pr, eps !prob(16)
         do i0=1, mpi
          oldpi(i0)=1.d0/mpi
         end do
         print*, 'Kraft et al. (2005)'
        eps=2.d0
        k=0
        do while((eps. gt. 0.00001d0).and.(k.lt.500))
        k=k+1
C        print*, 'oldpi= ', oldpi
        do i11=1,  mpi
          pi(i11)=0.d0
        end do
        do i0=1, nr
C         print*, 'check', 'i0= ', i0
         i=i0!2*(i0-1)+2
         if(i.eq.1) then 
          i1=1
         else 
          i1=ncount(i-1)+1
         endif
         i2=ncount(i)
         s(i0)=0.d0
          do j=i1, i2
            if(vecdip(j, 1).eq.vecdip(j, 2)) then
              prold(j)=oldpi(vecdip(j, 1))*oldpi(vecdip(j, 1))
              s(i0)=s(i0)+prold(j)
            else
              prold(j)=oldpi(vecdip(j, 1))*oldpi(vecdip(j, 2))
              s(i0)=s(i0)+prold(j)
            endif
          end do
          do j=i1, i2
            prold(j)=prold(j)/s(i0)
             pi(vecdip(j, 1))=pi(vecdip(j, 1))+prold(j)
             pi(vecdip(j, 2))=pi(vecdip(j, 2))+prold(j) 
          end do
C          print*,'i0= ', i0
          
        end do
C         print*, '************'
          eps=0.0d0
        do ij=1, mpi
          pi(ij)=pi(ij)/(2*nr)
          if(pi(ij). lt. 0.0000001d0) pi(ij)=0.d0
          if( pi(ij).ne.0.d0) eps=eps+abs(pi(ij)-oldpi(ij))/oldpi(ij)
          oldpi(ij)=pi(ij)
        end do
C          print*, 'pi= ', pi
        end do
C        print*, 'k= ', k, ' eps = ', eps, ' pi= ', pi
        return
        end
C ************************************************************
C *****************************************************************
C       This calculates the estimates of beta
        subroutine calbetakraft(nr, vecdip, ndip, ncount, prob, beta,
     *  l, mpi, nbet, nref)
        IMPLICIT NONE
        integer*4 nr, ndip, ns, mpi, i1, i2, j1, j2, k0, k, i10, i20, i
        integer*4 vecdip(ndip, 2), ncount(nr), nbet, nref
        double precision pi(mpi), beta(nbet), l, q1, q2, temp_q_1, 
     *  temp_1, prob(ndip), temp_q_2, temp_q_12
C       print*, 'check......'
        ns=nr/2
        l=0.d0
        do i=1, ns
           j1=2*(i-1)+1
           if( j1 .eq. 1) then 
            i1=1
           else
            i1=ncount(j1-1)+1
           endif
           i2=ncount(j1)
           temp_1=0.d0
          do k=i1, i2
             j2=2*(i-1)+2
             i10=ncount(j2-1)+1
             i20=ncount(j2)
             do k0=i10, i20
C             if(prob(k)*prob(k0).ne.0.d0) then 
              temp_q_1=0.d0
C              if(vecdip(k, 1).ne. 1) temp_q_1=beta(vecdip(k, 1)-1)
C              if(vecdip(k, 2).ne. 1) temp_q_1=temp_q_1+
C     *        beta(vecdip(k, 2)-1)
C              temp_q_2=0.d0
C              if(vecdip(k0, 1).ne. 1) temp_q_2=beta(vecdip(k0, 1)-1)
C              if(vecdip(k0, 2).ne. 1) temp_q_2=temp_q_2+
C     *        beta(vecdip(k0, 2)-1)
              if(vecdip(k, 1).ne. nref) temp_q_1=beta(vecdip(k, 1))
              if(vecdip(k, 2).ne. nref) temp_q_1=temp_q_1+
     *        beta(vecdip(k, 2))
              temp_q_2=0.d0
              if(vecdip(k0, 1).ne. nref) temp_q_2=beta(vecdip(k0, 1))
              if(vecdip(k0, 2).ne. nref) temp_q_2=temp_q_2+
     *        beta(vecdip(k0, 2))

              temp_q_12=temp_q_2-temp_q_1
              temp_1=temp_1+prob(k)*prob(k0)/
     *        (1.d0+exp(temp_q_12))
C              endif
             end do
          end do
          l=l+log(temp_1)
C          print*, 'i= ', i, ' l= ', l
        end do   
        return
        end
C **************************************************************
C       This calculates the estimates of beta
        subroutine calbetazhang(nr, vecdip, ndip, ncount, pi, beta, l,
     *   mpi, nbet, nref)
        integer*4 nr, ndip, ns, mpi, nbet, nref
        integer*4 vecdip(ndip, 2), ncount(nr)
        double precision pi(mpi), beta(nbet), l, q1, q2, temp_q_1, temp_1,
     *  temp_q_2, temp_2, temp_3, temp_4
C        print*, 'check......'
         temp_4=0.d0
         do i=1, mpi
         temp_2=0.d0
           if(i.ne.nref) temp_2=beta(i)!if(i.ne.1) temp_2=beta(i-1)
           do j=i, mpi
            if(j.ne.nref) then 
              temp_3=exp(temp_2+beta(j))
            else
              temp_3=exp(temp_2)
            endif
C            if(j.ne.1) then 
C              temp_3=exp(temp_2+beta(j-1))
C            else
C              temp_3=exp(temp_2)
C            endif
            if(j.eq.i) then 
             temp_4=temp_4+temp_3*pi(j)*pi(j)
            else
             temp_4=temp_4+temp_3*2.d0*pi(i)*pi(j)
            endif
           end do
         end do
C         print*, 'temp_4= ', temp_4
        ns=nr/2
        l=0.d0
        do i=1, ns
           j1=2*(i-1)+1
           if( j1 .eq. 1) then 
            i1=1
           else
            i1=ncount(j1-1)+1
           endif
           i2=ncount(j1)
           temp_1=0.d0
          do k=i1, i2
            temp_q_1=0.d0
C            if(vecdip(k, 1).ne. 1) temp_q_1=beta(vecdip(k, 1)-1)
C            if(vecdip(k, 2).ne. 1) temp_q_1=temp_q_1+beta(vecdip(k, 2)
C     *      -1)
           if(vecdip(k, 1).ne. nref) temp_q_1=beta(vecdip(k, 1))
          if(vecdip(k, 2).ne. nref) temp_q_1=temp_q_1+beta(vecdip(k, 2))
            if(vecdip(k, 1).eq. vecdip(k, 2)) then 
               temp_1=temp_1+pi(vecdip(k ,1))*pi(vecdip(k, 1))*
     *         exp(temp_q_1)
            else
               temp_1=temp_1+2.d0*pi(vecdip(k ,1))*pi(vecdip(k, 2))*
     *         exp(temp_q_1)
            endif
          end do
          if(temp_1.ne.0.d0)  l=l+log(temp_1)
        end do             
        l=l-ns*log(temp_4)
        return
        end
C **************************************************************
C       This subroutine does MCMC simulations
C       First we update the latent variables D_{ij}
C       pi=vector of 16 components represents the frequency of the haplotypes
C       beta=vector of the 3 components represents the 
        subroutine gibbs(nr, nstr, ncount, vecdip, ndip, pi, beta, 
     *  rn, isd)
        integer*4  nr, nstr, ndip
        integer*4 vecdip(ndip, 2), d(nr), rmultinom, ncount(nr), 
     *  rn(nr, 2) 
        double precision prob(8), beta, pi(nstr, 16), s, temp_q_1, 
     *  temp_q
        integer*4 is1, is2, seed1, seed2, isd, i, j, i1, i2, k, k1
        common /unif_seeds/ is1, is2
        save /unif_seeds/
        seed1=7321
        seed2=isd
        CALL set_uniform(seed1, seed2)
        
        do i=1, nstr
          j=2*(i-1)+1
           do l=1, 8
             prob(l)=0.d0
           end do

           if( j .eq. 1) then 
            i1=1
           else
            i1=ncount(j-1)+1
           endif
          i2=ncount(j)
          s=0.d0
          k1=0
          do k=i1, i2
           k1=k1+1
            temp_q_1=0.d0
            if(vecdip(k, 1).eq. 2) temp_q_1=beta
            if(vecdip(k, 2).eq. 2) temp_q_1=temp_q_1+beta
            temp_q=exp(temp_q_1)
            if(vecdip(k, 1).eq.vecdip(k, 2)) then 
             prob(k1)=pi(i, vecdip(k, 1))*pi(i, vecdip(k, 1))*temp_q
            else 
             prob(k1)=2*pi(i, vecdip(k, 1))*pi(i, vecdip(k, 2))*temp_q
            endif 
            s=s+prob(k1)
          end do
C          if( i.eq. 1) print*, 'i= ', i, 'j= ', j, ' k1= ', k1,
C     *     'i2= ', i2, ' i1= ', i1
C     *    i2-i1+1
           s=1.d0/s
          do k2=1, i2-i1+1
            prob(k2)=s*prob(k2)
          end do
C          if( i .eq. 1) print*, 'prob= ', (prob(iu), iu=1, i2-i1+1)
          d(j)=rmultinom(8, prob)
          rn(j, 1)=vecdip(i1+d(j)-1, 1)
          rn(j, 2)=vecdip(i1+d(j)-1, 2)
          

C          if(i .eq. 1) print*, 'd(j)=', d(j), ' rn(j, 1)= ',
C     *   rn(j, 1), 'rn(j, 2)= ', rn(j, 2) 
          j=2*(i-1)+2
           do l=1, 8
             prob(l)=0.d0
           end do
          i1=ncount(j-1)+1
          i2=ncount(j)
      
          s=0.d0
          k1=0
          do k=i1, i2
            k1=k1+1
            temp_q_0=0.d0
           if(vecdip(k, 1) .eq. vecdip(k, 2)) then 
               prob(k1)=pi(i, vecdip(k, 1))*pi(i, vecdip(k, 1))
            else 
              prob(k1)=2*pi(i, vecdip(k, 1))*pi(i, vecdip(k, 2))
            endif
            s=s+prob(k1)
          end do
           s=1.d0/s
          do k2=1, i2-i1+1
            prob(k2)=s*prob(k2)
          end do
C          print*, 'i2-i1+1= ', i2-i1+1
          d(j)=rmultinom(8, prob)
          rn(j, 1)=vecdip(i1+d(j)-1, 1)
          rn(j, 2)=vecdip(i1+d(j)-1, 2)
                    
C          print*, 'i1= ', i1, ' i2= ', i2, ' d(j)=', d(j), rn(j, 1), 
C     *    rn(j, 2) 
        end do 
C        print*, 'd= ', (d(i), i=1, 10)
        return
        end
C ----------------------------------
          INTEGER function rmultinom(n, p)
          integer*4 n
          double precision p(n), uniform, s
          integer*4 is1, is2, i
          common /unif_seeds/ is1, is2
          save /unif_seeds/
            s=uniform()
            s=s-p(1)
            i=1
           do while(s.gt.0.d0)
            i=i+1
            s=s-p(i)
           end do 
           rmultinom=i
           return 
           end 
C-----------------------------------
          subroutine dirichletran(n, p, r, isd)
          integer*4 n
          double precision p(n), r(n)
          double precision x(n), gam_ran, s
          integer*4 is1, is2, seed1, seed2, isd
          common /unif_seeds/ is1, is2
          save /unif_seeds/
          seed1=100
          seed2=isd
          CALL set_uniform(seed1, seed2)
           s=0.d0
           do i=1, n
            x(i)=gam_ran(p(i), 1.d0)
            s=s+x(i)
           end do 
            s=1.d0/s
           do i=1, n
             r(i)=x(i)*s
           end do
           return 
           end 
C ------------------------------------
         subroutine piupdt(conf_ind, psi, nstr, nr, k, dip, isd, 
     *   new_psi)
         integer*4 nstr, k, isd, nr
         integer*4 conf_ind(nstr), dip(nr, 2)
         double precision psi(nstr, 16), new_psi(nstr, 16)
         double precision lk, r(16), pr(16), uniform
         integer*4 is1, is2, seed1, seed2
          common /unif_seeds/ is1, is2
          save /unif_seeds/
          seed1=100
          seed2=isd
          CALL set_uniform(seed1, seed2)
          do i=1, k
             do i1=1, 16
              pr(i1)=0.02d0
             end do
            isd=isd+i
            do j=1, nstr
              j1=2*(j-1)+2
              if(conf_ind(j).eq. i) then 
                pr(dip(j1, 1))=pr(dip(j1, 1))+1.d0
                pr(dip(j1, 2))=pr(dip(j1, 2))+1.d0 
              endif
            end do
            call dirichletran(16, pr, r, isd) 
               do iu=1, 16
                psi(i, iu)=r(iu)
                new_psi(i, iu)=psi(i, iu)
               end do
          end do
          return 
          end
C ------------------------------------
C        This subroutine updates the configuration indicators
         subroutine cnfindupdt(conf_ind, clst_name, clst_size, 
     *   new_conf_ind, new_clst_name, new_clst_size, psi, alpha,
     *    nstr, k, dip, nr, newk, isd, psi_distinct)
C        input configuration indicator, name of the cluster 1 through k, 
C         presuming that k is the total number of clusters, and the
C        size of each cluster
C        conf_ind, clst_names, clst_size, k
         integer*4 nstr, rmultinom, c_star, nr, k
         integer*4 dip(nr, 2), newk
         integer*4 conf_ind(nstr), clst_name(nstr), clst_size(nstr), 
     *   new_conf_ind(nstr), new_clst_name(nstr), new_clst_size(nstr)
         double precision prob(nstr), uniform, gam_ran
         double precision  new_psi(16), psi(nstr, 16), const,
     *   alpha, psi_distinct(nstr, 16), pr(16), r(16), diff_l
         integer*4 is1, is2, seed1, seed2, isd, old_ind
          common /unif_seeds/ is1, is2
          save /unif_seeds/
          seed1=100
          seed2=isd
          CALL set_uniform(seed1, seed2)
        
C          print*, 'HIIIIIIIIII'
          do iu=1, 16
          pr(iu)=0.02d0
          end do
          const=1.d0/(nstr+alpha-1)
         do i=1, nstr
         old_ind=conf_ind(i)
         isd=isd+i
         ij=2*(i-1)+2
         s=0.d0 
         if(dip(ij, 1).eq. dip(ij, 2)) then
          
           do j=1, k
            prob(j)=clst_size(j)*const*psi(j, dip(ij, 1))*
     *      psi(j, dip(ij, 2))
            s=s+prob(j)
           end do
            s=s-prob(conf_ind(i))
            prob(conf_ind(i))=(clst_size(conf_ind(i))-1)*const*
     *      psi(conf_ind(i), dip(ij, 1))*psi(conf_ind(i), dip(ij, 1))
            prob(k+1)=alpha*const*0.048d0
            s=s+prob(conf_ind(i))+prob(k+1)            
         else  
           do j=1, k
            prob(j)=clst_size(j)*const*2*psi(j, dip(ij, 1))*
     *      psi(j, dip(ij, 2))
            s=s+prob(j)
           end do
            s=s-prob(conf_ind(i))
            prob(conf_ind(i))=(clst_size(conf_ind(i))-1)*const*
     *      2*psi(conf_ind(i), dip(ij, 1))*psi(conf_ind(i), dip(ij, 1))
            prob(k+1)=alpha*const*0.000946d0
            s=s+prob(conf_ind(i))+prob(k+1) 
         endif
         s=1.d0/s
         do iu=1, k+1 
         prob(iu)=s*prob(iu)
         end do
            c_star=rmultinom(nstr, prob)
            if(clst_size(conf_ind(i)). gt. 1) then
              if(c_star.le. k) then
               clst_size(c_star)=clst_size(c_star)+1
               clst_size(conf_ind(i))= clst_size(conf_ind(i))-1
               conf_ind(i)=c_star               
              else
                call dirichletran(16, pr, r, isd)
                 clst_size(k+1)=1
                 clst_name(k+1)=k+1
                 clst_size(conf_ind(i))= clst_size(conf_ind(i))-1
                 conf_ind(i)=k+1
                 k=k+1
                 do j1=1, 16
                  psi(k+1,j1)=r(j1)
                 end do
              endif
            else
              if(c_star.le. k) then
               clst_size(c_star)=clst_size(c_star)+1
               clst_size(conf_ind(i))= clst_size(conf_ind(i))-1
               conf_ind(i)=c_star
                  do j=1, nstr
                   if(conf_ind(j).gt.old_ind) conf_ind(j)=conf_ind(j)-1
                  end do
                  do j=old_ind+1, nstr
                   clst_size(j-1)=clst_size(j)
                  end do
                   clst_size(k)=0 
                   clst_name(k)=0
                   k=k-1          
              else
                call dirichletran(16, pr, r, isd)
                 do j1=1, 16
                  psi(conf_ind(i), j1)=r(j1)
                 end do
              endif
            endif
C          if( i.eq.1) then 
C           print*, 'i=', i, 'old_ind= ', old_ind, ' c_star=', c_star
C           print*, 'clst_size'
C           print*, clst_size
C           print*, 'clst_name='
C           print*, clst_name
C           print*, 'conf_ind='
C           print*, conf_ind
C          endif
C          print*, 'i=', i
         end do
         do i=1, nstr
               new_conf_ind(i)=conf_ind(i)
               new_clst_name(i)=clst_name(i)
               new_clst_size(i)=clst_size(i)
           do i1=1, 16
              psi_distinct(i, i1)=psi(i, i1)               
           end do
         end do
         newk=k
         print*, 'newk=', newk
         return
         end

C --------------------------------------
          double precision function gam_ran(aa, bb)
          double precision aa, bb, d, c, w, u, v, x, z, yy, e,uniform
          logical acc
          integer*4 is1, is2
          common /unif_seeds/ is1, is2
          save /unif_seeds/
C         integer*4 l

          if (aa.lt.1.d0) then
c               { ********************* Johnk's algorithm***}
         x = 2.d0
         yy=0.d0
         do while((x+yy).gt.1.d0)
           u = uniform()
           v = uniform()
           x = exp(log(u)/aa)
           yy = exp(log(v)/(1.d0-aa))
           enddo
         e = -log(uniform())
         z = e*x/(x + yy)
         gam_ran = z * bb

        else if (aa.eq.1.d0) then
        gam_ran = -bb * log(uniform())
        else
c               { ************************* Best's algorithm***}
        d = aa-1.d0
        c = 3.d0*aa - (.75d0)
        acc = .false.
        do while (.not.acc)
          u = uniform()
          v = uniform()
          w = u*(1.d0-u)
          yy = sqrt(c/w)*(u-0.5d0)
          x = d+yy
          if (x.gt.0.d0) then
            z = 64.d0*w**3*v**2
            acc = (z.le.(1.d0-2.d0*yy**2/x))
            if (.not.acc) acc = (log(z).le.(2.d0 * (d*log(x/d)-yy) ))
            endif
          enddo
        gam_ran = x * bb
       endif
       return
       end
c --------------------------------------------------------------------
c{*************Uniform random number generator***********************}
c
c
        double precision function uniform()
c
c       Generate uniformly distributed random numbers using the 32-bit
c       generator from figure 3 of:
c       L'Ecuyer, P. Efficient and portable combined random number
c       generators, C.A.C.M., vol. 31, 742-749 & 774-?, June 1988.
c
c       The cycle length is claimed to be 2.30584E+18
c
c       Seeds can be set by calling the routine set_uniform
c
c       It is assumed that the Fortran compiler supports long variable
c       names, and integer*4.
c
        integer*4 z, k, is1, is2
        common /unif_seeds/ is1, is2
        save /unif_seeds/
c
        k = is1 / 53668
        is1 = 40014 * (is1 - k * 53668) - k * 12211
        if (is1 .lt. 0) then
           is1 = is1 + 2147483563
        endif
c
        k = is2 / 52774
        is2 = 40692 * (is2 - k * 52774) - k * 3791
        if (is2 .lt. 0) then
           is2 = is2 + 2147483399
        endif
c
        z = is1 - is2
        if (z .lt. 1) then
           z = z + 2147483562
        endif
c
        uniform = z / 2147483563.d0
        return
        end

C----------------------------------------------------------------
        SUBROUTINE set_uniform(iseed1, iseed2)
c
c       Set seeds for the uniform random number generator.
c
        INTEGER*4 is1, is2, iseed1, iseed2
        common /unif_seeds/ is1, is2
        save /unif_seeds/

        is1 = iseed1
        is2 = iseed2
        return
        end
   


c{*************** normal random number generator *****************}
 
        DOUBLE PRECISION FUNCTION ran_nor(m,v)
        DOUBLE PRECISION m, v, uniform, pi,z,huge
        INTEGER*4 is1, is2
        common /unif_seeds/ is1, is2
        save /unif_seeds/
 
        pi = 3.1415926536d0
       huge = 1e100
 1     z=dsqrt(-2.d0*dlog(uniform()))*dcos(2.d0*pi*uniform())
       if (dabs(z).ge.huge) goto 1
       ran_nor = m + sqrt(v)*z
       return
       end

 

C-------------------------------------------------------------------
C                                                                       
        DOUBLE PRECISION FUNCTION random(L)                               
C                                                                       
C       algorithm AS 183 APPL STATIST. (1982) Vol 31 n 2                
C       modified according to AS R58                                    
C      returns a pseudo-random number rectangularly distributed         
C      between 0 and 1                                                  
C                                                                       
C      IXR, IYR, AND IZR should be set to integer values between        
C      1 and 30000before first entry                                    
C                                                                       
C      Integer arithmetic up to 30323 required                          
C                                                                       
       IMPLICIT NONE                                                     
       INTEGER IXR,IYR,IZR,L                                            
       COMMON /RAND/ IXR,IYR,IZR                                        
       IXR = 171 * MOD(IXR,177) -2 *(IXR/177)                           
       IYR = 172 * MOD(IYR,176) -35*(IYR/176)                           
       IZR = 170 * MOD(IZR,178) -63*(IZR/178)                           
C                                                                       
       IF (IXR.LT.0) IXR=IXR+30269                                      
       IF (IYR.LT.0) IYR=IYR+30307                                      
       IF (IZR.LT.0) IZR=IZR+30323                                      
C                                                                       
       RANDOM= AMOD(FLOAT(IXR)/30269.0 +FLOAT(IYR) /30308.0 +           
     +              FLOAT(IZR)/30323.0,1.0)                             
       IF (RANDOM.GT.0.0)  RETURN                                       
       RANDOM=   DMOD(DBLE(FLOAT(IXR))/30269.0D0 +    
     +  DBLE(FLOAT(IYR))/30307.0D0+DBLE(FLOAT(IZR))/30323.0D0,1.0D0)    
       IF (RANDOM.GE.1) RANDOM=0.999999                                 
       RETURN                                                           
       END                                                              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         double precision function rbinary(p)

C       GENERATE BERNOULII RANDOM VARIABLE WITH PROBABILITY p  
         DOUBLE PRECISION p, uniform ,r
         INTEGER*4 is1, is2
         COMMON /unif_seeds/ is1, is2
         SAVE /unif_seeds/
        
         r=uniform()
         
        if(r.le.p) then
         rbinary=1.0d0
        else
         rbinary=0.0d0
        endif
        return
        end
cccccccccccccccccccccccccccccccccccccccccc
 
