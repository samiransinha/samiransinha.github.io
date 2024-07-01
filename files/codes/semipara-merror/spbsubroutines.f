
C --------------------
        subroutine detbeta(outbeta, ntot, m, nc, nbet, combdata, 
     *   beta, isd, diffmat, del, dmatbeta, sigma2beta)
        IMPLICIT NONE
        integer*4 ntot, m, nc, nbet
        double precision outbeta(nbet), beta(nbet), priormn, sigma2beta,
     *   combdata(ntot, nc), diffmat(nbet, nbet), del, newvar(nbet)
        double precision loglk, newloglk, oldloglk, temp, oldbeta
        double precision ran_nor, uniform, dmatbeta(ntot, nbet)
        integer*4 i, i1, j, k, is1, is2, seed1, seed2, isd
        common /unif_seeds/ is1, is2
        save /unif_seeds/
        seed1=10005
        seed2=isd      
        CALL set_uniform(seed1, seed2)
C        print*, 'check         -I'
        newvar(1)=1.d0/(del+1.d0/sigma2beta)
        newvar(2)=1.d0/(5.d0*del+1.d0/sigma2beta)
        newvar(nbet-1)=newvar(2)
        newvar(nbet)=newvar(1)
        do k=3, nbet-2 
         newvar(k)=1.d0/(6.d0*del+1.d0/sigma2beta)
        end do

        loglk=0.d0
        do i=(m+1), ntot
          temp=0.d0
         do j= 1, nbet        
          temp=temp+dmatbeta(i, j)*beta(j)
         end do 
         loglk=loglk+combdata(i, 1)*temp-log(1.d0+exp(temp))
        end do
         oldloglk=loglk
C   For \beta_1
        do k=1, nbet
         oldbeta=beta(k)
         priormn=0.d0
         do i1=1, nbet
          priormn=priormn+beta(i1)*diffmat(k, i1)
         end do
         beta(k)=ran_nor(del*priormn*newvar(k), newvar(k))
         loglk=0.d0
         do i=(m+1), ntot
          temp=0.d0
          do j= 1, nbet        
           temp=temp+dmatbeta(i, j)*beta(j)
          end do 
          loglk=loglk+combdata(i, 1)*temp-log(1.d0+exp(temp))
         end do
         newloglk=loglk
         if(uniform().lt.exp(newloglk-oldloglk)) then 
          oldloglk=newloglk 
         else
          beta(k)=oldbeta
         endif
        end do
        do k=1, nbet
         outbeta(k)=beta(k)
        end do
C
         return 
         end
C -----------------------------------------------------------
C
        subroutine detalpha(outalpha, ntot, nc, nalph, combdata,  
     *   alpha, isd, diffmat, delalph, sigma2q, dmatalpha, sigma2alpha)
        IMPLICIT NONE
        integer*4 ntot, nc, nalph
        double precision outalpha(nalph), alpha(nalph), priormn, 
     *   prvar(nalph), combdata(ntot, nc), diffmat(nalph, nalph)
        double precision temp1, temp2, temp3, mn, var, delalph
        double precision sigma2q, uniform, sigma2alpha,   
     *   CDFNOR, PPND7, dmatalpha(ntot, nalph), sqrt_var, constant
        integer*4 i, i1, j, k, k1, is1, is2, seed1, seed2, isd, IFAULT
        common /unif_seeds/ is1, is2
        save /unif_seeds/
        seed1=10005
        seed2=isd      
        CALL set_uniform(seed1, seed2)
        sigma2alpha=10.d0**8
        do k=1, nalph
         outalpha(k)=alpha(k)
        end do
C For \alpha_1
        priormn=0.d0
        do i1=1, nalph
         priormn=priormn+alpha(i1)*diffmat(1, i1)
        end do
        temp1=0.d0
        temp2=0.d0
        do i=1, ntot
         temp2=temp2+dmatalpha(i, 1)*dmatalpha(i, 1)
         temp3=-dmatalpha(i, 1)*alpha(1)
         do j= 1, nalph        
          temp3=temp3+dmatalpha(i, j)*alpha(j)
         end do      
         temp1=temp1+(combdata(i, 2)-temp3)*dmatalpha(i, 1) 
        end do
        var=1.d0/(delalph+temp2/sigma2q+1.d0/sigma2alpha)
        mn=var*(delalph*priormn+temp1/sigma2q)        
        sqrt_var=sqrt(var)
        constant=uniform()*CDFNOR((alpha(2)-mn)/sqrt_var)
        alpha(1)=mn+sqrt_var* PPND7(constant)
        if( abs(alpha(1)-alpha(2)).gt.20) alpha(1)=alpha(2)-2.d0
C        print*, PPND7(constant)
C For \alpha_2        
        priormn=0.d0
        do i1=1, nalph
         priormn=priormn+alpha(i1)*diffmat(2, i1)
        end do
        temp1=0.d0
        temp2=0.d0
        do i=1, ntot
         temp2=temp2+dmatalpha(i, 2)*dmatalpha(i, 2)
         temp3=-dmatalpha(i, 2)*alpha(2)
         do j= 1, nalph        
          temp3=temp3+dmatalpha(i, j)*alpha(j)
         end do      
         temp1=temp1+(combdata(i, 2)-temp3)*dmatalpha(i, 2) 
        end do
        var=1.d0/(5*delalph+temp2/sigma2q+1.d0/sigma2alpha)
        mn=var*(delalph*priormn+temp1/sigma2q) 
        sqrt_var=sqrt(var)
        constant=CDFNOR((alpha(1)-mn)/sqrt_var)+ uniform()*(
     *  CDFNOR((alpha(3)-mn)/sqrt_var)-CDFNOR((alpha(1)-mn)/sqrt_var))
        alpha(2)=mn+sqrt_var*PPND7(constant)
C For \alpha_j j=3, \cdots, nalph-2
        do k=3, (nalph-2)
         priormn=0.d0
         do i1=1, nalph
          priormn=priormn+alpha(i1)*diffmat(k, i1)
         end do
         temp1=0.d0
         temp2=0.d0
         do i=1, ntot
          temp2=temp2+dmatalpha(i, k)*dmatalpha(i, k)
          temp3=-dmatalpha(i, k)*alpha(k)
          do j= 1, nalph        
           temp3=temp3+dmatalpha(i, j)*alpha(j)
          end do      
          temp1=temp1+(combdata(i, 2)-temp3)*dmatalpha(i, k) 
         end do
         var=1.d0/(6*delalph+temp2/sigma2q+1.d0/sigma2alpha)
         mn=var*(delalph*priormn+temp1/sigma2q) 
         sqrt_var=sqrt(var)
         constant=CDFNOR((alpha(k-1)-mn)/sqrt_var)+ uniform()*(
     *   CDFNOR((alpha(k+1)-mn)/sqrt_var)-CDFNOR((alpha(k-1)-mn)/
     *   sqrt_var))
         alpha(k)=mn+sqrt_var*PPND7(constant)
        end do
C For \alpha_{nalph-1} 
        priormn=0.d0
        do i1=1, nalph
         priormn=priormn+alpha(i1)*diffmat((nalph-1), i1)
        end do
        temp1=0.d0
        temp2=0.d0
        do i=1, ntot
         temp2=temp2+dmatalpha(i, (nalph-1))*dmatalpha(i, (nalph-1))
         temp3=-dmatalpha(i, (nalph-1))*alpha(nalph-1)
         do j= 1, nalph        
          temp3=temp3+dmatalpha(i, j)*alpha(j)
         end do      
         temp1=temp1+(combdata(i, 2)-temp3)*dmatalpha(i, (nalph-1)) 
        end do
        var=1.d0/(5*delalph+temp2/sigma2q+1.d0/sigma2alpha)
        mn=var*(delalph*priormn+temp1/sigma2q) 
        sqrt_var=sqrt(var)
        constant=CDFNOR((alpha(nalph-2)-mn)/sqrt_var)+ uniform()*(
     *   CDFNOR((alpha(nalph)-mn)/sqrt_var)-CDFNOR((alpha(nalph-2)-mn)/
     *   sqrt_var))
         alpha(nalph-1)=mn+sqrt_var*PPND7(constant)
C For \alpha_{nalph}
        priormn=0.d0
        do i1=1, nalph
         priormn=priormn+alpha(i1)*diffmat(nalph, i1)
        end do
        temp1=0.d0
        temp2=0.d0
        do i=1, ntot
         temp2=temp2+dmatalpha(i, nalph)*dmatalpha(i, nalph)
         temp3=-dmatalpha(i, nalph)*alpha(nalph)
         do j= 1, nalph        
          temp3=temp3+dmatalpha(i, j)*alpha(j)
         end do      
         temp1=temp1+(combdata(i, 2)-temp3)*dmatalpha(i, nalph) 
        end do
        var=1.d0/(delalph+temp2/sigma2q+1.d0/sigma2alpha)
        mn=var*(delalph*priormn+temp1/sigma2q) 
        sqrt_var=sqrt(var)
        constant=CDFNOR((alpha(nalph-1)-mn)/sqrt_var)+ uniform()*(1.d0-
     *  CDFNOR((alpha(nalph-1)-mn)/sqrt_var))
        alpha(nalph)=mn+sqrt_var*PPND7(constant)
        do k=1, nalph 
         outalpha(k)=alpha(k)
        end do
        return 
        end

C -----------------------------------------------------------
C       Resampling X
        subroutine samplex(m, ntot, nc, nalph, nbet, alpha, 
     *  beta, propx, x, combdata,  sigma2m, sigma2q, 
     *  outx, isd, temp1, temp2, dmatalpha, dmatbeta, phi, indicator, 
     *  propmean, propvar)
        IMPLICIT NONE
        integer*4 m, ntot, nc, nalph, nbet
        integer*4 indicator(ntot)
        double precision uniform, x(ntot), propx(ntot), temp1(ntot, 
     *  nalph), temp2(ntot, nbet), combdata(ntot, nc), dmatalpha(ntot,
     *  nalph), dmatbeta(ntot, nbet), g, qnty1, qnty2, qnty3, qnty4, 
     *  CDFNOR, lambdaq, lambdam, outx(ntot), alpha(nalph), 
     *  beta(nbet), sigma2m, sigma2q, phi(2, ntot), propmean(ntot), 
     *  propvar(ntot), var,  mean
        integer*4 i, is1, is2, seed1, seed2, isd, j
        common /unif_seeds/ is1, is2
        save /unif_seeds/
        seed1=10005
        seed2=isd      
        CALL set_uniform(seed1, seed2)
        
        do i=1, m
         qnty1=0.d0
         qnty2=0.d0
         do j=1, nalph
          qnty1=qnty1+temp1(i, j)*alpha(j)
          qnty2=qnty2+dmatalpha(i, j)*alpha(j)
         end do
          var=1/(1/phi(2, indicator(i)) +2/sigma2m)
          mean=var*(phi(1, indicator(i))/phi(2, indicator(i))+
     *    (combdata(i, 3)+combdata(i, 4))/sigma2m)
         g=exp((-(combdata(i, 2)-qnty1)**2+(combdata(i, 2)-qnty2)**2
     *    )/(2*sigma2q)-
     *    (propx(i)-mean)**2/(2*var)+(x(i)-mean)**2/(2*var)-
     *   (x(i)-propmean(i))**2/(2*propvar(i))+
     *    (propx(i)-propmean(i))**2/(2*propvar(i))
     *    )
         if(uniform().lt.g) x(i)=propx(i)
        end do
        do i=(m+1), ntot
         qnty1=0.d0
         qnty2=0.d0
         do j=1, nalph
          qnty1=qnty1+temp1(i, j)*alpha(j)
          qnty2=qnty2+dmatalpha(i, j)*alpha(j)
         end do
         qnty3=0.d0
         qnty4=0.d0
         do j=1, nbet
          qnty3=qnty3+temp2(i, j)*beta(j)
          qnty4=qnty4+dmatbeta(i, j)*beta(j)
         end do
         g=exp((-(combdata(i, 2)-qnty1)**2+(combdata(i, 2)-qnty2)**2
     *    )/(2*sigma2q) +(qnty3-qnty4)*combdata(i, 1) +
     *   log(1+exp(qnty4))-log(1+exp(qnty3)) -
     *   (propx(i)-phi(1, indicator(i)))**2/(2*phi(2, indicator(i)))+
     *   (x(i)-phi(1, indicator(i)))**2/(2*phi(2, indicator(i)))-
     *   (x(i)-propmean(i))**2/(2*propvar(i))+
     *   (propx(i)-propmean(i))**2/(2*propvar(i)))
         if(uniform().lt.g) x(i)=propx(i)
        end do
        do i=1, ntot
         outx(i)=x(i)
        end do
        return
        end

C ----------------------------------------------------------
C this subroutine is used to draw phi's for DPP
          subroutine dpp(ntot, ncount, phi, a0, b0, alpha0, tau, 
     *    ncountout, phiout, indicatorout, ndistinctout, isd, indicator,
     *    ndistinct, x, const1)
          IMPLICIT NONE
          integer*4 ntot, ndistinct, ndistinctout
          integer*4 ncount(ntot), ncountout(ntot), indicator(ntot)
          double precision phi(2, ntot), phiout(2, ntot), alpha0
          double precision uniform, ran_nor, x(ntot), pi, q0, b, 
     *    temp1(ntot), r, prob(ntot), var, scale, shape,  sum, const1,
     *    tau, oldphi(2), mn, gam_ran, tempsum, a0, b0, temp2, temp3
          integer*4 i, j, j1, is1, is2, seed1, seed2, isd
          integer*4 indicatorout(ntot), i1, j2
          common /unif_seeds/ is1, is2
          save /unif_seeds/
          seed1=10005
          seed2=isd      
          CALL set_uniform(seed1, seed2)
          pi=3.141593d0          
C Step 1, updating configuration indicator function
          do i=1, ntot
           if(ncount(indicator(i)).gt.1) then
            ncount(indicator(i))=ncount(indicator(i))-1
            q0=const1/((1+b0*x(i)**2/(2*(1+tau))
     *      )**(a0+0.5))
            b=alpha0*q0
            do j=1, ndistinct
             temp1(j)=ncount(j)*exp(-(x(i)-phi(1, j))**2/(2*phi(2, j)))/
     *       sqrt(2*pi*phi(2, j))
             b=b+temp1(j)
            end do
            b=1.d0/b
            temp1(ndistinct+1)=q0*alpha0
            r=uniform()
            sum=0.d0
            j1=0
            do while(r.gt.sum)
             j1=j1+1
             prob(j1)=b*temp1(j1)
             sum=sum+prob(j1)
            end do
            if(j1.lt.(ndistinct+1)) then 
             indicator(i)=j1
             ncount(j1)=ncount(j1)+1
            else
             indicator(i)=ndistinct+1
             ncount(ndistinct+1)=1
             ndistinct=ndistinct+1 
             scale=(1/b0+(x(i)**2)/(2*(1+tau)))
             scale=1/scale
             phi(2, indicator(i))=1.d0/gam_ran((a0+0.5d0),  scale)
             var=phi(2, indicator(i))*tau/(1+tau)
             mn=x(i)*tau/(1+tau)
             phi(1, indicator(i))=ran_nor(mn, var)              
            endif
           else 
C Redefining the clusters
            oldphi(1)=phi(1, indicator(i))
            oldphi(2)=phi(2, indicator(i))
            do j=indicator(i), ndistinct
             ncount(j)=ncount(j+1)
             phi(1, j)=phi(1, j+1)
             phi(2, j)=phi(2, j+1) 
            end do       
C Redefining the configuration indicators     
            do j=1, ntot
             if(indicator(j).gt.indicator(i)) indicator(j)=
     *       indicator(j)-1
            end do
            q0=const1/((1+b0*x(i)**2/(2*(1+tau))
     *      )**(a0+0.5))
            b=alpha0*q0
            do j=1, (ndistinct-1)
             temp1(j)=ncount(j)*exp(-(x(i)-phi(1, j))**2/(2*phi(2, j)))/
     *       sqrt(2*pi*phi(2, j))
             b=b+temp1(j)
            end do
            b=1/b
            temp1(ndistinct)=q0*alpha0
            r=uniform()
            sum=0.d0
            j1=0
            do while(r.gt.sum)
             j1=j1+1
             prob(j1)=b*temp1(j1)
             sum=sum+prob(j1)
            end do
            if(j1.lt.(ndistinct)) then 
             indicator(i)=j1
             ncount(j1)=ncount(j1)+1
             ndistinct=ndistinct-1
            else 
             indicator(i)=ndistinct
             ncount(ndistinct)=1
             scale=(1/b0+(x(i)**2)/(2*(1+tau)))
             scale=1/scale
             phi(2, ndistinct)=1.d0/gam_ran((a0+0.5d0), scale)
             var= phi(2, ndistinct)*tau/(1+tau)
             mn=x(i)*tau/(1+tau)
             phi(1, ndistinct)=ran_nor(mn, var) 
            endif
           endif
          end do
C End of Step 1
C Step 2, resampling  the distinct elements.
           do j=1, ndistinct 
            temp2=0.d0
            temp3=0.d0
            do i=1, ntot
             if(indicator(i).eq.j) then
              temp2=temp2+(x(i)-phi(1, j))*(x(i)-phi(1, j))
              temp3=temp3+x(i)
             endif
            end do
            shape=0.5d0*ncount(j)+a0+0.5d0
            scale=1.d0/(1.d0/b0+(0.5d0*phi(1, j)*phi(1, j)/tau)+
     *      0.5d0*temp2)
            phi(2, j)=1.d0/gam_ran(shape, scale)
            var=phi(2, j)/(ncount(j)+1.d0/tau)
            mn=var*temp3/phi(2, j)
            phi(1, j)=ran_nor(mn, var)
           end do
C end of Step 2
          ndistinctout=ndistinct
          do i=1, ndistinct
           phiout(1, i)=phi(1, i)
           phiout(2, i)=phi(2, i)
           ncountout(i)=ncount(i)
          end do
          do i=1, ntot
           indicatorout(i)=indicator(i)
          end do
         return
         end      
C ------------------------------------------------------------
          subroutine detrminealpha0(mnofk, ntot, alpha0, alpha0out)          
          IMPLICIT NONE
          integer*4 ntot, i, ncnt
          double precision mnofk, alpha0out, sum, eps, alpha, 
     *    alpha1, alpha2, alpha0
          alpha1=0.001d0 
          alpha2=10000.0d0
          eps=2.d0
          ncnt=0
          do while((eps.gt.0.01d0). and.(ncnt.lt.100))
           ncnt=ncnt+1
           alpha=(alpha1+alpha2)/2
           sum=-mnofk
           do i=1, ntot
            sum=sum+alpha/(alpha+i-1)
           end do
           if(sum.gt.0) then 
            alpha2=alpha
           else
            alpha1=alpha
           endif
           eps=abs(sum)
           
          end do  
          if(ncnt.lt.90) then 
           alpha0out=alpha
          else
           alpha0out=alpha0
          endif

          return
          end
C -----------------------------------------------------------
          subroutine gamoutput(isd, u, nr)
          IMPLICIT NONE
          integer*4 isd, nr
          double precision u(nr), gam_ran
          integer*4 i, is1, is2, seed1, seed2
          common /unif_seeds/ is1, is2
          save /unif_seeds/
          seed1=10005
          seed2=isd      
          CALL set_uniform(seed1, seed2)
          do i=1, nr
           u(i)=gam_ran(2.d0, 2.d0)
          end do
          return
          end 
          
C------------------------------------------------------------
          double precision function gam_ran(aa, bb)
          IMPLICIT NONE
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
        IMPLICIT NONE
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
        IMPLICIT NONE
        INTEGER*4 is1, is2, iseed1, iseed2
        common /unif_seeds/ is1, is2
        save /unif_seeds/

        is1 = iseed1
        is2 = iseed2
        return
        end
   


c{*************** normal random number generator *****************}
 
        DOUBLE PRECISION FUNCTION ran_nor(m,v)
        IMPLICIT NONE
        DOUBLE PRECISION m, v, uniform, pi,z,huge
        INTEGER*4 is1, is2
        common /unif_seeds/ is1, is2
        save /unif_seeds/
 
        pi = 3.1415926536d0
       huge = 1e35
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
C
      DOUBLE PRECISION FUNCTION DGAMLN (X)
C
C-----------------------------------------------------------------------
C   DGAMLN   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
C            MARYLAND  20899
C
C   FOR: COMPUTING THE DOUBLE PRECISION LOG OF THE GAMMA FUNCTION WITH
C        SINGLE PRECISION PARAMETER X>0.  THE MAXIMUM TRUNCATION ERROR
C        IN THE INFINITE SERIES (SEE REFERENCE 1) IS DETERMINED BY THE
C        CONSTANT XMIN.  WHEN X<XMIN A RECURRENCE RELATION IS USED IN 
C        ORDER TO ACHIEVE THE REQUIRED ABSOLUTE ACCURACY.  THE TABLE
C        BELOW GIVES THE MINIMUM VALUE OF X WHICH YIELDS THE CORRES-
C        PONDING ABSOLUTE ACCURACY IN DGAMLN(X) ASSUMING THE MACHINE
C        CARRIES ENOUGH DIGITS WHEN THOSE TO THE LEFT OF THE DECIMAL
C        ARE CONSIDERED (SEE REFERENCE 2 FOR FURTHER DISCUSSION).  IF 
C        THE LATTER CONDITION IS NOT MET, AN ERROR MESSAGE IS PRINTED.
C
C        THE CYBER 180/855 AT NBS CARRIES ABOUT 15 DIGITS IN SINGLE
C        PRECISION, THEREFORE THE PRE-SET VALUE OF ABSACC IS 10**(-15)
C        AND THE CORRESPONDING VALUE OF XMIN IS 6.894.  ON A DIFFERENT
C        MACHINE THESE CONSTANTS SHOULD BE CHANGED ACCORDINGLY.
C
C           XMIN    ACCURACY    XMIN    ACCURACY    XMIN    ACCURACY
C          ------   --------   ------   --------   ------   --------
C           1.357     1E-3      4.592     1E-12    15.539     1E-21
C           2.037     1E-6      6.894     1E-15    23.330     1E-24
C           3.059     1E-9     10.351     1E-18    35.025     1E-27
C
C   NOTE: THIS IS EXACTLY THE SAME SOFTWARE AS SUBROUTINE DGAMLN
C
C   SUBPROGRAMS CALLED: -NONE-
C
C   CURRENT VERSION COMPLETED MAY 1, 1989
C
C   REFERENCES: 
C
C   1) ABRAMOWITZ, MILTON AND STEGUN, IRENE, 'HANDBOOK OF MATHEMATICAL
C      FUNCTIONS', NBS APPLIED MATHEMATICS SERIES 55, NOV. 1970,
C      EQ. 6.1.40, P 257.
C
C   2) REEVE, CHARLES P., 'ACCURATE COMPUTATION OF THE LOG OF THE GAMMA
C      FUNCTION', STATISTICAL ENGINEERING DIVISION NOTE 86-1, OCTOBER 
C      1986.
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION ABSACC,B1,B2,B3,B4,B5,B6,B7,
     * B8,C,DX,Q,R,XMIN,XN, X
      DATA XMIN,ABSACC / 6.894D0,1.0D-15 /
      DATA C / 0.918938533204672741780329736D0 /
      DATA B1 / 0.833333333333333333333333333D-1 /
      DATA B2 / -0.277777777777777777777777778D-2 /
      DATA B3 / 0.793650793650793650793650794D-3 /
      DATA B4 / -0.595238095238095238095238095D-3 /
      DATA B5 / 0.841750841750841750841750842D-3 /
      DATA B6 / -0.191752691752691752691752692D-2 /
      DATA B7 / 0.641025641025641025641025641D-2 /
      DATA B8 / -0.295506535947712418300653595D-1 /
C
C--- TERMINATE EXECUTION IF X<=0.0
C
      IF (X.LE.0.0) STOP  '*** X<=0.0 IN FUNCTION DGAMLN ***'
      DX = DBLE(X)
      N = MAX0(0,INT(XMIN-DX+1.0D0))
      XN = DX+DBLE(N)
      R = 1.0D0/XN
      Q = R*R
      DGAMLN = R*(B1+Q*(B2+Q*(B3+Q*(B4+Q*(B5+Q*(B6+Q*(B7+Q*B8)))))))+C+
     *   (XN-0.5D0)*DLOG(XN)-XN
C
C--- USE RECURRENCE RELATION WHEN N>0 (X<XMIN)
C
      IF (N.GT.0) THEN
         Q = 1.0D0
         DO 20 I  =  0, N-1
            Q = Q*(DX+DBLE(I))
   20    CONTINUE
         DGAMLN=DGAMLN-DLOG(Q)
      ENDIF
C
C--- PRINT WARNING IF ABSOLUTE ACCURACY HAS NOT BEEN ATTAINED
C
      IF (DGAMLN+ABSACC.EQ.DGAMLN) THEN 
CCCCC    PRINT *,' ********* WARNING FROM FUNCTION DGAMLN *********'
CCCCC    PRINT *,' REQUIRED ABSOLUTE ACCURACY NOT ATTAINED FOR X = ',X
      ENDIF
      RETURN
      END 
            SUBROUTINE CDFGAM (X,ALPHA,EPS,IFLAG,CDFX)
C
C-----------------------------------------------------------------------
C   CDFGAM   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
C            MD  20899
C
C   FOR: COMPUTING THE CUMULATIVE DISTRIBUTION FUNCTION OF THE GAMMA
C        DISTRIBUTION (ALSO KNOWN AS THE INCOMPLETE GAMMA RATIO) TO A 
C        SPECIFIED ACCURACY (TRUNCATION ERROR IN THE INFINITE SERIES).
C        THE ALGORITHM, DESCRIBED IN REFERENCE 2, IS A MODIFICATION OF
C        THE ALGORITHM OF REFERENCE 1.  THREE FEATURES HAVE BEEN ADDED: 
C
C        1) A PRECISE METHOD OF MEETING THE TRUNCATION ACCURACY,
C        2) COMPUTATION OF THE UPPER TAIL AREA BY DECREMENTING ALPHA
C           WHEN THAT METHOD IS MORE EFFICIENT, AND
C        3) A CONSTANT UFLO >= THE UNDERFLOW LIMIT ON THE COMPUTER.
C
C   SUBPROGRAMS CALLED: DGAMLN (LOG OF GAMMA FUNCTION)
C
C   CURRENT VERSION COMPLETED OCTOBER 29, 1986
C
C   REFERENCES: 
C
C   1) LAU, CHI-LEUNG, 'A SIMPLE SERIES FOR THE INCOMPLETE GAMMA
C      INTEGRAL', ALGORITHM AS 147, APPLIED STATISTICS, VOL. 29,
C      NO. 1, 1980, PP. 113-114.
C
C   2) REEVE, CHARLES P., 'AN ALGORITHM FOR COMPUTING THE GAMMA C.D.F.
C      TO A SPECIFIED ACCURACY', STATISTICAL ENGINEERING DIVISION
C      NOTE 86-2, OCTOBER 1986.
C-----------------------------------------------------------------------
C   DEFINITION OF PASSED PARAMETERS: 
C
C       * X = VALUE AT WHICH THE C.D.F IS TO BE COMPUTED (REAL)
C
C   * ALPHA = PARAMETER OF THE GAMMA FUNCTION (>0) (REAL)
C
C     * EPS = THE DESIRED ABSOLUTE ACCURACY OF THE C.D.F (>0) (REAL)
C
C     IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION: 
C             0 -> NO ERRORS DETECTED
C             1 -> EITHER ALPHA OR EPS IS <= UFLO 
C             2 -> NUMBER OF TERMS EVALUATED IN THE INFINITE SERIES
C                  EXCEEDS IMAX.
C
C      CDFX = THE C.D.F. EVALUATED AT X (REAL)
C
C   * INDICATES PARAMETERS REQUIRING INPUT VALUES 
C-----------------------------------------------------------------------
C
      
     
      DOUBLE PRECISION DX, DGAMLN, X, ALPHA, EPS, CDFX
      DOUBLE PRECISION UFLO, PDFL, P, U, ETA, BL
      INTEGER*4 IFLAG, K
       LOGICAL LL
      DATA IMAX,UFLO / 5000,1.0E-35 /
      CDFX = 0.0
C
C--- CHECK FOR VALIDITY OF ARGUMENTS ALPHA AND EPS
C
      IF (ALPHA.LE.UFLO.OR.EPS.LE.UFLO) THEN
         IFLAG = 1
         RETURN
      ENDIF
      IFLAG = 0
C
C--- CHECK FOR SPECIAL CASE OF X
C
      IF (X.LE.0.0) RETURN
C
C--- EVALUATE THE GAMMA P.D.F. AND CHECK FOR UNDERFLOW
C
      DX = DBLE(X)
      PDFL = SNGL(DBLE(ALPHA-1.0)*DLOG(DX)-DX-DGAMLN(ALPHA))
      IF (PDFL.LT.LOG(UFLO)) THEN
         IF (X.GE.ALPHA) CDFX = 1.0
      ELSE
         P = ALPHA
         U = EXP(PDFL)
C
C--- DETERMINE WHETHER TO INCREMENT OR DECREMENT ALPHA (A.K.A. P)
C
         LL = .TRUE.
         IF (X.GE.P) THEN
            K = INT(P)
            IF (P.LE.REAL(K)) K = K-1
            ETA = P-REAL(K)
            BL = SNGL(DBLE(ETA-1.0)*DLOG(DX)-DX-DGAMLN(ETA))
            LL = BL.GT.LOG(EPS)
         ENDIF
         EPSX=EPS/X 
         IF (LL) THEN
C
C--- INCREMENT P
C
            DO 10 I = 0, IMAX 
               IF (U.LE.EPSX*(P-X)) RETURN
               U = X*U/P
               CDFX = CDFX+U
               P = P+1.0
   10       CONTINUE
            IFLAG = 2
         ELSE
C
C--- DECREMENT P
C
            DO 20 J = 1, K
               P = P-1.0
               IF (U.LE.EPSX*(X-P)) GO TO 30
               CDFX = CDFX+U
               U = P*U/X
   20       CONTINUE
   30       CDFX = 1.0-CDFX
         ENDIF
      ENDIF
      RETURN
      END 
            DOUBLE PRECISION FUNCTION CDFNOR (Z)
C
C-----------------------------------------------------------------------
C   CDFNOR   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C            DIVISION, NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY,
C            GAITHERSBURG, MARYLAND  20899
C
C   FOR: COMPUTING THE CUMULATIVE DISTRIBUTION FUNCTION OF THE STANDARD
C        NORMAL DISTRIBUTION TO A SPECIFIED ACCURACY.  THE C.D.F. OF
C        THE GAMMA DISTRIBUTION IS FIRST COMPUTED, THEN TRANSFORMED TO
C        THE C.D.F. OF THE NORMAL DISTRIBUTION.
C
C   NOTE: TIMING TESTS ON THE CYBER 180/855 COMPUTER AT N.I.S.T.
C         INDICATE THAT THE AVERAGE CPU TIME FOR COMPUTING ONE C.D.F. 
C         IS ABOUT 0.0004 SECOND.
C
C   SUBPROGRAMS CALLED:  (C.D.F. OF THE GAMMA DISTRIBUTION)
C
C   CURRENT VERSION COMPLETED APRIL 7, 1989
C-----------------------------------------------------------------------
C   DEFINITION OF PASSED PARAMETERS: 
C
C     * Z = THE VALUE FOR WHICH THE NORMAL C.D.F. IS TO BE COMPUTED
C           [REAL]
C
C   * EPS = THE ABSOLUTE ACCURACY REQUIREMENT FOR THE C.D.F. [REAL]
C
C   IFLAG = ERROR INDICATOR ON OUTPUT [INTEGER]   INTERPRETATION: 
C           0 -> NO ERRORS DETECTED.
C           1 -> ERROR FLAG FROM CDFGAM 
C           2 -> ERROR FLAG FROM CDFGAM 
C
C    CDFZ = THE C.D.F. OF THE STANDARD NORMAL DISTRIBUTION EVALUATED
C           AT X [REAL]
C
C   * INDICATES PARAMETERS REQUIRING INPUT VALUES 
C-----------------------------------------------------------------------
C
       IMPLICIT NONE
       DOUBLE PRECISION Z, EPS, CDFZ, DEL, X, CDFX
       INTEGER*4 IFLAG
       EPS=0.0000001d0
       DEL = 2.0*EPS 
       IF (Z.EQ.0.0) THEN
        CDFZ = 0.50 
       ELSE
        X = 0.5*Z*Z
        CALL CDFGAM (X,0.5d0,DEL,IFLAG,CDFX)
        IF (IFLAG.NE.0) RETURN
        IF (Z.GT.0.0) THEN
         CDFZ = 0.5+0.5*CDFX
        ELSE
         CDFZ = 0.5-0.5*CDFX
        ENDIF
       ENDIF
        CDFNOR=CDFZ
      RETURN
C
      END 

C---------------------------------------------
	DOUBLE PRECISION FUNCTION PPND7 (P)
C
C	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-
C	484.
C
C	Produces the normal deviate Z corresponding to a given lower
C	tail area of P; Z is accurate to about 1 part in 10**7.
C
C	The hash sums below are the sums of the mantissas of the
C	coefficients.   They are included for use in checking
C	transcription.
C
	DOUBLE PRECISION ZERO, ONE, HALF, SPLIT1, SPLIT2, CONST1, 
     *  CONST2, A0, A1,
     *	A2, A3, B1, B2, B3, C0, C1, C2, C3, D1, D2, E0, E1, E2,
     *	E3, F1, F2, P, Q, R
	PARAMETER (ZERO = 0.0, ONE = 1.0, HALF = 0.5,
     *		SPLIT1 = 0.425, SPLIT2 = 5.0,
     *		CONST1 = 0.180625, CONST2 = 1.6)
        INTEGER IFAULT
C
C	Coefficients for P close to 0.5
C
	PARAMETER (A0 = 3.38713 27179E+00, A1 = 5.04342 71938E+01,
     *		   A2 = 1.59291 13202E+02, A3 = 5.91093 74720E+01,
     *		   B1 = 1.78951 69469E+01, B2 = 7.87577 57664E+01,
     *		   B3 = 6.71875 63600E+01)
C	HASH SUM AB    32.31845 77772
C
C	Coefficients for P not close to 0, 0.5 or 1.
C
	PARAMETER (C0 = 1.42343 72777E+00, C1 = 2.75681 53900E+00,
     *		   C2 = 1.30672 84816E+00, C3 = 1.70238 21103E-01,
     *		   D1 = 7.37001 64250E-01, D2 = 1.20211 32975E-01)
C	HASH SUM CD    15.76149 29821
C
C	Coefficients for P near 0 or 1.
C
	PARAMETER (E0 = 6.65790 51150E+00, E1 = 3.08122 63860E+00,
     *		   E2 = 4.28682 94337E-01, E3 = 1.73372 03997E-02,
     *		   F1 = 2.41978 94225E-01, F2 = 1.22582 02635E-02)
C	HASH SUM EF    19.40529 10204
C
	IFAULT = 0
	Q = P - HALF
	IF (ABS(Q) .LE. SPLIT1) THEN
	  R = CONST1 - Q * Q
	  PPND7 = Q * (((A3 * R + A2) * R + A1) * R + A0) /
     *		      (((B3 * R + B2) * R + B1) * R + ONE)
	  RETURN
	ELSE
	  IF (Q .LT. ZERO) THEN
	    R = P
	  ELSE
	    R = ONE - P
	  END IF
	  IF (R .LE. ZERO) THEN
	    IFAULT = 1
	    PPND7 = ZERO
	    RETURN
	  END IF
	  R = SQRT(-LOG(R))
	  IF (R .LE. SPLIT2) THEN
	    R = R - CONST2
	    PPND7 = (((C3 * R + C2) * R + C1) * R + C0) /
     *		     ((D2 * R + D1) * R + ONE)
	  ELSE
	    R = R - SPLIT2
	    PPND7 = (((E3 * R + E2) * R + E1) * R + E0) /
     *		     ((F2 * R + F1) * R + ONE)
	  END IF
	  IF (Q .LT. ZERO) PPND7 = - PPND7
	  RETURN
	END IF
	END

