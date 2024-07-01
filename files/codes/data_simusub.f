      double precision function pr(x)
      double precision x
      if(x.gt.18.d0) then
       pr=1.d0
      else
       pr=exp(x)/(1+exp(x))
      endif
      return
      end
C ***********************************
       subroutine paikstage1(matchcc, nrow, ncol, theta, theta_old, 
     * ndthet, lglk, thetax)!data=data.frame(y, s, z, x) for binary exposure variable
       integer*4 nrow, ncol, ndthet
       double precision matchcc(nrow, ncol), theta(ndthet), 
     * theta_old(ndthet), lglk, thetax
       double precision eta_delta, eta_delta_0, eta_delta_1, eta_x, 
     * eta_x_old, eta_delta_0old, eta_delta_1old, compo1, compo2, 
     * compo3, pr_delta, pr_delta_0, pr_delta_1, pr_x, pr_x_old, pr,
     * temp

C      theta=c(delta0, deltad, deltas, deltaz,
C      gamma0, gammad, gammas, gammaz) 
        compo1=0.d0
        compo2=0.d0
        compo3=0.d0       

       do i=1, nrow
       eta_delta=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*matchcc(i, 2)
     * +thetax*matchcc(i, 4)+theta(4)*matchcc(i, 3)
       eta_delta_1=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*
     * matchcc(i, 2)+thetax*1.d0+theta(4)*matchcc(i, 3)
       eta_delta_0=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*
     * matchcc(i, 2)+thetax*0.d0+theta(4)*matchcc(i, 3)
       eta_x=theta(5)+theta(6)*matchcc(i, 1)+theta(7)*matchcc(i, 2)+
     * theta(8)*matchcc(i, 3)
C -------------------------------
       eta_x_old=theta_old(5)+theta_old(6)*matchcc(i, 1)+theta_old(7)*
     * matchcc(i, 2)+theta_old(8)*matchcc(i, 3)
       eta_delta_1old=theta_old(1)+theta_old(2)*matchcc(i, 1)+
     * theta_old(3)*matchcc(i, 2)+thetax*1.d0+theta_old(4)*
     * matchcc(i, 3)
       eta_delta_0old=theta_old(1)+theta_old(2)*matchcc(i, 1)+
     * theta_old(3)*matchcc(i, 2)+thetax*0.d0+theta_old(4)*
     * matchcc(i, 3)
C -------------------------------
       pr_x=max(0.000001d0, pr(eta_x))
       
       pr_delta=max(0.000001d0, pr(eta_delta))
       pr_delta_0=max(0.000001d0, pr(eta_delta_0))
       pr_delta_1=max(0.000001d0, pr(eta_delta_1))
       temp=eta_x_old+log(max(0.000001d0, pr(eta_delta_1old))/max(
     * 0.000001d0, pr(eta_delta_0old)))
       pr_x_old=max(0.000001d0, pr(temp))

       compo1=compo1+log(pr_delta)*matchcc(i, 5)
     *  + log(pr_x)*matchcc(i, 4)*matchcc(i, 5)+
     *  log(1-pr_x)*(1-matchcc(i, 4))*matchcc(i, 5)
        compo2=compo2+log(1-pr_delta_1)*
     *  (1-matchcc(i, 5))*pr_x_old+log(1-pr_delta_0)*
     *  (1-matchcc(i, 5))*(1-pr_x_old)
        compo3=compo3+(1-matchcc(i, 5))*pr_x_old*log(pr_x)*
     *  (1-matchcc(i, 5))*(1-pr_x_old)*
     *  log( 1-pr_x) 
        end do
        lglk=compo1+compo2+compo3
        return
        end
C *************************************************
       subroutine newpaikstage1(matchcc, nrow, ncol, theta, theta_old, 
     * ndthet, lglk)!data=data.frame(y, s, z, x) for binary exposure variable
       integer*4 nrow, ncol, ndthet
       double precision matchcc(nrow, ncol), theta(ndthet), 
     * theta_old(ndthet), lglk
       double precision eta_delta, eta_delta_0, eta_delta_1, eta_x, 
     * eta_x_old, eta_delta_0old, eta_delta_1old, compo1, compo2, 
     * compo3, pr_delta, pr_delta_0, pr_delta_1, pr_x, pr_x_old, pr,
     * temp

C      theta=c(delta0, deltad, deltas, deltax, deltaz,
C      gamma0, gammad, gammas, gammaz) 
        compo1=0.d0
        compo2=0.d0
        compo3=0.d0       

       do i=1, nrow
       eta_delta=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*matchcc(i, 2)
     * +theta(4)*matchcc(i, 4)+theta(5)*matchcc(i, 3)
       eta_delta_1=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*
     * matchcc(i, 2)+theta(4)*1.d0+theta(5)*matchcc(i, 3)
       eta_delta_0=theta(1)+theta(2)*matchcc(i, 1)+theta(3)*
     * matchcc(i, 2)+theta(4)*0.d0+theta(5)*matchcc(i, 3)
       eta_x=theta(6)+theta(7)*matchcc(i, 1)+theta(8)*matchcc(i, 2)+
     * theta(9)*matchcc(i, 3)
C -------------------------------
       eta_x_old=theta_old(6)+theta_old(7)*matchcc(i, 1)+theta_old(8)*
     * matchcc(i, 2)+theta_old(9)*matchcc(i, 3)
       eta_delta_1old=theta_old(1)+theta_old(2)*matchcc(i, 1)+
     * theta_old(3)*matchcc(i, 2)+theta_old(4)*1.d0+theta_old(5)*
     * matchcc(i, 3)
       eta_delta_0old=theta_old(1)+theta_old(2)*matchcc(i, 1)+
     * theta_old(3)*matchcc(i, 2)+theta_old(4)*0.d0+theta_old(5)*
     * matchcc(i, 3)
C -------------------------------
       pr_x=max(0.000001d0, pr(eta_x))
       
       pr_delta=max(0.000001d0, pr(eta_delta))
       pr_delta_0=max(0.000001d0, pr(eta_delta_0))
       pr_delta_1=max(0.000001d0, pr(eta_delta_1))
       temp=eta_x_old+log(max(0.000001d0, pr(eta_delta_1old))/max(
     * 0.000001d0, pr(eta_delta_0old)))
       pr_x_old=max(0.000001d0, pr(temp))

       compo1=compo1+log(pr_delta)*matchcc(i, 5)
     *  + log(pr_x)*matchcc(i, 4)*matchcc(i, 5)+
     *  log(1-pr_x)*(1-matchcc(i, 4))*matchcc(i, 5)
        compo2=compo2+log(1-pr_delta_1)*
     *  (1-matchcc(i, 5))*pr_x_old+log(1-pr_delta_0)*
     *  (1-matchcc(i, 5))*(1-pr_x_old)
        compo3=compo3+(1-matchcc(i, 5))*pr_x_old*log(pr_x)*
     *  (1-matchcc(i, 5))*(1-pr_x_old)*
     *  log( 1-pr_x) 
        end do
        lglk=compo1+compo2+compo3
        return
        end
C *************************************************

C *************************************************
         subroutine sm(matchcc, nrow, ncol, alpha, alpha_old, ndalpa,
     *   lglk, alphax)!data=data.frame(y, s, z, x) for binary exposure variable
         integer*4 nrow, ncol, ndalpa
         double precision matchcc(nrow, ncol), alpha(ndalpa), 
     *   alpha_old(ndalpa), lglk, alphax
C  alpha=(delta0, deltad, deltas,  deltaz,
C      gamma0, gammas, gammaz, beta1, beta2)
         integer*4 nrowby2, i1, i2
         double precision compo1, compo2, compo3
         double precision eta_x, eta_delta, eta_delta_1, eta_delta_0
         double precision eta_x_old, eta_delta_1old, eta_delta_0old
         double precision pr_x, pr_delta, pr_delta_0, pr_delta_1
         double precision temp, q1, q2, pr_x_old, pr
         double precision eta_xwd, eta_xwotd
         
          compo1=0.d0
          compo2=0.d0
          compo3=0.d0 
          do i=1, nrow
          eta_x=alpha(5)+alpha(6)*matchcc(i, 2)+alpha(7)*matchcc(i, 3)+
     *    alpha(9)*matchcc(i, 1)! it takes into account eta_x and eta_star_x
          eta_delta=alpha(1)+alpha(2)*matchcc(i, 1)+
     *    alpha(3)*matchcc(i, 2)+alphax*matchcc(i, 4)+
     *    alpha(4)*matchcc(i, 3)
          eta_delta_1=alpha(1)+alpha(2)*matchcc(i, 1)+
     *    alpha(3)*matchcc(i, 2)+alphax*1.d0+
     *    alpha(4)*matchcc(i, 3)
          eta_delta_0=alpha(1)+alpha(2)*matchcc(i, 1)+
     *    alpha(3)*matchcc(i, 2)+alphax*0.d0+
     *    alpha(4)*matchcc(i, 3)
C ----------------------------
          pr_x=max(0.000001d0, pr(eta_x))
          pr_delta=max(0.000001d0, pr(eta_delta))
          pr_delta_0=max(0.000001d0, pr(eta_delta_0))
          pr_delta_1=max(0.000001d0, pr(eta_delta_1))
C ----------------------------
          eta_x_old=alpha_old(5)+alpha_old(6)*matchcc(i, 2)+
     *    alpha_old(7)*matchcc(i, 3)+alpha_old(9)*matchcc(i, 1)! it 
C         takes into account eta_x and eta_star_x
          eta_delta_1old=alpha_old(1)+alpha_old(2)*matchcc(i, 1)+
     *    alpha_old(3)*matchcc(i, 2)+alphax*1.d0+
     *    alpha_old(4)*matchcc(i, 3)
          eta_delta_0old=alpha_old(1)+alpha_old(2)*matchcc(i, 1)+
     *    alpha_old(3)*matchcc(i, 2)+alphax*0.d0+
     *    alpha_old(4)*matchcc(i, 3)
C ----------------------------
          temp=eta_x_old+log(max(0.000001d0, pr(eta_delta_1old))/
     *    max(0.000001d0, pr(eta_delta_0old)))
          pr_x_old=max(0.000001d0, pr(temp))
C ----------------------------
          compo2=compo2+matchcc(i, 5)*(matchcc(i, 4)*log(pr_x)+
     *    (1-matchcc(i, 4))*log(1-pr_x))+
     *    matchcc(i, 5)*log(pr_delta)

          compo3=compo3+(1-matchcc(i, 5))*(pr_x_old*log(pr_x)+
     *    (1-pr_x_old)*log(1-pr_x))+(1-matchcc(i, 5))*(pr_x_old*
     *    log(1-pr_delta_1) +(1-pr_x_old)*log(1-pr_delta_0))
          end do
          nrowby2=nrow/5
          do j=1, nrowby2
          i1=5*(j-1)+1          
          eta_xwd=alpha(5)+alpha(6)*matchcc(i1, 2)+
     *    alpha(7)*matchcc(i1, 3)+alpha(9)*1.d0! it 
          eta_xwotd=alpha(5)+alpha(6)*matchcc(i1, 2)+
     *    alpha(7)*matchcc(i1, 3) 
          q1=alpha(8)*matchcc(i1, 3)+log(1+exp(eta_xwd))-
     *    log(1+exp(eta_xwotd))
          q2=exp(q1)

          do il=2, 5
           i2=5*(j-1)+il
           eta_xwd=alpha(5)+alpha(6)*matchcc(i2, 2)+
     *     alpha(7)*matchcc(i2, 3)+alpha(9)*1.d0! it 
           eta_xwotd=alpha(5)+alpha(6)*matchcc(i2, 2)+
     *     alpha(7)*matchcc(i2, 3) 
           q2=q2+exp(alpha(8)*matchcc(i2, 3)+log(1+exp(eta_xwd))-
     *     log(1+exp(eta_xwotd)))
          end do
          compo1=compo1+q1-log(q2)
          end do
          lglk=compo1+compo2+compo3
          return
          end 
C **********************************************************************
         subroutine mmarsm(matchcc, nrow, ncol, phi, phi_old, 
     *   ndphi, lglk)!data=data.frame(y, s, z, x) for binary exposure variable
         integer*4 nrow, ncol, ndphi
         double precision matchcc(nrow, ncol), phi(ndphi), 
     *   phi_old(ndphi), lglk
C  phi=(gamma0, gammas, gammaz, beta1, beta2)
         integer*4 nrowby2, i1, i2
         double precision compo1, compo2, compo3
         double precision eta_x, eta_x_old
         double precision pr_x, pr_x_old, pr
         double precision q1, q2 
         double precision eta_xwd, eta_xwotd
          compo1=0.d0
          compo2=0.d0
          compo3=0.d0 
          do i=1, nrow
          eta_x=phi(1)+phi(2)*matchcc(i, 2)+phi(3)*matchcc(i, 3)+
     *    phi(5)*matchcc(i, 1)! it takes into account eta_x and eta_star_x
C ----------------------------
          pr_x=max(0.000001d0, pr(eta_x))
C ----------------------------
          eta_x_old=phi_old(1)+phi_old(2)*matchcc(i, 2)+
     *    phi_old(3)*matchcc(i, 3)+phi_old(5)*matchcc(i, 1)! it 
          pr_x_old=max(0.000001d0, pr(eta_x_old))
C ----------------------------
          compo2=compo2+matchcc(i, 5)*(matchcc(i, 4)*log(pr_x)+
     *    (1-matchcc(i, 4))*log(1-pr_x))

          compo3=compo3+(1-matchcc(i, 5))*(pr_x_old*log(pr_x)+
     *    (1-pr_x_old)*log(1-pr_x))
          end do
          nrowby2=nrow/2
          do j=1, nrowby2
          i1=2*(j-1)+1
          i2=2*(j-1)+2
          eta_xwd=phi(1)+phi(2)*matchcc(i1, 2)+
     *    phi(3)*matchcc(i1, 3)+phi(5)*1.d0! it 
          eta_xwotd=phi(1)+phi(2)*matchcc(i1, 2)+
     *    phi(3)*matchcc(i1, 3) 
          q1=phi(4)*matchcc(i1, 3)+log(1+exp(eta_xwd))-
     *    log(1+exp(eta_xwotd))

          eta_xwd=phi(1)+phi(2)*matchcc(i2, 2)+
     *    phi(3)*matchcc(i2, 3)+phi(5)*1.d0! it 
          eta_xwotd=phi(1)+phi(2)*matchcc(i2, 2)+
     *    phi(3)*matchcc(i2, 3) 
          q2=phi(4)*matchcc(i2, 3)+log(1+exp(eta_xwd))-
     *    log(1+exp(eta_xwotd))
          
          compo1=compo1+q1-log(exp(q1)+exp(q2))
          end do
          lglk=compo1+compo2+compo3
          return
          end 

