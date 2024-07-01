
#### The function is defined below. 
msw.simu=function(matched.case.control=matcdata, no.cases=n, 
                  no.controls.per.set=capm){
  ##########################
  #### This function conducts the analysis using the true treatment/exposure,
  #### the mismeasured treatment only, and using the instrumental variable 
  #### and the mismeasured  treatment. Note that for a real dataset, one cannot 
  #### apply the first method because the true treatment/exposure is not available
  
  myv=rep(rnorm(no.cases), each=(no.controls.per.set+1))
  myst= rep(1:n, each=(no.controls.per.set+1))
  outtr=coxph(Surv(myv, y)~x+z+strata(myst), matched.case.control)
  ##
  ##
  outw=coxph(Surv(myv, y)~w+z+strata(myst), matched.case.control)
  ##
  ##
  ##
  ##
  ### Our proposed method 
  ### of estimation of the first likelihood
  ### This function calculates the negative log likelihood used to 
  ### find the parameters of gamma and eta
  ### gamma parameters are par[1],par[2],par[3]
  ### eta parameters are par[4], par[5]
  
  nlglik1=function(par){
    # prpi=pr(T=1|T^*, Z) # a logistic model is assumed 
    prt1=1/(1+exp(-par[1]-par[2]*matched.case.control$xstar-par[3]*matched.case.control$sv))
    ## alpha.0=pr(W=1|X=0, Z)     
    alpha.0=1/(1+exp(-par[4])+exp(par[5]-par[4]) )
    ##   alpha.1=pr(W=0|X=1, Z)
    alpha.1=1/(1+exp(-par[5])+exp(par[4]-par[5]) )
    
    #In the event that we get an extremely high or extremely low
    #probability of treatment, set it to a specific value to avoid possible 
    #convergence issue.
    prt1[prt1>0.99999]=0.99999
    prt1[prt1<0.00001]=0.00001
    
    # prw.1=pr(W=1|X^*, Z)
    prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
    prw1[prw1>0.99999]=0.99999
    prw1[prw1<0.00001]=0.00001
    ########################
    #This section is used to calculate the analytic derivative used
    #in the penalty
    
    analytic.deriv=0*diag(5)
    
    term1=(
      (1-matched.case.control$y)*(1-alpha.0-alpha.1)^2*prt1^2*
        (1-prt1)^2/(prw1*(1-prw1))
    )
    #
    term2=(
      (1-matched.case.control$y)*(1-alpha.0-alpha.1)*prt1*(1-prt1)/(prw1*(1-prw1))
    )
    # 
    term3=(
      (1-matched.case.control$y)/(prw1*(1-prw1))
    )
    #
    qnty1=alpha.0*(1-alpha.0)+ (-alpha.0*(1-alpha.0)+alpha.0*alpha.1)*prt1
    #
    qnty2=-alpha.0*alpha.1+(alpha.0*alpha.1-alpha.1*(1-alpha.1))*prt1
    #
    analytic.deriv[1, 1]=analytic.deriv[1, 1]+sum(term1)
    analytic.deriv[1, 2]=analytic.deriv[1, 2]+sum(term1*matched.case.control$xstar)
    analytic.deriv[1, 3]=analytic.deriv[1, 3]+sum(term1*matched.case.control$sv)
    analytic.deriv[1, 4]=analytic.deriv[1, 4]+sum(term2*qnty1)
    analytic.deriv[1, 5]=analytic.deriv[1, 5]+sum(term2*qnty2)
    #
    analytic.deriv[2, 2]=analytic.deriv[2, 2]+sum(term1*matched.case.control$xstar^2)
    analytic.deriv[2, 3]=analytic.deriv[2, 3]+sum(term1*matched.case.control$xstar*matched.case.control$sv)
    analytic.deriv[2, 4]=analytic.deriv[2, 4]+sum(term2*qnty1*matched.case.control$xstar)
    analytic.deriv[2, 5]=analytic.deriv[2, 5]+sum(term2*qnty2*matched.case.control$xstar)
    #
    analytic.deriv[3, 3]=analytic.deriv[3, 3]+sum(term1*matched.case.control$sv^2)
    analytic.deriv[3, 4]=analytic.deriv[3, 4]+sum(term2*qnty1*matched.case.control$sv)
    analytic.deriv[3, 5]=analytic.deriv[3, 5]+sum(term2*qnty2*matched.case.control$sv)
    #
    #
    analytic.deriv[4, 4]=analytic.deriv[4, 4]+sum(term3*qnty1*qnty1)
    analytic.deriv[4, 5]=analytic.deriv[4, 5]+sum(term3*qnty1*qnty2)
    #
    analytic.deriv[5, 5]=analytic.deriv[5, 5]+sum(term3*qnty2*qnty2)
    
    
    analytic.deriv[2:5, 1]=analytic.deriv[1, 2:5]
    analytic.deriv[3:5, 2]=analytic.deriv[2, 3:5]
    analytic.deriv[4:5, 3]=analytic.deriv[3, 4:5]
    analytic.deriv[5, 4]=analytic.deriv[4, 5]
    #
    #
    #
    #########################
    
    neglk=- sum ((1-matched.case.control$y)*(log(prw1)*matched.case.control$w+(1-matched.case.control$w)*log(1-prw1)))-
      0.5*log(det(analytic.deriv))
    
    return(neglk)}
  
  ### We use the optim function to estimate the optimal values for the parameters from our first likelihood
  out2=optim(runif(5, -0.75, 0.75), nlglik1, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
  
  #Now store the optimal values as par tocalculate alpha.0,alpha.1,prt1,prw1
  par=out2$par
  alpha.0=1/(1+exp(-par[4])+exp(par[5]-par[4]) )# alpha.0=pr(W=1|X=0, Z)
  alpha.1=1/(1+exp(-par[5])+exp(par[4]-par[5]) )# alpha.1=pr(W=0|X=1, Z)
  prt1=1/(1+exp(-par[1]-par[2]*matched.case.control$xstar-par[3]*matched.case.control$sv)) #=pr(T=1|T^*, SV)
  prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
  
  ## pr(x=1|xstar, w=1, z)=pr(w=1|x=1)pr(x=1|xstar, z)/pr(w=1|xstar, z), and 
  ## pr(w=1|xstar, z)=pr(w=1|x=1, xstar, z)pr(x=1|xstar, z)+pr(w=1|x=0, xstar, z)pr(x=0|xstar, z)
  #
  #
  ## pr(x=1|xstar, w=0, z)=pr(w=0|x=1)pr(x=1|xstar, z)/pr(w=0|xstar, z), and 
  ## pr(w=0|xstar, z)=pr(w=0|x=1, xstar, z)pr(x=1|xstar, z)+pr(w=0|x=0, xstar, z)pr(x=0|xstar, z)
  ##
  ##
  #
  ###
  ###
  ###
  ##  Estimation of the beta parameters
  ###
  ## We use the nlglik2 function to create the log likelihood for the beta parameters, then use 
  ## the optim function to determine the optimal parameters for this function
  
  nlglik2=function(beta){
    exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
    exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)
    
    
    deriv.beta1.g.w1=(exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))
    
    deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))
    
    deriv.beta1.g=deriv.beta1.g.w1*matched.case.control$w+(1-matched.case.control$w)*deriv.beta1.g.w0
    
    exp.g=exp.g.w1*matched.case.control$w+(1-matched.case.control$w)*exp.g.w0
    
    ################ conditional probability calculation 
    term1=exp.g*exp(beta[2]*matched.case.control$z) 
    term2=matrix(term1, ncol=(no.controls.per.set+1), byrow=T)
    term3=rep(apply(term2, 1, sum), each=(no.controls.per.set+1))
    cond.prob=term1/term3
    
    #
    mat2=0*diag(2)
    newterm1=cond.prob*deriv.beta1.g
    newterm2=cond.prob*matcdata$z
    am1=matrix(newterm1, ncol=(no.controls.per.set+1), byrow=T)
    am2=matrix(newterm2, ncol=(no.controls.per.set+1), byrow=T)
    newterm3=apply(am1, 1, sum)
    newterm4=apply(am2, 1, sum)
    
    mat2[1,1]=sum(cond.prob*deriv.beta1.g^2)-sum(newterm3*newterm3)
    mat2[2, 2]=sum(cond.prob*matched.case.control$z^2)-sum(newterm4*newterm4)
    mat2[1, 2]=sum(cond.prob*deriv.beta1.g*matched.case.control$z)-sum(newterm3*newterm4)
    mat2[2, 1]=mat2[1, 2]
    
    ee=-sum(matched.case.control$y*log(cond.prob))-0.5*log(det(mat2)) 
    return(ee)
  }
  
  #As a starting point for possible parameters of beta, we use the
  #estimates from the mismeasured treatment analysis 
  #not necessary
  
  out3=optim(as.numeric(outw$coef), nlglik2, method="L-BFGS-B", control=list(maxit=1500))
  #  
  theta=out3$par
  
  #### Standard error calculation 
  #### Next we calculate the standard errors of all of the parameters
  #### First we need to calculate the "middle term" of the asymptotic variance
  #### It is composed of the U equations
  #### We need the output from the first analysis
  
  allpara=c(par, theta)
  npara=length(par)
  ntheta=length(theta)
  nallpara=length(allpara)
  
  par=allpara[1:npara]
  beta= allpara[(npara+1):nallpara]
  #
  alpha.0=1/(1+exp(-par[4])+exp(par[5]-par[4]) )# alpha.0=pr(W=1|X=0, Z)
  alpha.1=1/(1+exp(-par[5])+exp(par[4]-par[5]) )# alpha.1=pr(W=0|X=1, Z)
  prt1=1/(1+exp(-par[1]-par[2]*matched.case.control$xstar-par[3]*matched.case.control$sv)) #=pr(x=1|xstar, z)
  prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
  
  #
  exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
  exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)
  
  #der.g.w1=deriv of log(exp.g.w1)= (exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0(1-prt1))
  deriv.beta1.g.w1=(exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))
  deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))
  deriv.beta1.g=deriv.beta1.g.w1*matched.case.control$w+(1-matched.case.control$w)*deriv.beta1.g.w0
  exp.g=exp.g.w1*matched.case.control$w+(1-matched.case.control$w)*exp.g.w0
  
  ################ conditional probability calculation 
  term1=exp.g*exp(beta[2]*matched.case.control$z) 
  term2=matrix(term1, ncol=(no.controls.per.set+1), byrow=T)
  term3=rep(apply(term2, 1, sum), each=(no.controls.per.set+1))
  cond.prob=term1/term3
  #
  mid.term1=0
  term100=0
  #
  ### Here the estimates for the U equations are estimated
  uvec=rep(0, 7);
  uvec=matrix(NA,nrow=no.cases,ncol=7,byrow=T)
  
  term11=(
    (1-matched.case.control$y)*(1-alpha.0-alpha.1)*prt1*
      (1-prt1)*(matched.case.control$w-prw1)/(prw1*(1-prw1))
  )
  
  term1.mat=matrix(term11, ncol=(no.controls.per.set+1), byrow=T)
  
  uvec[,1]=apply(term1.mat,1,sum)
  
  uvec21=term11*matched.case.control$xstar
  uvec21.mat=matrix(uvec21,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,2]=apply(uvec21.mat,1,sum)
  uvec31=term11*matched.case.control$sv
  uvec31.mat=matrix(uvec31,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,3]=apply(uvec31.mat,1,sum)
  
  term2= (1-matched.case.control$y)*(matched.case.control$w-prw1)/(prw1*(1-prw1))
  
  uvec41=term2*( alpha.0*(1-alpha.0)-alpha.0*(1-alpha.0-alpha.1)*prt1 )
  uvec41.mat=matrix(uvec41,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,4]=apply(uvec41.mat,1,sum)
  uvec51= term2*( -alpha.0*alpha.1-alpha.1*(1-alpha.0-alpha.1)*prt1) 
  uvec51.mat=matrix(uvec51,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,5]=apply(uvec51.mat,1,sum)
  
  
  newterm1=matched.case.control$y-cond.prob
  
  mult=(matched.case.control$y-cond.prob)*deriv.beta1.g
  mult.mat=matrix(mult,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,6]=apply(mult.mat,1,sum)
  
  uvec71=newterm1*matched.case.control$z
  uvec71.mat=matrix(uvec71,ncol=(no.controls.per.set+1),byrow=T)
  uvec[,7]=apply(uvec71.mat,1,sum)
  
  ##After calculating the uvec, we can now calculate the middle term 
  ##mid.term
  
  mid.term1=t(uvec)%*%uvec
  mid.term1=mid.term1/no.cases
  mid.term=mid.term1
  
  #### Derivative calculation using numerical approach #####
  #### After calculating the middle term, we now need to calculate A-hat
  #### 
  forderiveq2=function(allpara){
    par=allpara[1:5]
    beta= allpara[6:7]
    #
    alpha.0=1/(1+exp(-par[4])+exp(par[5]-par[4]) )# alpha.0=pr(W=1|T=0, Z)
    alpha.1=1/(1+exp(-par[5])+exp(par[4]-par[5]) )# alpha.1=pr(W=0|T=1, Z)
    prt1=1/(1+exp(-par[1]-par[2]*matched.case.control$xstar-par[3]*matched.case.control$sv)) #=pr(X=1|X^*, Z)
    prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
    # 
    exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
    exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)
    #
    #
    deriv.beta1.g.w1=(exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))
    deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))
    deriv.beta1.g=deriv.beta1.g.w1*matched.case.control$w+(1-matched.case.control$w)*deriv.beta1.g.w0
    exp.g=exp.g.w1*matched.case.control$w+(1-matched.case.control$w)*exp.g.w0
    ##
    ##
    #### Conditional probability calculation 
    term1=exp.g*exp(beta[2]*matched.case.control$z) 
    term2=matrix(term1, ncol=(no.controls.per.set+1), byrow=T)
    term3=rep(apply(term2, 1, sum), each=(no.controls.per.set+1))
    cond.prob=term1/term3
    
    mat2=0*diag(2)
    newterm1=cond.prob*deriv.beta1.g
    newterm2=cond.prob*matched.case.control$z
    am1=matrix(newterm1, ncol=(no.controls.per.set+1), byrow=T)
    am2=matrix(newterm2, ncol=(no.controls.per.set+1), byrow=T)
    newterm3=apply(am1, 1, sum)
    newterm4=apply(am2, 1, sum)
    
    
    mat2[1,1]=sum(cond.prob*deriv.beta1.g^2)-sum(newterm3*newterm3)
    mat2[2, 2]=sum(cond.prob*matched.case.control$z^2)-sum(newterm4*newterm4)
    mat2[1, 2]=sum(cond.prob*deriv.beta1.g*matched.case.control$z)-sum(newterm3*newterm4)
    mat2[2, 1]=mat2[1, 2]
    
    ee=-sum(matched.case.control$y*log(cond.prob))-0.5*log(det(mat2)) 
    return(ee)
  }
  #
  #
  #
  deriv=0*diag(7)#jacobian(forderiv, allpara)/no.cases
  deriv[1:5, 1:5]=-out2$hessian/no.cases
  out5=hessian(forderiveq2, allpara)
  deriv[6:7, ]=-out5[6:7, ]/no.cases
  
  ##########
  ##########
  ##Finally calculate the A hat, then calculate the asymptotic variance
  inv.deriv=solve(deriv)
  var.cov= inv.deriv%*%mid.term%*%t(inv.deriv)/no.cases
  theta.sd=sqrt(diag(var.cov))
  
  ########################################################################### 
  #Finally, we save our results 
  M1= list(as.numeric(summary(outtr)$coeff[, 1]), as.numeric(summary(outtr)$coeff[, 3]))
  names(M1)=c("Est", "SE")
  M2= list(as.numeric(summary(outw)$coeff[, 1]), as.numeric(summary(outw)$coeff[, 3]))
  names(M2)=c("Est", "SE")
  M3=list(theta, theta.sd[6:7])
  names(M3)=c("Est", "SE")
  result=list(M1, M2, M3)
  names(result)<-c("M1", "M2", "M3")
  #
  #
  #
  return(result)}

