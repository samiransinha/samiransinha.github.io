We describe how to use the SPB function to get an estimate of 
$\theta_{risk}$ and $\theta_{cal}$ functions described in the paper 
"Semiparametric Bayesian analysis of Nutritional Epidemiology Data in 
the Presence of Measurement Error" by Sinha, S., Mallick, B. K., 
Kipnis, V., and Carroll, R. J. 


The subroutines for the program is written in Fortran. Therefore, before 
the use of SPB one should first compile and load the subroutine in R 
shared library. Here are the necessary commands to use in a linux system. 

[user@machine:~]% gforan -O2 -c spbsubroutines.f

[user@machine:~]% R CMD SHLIB spbsubroutines.o


Here is an example of how to use SPB function based on a simulated data set.  

rm(list=ls(all=TRUE))
# Here both  both theta_cal and theta_risk are treated nonparametrically.
# First we generate the data mimicking the original data, then we 
# analyze it by the proposed method.
# Data simulation
n=10000; # cohort size
m=1000 # calibration data size
sigmaq=0.2; 
sigmam=0.5; 
mnofx=3.5; sigmax=0.5; # We need this information if we simulated X from a 
# normal distribution. 
x=rep(0, n)
x=rnorm(n, mnofx, sigmax);
#for( i in 1:n){ if(rbinom(1, 1, 0.5)==1) x[i]=rgamma(1, 65, (1/0.0455)) else x[i]=rnorm(1, 3.8, 0.10)} # if 
# we generate X from some non-normal distribution. 
errm1=rnorm(n, 0, sigmam);
errm2=rnorm(n, 0, sigmam);
w1=x+errm1; w2=x+errm2;
errq=rnorm(n, 0, sigmaq);
#theta.cal=2+3*exp(4*(x-mnofx))/(1+exp(4*(x-mnofx))); # this shows a nonlinear association between x and qx
theta.cal=1.25+0.6*x# this shows a linear relationship between x and qx
q=theta.cal+errq;
y=rep(0, n);
theta.risk=-2.8+1.5*exp(10*(x-mnofx))/(1+exp(10*(x-mnofx)));
prob=exp(theta.risk)/(1+exp(theta.risk))
for( i in 1:n)y[i]=rbinom(1, 1, prob[i])
cohort.data=cbind(y, q)
calibration.data=cbind(y[1:m], q[1:m], w1[1:m], w2[1:m])
# End of data simulation



#
# Beginning of the data analysis


output
alpha: MCMC sample of \alpha_1, \cdots, \alpha_q

beta: MCMC sample of \alpha_1, \cdots, \alpha_p

sigma2q: MCMC sample of sigma2q

sigma2m: MCMC sample of sigma2m

k: number of clusters of the Dirichlet process prior in each MCMC sample

alpha0: MCMC sample of \alpha_0 of the precision parameter of the Dirichlet process prior,
        which is determined emprically

# In the following we initialize the parameters:
>a0=2.1; b0=0.25
>totiter=30000

>out=spb(cohort.data, calibration.data, totiter, a0, b0)


