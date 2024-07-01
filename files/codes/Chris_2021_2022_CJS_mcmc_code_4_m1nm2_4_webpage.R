# First, data are generated, and then 
# analyzed using methods M1 or M2. In this analysis 
# MCMC method is used. 

# The following library is needed
library(MCMCpack)

# Data simulation 
# set the seed 
set.seed(5)
# Set the sample size
n=10000
# generation of error free covariates
z1=rbinom(n, 1, 0.23)
z2=rbinom(n, 1, 0.31)
# Generation of an instrumental variable with three categories. 
# Note that the instrumental variable need not be a categorical variable. 
# Also, there could be multiple instrumental variables. 
x0star=t(rmultinom(n, 1, prob=c(0.74, 0.22, 0.04)))%*%c(1, 2, 3)
# Two dummy variables corresponding to the three-category
# nominal instrumental variable.  
x1star=ifelse(x0star==2, 1, 0)
x2star=ifelse(x0star==3, 1, 0)
# Generation of X
## If needed change the coefficients of f(X|X^*, Z) of the following two lines
##
pr2= exp(0.47+0.8*x1star+0.36*x2star-0.99*z1-2*z2)
pr3= exp(-0.09 - 4.12*x1star+0.87*x2star-0.45*z1+1.75*z2)
deno=1+pr2+pr3
pr1=1/deno
pr2=pr2/deno
pr3=pr3/deno
prob=cbind(pr1, pr2, pr3)
f1=function(i) sum(rmultinom(1, 1, prob[i, ])*(1:3))
myx=unlist(lapply(seq(1:n), f1))
# Generation of the response
## If needed change the coefficients of f(Y|X, Z) of the following line
eta=1.13+ 0*as.numeric(myx==2) +1.98*as.numeric(myx==3) - 0.5*z1-2*z2
pry= 1/(1+exp(-eta))
y=rbinom(n, 1, pry)
### data where x is being recorded  
originaldata=data.frame(myx, y, z1, z2)
head(originaldata)
### generation of the surrogate W
### Generation of the misclassification matrix 
misclassmat=matrix(c(0.96, 0.025, 0.025, 0.1, 0.85, 0.05, 0.2, 0.1, 0.7), byrow=F, ncol=3)
### For a different misclassification, change the entries of the misclassification matrix. 
f2=function(i) sum(rmultinom(1, 1, misclassmat[,i])*(1:3))
myw=unlist(lapply(myx, f2))
# For checking you may print 
table(myw, myx)
###
###
## setting the standard deviation of the prior distribution 
## of the parameters
 priorsd=2 # later change it to 5
# DO NOT MAKE ANY CHANGES AFTER THIS LINE
##
##
## data without x, rather the data contain the binary response, w, 
## instrumental variable(s) and the error-free covariates 
mydata=data.frame(myw, x1star, x2star, y, z1, z2)

## The frequentist analysis of the data that does not include the 
## true exposure, myx, rather include the surrogate exposure 
## measurements. This analysis is referred to as M2 in the article. 
## For obtaining results under M1, one should replace mydata by 
## originaldata in the following lines of code. 
## 
out1=glm(y~as.factor(myw)+z1+z2, family=binomial, data=mydata)
## The Bayesian analysis of the data 
## Setting the mcmc iterations. 
nmcmc=15000
## Setting the number of burn-in samples
burnin=5000
myk=ncol(model.matrix(out1))

logpriorfun <- function(beta){
   sum(dnorm(beta, mean=0, sd=priorsd, log=TRUE))
 }
posterior.m2 <- MCMClogit(y ~ as.factor(myw) + z1 + z2,
                          data = mydata, 
                          user.prior.density=logpriorfun, mcmc=nmcmc, 
                          logfun=TRUE)
## Collecting MCMC samples after the burn-in samples
mcmcsamplesm2=posterior.m2 [-c(1:burnin), ]  
qntfunm2=function(i) quantile(mcmcsamplesm2[, i], prob=c(0.025, 0.975))
## The following line of code returns the posterior mean, median, standard 
## deviation, and the 95% credible interval for each of the model parameters. 
myresultm2=list(posterior_mean=apply(mcmcsamplesm2, 2, mean), 
                posterior_median=apply(mcmcsamplesm2, 2, median), 
                posterior_sd=apply(mcmcsamplesm2, 2, sd), 
                credible_interval=as.numeric(unlist(lapply(1:myk, qntfunm2))))

               


