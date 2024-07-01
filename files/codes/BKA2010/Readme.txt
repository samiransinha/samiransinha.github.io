This program can be used to obtain the parameter estimates and the corresponding 
standard error of the model used in Chatterjee  et al. (2010). 

1. In order to use it first one needs to save the two files in the same 
directory. 

2. In a linux system use the following command to load  the subroutines in 
R shared library.

   R CMD SHLIB subroutineBKA2010.f


3. In an R environment type in 
> source("ee_code_online.txt")

Now you ready to use the R function "ee".

#########
Following are the input variables for the function. 
1. Scalar covariate X
2. V=min(T, C), T: Failure time and C: Censoring time 
3. delta=I(T<=C)
4. ndischar: Number of disease characteristics
5. numlevels: a vector of integers which represent the number of levels of each characteristics. For 
  example, if the first, second and the third characteristics have 2, 4, and 2 nominal levels, respectively, 
  numlevels=c(2, 4, 2).  
6. disease.char: Matrix of disease characteristics for all the subjects. This is an n x ndischar matrix.
   Missing disease characteristics will be denoted by 0. Do not use NA values.
#########
Following are the output variables.
psi.est.ee:  Estimate of the main parameters (psi)
psi.se.ee:   Standard error of the estimates
alpha.est.ee: Estimate of alpha, the nuisance parameters involved in the baseline modelling
alpha.se.ee:  Standard error of the estimate of alpha
tolerance:    Tolerance acheived for the estimating equation approach
require.no.of.iteration:  Number of iteration taken for the proposed estimating equation appraoch
value.of.estimating.equations: Value of the estimating equations at when the estimates converged.
psi.est.pl:   Estimate of the main parameters using completely observed data and the method is based on the 
              partial likelihood technique
psi.se.pl:    Standard error of psi.est.pl
##########


#### Simulating survival data with multiple disease characteristics
n=8000;
x=rnorm(n, 0, 1)
ndischar=2;
numlevels=c(2, 2)
r1=1+rbinom(n, 1, .5)
r2=1+rbinom(n, 1, 0.5)
disease.char=cbind(r1, r2) # two disease characteristics
t=rexp(n, exp(.4+ .2*(r1-1)+0.5*(r2-1)) ) # true failure times
c=rnorm(n, 1, 0.7)
v=apply(cbind(t, c), 1, min)
delta=rep(0, n)
for( i in 1:n){ if(t[i]<c[i]) delta[i]=1}


#### Creating missing disease characteristics
m.r1=rbinom(n, 1, 0.25)
m.r2=rbinom(n, 1, 0.25)
disease.char[, 1][m.r1==0]<-0
disease.char[, 2][m.r2==0]<-0


source("ee_code_online_R.txt")
out=ee(x, v, delta, ndischar, numlevels, disease.char)  


   
