#########################################################################
# This R program calculates required sample size for 1:3 matched case-control
# study with a categorical exposure variable with 3 nominal levels.
# p0i=Prevalence of level i of the exposure variable among the control polulation, i=0, 1,2
# p1i=Prevalence of level i of the exposure variable among the case polulation, i=0, 1,2
# psi_1, psi_2 are the two odds ratios. 
# M: number of controls in every matched set
# Function "sample.size" calculates the required sample sizes using other functions such as
# "possible.combination", "power", "pbv", "ptv", "cnd.prb", and "cond.prob".
# Input variable(s): p01, p02, psi1, psi2, M
# Output of the function "sample.size": Number of matched sets 
##########################################################################
 library(nnet)
 possible.combination=function(M){
 M1=M*(M-1)/2
 amatrix=matrix(0, ncol=2, nrow=M1)
 x1=1; x2=0
 i=0
 while(x1+x2<M){
 while(x1+x2<M){
 i=i+1
 x2=x2+1
 amatrix[i,]=c(x1, x2)}
 x1=x1+1
 x2=0}
# nrw=which.max(amatrix)
# vec=amatrix[1:nrw,]
 return(amatrix)}



### Following function calculates the probability that
### a matched set has m subjects exposed at the level r1
### and the rest are exposed at the level r2; m=1,..M.
 pbv=function(r1,r2,p01,p02,psi1,psi2,m,M){
 p00=1-p01-p02;

 p10=1/(1+(psi1*p01/p00)+(psi2*p02/p00))
 p11=p10*psi1*p01/p00
 p12=p10*psi2*p02/p00

      if((r1==1)&(r2==2)) {pbv=p11*choose(M, m-1)*(p01^(m-1))*(p02^(M-m+1)) +p12*choose(M, m)*(p01^m)*(p02^(M-m)) }
 else if((r1==1)&(r2=0))  {pbv=p11*choose(M, m-1)*(p01^(m-1))*(p00^(M-m+1)) +p10*choose(M, m)*(p01^m)*(p00^(M-m)) }
 else                     {pbv=p12*choose(M, m-1)*(p02^(m-1))*(p00^(M-m+1)) +p10*choose(M, m)*(p02^m)*(p00^(M-m)) }
 return(pbv)}



 ## Following function calculates the probability that a
 ## matched set has m1 subjects exposed at the level 1 and
 ## m2 are exposed at the level 2; m1,m2=1,..M, and m1+m2<=M.
 ptv=function(p01,p02,psi1,psi2, m1,m2,M){
 p00=1-p01-p02;

 p10=1/(1+(psi1*p01/p00)+(psi2*p02/p00))
 p11=p10*psi1*p01/p00
 p12=p10*psi2*p02/p00

 ## Probability that a matched set has m1 of category 1 and m2 of category 2
 if(M-m1-m2<0)q1=0 else{
 q1=p10*(factorial(M)/(factorial(m1)*factorial(m2)*factorial(M-m1-m2)))*(p01^m1)*(p02^m2)*
 (1-p01-p02)^(M-m1-m2)}
 if(m1-1<0) q2=0 else{
 q2=p11*(factorial(M)/(factorial(m1-1)*factorial(m2)*factorial(M-m1-m2+1)))*(p01^(m1-1))*(p02^m2)*
 (1-p01-p02)^(M-m1-m2+1)}
 if(m2-1<0)q3=0 else{
 q3=p12*(factorial(M)/(factorial(m1)*factorial(m2-1)*factorial(M-m1-m2+1)))*(p01^(m1))*
 (p02^(m2-1))*
 (1-p01-p02)^(M-m1-m2+1)}
 ptv=q1+q2+q3
 return(ptv)
 }

## Function
# Conditional Proability
 cnd.prb=function(r,psi1,psi2,m1,m2,M){
 if(r==1) cnd.prb=m1*psi1/(m1*psi1+m2*psi2+(M-m1-m2+1)) else{
 cnd.prb=m2*psi2/(m1*psi1+m2*psi2+(M-m1-m2+1))}
 return(cnd.prb)}
###################################################
 cond.prob=function(r1,r2,psi1,psi2,m,M){
 if((r1==1)&(r2==2)){cond.prob=m*psi1/(m*psi1+(M-m+1)*psi2)
 } else if((r1==1)&(r2==0)) {cond.prob=m*psi1/(m*psi1+(M-m+1))}
 else {cond.prob=m*psi2/(m*psi2+(M-m+1))}
 return(cond.prob)}
##################################################
 power=function(p01,p02,psi1,psi2,N,M){
 vec.m=possible.combination(M)
 
 d11.c1=0; d12.c1=0;d22.c1=0;mu1.c1=0;mu2.c1=0;
 d11h1.c1=0;d12h1.c1=0;d22h1.c1=0;
 for(i in 1:nrow(vec.m)){
 d11.c1=d11.c1+N*ptv(p01,p02,1,1, vec.m[i,1],vec.m[i,2],M)*
 vec.m[i,1]*(M+1-vec.m[i,1])/((M+1)^2)
 d12.c1=d12.c1-N*ptv(p01,p02,1,1, vec.m[i,1],vec.m[i,2],M)*vec.m[i,1]*
 vec.m[i,2]/((M+1)^2)
 d22.c1=d22.c1+N*ptv(p01,p02,1,1, vec.m[i,1],vec.m[i,2],M)*
 vec.m[i,2]*(M+1-vec.m[i,2])/((M+1)^2)

 d11h1.c1=d11h1.c1+N*ptv(p01,p02,psi1,psi2, vec.m[i,1],vec.m[i,2],M)*
 cnd.prb(1,psi1,psi2,vec.m[i,1],vec.m[i,2],M)*(1-
 cnd.prb(1,psi1,psi2,vec.m[i,1],vec.m[i,2],M))

 d12h1.c1=d12h1.c1-N*ptv(p01,p02,psi1,psi2, vec.m[i,1],vec.m[i,2],M)*
 cnd.prb(1,psi1,psi2,vec.m[i,1],vec.m[i,2],M)*
 cnd.prb(2,psi1,psi2,vec.m[i,1],vec.m[i,2],M)

 d22h1.c1=d22h1.c1+N*ptv(p01,p02,psi1,psi2, vec.m[i,1],vec.m[i,2],M)*
 cnd.prb(2,psi1,psi2,vec.m[i,1],vec.m[i,2],M)*(1-
 cnd.prb(2,psi1,psi2,vec.m[i,1],vec.m[i,2],M))

 mu1.c1=mu1.c1+N*ptv(p01,p02,psi1,psi2, vec.m[i,1],vec.m[i,2],M)*
 cnd.prb(1,psi1,psi2,vec.m[i,1],vec.m[i,2],M)-N*ptv(p01,p02,1,1, vec.m[i,1],vec.m[i,2],M)*(vec.m[i,1]/(M+1))
 mu2.c1=mu2.c1+N*ptv(p01,p02,psi1,psi2, vec.m[i,1],vec.m[i,2],M)*cnd.prb(2,
 psi1,psi2,vec.m[i,1],vec.m[i,2],M)-N*ptv(p01,p02,1,1, vec.m[i,1],vec.m[i,2],M)*(vec.m[i,2]/(M+1))}

 d11.c2=0; d12.c2=0;d22.c2=0;d11h1.c2=0;d12h1.c2=0;d22h1.c2=0;
 mu1.c2=0;mu2.c2=0;

 for( i in 1:M){
 d11.c2=d11.c2+N*(pbv(1,2,p01,p02,1,1,i,M)+
 pbv(1,0,p01,p02,1,1,i,M))*i*(M-i+1)/((M+1)^2)
 d22.c2=d22.c2+N*(pbv(1,2,p01,p02,1,1,i,M)
 +pbv(2,0,p01,p02,1,1,i,M))*i*(M-i+1)/((M+1)^2)
 d12.c2=d12.c2-N*pbv(1,2,p01,p02,1,1,i,M)*i*(M-i+1)/((M+1)^2)
 d11h1.c2=d11h1.c2+N*pbv(1,2,p01,p02,psi1,psi2,i,M)*
 cond.prob(1,2,psi1,psi2,i,M)*(1-cond.prob(1,2,psi1,psi2,i,M))
 +N*pbv(1,0,p01,p02,psi1,psi2,i,M)*
 cond.prob(1,0,psi1,psi2,i,M)*(1-cond.prob(1,0,psi1,psi2,i,M))

 d22h1.c2=d22h1.c2+N*pbv(1,2,p01,p02,psi1,psi2,i,M)*
 cond.prob(1,2,psi1,psi2,i,M)*(1-cond.prob(1,2,psi1,psi2,i,M))
 +N*pbv(2,0,p01,p02,psi1,psi2,i,M)*
 cond.prob(2,0,psi1,psi2,i,M)*(1-cond.prob(2,0,psi1,psi2,i,M))

 d12h1.c2=d12h1.c2-N*pbv(1,2,p01,p02,psi1,psi2,i,M)*
 cond.prob(1,2,psi1,psi2,i,M)*(1-cond.prob(1,2,psi1,psi2,i,M))

 mu1.c2=mu1.c2+N*pbv(1,2,p01,p02,psi1,psi2,i,M)*cond.prob(1,2,psi1,psi2,i,M)-
N*pbv(1,2,p01,p02,1,1,i,M)*(i/(M+1))+N*pbv(1,0,p01,p02,psi1,psi2,i,M)*cond.prob(1,0,psi1,psi2,i,M)-N*pbv(1,0,p01,p02,1,1,i,M)*(i/(M+1))

 mu2.c2=mu2.c2+N*pbv(1,2,p01,p02,psi1,psi2,i,M)*(1-cond.prob(1,2,psi1,psi2,i,M))-
N*pbv(1,2,p01,p02,1,1,i,M)*((M-i+1)/(M+1))+N*pbv(2,0,p01,p02,psi1,psi2,i,M)*cond.prob(2,0,psi1,psi2,i,M)-
N*pbv(2,0,p01,p02,1,1,i,M)*(i/(M+1))
 }
 d11=d11.c1+d11.c2; d12=d12.c1+d12.c2; d22=d22.c1+d22.c2
 d11h1=d11h1.c1+d11h1.c2; d12h1=d12h1.c1+d12h1.c2; d22h1=d22h1.c1+d22h1.c2
 mu1=mu1.c1+mu1.c2
 mu2=mu2.c1+mu2.c2
 a=solve(matrix(c(d11,d12,d12,d22),nrow=2))
 ny=eigen(a,symmetric=TRUE, only.values=FALSE)
 mu=c(mu1,mu2)
 dh1=matrix(c(d11h1,d12h1,d12h1,d22h1),ncol=2)
 delta=rep(0,2)
# delta[1]=((ny$vectors[,1]%*%mu)^2)*(ny$vectors[,1]%*%dh1%*%ny$vectors[,1])
# delta[2]=((ny$vectors[,2]%*%mu)^2)*(ny$vectors[,2]%*%dh1%*%ny$vectors[,2])
delta[1]=((ny$vectors[,1]%*%mu)^2)
delta[2]=((ny$vectors[,2]%*%mu)^2)
const=rep(0,2)
const[1]=(ny$vectors[,1]%*%dh1%*%ny$vectors[,1])
const[2]=(ny$vectors[,2]%*%dh1%*%ny$vectors[,2])
critical.value=5.9914
nu=max(1,(2*sum(ny$values*const*(1+delta))-1)/sum((ny$values^2)*(const^2)*(1+2*delta)))
#nu=(1+2*max(sum(ny$values*(1+delta))-1,0)/sum((ny$values^2)*(1+2*delta)))
ncp=max(nu*(sum(ny$values*const*(1+delta))-1),0)
power=pchisq((critical.value*nu),nu,ncp,lower.tail=FALSE)
return(power)}
########################################################
sample.size=function(p01, p02, psi1, psi2, M){
# set an initial value of N=10
N=10
while(power(p01, p02, psi1, psi2, N, M)<0.80){ N=N+1}
sample.size=N
return(sample.size)}

