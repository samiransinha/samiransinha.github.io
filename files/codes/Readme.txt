Hap.mathced.case.control is an R function that calculates the relative risk 
parameters corresponding to the Haplotype of interest based on a matched 
case-control study. This program fits only additive model. However, with 
some simple modifications this program can also be used to fit dominant
and recessive model. 

The arguments of the function are 
cc.data: matched case-control data
ll: the column of cc.data matrix which provides geneotype information 
for SNP-I.
m: number of control for each case
nbet: number of relative risk parameters
nref: the reference haplotype, it should be the last haplotype. For example 
if we have 8 haplotypes, then the nrefshould equal to 8.
p: number of SNPs
nstr: number of strata
tol: tolerance, the default value is 0.00001
maxit: maximum number of iterations used for the ECM algorithm
If the value of maxit is not provided, then the program uses the default value.  

In the data missing genotype must be repalced by 9. 

The function returns the following values.
estimate: estimate of the relative risk parameters.
Hap.freq.control: Haplotype frequency in the control population.
eps: the achieved level of tolerance.
Num.iter: Number iterations required for the convergence.


The subroutines are written in Fortran and they are saved in 
new-hapdataanalysissub1.f. In order to run the code one first need to compile
the Fortran program with the following command in a Unix matching.

f77 -O2 -c new-hapdataanalysissub1.f

Then load the code in R's shared library  by the following command 

R CMD SHLIB new-hapdataanalysissub1.so 

Make sure to run the R code in the same directory where 
new-hapdataanalysissub1.so is stored. 





Following is an example on how to use the code. 

Example: We create a 1:2 matched case-control data with 
         2 SNPs. Suppose the number of strata is 200.   
p=2; ncase=200; m=2
numofsnp=p
numofhap=2^p
new.h=matrix(0, nrow=numofhap, ncol=numofsnp)
for( i in 1: numofhap){
j=i-1
for(k in 1: numofsnp){
k1=numofsnp-k+1
new.h[i, k1]=as.integer(j/2^(k1-1))
j=j-new.h[i, k1]*2^(k1-1)
                }    }
pr_multi=c(0.12, 0.15, 0.2,  0.53)
n=15000;
s=rep(0, n)
dip=matrix(0, ncol=2, nrow=n)
h=new.h; h.val=h%*%2^((1:p)-1)+1
vec=1:numofhap
g=matrix(0, ncol=p, nrow=n); t=rep(0, n); ind=rep(0, n); 
hapc=rep(0, n);
s=rgamma(n, 0.5, 1)
for( i in 1: n){
r1=rmultinom(1, 1, pr_multi); r2=rmultinom(1, 1, pr_multi)
dip[i, ]=c(vec[r1!=0], vec[r2!=0]);
g[i, ]=h[vec[r1!=0], ]+h[vec[r2!=0], ]; 
if(dip[i, 1]==3) hapc[i]=1
if(dip[i, 2]==3) hapc[i]=hapc[i]+1 
prob=exp(-3.2+0.51*s[i]+0.5*hapc[i])/(1+exp(-3.2+0.51*s[i]+0.5*hapc[i]))## intercept: -3.2
ind[i]=rbinom(1, 1, prob);
}
cc.data=NULL; 
i5=0; ncs=0 
while((i5<n-1)&(ncs<ncase)){
i5=i5+1
if(ind[i5]==1) {
n1=(1:n)[ind==0& abs(s-s[i5])<0.05];  
if(length(n1)>(m-1))  {ncs=ncs+1; i6=sample(n1, size=m);
cc.data=rbind(cc.data,c(g[i5, ], ind[i5]));
cc.data=rbind(cc.data,cbind(g[i6, ], ind[i6]));
                   }
                }
}


out=Hap.mathced.case.control(cc.data, ll=1, m=2, nbet=3, nref=4, p=2, nstr=200)
 out
$estimate
[1] -0.3504115 -0.1586919  0.2346554

$Hap.freq.control
[1] 0.1820996 0.1181753 0.1565780 0.5431471

$eps
[1] 0.002852745

$Num.iter
[1] 200

                    





