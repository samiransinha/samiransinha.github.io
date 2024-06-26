Hap.mathced.case.control<-function(cc.data, ll, m, nbet, nref, p, nstr, tol, maxit){
if(missing(tol)) tol=0.00001
if(missing(maxit)) maxit=200
dyn.load("new-hapdataanalysissub1.so")
library(survival)
####### Input data##############
# Note missing geneotypes are replaced by 9
# Following program is for an additive model
# the value of nstr must be provided by the user
# the value of m must be provided by the user 
# here m=number of controls for each case
# cc.data=is the dataset
################################
################################
# This calculates possible haplotypes
# The value of p must be provided 
 nr=nrow(cc.data); 
 print('nr= '); print(nr)
 numofsnp=p
 numofhap=2^p

new.h=matrix(0, nrow=numofhap, ncol=numofsnp)
for( i in 1: numofhap){
j=i-1
for(k in 1: numofsnp){
k1=numofsnp-k+1
new.h[i, k1]=as.integer(j/2^(k1-1))
j=j-new.h[i, k1]*2^(k1-1)
               }
                }
#----------------------------------------------------------
#------ This part determines possible haplotypes for ------
#------ unphased genotypes --------------------------------
uvecdip=NULL; ncount=NULL;
# m=number of controls
# ll=column number of the first SNP must be provided by the user
ul=ll+p-1 # column number of the second SNP must be provided by the user
const=2*p+1 
for( it in 1: nrow(cc.data)){
if( sum(cc.data[it, ll:ul])>const){
a1=cc.data[it, ll:ul];
a2=(1:p)[a1!=9];
n.dip=NULL; 
for( i1 in 1: numofhap){
    for( j1 in i1 : numofhap){
if(new.h[i1, a2]+new.h[j1, a2]==a1[a2]) {t1=new.h[i1, ]%*%2^((1:p)-1)+1; 
t2=new.h[j1, ]%*%2^((1:p)-1)+1; n.dip=rbind(n.dip, c(t1, t2))}
                      } 
                 }
uvecdip=rbind(uvecdip, n.dip); ncount=append(ncount, nrow(n.dip))
}  else{new.g1=cc.data[it, ll:ul ]/2
vec=c(1:2)[new.g1!=.5]
veclen=length(vec)
set=NULL
for( i in 1: numofhap){
if(identical(as.numeric(new.h[i, vec]), as.numeric(new.g1[vec]))) set=rbind(set, new.h[i,] )
};ns=nrow(set); n.dip=NULL;
t=2^((1:p)-1)
for( i in 1: ns){
 for( j in i: ns){
 if(identical(as.numeric(set[i, ]+set[j, ]), as.numeric(cc.data[it, ll:ul ])))
 {t1=t%*%set[i, ]+1; t2=t%*%set[j, ]+1;  n.dip=rbind(n.dip, c(t1, t2))}
                 } 
                }
uvecdip=rbind(uvecdip, n.dip); ncount=append(ncount, nrow(n.dip))
}}
# -----------------------

if(nref!=numofhap){
uvecdip[uvecdip==numofhap]<-100
uvecdip[uvecdip==nref]<-numofhap
uvecdip[uvecdip==100]<-nref
                  }


ndip=nrow(uvecdip);
beta=rep(0.1, nbet)
ncnt=cumsum(ncount)
storage.mode(uvecdip)<-"integer"
# --------------- EM based methods ----------
upi=rep((1/numofhap), numofhap); upi=as.double(upi);
gn1=.Fortran("expem", nr=as.integer(nr), ns=as.integer(nstr), 
ncnt=as.integer(ncnt), uvecdip, ndip=as.integer(ndip), output=upi, 
mupi=as.integer(numofhap), numofhap=as.integer(numofhap), m=as.integer(m))
upi=gn1$output
loglk=as.double(0); upi=as.double(gn1$output);




ff1=function(beta){
fn1=.Fortran("funcbetacc", nr=as.integer(nr), nstr=as.integer(nstr), 
   ndip=as.integer(ndip), uvecdip, ncnt=as.integer(ncnt), 
    beta=as.double(beta), pi=as.double(upi), output=loglk,
 nref=as.integer(numofhap), nbet=as.integer(nbet), mupi=as.integer(numofhap), 
m=as.integer(m)) 
fn1$output
}
b=nlm(ff1, beta)
beta=as.numeric(b$estimate)
#--------- Proposed method -------------------
#----------ECM -------------------------------
nonzeroupi=(1:numofhap)[upi!=0]
nl=length(nonzeroupi)
#---------- Function1 ------------------------
crs=as.double(rep(0, sum(ncount)))
hn2=function(beta, upi){
fn2=.Fortran("estepofem", nr=as.integer(nr), nstr=as.integer(nstr), 
      ndip=as.integer(ndip), ncnt=as.integer(ncnt), uvecdip, 
output=crs, upi=as.double(upi), beta=as.double(beta),
 mupi=as.integer(numofhap), nbet=as.integer(nbet), nref=as.integer(numofhap),
m=as.integer(m))
fn2$output
}
crs=hn2(beta, upi)
#--------- Function2 ----------------------
loglk=as.double(0)
hn1=function(beta){
hn11=.Fortran("funcbeta", nr=as.integer(nr), nstr=as.integer(nstr), 
ndip=as.integer(ndip), crs=as.double(crs), uvecdip, ncnt=as.integer(ncnt),
 beta=as.double(beta), upi=as.double(upi), output=loglk, mpi=as.integer(numofhap), 
nbet=as.integer(nbet), nref=as.integer(numofhap), m=as.integer(m)) 
hn11$output
}
#--------- Beginning of the E-step ---------
eps=5.0
rpeat=0
while((eps>tol) &(rpeat<maxit)){
rpeat=rpeat+1
old_beta=beta
crs=hn2(beta, upi)
#--------- Maximization of beta  --------------
out=nlm(hn1, old_beta)
beta=out$estimate
#--------- Determination of upi ----------------
old_pi=upi
new_upi=as.double(rep(0, numofhap)); 
hn4=.Fortran("determinepi", nr=as.integer(nr), nstr=as.integer(nstr), 
 ncnt=as.integer(ncnt), nl=as.integer(nl), 
nonzeroupi=as.integer(nonzeroupi), pcrs=as.double(crs),  
    ndip=as.integer(ndip), uvecdip, beta=as.double(beta), 
upi=as.double(upi),  output=new_upi, mupi=as.integer(numofhap), 
nbet=as.integer(nbet), nref=as.integer(numofhap), m=as.integer(m)) 
upi=hn4$output
sm=(1:numofhap)[old_pi!=0]
temp.vec=old_pi[sm]
temp.vec[temp.vec<0.08]<-0.08
eps=sum(abs((old_beta-beta)/old_beta))+sum(abs((old_pi[sm]-upi[sm])/temp.vec))# checking convergence
#print('beta= '); print(beta)
}
out.beta=beta;
out.upi=upi;
out.eps=eps;
out.rpeat=rpeat;

print("Please wait while we calculate the standard error......")
# -----------------
oldncount=ncount
olduvecdip=uvecdip


method=NULL;
numofrep=as.integer(nstr/2);
tempo.1=as.integer(nstr/numofrep)
cc.org=cc.data
for(urep in 1: numofrep){
sv1=tempo.1*(m+1)*(urep-1)+1
sv2=sv1+tempo.1*(m+1)-1

cc.data=cc.org[-c(sv1:sv2), ]
new.nstr=nrow(cc.data)/(m+1)
ncount=oldncount[-c(sv1:sv2)]
balti=cumsum(oldncount)
if(sv1==1)glass=(1: balti[sv2]) else glass=((balti[sv1-1]+1):balti[sv2])
uvecdip=olduvecdip[-glass, ]

nr=nrow(cc.data); 
ndip=nrow(uvecdip);
beta=rep(0.1, nbet)
ncnt=cumsum(ncount)
storage.mode(uvecdip)<-"integer"
# --------------- EM based methods ----------
upi=rep((1/numofhap), numofhap); upi=as.double(upi);
gn1=.Fortran("expem", nr=as.integer(nr), ns=as.integer(new.nstr), 
ncnt=as.integer(ncnt), uvecdip, ndip=as.integer(ndip), output=upi, 
mupi=as.integer(numofhap), numofhap=as.integer(numofhap), m=as.integer(m))
upi=gn1$output
loglk=as.double(0); upi=as.double(gn1$output);




ff1=function(beta){
fn1=.Fortran("funcbetacc", nr=as.integer(nr), nstr=as.integer(new.nstr), 
   ndip=as.integer(ndip), uvecdip, ncnt=as.integer(ncnt), 
    beta=as.double(beta), pi=as.double(upi), output=loglk,
 nref=as.integer(numofhap), nbet=as.integer(nbet), mupi=as.integer(numofhap), 
m=as.integer(m)) 
fn1$output
}
b=nlm(ff1, beta)
beta=as.numeric(b$estimate)
#--------- Proposed method -------------------
#----------ECM -------------------------------
nonzeroupi=(1:numofhap)[upi!=0]
nl=length(nonzeroupi)
#---------- Function1 ------------------------
crs=as.double(rep(0, sum(ncount)))
hn2=function(beta, upi){
fn2=.Fortran("estepofem", nr=as.integer(nr), nstr=as.integer(new.nstr), 
      ndip=as.integer(ndip), ncnt=as.integer(ncnt), uvecdip, 
output=crs, upi=as.double(upi), beta=as.double(beta),
 mupi=as.integer(numofhap), nbet=as.integer(nbet), nref=as.integer(numofhap),
m=as.integer(m))
fn2$output
}
crs=hn2(beta, upi)
#--------- Function2 ----------------------
loglk=as.double(0)
hn1=function(beta){
hn11=.Fortran("funcbeta", nr=as.integer(nr), nstr=as.integer(new.nstr), 
ndip=as.integer(ndip), crs=as.double(crs), uvecdip, ncnt=as.integer(ncnt),
 beta=as.double(beta), upi=as.double(upi), output=loglk, mpi=as.integer(numofhap), 
nbet=as.integer(nbet), nref=as.integer(numofhap), m=as.integer(m)) 
hn11$output
}
#--------- Beginning of the E-step ---------
eps=5.0
rpeat=0
while((eps>tol) &(rpeat<maxit)){
rpeat=rpeat+1
old_beta=beta
crs=hn2(beta, upi)
#--------- Maximization of beta  --------------
out=nlm(hn1, old_beta)
beta=out$estimate
#--------- Determination of upi ----------------
old_pi=upi
new_upi=as.double(rep(0, numofhap)); 
hn4=.Fortran("determinepi", nr=as.integer(nr), nstr=as.integer(new.nstr), 
 ncnt=as.integer(ncnt), nl=as.integer(nl), 
nonzeroupi=as.integer(nonzeroupi), pcrs=as.double(crs),  
    ndip=as.integer(ndip), uvecdip, beta=as.double(beta), 
upi=as.double(upi),  output=new_upi, mupi=as.integer(numofhap), 
nbet=as.integer(nbet), nref=as.integer(numofhap), m=as.integer(m)) 
upi=hn4$output
sm=(1:numofhap)[old_pi!=0]
temp.vec=old_pi[sm]
temp.vec[temp.vec<0.08]<-0.08
eps=sum(abs((old_beta-beta)/old_beta))+sum(abs((old_pi[sm]-upi[sm])/temp.vec))# checking convergence
#print('beta= '); print(beta)
}
method=rbind(method, beta)
} # end of the loop of urep 

if(nref!=numofhap){
temporary.1=out.beta[nref]
temporary.2=method[, nref]
for( i2 in nref:(nbet-1)) {out.beta[i2]=out.beta[i2+1]
method[, i2]=method[, (i2+1)]
}
out.beta[nbet]=temporary.1
method[, nbet]=temporary.2

temporary.3=out.upi[nref]
out.upi[nref]=out.upi[numofhap]
out.upi[numofhap]=temporary.3
}


cov.mat=matrix(0, nrow=nbet, ncol=nbet)
for( i in 1: nbet){
 for( j in 1: nbet){
cov.mat[i, j]=((numofrep-1)/numofrep)*sum((method[, i]-out.beta[i])*(method[, j]-out.beta[j]))
}}




output=list(out.beta, cov.mat,  out.upi,  out.eps, out.rpeat)
names(output)<-c("estimate", "cov.matrix", "Hap.freq.control", "eps", "Num.iter")
return (output)
}
