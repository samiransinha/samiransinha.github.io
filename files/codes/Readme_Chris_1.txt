### Code for analyzing simulated datasets 
#### Data simulation 

nt=1
print(nt) 
set.seed(nt) 
####
#### Simulation of a population 
capn=80000 
capm=2; # the number of controls to be matched with each case
n=1000 # the number of cases 
####
#### Simulate variables used in the analysis
sv=runif(capn, -1, 1) ## stratification (confounding) variable 
z=rnorm(capn, 0, 0.5) ##     
#### simulation of the instrumental variable xstar
xstar=rnorm(capn, 0, 0.5)
#### simulation of the actual treatment/exposure variable 
#### True treatment variable X
#### Set probability for exposure for strong association
probx=1/(1+exp(-(-1+2*xstar+1*(sv)))) 
x=rbinom(capn, 1, probx)

## w is the mismeasured form of x
# Misclassification alpha1=alpha0=0.2
capbw=rbinom(capn, 1, 0.8) #capbw=pr(W=1|X=0) =1-alpha1
capbws=rbinom(capn, 1, 0.8)#capbws=pr(W=0|X=1)=1-alpha0
w=capbw*x+(1-capbws)*(1-x)

## simulation of y 
proby=1/(1+exp(-(-2-2*(sv)+1*x+0.5*z)))
y=rbinom(capn,1,proby)

## Combine all the variables of the population 
mydata=cbind(sv, x, xstar, w, y, z)

############# Forming matched data ####
## First randomly select the id of cases from the data
index=rep(0, ((capm+1)*n))
idc=which(y==1)
samplec=sample(idc, n, replace=F)

### Next, for each case, find a matching control by comparing difference 
### of the stratification variable for case and control, using tolerances of 
### 0.01, 0.05, or 0.1, whichever first gives non-empty set with cardinality 
### capm 

for( i1 in 1: n){
k1=(capm+1)*(i1-1)+1
ids=which(abs(sv-sv[samplec[i1]])<0.01 & y==0)
if(length(ids)>capm) {temp1=sample(ids, capm, replace=F)} else {
 ids=which(abs(sv-sv[samplec[i1]])<0.05 & y==0) 
 if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
ids=which(abs(sv-sv[samplec[i1]])<0.1 & y==0) 
 temp1=sample(ids, capm, replace=F)
}
} 
index[k1]=samplec[i1]
index[k1+(1:capm)]=temp1
} 
#### Presenting the matched data as a matrix where  
#### the columns of the matched data 
####  are sv, x, xstar, w, y, z
matcdata=mydata[index, ] 
matcdata=data.frame(matcdata) 
### End of the simulation of a matched data
###
###
#### Beginning of the analysis 
####
####  Install necessary packages for the analysis
require(MASS)
require(survival)
require(numDeriv)
#### Before running the following code, the function msw.simu function needs to be run 
####
mynout=msw.simu(matched.case.control=matcdata, no.cases=n, no.controls.per.set=capm)

#### The above function(msw.simu) calculates the parameter estimates using 
#### true X and then surrogate W in the conditional logistic regression model,
#### then and using the proposed method.
#### Input variables: matched case control data, the number of cases (no.cases), and 
#### the number of controls in each matched set (no.controls.per.set)

#### There will be three main outputs, M1, M2, and M3 for three different methods. Under 
#### each method, parameter estimates (Est) and the corresponding standard errors (SE)
#### are returned.
#
#
