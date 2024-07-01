The code is for calculating bias corrected estimates of 
the 
parameter based on the modified score approach (Sun et al., 2008).



The following abbreviations have used in the subroutines.


Input variables:

y:  response variable which takes only 0 or 1 (depending on case or control)
covariates: values of the covariate (it can take maximum two covariates)
ncase: number of cases in the dataset
M:  the number of controls in each stratum
eps: tolerance 
beta0: the initial value of the parameter

 

Output variables: 

estimae: parameter estimate
std: standard error of the parameter





Reference:
 
Sun, J., Sinha, S., Wang, S. J., and Maity, T. (2008).
 Bias corrected 
inference for the conditional logistic regression. Submitted.  


Example: 
# Suppose we have one 1:3 matched case-control data with two covariates, 
# and given  as follows. The number of strata is n=29. 

>data
y	cov1	cov2
1	1.3	0
0	1.12	0
0	1.35	0
0	0.95	0
1	1.3	1
0	1.03	0
0	1.22	0
0	1.13	0
1	1.2	0
0	1.13	0
0	1.19	0
0	1.19	0
1	1.48	0
0	1	0
0	0.9	0
0	2.29	0
1	1.1	1
0	1.07	0
0	1	0
0	0.9	0
1	0.91	1
0	1.38	0
0	1.89	0
0	1.47	0
1	1.02	0
0	1.5	0
0	2.35	0
0	1.84	0
1	1.12	0
0	1.82	0
0	0.95	0
0	1.32	0
1	1.5	0
0	1.2	0
0	1.05	0
0	1.41	1
1	1.2	0
0	1.03	0
0	1.27	0
0	1.7	0
1	1.21	1
0	1.69	1
0	1.21	0
0	1.2	0
1	2	0
0	1.08	0
0	1.24	0
0	1.85	0
1	1	1
0	1.6	0
0	1.1	0
0	1.15	0
1	1.3	1
0	0.85	0
0	1.3	0
0	1.25	0
1	1.3	0
0	1.2	0
0	1.12	1
0	1.69	0
1	0.97	0
0	1.3	0
0	1.19	0
0	1.23	0
1	1.1	1
0	1.28	0
0	1.9	0
0	1.1	0
1	1.32	0
0	1.15	0
0	1.15	0
0	1.1	0
1	1.38	0
0	0.9	1
0	1.33	0
0	1.16	0
1	0.85	0
0	1.18	0
0	1.25	0
0	1.2	0
1	0.92	0
0	1.2	0
0	1.4	0
0	2.41	0
1	1.05	1
0	1.55	0
0	0.95	1
0	1.3	0
1	1.9	0
0	1.13	0
0	1.68	0
0	1.6	0
1	1.2	1
0	1.4	0
0	2.5	0
0	1.34	0
1	0.95	0
0	1.2	0
0	1.2	0
0	1.3	0
1	1.3	0
0	1.5	0
0	1.35	0
0	1.3	0
1	1.42	1
0	1.07	1
0	1.53	0
0	1.37	0
1	1.02	1
0	1	0
0	1.5	0
0	1.2	0
1	1.05	0
0	1.21	0
0	1.32	0
0	1.34	1


>y <- data[,1]
>covariate <- data[,c(2,3)]
>out <- MDS(y, covariate, ncase=29, M=3, eps=0.0001, beta0=c(0.2,0.4))

> out
       betaP       sdP
1 -0.6653633 0.7779791
2  1.7590951 0.5502278
