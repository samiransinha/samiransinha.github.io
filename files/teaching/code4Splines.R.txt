newbs=function(x, degree, inner.knots, Boundary.knots) {
Boundary.knots=sort(Boundary.knots);
knots=c(rep(Boundary.knots[1], (degree+1)), sort(inner.knots),
rep(Boundary.knots[2], (degree+1)));
np=degree+length(inner.knots)+1
s=rep(0, np)

if(x==Boundary.knots[2]) {s[np]=1} else {for( i in 1: np)s[i]=basis(x, degree, i, knots)}

return(s)}
basis=function(x, degree, i, knots)
{ if(degree==0){ if((x<knots[i+1])&(x>=knots[i])) y=1 else
y=0}else{
if((knots[degree+i]-knots[i])==0) {temp1=0} else {temp1=
(x-knots[i])/(knots[degree+i]-knots[i])};
if((knots[i+degree+1]-knots[i+1])==0) {temp2=0} else {temp2=
(knots[i+degree+1]-x)/(knots[i+degree+1]-knots[i+1])}
y= temp1*basis(x, (degree-1), i, knots) +temp2*basis(x, (degree-1),
(i+1), knots)}
return(y)}

 newbs(2, degree=3, inner.knots=c(-0.25, -0.5, 0, 0.25, 0.5),
 Boundary.knots=c(-4, 4))



#
#
#
#
#
#
#
#
# Use of the Fortran code
dyn.load("spline.so")
 n=1;
 m1=3;
 d=2;
 m=m1+2*(d+1);
 knots=c(-0.25, 0.0, 0.25)
 boundaryknots=c(-3, 3)
 x=rnorm(n)
 k=d+m1+1;

 basis=matrix(0, nrow=n, ncol=k);
 storage.mode(basis)<-"double"
 f1=.Fortran("splinebasis", d=as.integer(d), n=as.integer(n),
 m=as.integer(m), m1=as.integer(m1),k=as.integer(k),x=as.double(x), 
knots=as.double(knots), boundaryknots=as.double(boundaryknots), output=basis)
