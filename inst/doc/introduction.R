## ----include=FALSE------------------------------------------------------------
library(SA23204189)

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector Cauchy_paraest_mu(NumericVector x, int a) {
#     int mu_hat;
#  
#     if(a==1){
#       mu_hat = quantile(x,0.5);
#     }
#     else if(a==2){
#       mu_hat = 0.5*(quantile(x,0.25)+quantile(x,0.75));
#     }
#  
#     return mu_hat;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector Cauchy_paraest_lambda(NumericVector x, int b) {
#     int lambda_hat;
#  
#     if(b==1){
#       lambda_hat = quantile(x,0.75)-quantile(x,0.5);
#     }
#     else if(b==2){
#       lambda_hat = quantile(x,0.5)-quantile(x,0.25);
#     }
#     else if(b==3){
#       lambda_hat = 0.5*(quantile(x,0.75)-quantile(x,0.25));
#     }
#  
#     return lambda_hat;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Cauchy_paraest_IME_mu<-function(x,lambda){
#    f<-function(mu){
#      sum(atan((x-mu)/lambda))
#    }
#  
#    xmin<-min(x)
#    xmax<-max(x)
#    uniroot(f,c(xmin,xmax))$root
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Cauchy_intervalparaest_IME_mu<-function(x,lambda,alpha){
#    f<-function(mu){
#      -2*sum(log(0.5+1/pi*atan((x-mu)/lambda)))
#    }
#  
#    xmin<-min(x)
#    xmax<-max(x)
#    n<-length(x)
#  
#    f1<-function(mu){
#      f(mu)-qchisq(alpha/2,2*n)
#    }
#  
#    f2<-function(mu){
#      f(mu)-qchisq(1-alpha/2,2*n)
#    }
#  
#    mu1_hat<-uniroot(f1,c(xmin,xmax))$root
#    mu2_hat<-uniroot(f2,c(xmin,xmax))$root
#    c(mu1_hat,mu2_hat)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Cauchy_paraest_IME_lambda<-function(x,mu,C){
#    y<-as.numeric(abs(x-mu)<=C)
#    stopifnot(sum(y)>0&&sum(y)<length(x))
#    n<-length(x)
#    C/tan(pi*sum(y)/(2*n))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Cauchy_intervalparaest_IME_lambda<-function(x,mu,C,alpha1,alpha2){
#    y<-as.numeric(abs(x-mu)<=C)
#    stopifnot(sum(y)>0&&sum(y)<length(x))
#  
#    n<-length(x)
#    y0<-sum(y)
#  
#    p1<-1/(1+(n-y0+1)/y0*qf(1-alpha1,2*(n-y0+1),2*y0))
#    p2<-1-1/(1+(y0+1)/(n-y0)*qf(1-alpha2,2*(y0+1),2*(n-y0)))
#  
#    lambda1_hat<-C/tan(pi*p2/2)
#    lambda2_hat<-C/tan(pi*p1/2)
#    c(lambda1_hat,lambda2_hat)
#  }

## -----------------------------------------------------------------------------
set.seed(5)
x1<-rcauchy(10)
x2<-rcauchy(50)
x3<-rcauchy(200)
x4<-rcauchy(1000)

mu<-0
lambda<-1
C<-1.5
alpha<-0.05
alpha1<-alpha2<-0.025

order_mu_est<-c(Cauchy_paraest_mu(sort(x1),2),Cauchy_paraest_mu(sort(x2),2),Cauchy_paraest_mu(sort(x3),2),Cauchy_paraest_mu(sort(x4),2))

order_lambda_est<-c(Cauchy_paraest_lambda(sort(x1),3),Cauchy_paraest_lambda(sort(x2),3),Cauchy_paraest_lambda(sort(x3),3),Cauchy_paraest_lambda(sort(x4),3))

IME_mu_est<-c(Cauchy_paraest_IME_mu(x1,lambda),Cauchy_paraest_IME_mu(x2,lambda),Cauchy_paraest_IME_mu(x3,lambda),Cauchy_paraest_IME_mu(x4,lambda))

IME_lambda_est<-c(Cauchy_paraest_IME_lambda(x1,mu,C),Cauchy_paraest_IME_lambda(x2,mu,C),Cauchy_paraest_IME_lambda(x3,mu,C),Cauchy_paraest_IME_lambda(x4,mu,C))

IME_mu_intervalest<-matrix(c(Cauchy_intervalparaest_IME_mu(x1,lambda,alpha),Cauchy_intervalparaest_IME_mu(x2,lambda,alpha),Cauchy_intervalparaest_IME_mu(x3,lambda,alpha),Cauchy_intervalparaest_IME_mu(x4,lambda,alpha)),nrow=2)

IME_lambda_intervalest<-matrix(c(Cauchy_intervalparaest_IME_lambda(x1,mu,C,alpha1,alpha2),Cauchy_intervalparaest_IME_lambda(x2,mu,C,alpha1,alpha2),Cauchy_intervalparaest_IME_lambda(x3,mu,C,alpha1,alpha2),Cauchy_intervalparaest_IME_lambda(x4,mu,C,alpha1,alpha2)),nrow=2)

## -----------------------------------------------------------------------------
order_mu_est

## -----------------------------------------------------------------------------
IME_mu_est

## -----------------------------------------------------------------------------
IME_mu_intervalest

## -----------------------------------------------------------------------------
order_lambda_est

## -----------------------------------------------------------------------------
IME_lambda_est

## -----------------------------------------------------------------------------
IME_lambda_intervalest

