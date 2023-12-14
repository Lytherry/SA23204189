#' @title Predict the positional parameter of a Cauchy distribution using R.
#' @description Inverse moment point estimation of the positional parameter of a Cauchy distribution if the scale parameter is known.
#' @param x the i.i.d Cauchy variables (numeric)
#' @param lambda the scale parameter of a Cauchy distribution (numeric)
#' @return the point estimation of the positional parameter \code{n}
#' @examples
#' \dontrun{
#' x<-rcauchy(50,0,1)
#' lambda<-1
#' Cauchy_paraest_IME_mu(x,lambda)
#' }
#' @importFrom stats uniroot
#' @export
Cauchy_paraest_IME_mu<-function(x,lambda){
  f<-function(mu){
    sum(atan((x-mu)/lambda))
  }
  
  xmin<-min(x)
  xmax<-max(x)
  uniroot(f,c(2*xmin,2*xmax))$root
}

#' @title Predict the scale parameter of a Cauchy distribution using R.
#' @description Inverse moment point estimation of the scale parameter of a Cauchy distribution if the positional parameter is known.
#' @param x the i.i.d Cauchy variables (numeric)
#' @param mu the positional parameter of a Cauchy distribution (numeric)
#' @param C a suitable positive constant
#' @return the point estimation of the scale parameter \code{n}
#' @examples
#' \dontrun{
#' x<-rcauchy(50,0,1)
#' mu<-0
#' C<-1
#' Cauchy_paraest_IME_lambda(x,mu,C)
#' }
#' @export
Cauchy_paraest_IME_lambda<-function(x,mu,C){
  y<-as.numeric(abs(x-mu)<=C)
  stopifnot(sum(y)>0&&sum(y)<length(x))
  n<-length(x)
  C/tan(pi*sum(y)/(2*n))
}

#' @title Interval estimation of the positional parameter of a Cauchy distribution using R.
#' @description Inverse moment interval estimation of the positional parameter of a Cauchy distribution if the scale parameter is known.
#' @param x the i.i.d Cauchy variables (numeric)
#' @param lambda the scale parameter of a Cauchy distribution (numeric)
#' @param alpha the confidence coefficient of the confidence interval (numeric)
#' @return the interval estimation of the positional parameter \code{n}
#' @examples
#' \dontrun{
#' x<-rcauchy(50,0,1)
#' lambda<-1
#' alpha<-0.05
#' Cauchy_intervalparaest_IME_mu(x,lambda,alpha)
#' }
#' @importFrom stats qchisq
#' @export
Cauchy_intervalparaest_IME_mu<-function(x,lambda,alpha){
  f<-function(mu){
    -2*sum(log(0.5+1/pi*atan((x-mu)/lambda)))
  }
  
  xmin<-min(x)
  xmax<-max(x)
  n<-length(x)
  
  f1<-function(mu){
    f(mu)-qchisq(alpha/2,2*n)
  }
  
  f2<-function(mu){
    f(mu)-qchisq(1-alpha/2,2*n)
  }
  
  mu1_hat<-uniroot(f1,c(2*xmin,2*xmax))$root
  mu2_hat<-uniroot(f2,c(2*xmin,2*xmax))$root
  c(mu1_hat,mu2_hat)
}

#' @title Interval estimation of the positional parameter of a Cauchy distribution using R.
#' @description Inverse moment interval estimation of the positional parameter of a Cauchy distribution if the scale parameter is known.
#' @param x the i.i.d Cauchy variables (numeric)
#' @param mu the positional parameter of a Cauchy distribution (numeric)
#' @param C a suitable positive constant
#' @param alpha1 alpha1+alpha2=alpha, in which alpha is the confidence coefficient of the confidence interval (numeric)
#' @param alpha2 alpha1+alpha2=alpha, in which alpha is the confidence coefficient of the confidence interval (numeric)
#' @return the interval estimation of the scale parameter \code{n}
#' @examples
#' \dontrun{
#' x<-rcauchy(50,0,1)
#' mu<-0
#' C<-1
#' alpha1<-0.025
#' alpha2<-0.025
#' Cauchy_intervalparaest_IME_lambda(x,mu,C,alpha1,alpha2)
#' }
#' @importFrom stats qf
#' @export
Cauchy_intervalparaest_IME_lambda<-function(x,mu,C,alpha1,alpha2){
  y<-as.numeric(abs(x-mu)<=C)
  stopifnot(sum(y)>0&&sum(y)<length(x))
  
  n<-length(x)
  y0<-sum(y)
  
  p1<-1/(1+(n-y0+1)/y0*qf(1-alpha1,2*(n-y0+1),2*y0))
  p2<-1-1/(1+(y0+1)/(n-y0)*qf(1-alpha2,2*(y0+1),2*(n-y0)))
  
  lambda1_hat<-C/tan(pi*p2/2)
  lambda2_hat<-C/tan(pi*p1/2)
  c(lambda1_hat,lambda2_hat)
}

