#include <Rcpp.h>
using namespace Rcpp;

//' @title Produce a quantile number using Rcpp.
//' @description Produce the quantile of order q of a numeric vector.
//' @param x a numeric vector
//' @param q the order of the quantile
//' @return the quantile of order q of vector x
//' @import Rcpp
//' @examples
//' \dontrun{
//' x<-rcauchy(5,0,1)
//' quantile(x,0.5)
//' }
//' @useDynLib SA23204189
//' @export
// [[Rcpp::export]]
double quantile0(NumericVector x, double q)
{
  assert(q >= 0.0 && q <= 1.0);
  const int n = x.size();
  double id = (n-1)*q;
  int lo = floor(id);
  int hi = ceil(id);
  double qs = x[lo];
  double h = (id-lo);
  return (1.0 - h) * qs + h * x[hi];
}

//' @title Predict the positional parameter of a Cauchy distribution using Rcpp.
//' @description Predict the positional parameter of a Cauchy distribution using Rcpp.
//' @param x the sorted i.i.d Cauchy variables (NumericVector)
//' @param a the method (1 or 2) (int)
//' @return the point estimation of the positional parameter
//' @examples
//' \dontrun{
//' x<-rcauchy(50,0,1)
//' Cauchy_paraest_mu(x,1)
//' }
//' @useDynLib SA23204189
//' @export
// [[Rcpp::export]]
double Cauchy_paraest_mu(NumericVector x, int a) {
   double mu_hat;
  
   if(a==1){
     mu_hat = quantile0(x,0.5);
   }
   else if(a==2){
     mu_hat = 0.5*(quantile0(x,0.25)+quantile0(x,0.75));
   }
   
   return mu_hat;
}

//' @title Predict the scale parameter of a Cauchy distribution using Rcpp.
//' @description Predict the scale parameter of a Cauchy distribution using Rcpp.
//' @param x the sorted i.i.d Cauchy variables (NumericVector)
//' @param b the method (1, 2 or 3) (int)
//' @return the point estimation of the scale parameter
//' @examples
//' \dontrun{
//' x<-rcauchy(50,0,1)
//' Cauchy_paraest_lambda(x,1)
//' }
//' @useDynLib SA23204189
//' @export
// [[Rcpp::export]]
double Cauchy_paraest_lambda(NumericVector x, int b) {
   double lambda_hat;
  
   if(b==1){
     lambda_hat = quantile0(x,0.75)-quantile0(x,0.5);
   }
   else if(b==2){
     lambda_hat = quantile0(x,0.5)-quantile0(x,0.25);
   }
   else if(b==3){
     lambda_hat = 0.5*(quantile0(x,0.75)-quantile0(x,0.25));
   }
   
   return lambda_hat;
}


