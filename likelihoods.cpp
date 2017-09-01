#include <Rcpp.h>
#include <vector>
#include "stan/math.hpp"
#include "stan/math/fwd/mat.hpp"
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]

List squared_loss(NumericVector mu_, NumericVector x) {
  using namespace Eigen;
  using namespace stan;
  int N = mu_.size();
  std::vector<math::fvar<math::var> > mu(N);
  NumericVector jacobian(N);
  NumericMatrix hessian(N, N);
  
  List retval;
  
  for(size_t j = 0; j < N; j++) {
    math::fvar<math::var> sum(0.0);
    std::vector<double> tmp(N);
    for(size_t i = 0; i < N; i++) {
      mu[i] = mu_[i];
      
      if(i == j)
        mu[i].d_ = 1.0;
      else
        mu[i].d_ = 0.0;
    }
  
    for(size_t i = 0; i < N; i++) {
      sum += (mu[i] - x[i]) * (mu[i] - x[i]);
    }
    
    retval["nlp"] = sum.val().val();
    jacobian[j] = sum.tangent().val();
    math::set_zero_all_adjoints();
    sum.tangent().grad();
    for(size_t i = 0; i < N; i++) {
      hessian(j, i) = mu[i].val().adj();
    }
  }
  
  retval["jacobian"] = jacobian;
  retval["hessian"] = hessian;
  
  return retval;
}