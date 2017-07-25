#include <Rcpp.h>
#include <random>
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

// [[Rcpp::export]]
double solo(NumericMatrix x) {
  double total = 0.0;
  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < x.ncol(); j++) {
      total += x(i, j);
    }
  }
  return total;
}

// [[Rcpp::export]]
double pairs(NumericMatrix x) {
  double total = 0.0;
  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < x.ncol(); j++) {
      if(i < x.nrow() - 1) {
        total += x(i, j) * x(i + 1, j);
      }
          
      if(j < x.ncol() - 1) {
        total += x(i, j) * x(i, j + 1);
      }
    }
  }
  return total;
}

// [[Rcpp::export]]
List ising(NumericMatrix x, double mu, double beta, double kT, int S) {
  NumericMatrix out(S, 2);
  colnames(out) = CharacterVector::create("solo", "pairs");
  
  int N = x.nrow();
  int accept = 0;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  std::uniform_int_distribution<> dis(0, N - 1);
  double solo_ = solo(x);
  double pairs_ = pairs(x);
  double E = mu * solo_ + beta * pairs_;
  for(int s = 0; s < S; s++) {
    int i = dis(gen);
    int j = dis(gen);
    
    double dsolo_ = -2 * x(i, j);
    double dpairs_ = ((i < x.nrow() - 1) ? -2.0 * x(i, j) * x(i + 1, j) : 0) +
      ((i > 0) ? -2.0 * x(i - 1, j) * x(i, j) : 0) +
      ((j < x.ncol() - 1) ? -2.0 * x(i, j) * x(i, j + 1) : 0) +
      ((j > 0) ? -2.0 * x(i, j - 1) * x(i, j) : 0);
     
    double dE = mu * dsolo_ + beta * dpairs_;
    double r = runif(gen);
    
    if(dE < 0.0 | r < std::exp(-dE / kT)) {
      solo_ = solo_ + dsolo_;
      pairs_ = pairs_ + dpairs_;
      x(i, j) = x(i, j) * -1;
      accept += 1;
    }
    
    out(s, 0) = solo_;
    out(s, 1) = pairs_;
  }
  
  List ret;
  ret["x"] = x;
  ret["states"] = out;
  ret["accept"] = accept / double(S);
  return ret;
}

//// [[Rcpp::export]]
/*List ising_gibbs(NumericMatrix x, double mu, double beta, double kT, int S) {
  NumericMatrix out(S, 2);
  colnames(out) = CharacterVector::create("solo", "pairs");
  
  int N = x.nrow();
  int accept = 0;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  std::uniform_int_distribution<> dis(0, N - 1);
  double solo_ = solo(x);
  double pairs_ = pairs(x);
  double E = mu * solo_ + beta * pairs_;
  for(int s = 0; s < S / (N * N); s++) {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        double dsolo_ = -2 * x(i, j);
        double dpairs_ = ((i < x.nrow() - 1) ? -2.0 * x(i, j) * x(i + 1, j) : 0) +
          ((i > 0) ? -2.0 * x(i - 1, j) * x(i, j) : 0) +
          ((j < x.ncol() - 1) ? -2.0 * x(i, j) * x(i, j + 1) : 0) +
          ((j > 0) ? -2.0 * x(i, j - 1) * x(i, j) : 0);
        
        double dE = mu * dsolo_ + beta * dpairs_;
        double r = runif(gen);
        
        if(dE < 0.0 | r < std::exp(-dE / kT)) {
          solo_ = solo_ + dsolo_;
          pairs_ = pairs_ + dpairs_;
          E = mu * solo_ + beta * pairs_;
          x(i, j) = x(i, j) * -1;
          accept += 1;
        }
        
        out(s, 0) = solo_;
        out(s, 1) = pairs_;
      }
    }
  }
  
  List ret;
  ret["x"] = x;
  ret["states"] = out;
  ret["accept"] = accept / double(S);
  return ret;
}*/
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
