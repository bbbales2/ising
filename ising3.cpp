#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double solo(const NumericMatrix &x) {
  double total = 0.0;
  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < x.ncol(); j++) {
      total += x(i, j);
    }
  }
  return total;
}

// [[Rcpp::export]]
double pairs(const NumericMatrix &x, int oi, int oj) {
  double total = 0.0;
  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < x.ncol(); j++) {
      total += x(i, j) * x((i + oi + x.nrow()) % x.nrow(), (j + oj + x.ncol()) % x.ncol());
      total += x(i, j) * x((i + oj + x.nrow()) % x.nrow(), (j - oi + x.ncol()) % x.ncol());
    }
  }
  
  if(abs(oi) != abs(oj)) {
    for(int i = 0; i < x.nrow(); i++) {
      for(int j = 0; j < x.ncol(); j++) {
        total += x(i, j) * x((i + oi + x.nrow()) % x.nrow(), (j - oj + x.ncol()) % x.ncol());
        total += x(i, j) * x((i - oj + x.nrow()) % x.nrow(), (j - oi + x.ncol()) % x.ncol());
      }
    }
  }
  
  return total;
}

// There's twice as many additions here as in the pairs function cause the pairs functiona
// needs to avoid double counting and this needs to do that counting
// [[Rcpp::export]]
double dpairs(const NumericMatrix &x, int i, int j, int oi, int oj) {
  double total = 0.0;
  total += x(i, j) * x((i + oi + x.nrow()) % x.nrow(), (j + oj + x.ncol()) % x.ncol());
  total += x(i, j) * x((i + oi + x.nrow()) % x.nrow(), (j - oj + x.ncol()) % x.ncol());
  total += x(i, j) * x((i - oi + x.nrow()) % x.nrow(), (j + oj + x.ncol()) % x.ncol());
  total += x(i, j) * x((i - oi + x.nrow()) % x.nrow(), (j - oj + x.ncol()) % x.ncol());

  if(abs(oi) != abs(oj)) {
    total += x(i, j) * x((i - oj + x.nrow()) % x.nrow(), (j + oi + x.ncol()) % x.ncol());
    total += x(i, j) * x((i - oj + x.nrow()) % x.nrow(), (j - oi + x.ncol()) % x.ncol());
    total += x(i, j) * x((i + oj + x.nrow()) % x.nrow(), (j + oi + x.ncol()) % x.ncol());
    total += x(i, j) * x((i + oj + x.nrow()) % x.nrow(), (j - oi + x.ncol()) % x.ncol());
  }

  return total;
}

// [[Rcpp::export]]
double triplets(const NumericMatrix &x, int w) {
  double total = 0.0;
  
  if(w == 0) {
    for(int i = 0; i < x.nrow(); i++) {
      for(int j = 0; j < x.ncol(); j++) {
        // ##
        //  x
        total += x(i, j) *
          x((i + 1) % x.nrow(), (j - 1 + x.ncol()) % x.ncol()) *
          x((i + 1) % x.nrow(), j);
        //   #
        // x #
        total += x(i, j) *
          x(i, (j + 1) % x.ncol()) *
          x((i + 1) % x.nrow(), (j + 1) % x.ncol());
        // x
        // # #
        total += x(i, j) *
          x((i - 1 + x.nrow()) % x.nrow(), (j + 1) % x.ncol()) *
          x((i - 1 + x.nrow()) % x.nrow(), j);
        // # x
        // #  
        total += x(i, j) *
          x(i, (j - 1 + x.ncol()) % x.ncol()) *
          x((i - 1 + x.nrow()) % x.nrow(), (j - 1 + x.ncol()) % x.ncol());
      }
    }
  }
  
  return total;
}

// [[Rcpp::export]]
double dtriplets(const NumericMatrix &x, int i, int j, int w) {
  double total = 0.0;

  if(w == 0) {
    // # #
    //   x
    total += x(i, j) *
      x((i + 1) % x.nrow(), (j - 1 + x.ncol()) % x.ncol()) *
      x((i + 1) % x.nrow(), j);
    //   #
    // x #
    total += x(i, j) *
      x(i, (j + 1) % x.ncol()) *
      x((i + 1) % x.nrow(), (j + 1) % x.ncol());
    // x
    // # #
    total += x(i, j) *
      x((i - 1 + x.nrow()) % x.nrow(), (j + 1) % x.ncol()) *
      x((i - 1 + x.nrow()) % x.nrow(), j);
    // # x 
    // #  
    total += x(i, j) *
      x(i, (j - 1 + x.ncol()) % x.ncol()) *
      x((i - 1 + x.nrow()) % x.nrow(), (j - 1 + x.ncol()) % x.ncol());
    
    // # #
    // x
    total += x(i, j) *
      x((i + 1) % x.nrow(), j) *
      x((i + 1) % x.nrow(), (j + 1) % x.ncol());
    // x #
    //   #
    total += x(i, j) *
      x(i, (j + 1) % x.ncol()) *
      x((i - 1 + x.nrow()) % x.nrow(), (j + 1) % x.ncol());
    //   x
    // # #
    total += x(i, j) *
      x((i - 1 + x.nrow()) % x.nrow(), (j - 1 + x.ncol()) % x.ncol()) *
      x((i - 1 + x.nrow()) % x.nrow(), j);
    // #  
    // # x
    total += x(i, j) *
      x(i, (j - 1 + x.ncol()) % x.ncol()) *
      x((i + 1) % x.nrow(), (j - 1 + x.ncol()) % x.ncol());

    // x #
    // #
    total += x(i, j) *
      x((i - 1 + x.nrow()) % x.nrow(), j) *
      x(i, (j + 1) % x.ncol());
    // # x
    //   #
    total += x(i, j) *
      x(i, (j - 1 + x.ncol()) % x.ncol()) *
      x((i - 1 + x.nrow()) % x.nrow(), j);
    //   #
    // # x
    total += x(i, j) *
      x(i, (j - 1 + x.ncol()) % x.ncol()) *
      x((i + 1) % x.nrow(), j);
    // #
    // x #
    total += x(i, j) *
      x(i, (j + 1) % x.ncol()) *
      x((i + 1) % x.nrow(), j);
  }
  
  return total;
}

// [[Rcpp::export]]
List ising_gibbs(NumericMatrix x_, double mu, NumericVector beta, 
                 NumericVector gamma, int S, int seed = 0) {
  NumericMatrix x(clone(x_));
  NumericMatrix out(S, beta.size() + gamma.size() + 1);
  colnames(out) = CharacterVector::create("X0", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7");
  
  const int ois[] = { -1, -1, -2, -2, -2 };
  const int ojs[] = { 0, 1, 0, 1, 2 };
  
  std::vector<int> is(x.size());
  std::vector<int> js(x.size());
  
  std::vector<int> idxs(x.size());
  
  for(int j = 0; j < x.ncol(); j++) {
    for(int i = 0; i < x.nrow(); i++) {
      idxs[i + j * x.nrow()] = i + j * x.nrow();
      is[i + j * x.nrow()] = i;
      js[i + j * x.nrow()] = j;
    }
  }
  
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> runif(0.0, 1.0);
  std::uniform_int_distribution<> dis(0, x.nrow() - 1);
  
  double X0 = solo(x);
  NumericVector Q(beta.size());
  for(int k = 0; k < beta.size(); k++)
    Q(k) = pairs(x, ois[k], ojs[k]);
  
  for(int s = 0; s < S; s++) {
    std::shuffle(idxs.begin(), idxs.end(), gen);
    
    for(int idx : idxs) {
      int i = is[idx];
      int j = js[idx];
      
      double dX0 = -2 * x(i, j);
      NumericVector dQ(beta.size());
      for(int k = 0; k < beta.size(); k++)
        dQ(k) = -2 * dpairs(x, i, j, ois[k], ojs[k]);
      double dE = mu * dX0 + std::inner_product(beta.begin(), beta.end(),
                                                dQ.begin(), 0.0);
      
      double r = runif(gen);
      if(r > 1.0 / (1.0 + std::exp(-dE))) {
        X0 += dX0;
        Q += dQ;
        x(i, j) = x(i, j) * -1;
      }
    }
    
    out(s, 0) = X0;
    out(s, 1) = Q(0);
    out(s, 2) = Q(1);
    out(s, 3) = Q(2);
    out(s, 4) = Q(3);
    out(s, 5) = Q(4);
  }
  
  List ret;
  ret["x"] = x;
  ret["states"] = out;
  return ret;
}
