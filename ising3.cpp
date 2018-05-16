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
  } else if(w == 1) {
    // ##x
    total += x(i, j) *
      x(i, (j - 2 + x.ncol()) % x.ncol()) *
      x(i, (j - 1 + x.ncol()) % x.ncol());
    // #x#
    total += x(i, j) *
      x(i, (j - 1 + x.ncol()) % x.ncol()) *
      x(i, (j + 1) % x.ncol());
    // x##
    total += x(i, j) *
      x(i, (j + 2) % x.ncol()) *
      x(i, (j + 1) % x.ncol());
    
    // #
    // #
    // x
    total += x(i, j) *
      x((i + 2) % x.nrow(), j) *
      x((i + 1) % x.nrow(), j);
    // #
    // x
    // #
    total += x(i, j) *
      x((i + 1) % x.nrow(), j) *
      x((i - 1 + x.nrow()) % x.nrow(), j);
    // x
    // #
    // #
    total += x(i, j) *
      x((i - 2 + x.nrow()) % x.nrow(), j) *
      x((i - 1 + x.nrow()) % x.nrow(), j);
  } else {
    throw std::logic_error("w value not implemented"); 
  }
  
  return total;
}

// [[Rcpp::export]]
double triplets(const NumericMatrix &x, int w) {
  double total = 0.0;
  
  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < x.ncol(); j++) {
      total += dtriplets(x, i, j, w) / 3.0;        
    }
  }

  return total;
}

// [[Rcpp::export]]
List ising_gibbs(NumericMatrix x_, double mu, NumericVector beta, 
                 NumericVector gamma, int S, int seed = 0, double kT = 1.0) {
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
  NumericVector Q(beta.size() + gamma.size());
  for(int k = 0; k < beta.size(); k++)
    Q(k) = pairs(x, ois[k], ojs[k]);
  for(int k = 0; k < gamma.size(); k++)
    Q(beta.size() + k) = triplets(x, k);

  for(int s = 0; s < S; s++) {
    std::shuffle(idxs.begin(), idxs.end(), gen);

    for(int idx : idxs) {
      int i = is[idx];
      int j = js[idx];

      double dX0 = -2 * x(i, j);
      NumericVector dQ(beta.size() + gamma.size());
      for(int k = 0; k < beta.size(); k++)
        dQ(k) = -2 * dpairs(x, i, j, ois[k], ojs[k]);
      for(int k = 0; k < gamma.size(); k++)
        dQ(beta.size() + k) = -2 * dtriplets(x, i, j, k);
      double dE = mu * dX0;
      for(int k = 0; k < beta.size(); k++)
        dE += beta[k] * dQ(k);
      for(int k = 0; k < gamma.size(); k++)
        dE += gamma[k] * dQ(beta.size() + k);
      
      double r = runif(gen);
      if(r > 1.0 / (1.0 + std::exp(-dE / kT))) {
        X0 += dX0;
        Q += dQ;
        x(i, j) = x(i, j) * -1;
      }
    }
    
    out(s, 0) = X0;
    for(int k = 0; k < beta.size(); k++)
      out(s, k + 1) = Q(k);
    for(int k = 0; k < gamma.size(); k++)
      out(s, beta.size() + k + 1) = Q(beta.size() + k);
  }
  
  List ret;
  ret["x"] = x;
  ret["states"] = out;
  return ret;
}
