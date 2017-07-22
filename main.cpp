#include <Eigen/Core>
#include <iostream>
#include <cstdlib>

using namespace Eigen;
using namespace std;

double E(const MatrixXd& X) {
  double total = 0.0;

  for(int i = 0; i < X.rows(); i++)
    for(int j = 0; j < X.cols(); j++) {
      total += mu * X(i, j);
      
      if(i < X.rows() - 1)
	total += beta * X(i, j) * X(i + 1, j);

      if(j < X.cols() - 1)
	total += beta * X(i, j) * X(i, j + 1);
    }

  return total;
}

int N = 5;
double mu = 0.0;
double beta = 0.1;

int main(int argc, char **argv) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> unif(0.0, 1.0);
  std::uniform_int_distribution<> iunif(0, N - 1);

  MatrixXd X(N, N);

  srand(0);
  
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      X(i, j) = (rand() % 2) * 2 - 1.0;

  int S = 1000;
  double kT = 1.0;

  for(int s = 0; s < S; s++) {
    double Ec = E(X);

    int i = iunif(gen);
    int j = iunif(gen);

    MatrixXd Xn = X.copy();

    Xn(i, j) *= -1;

    double En = E(Xn);

    if(exp((En - Ec) / kT) < unif(gen)) {
      X(i, j) *= -1;
    }
  }
}
