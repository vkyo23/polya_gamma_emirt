#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector draw_theta(NumericMatrix Y, NumericVector alpha, NumericVector beta, 
                         NumericMatrix Omega, int constraint) {
  int I = Y.nrow();
  int J = Y.ncol();
  
  NumericVector theta(I);
  for (int i = 0; i < I; i++) {
    double left_part = 0;
    double right_part = 0;
    for (int j = 0; j < J; j++) {
      if (NumericVector::is_na(Y(i, j))) {
        left_part += 0;
        right_part += 0;
      } else {
        double omega_ij = Omega(i, j);
        double k_ij = Y(i, j) - 0.5;
        left_part += std::pow(beta[j], 2.0) * omega_ij;
        right_part += beta[j] * k_ij - alpha[j] * beta[j] * omega_ij;
      }
    }
    left_part = left_part + 1;
    theta[i] = (1 / left_part) * right_part;
  }
  
  if (theta[constraint - 1] < 0) {
    theta = -theta;
  }
  
  return(theta);
}