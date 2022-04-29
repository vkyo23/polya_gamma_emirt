#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector draw_beta(NumericMatrix Y, NumericVector alpha, NumericVector theta, 
                        NumericMatrix Omega, int b0, int B0) {
  int I = Y.nrow();
  int J = Y.ncol();
  
  NumericVector beta(J);
  for (int j = 0; j < J; j++) {
    double left_part = 0;
    double right_part_1 = 0;
    double right_part_2 = 0;
    for (int i = 0; i < I; i++) {
      if (NumericVector::is_na(Y(i, j))) {
        left_part += 0;
        right_part_1 += 0;
        right_part_2 += 0;
      } else {
        double omega_ij = Omega(i, j);
        double k_ij = Y(i, j) - 0.5;
        left_part += std::pow(theta[i], 2.0) * omega_ij;
        right_part_1 += k_ij * theta[i];
        right_part_2 += theta[i] * omega_ij;
      }
    }
    left_part = left_part + 1/B0;
    double right_part = b0/B0 + right_part_1 - alpha[j] * right_part_2;
    beta[j] = (1 / left_part) * right_part;
  }
  
  return(beta);
}