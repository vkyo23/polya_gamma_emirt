#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector draw_alpha(NumericMatrix Y, NumericVector beta, NumericVector theta, 
                         NumericMatrix Omega, int a0, int A0) {
  int I = Y.nrow();
  int J = Y.ncol();
  
  NumericVector alpha(J);
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
        left_part += omega_ij;
        right_part_1 += k_ij;
        right_part_2 += theta[i] * omega_ij;
      }
    }
    left_part = left_part + 1/A0;
    double right_part = a0/A0 + right_part_1 - beta[j] * right_part_2;
    alpha[j] = (1 / left_part) * right_part;
  }
  
  return(alpha);
}