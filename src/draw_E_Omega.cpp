//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat draw_E_Omega(arma::vec alpha, arma::vec beta, arma::vec theta) {
  // Generate I length row vector filled with 1
  arma::vec intercept(theta.n_rows);
  intercept = intercept.ones();
  
  
  // 2 by J item parameter matrix 
  arma::mat beta_tilde = join_horiz(alpha, beta).t();
  
  // I by 2 individual parameter matrix
  // 1st column is intercept
  arma::mat Theta =  join_horiz(intercept, theta);
   
  // Calculate psi
  arma::mat psi = Theta * beta_tilde;
  
  // Calculate E(omega)
  arma::mat Omega = arma::tanh(psi / 2) / (2 * psi);
  
  return(Omega);
}