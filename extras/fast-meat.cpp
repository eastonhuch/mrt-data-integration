#include<RcppArmadillo.h>
#include <stdexcept>

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat getMeat(arma::mat D, arma::colvec v, arma::colvec r, int t) {
  // Check inputs
  int nt = D.n_rows;
  if (nt % t != 0) {
    throw std::invalid_argument("D.n_rows is not a multiple of t");
  }
  
  if (v.n_rows != nt) {
    throw std::invalid_argument("D and v (as a column vector) must have the same number of rows");
  }
  
  if (r.n_rows != nt) {
    throw std::invalid_argument("D and r (as a column vector) must have the same number of rows");
  }
  
  // Create meat matrix
  int p = D.n_cols;
  arma::mat meat = arma::zeros(p, p);
  
  // Loop through participants
  int idx_upper;
  arma::mat D_subset;
  arma::colvec v_subset;
  arma::colvec r_subset;
  arma::rowvec z;
  for (int idx_lower = 0; idx_lower < nt; idx_lower += t) {
    idx_upper = idx_lower + t - 1;
    D_subset = D.rows(idx_lower, idx_upper);
    v_subset = v.rows(idx_lower, idx_upper);
    r_subset = r.rows(idx_lower, idx_upper);
    z = (r_subset / v_subset).t() * D_subset;
    meat += z.t() * z;
  }

  return meat;
}

/*** R
require(microbenchmark)
P_val <- 5
K_val <- 34
T_val <- 11
D <- matrix(rnorm(P_val*K_val*T_val), nrow=K_val*T_val, ncol=P_val)
v <- exp(rnorm(K_val * T_val))
r <- rnorm(K_val * T_val)

meat_cpp <- getMeat(D, v, r, T_val)

getMeatR <- function(D, v, r, t) {
  meat_R <- matrix(0, nrow=P_val, ncol=P_val)
  for (k in seq(0, K_val-1)) {
    idx_lower <- k * T_val
    idx_upper <- (k+1) * T_val - 1
    idx <- (idx_lower:idx_upper) + 1
    z <- t(r[idx] / v[idx]) %*% D[idx,]
    meat_inc <- crossprod(z)
    meat_R <- meat_R + meat_inc
  }
  meat_R
}
meat_R <- getMeatR(D, v, r, T_val)

all(meat_cpp == meat_R)
microbenchmark("R" = {getMeatR(D, v, r, T_val)}, "Rcpp" = {getMeat(D, v, r, T_val)}, times=40)
*/
