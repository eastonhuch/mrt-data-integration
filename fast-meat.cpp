#include<RcppArmadillo.h>
#include <stdexcept>

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat getMeat(arma::mat bread, arma::cube D, arma::mat v, arma::mat r) {
  // Check inputs
  int P = bread.n_rows;
  if (bread.n_cols != P) {
    throw std::invalid_argument("bread must be symmetric");
  }
    
  int K = D.n_rows;
  int T = D.n_cols;
  if (D.n_slices != P) {
    throw std::invalid_argument("Last dimension of D must == bread.n_cols");
  }
  
  if (v.n_rows != K) {
    throw std::invalid_argument("First dimension of v must == D.n_slices");
  }
  
  if (v.n_cols != T) {
    throw std::invalid_argument("Second dimension of v must == D.n_rows");
  }
  
  if (r.n_rows != K) {
    throw std::invalid_argument("First dimension of r must == D.n_slices");
  }
  
  if (r.n_cols != T) {
    throw std::invalid_argument("Second dimension of r must == D.n_rows");
  }
  
  // Create meat matrix
  arma::mat meat = arma::zeros(P, P);
  
  // Loop through participants
  arma::mat D_k;
  arma::mat Dt_vi_k;
  arma::rowvec v_k;
  arma::colvec z;
  arma::mat H;
  for (int k = 0; k < K; k++) {
    D_k = D.row(k);
    v_k = v.row(k);
    Dt_vi_k = D_k.t();
    Dt_vi_k.each_row() /= v_k;
    H = D_k * bread * Dt_vi_k;
    z = Dt_vi_k * arma::solve(arma::eye(T, T) - H, r.row(k).t());
    meat += z * z.t();
  }

  return meat;
}

/*** R
P_val <- 10
K_val <- 1000
T_val <- 20
x <- matrix(rnorm(P_val^2), nrow=P_val, ncol=P_val)
bread <- crossprod(x)
D <- array(rnorm(P_val*K_val*T_val), dim=c(K_val, T_val, P_val))
v <- matrix(exp(rnorm(K_val * T_val)), nrow=K_val, ncol=T_val)
r <- matrix(rnorm(K_val * T_val), nrow=K_val, ncol=T_val)

meat_cpp <- getMeat(bread, D, v, r)

getMeatR <- function(bread, D, v, r) {
  meat_R <- matrix(0, nrow=P_val, ncol=P_val)
  for (k in seq(K_val)) {
    D_k <- D[k,,]
    DtVi <- t(D_k / v[k,])
    H <- D_k %*% bread %*% DtVi
    z <- DtVi %*% solve(diag(T_val) - H, c(r[k,]))
    meat_inc <- tcrossprod(z)
    meat_R <- meat_R + meat_inc
  }
  meat_R
}
meat_R <- getMeatR(bread, D, v, r)

all(meat_cpp == meat_R)
microbenchmark("R" = {getMeatR(bread, D, v, r)}, "Rcpp" = {getMeat(bread, D, v, r)}, times=10)
*/
