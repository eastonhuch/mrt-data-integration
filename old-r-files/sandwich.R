# The only part of this function I'm not super confident is 
eastons_sandwich <- function(data, models, beta_h_formula, beta_s_formula) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  is_internal <- data$is_internal
  x2_internal <- data$x2[is_internal]
  
  # Construct design matrices
  X_alpha_s <- model.matrix(models$p_s$formula, data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_s <- model.matrix(beta_s_formula, data=data)
  X_beta_s_raw <- X_beta_s / a_centered
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  X_gamma_x2_internal <- model.matrix(models$s$formula, data=data[data$is_internal,])
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_s <- ncol(X_alpha_s)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_r <- ncol(X_gamma_x2_internal)
  d <- d_alpha_s + d_h + d_s + d_r
  
  # Extract coefficients
  alpha_s <- coef(models$p_s)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  gamma_x2 <- coef(models$s)
  
  # Construct position vectors
  pos_alpha_s <- seq(d_alpha_s)
  pos_beta_h <- d_alpha_s + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_gamma_x2 <- max(pos_beta_s) + seq(d_r)
  
  # p_s score/hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  p_s <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))  # Assume logit link
  scores[, pos_alpha_s] <- (a - p_s) * X_alpha_s
  sd_p_s <- sqrt(p_s * (1-p_s))
  X_alpha_s_scaled <- sd_p_s * X_alpha_s
  hessian[pos_alpha_s, pos_alpha_s] <- crossprod(X_alpha_s_scaled)
  
  # WCLS scores and Hessian
  p_s_a <- a*p_s + (1-a)*(1-p_s)
  w <- p_s_a / p_h_a
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_s_fitted_values <- c(X_beta_s %*% beta_s)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_s_fitted_values)
  wcls_weighted_resids <- w * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_s] <- wcls_weighted_resids * X_beta_s
  
  X_beta_hs_scaled <- sqrt(w) * X_beta_hs
  hessian[pos_beta_hs, pos_beta_hs] <- crossprod(X_beta_hs_scaled)
  
  # This is the only part I couldn't check with geeglm
  p_s_a_deriv <- -(a/(1-p_s) + (1-a)*p_s)
  p_s_deriv <- -1/(1-p_s)
  p_s_X_beta_s <- (-p_s) * X_beta_s_raw
  hessian[pos_beta_hs, pos_alpha_s] <- 
    t(X_beta_hs * wcls_weighted_resids) %*% p_s_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), p_s_X_beta_s) * wcls_weighted_resids) %*% p_s_deriv +
    t(X_beta_hs * (p_s * wcls_s_fitted_values * w)) %*% p_s_deriv
  hessian[pos_alpha_s, pos_beta_hs] <- t(hessian[pos_beta_hs, pos_alpha_s])
  
  # Gamma score
  scores[is_internal, pos_gamma_x2] <- (x2_internal - c(X_gamma_x2_internal %*% gamma_x2)) * X_gamma_x2_internal
  hessian[pos_gamma_x2, pos_gamma_x2] <- crossprod(X_gamma_x2_internal)
  
  meat <- crossprod(scores)
  half_sandwich <- solve(hessian, t(chol(meat)))
  sandwich <- tcrossprod(half_sandwich)
  sandwich
}
