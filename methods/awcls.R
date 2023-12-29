# NOTE: This code has not been checked thoroughly; use with caution
awcls_sandwich <- function(data, models, beta_h_formula, beta_s_formula) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  is_internal <- data$is_internal
  x2_internal <- data$x2[is_internal]
  
  # Construct design matrices
  X_alpha_s <- model.matrix(formula(models$p_s), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_s <- model.matrix(beta_s_formula, data=data)
  X_beta_s_raw <- X_beta_s / a_centered
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  X_gamma_x2_internal <- model.matrix(formula(models$s), data=data[data$is_internal,])
  
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
  
  # p_s_hat score/hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  p_s_hat <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))  # Assume logit link
  scores[, pos_alpha_s] <- (a - p_s_hat) * X_alpha_s
  sd_p_s_hat <- sqrt(p_s_hat * (1-p_s_hat))
  X_alpha_s_scaled <- sd_p_s_hat * X_alpha_s
  hessian[pos_alpha_s, pos_alpha_s] <- crossprod(X_alpha_s_scaled)
  
  # WCLS scores and Hessian
  p_s_hat_a <- a*p_s_hat + (1-a)*(1-p_s_hat)
  w <- p_s_hat_a / p_h_a
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_s_fitted_values <- c(X_beta_s %*% beta_s)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_s_fitted_values)
  wcls_weighted_resids <- w * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_s] <- wcls_weighted_resids * X_beta_s
  
  X_beta_hs_scaled <- sqrt(w) * X_beta_hs
  hessian[pos_beta_hs, pos_beta_hs] <- crossprod(X_beta_hs_scaled)
  
  p_s_hat_a_deriv <- -((2*a-1) * p_s_hat * (1 - p_s_hat) / p_s_hat_a)
  p_s_deriv <- -(1-p_s_hat)
  p_s_X_beta_s <- p_s_hat * X_beta_s_raw
  hessian[pos_beta_hs, pos_alpha_s] <- 
    t(X_beta_hs * wcls_weighted_resids) %*% p_s_hat_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_s) * wcls_weighted_resids) %*% p_s_deriv +
    t(X_beta_hs * (p_s_hat * wcls_s_fitted_values / a_centered * w)) %*% p_s_deriv
  
  # Gamma score
  scores[is_internal, pos_gamma_x2] <- (x2_internal - c(X_gamma_x2_internal %*% gamma_x2)) * X_gamma_x2_internal
  hessian[pos_gamma_x2, pos_gamma_x2] <- crossprod(X_gamma_x2_internal)
  
  n_users <- max(data$user_id)
  t_max <- floor(n / n_users)
  scores_agg <- apply(
    aperm(
      array(scores, dim = c(t_max, n_users, d)),
      c(2,1,3)),
    MARGIN=c(1,3), FUN=sum)
  meat <- crossprod(scores_agg)
  half_sandwich <- solve(hessian, t(chol(meat)), tol=1e-50)
  sandwich <- tcrossprod(half_sandwich) * n / (n-d)
  list(
    sandwich=sandwich,
    bread=hessian,
    meat=meat
  )
}

awcls <- function(data) {
  make_Gamma <- function(gamma_x2) {
    cbind(
      c(1, 0, 0, 0),
      c(0, 1, 0, 0),
      gamma_x2)
  }
  beta_s_true <- c(1, 2, -3)
  gamma_x2_true <- c(2, 1, -0.3, -0.1)
  beta_r_true <- c(-5, -1, 0.9, 0.3)
  theta_true <- c(beta_s_true, gamma_x2_true)
  
  # Get point estimates
  # Propensity score model
  p_s_mod <- glm(a ~ 1, data=data, family=binomial())
  alpha_s <- coef(p_s_mod)
  data$p_s_hat <- predict(p_s_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_s_hat
  data$p_s_hat_a <- data$a * data$p_s_hat + (1 - data$a) * (1 - data$p_s_hat)
  data$w <- data$p_s_hat_a / data$p_h_a
  
  # WCLS
  beta_h_formula <- y ~ x1 + x2 + x3
  beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  d_s <- length(beta_s)
  
  # Gamma
  s_formula <- x2 ~ x1 + I(x1^2) + I(x1^3)
  s_mod <- glm(s_formula, data=data[data$is_internal,])
  gamma_x2 <- coef(s_mod)
  Gamma <- make_Gamma(gamma_x2)
  row.names(Gamma) <- names(gamma_x2)
  
  # beta_r
  beta_r <- c(Gamma %*% beta_s)
  d_r <- length(beta_r)
  names(beta_r) <- row.names(Gamma)
  
  # Models list
  models <- list(
    p_s=p_s_mod,
    wcls=wcls_mod,
    s=s_mod
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_s,
    beta_h,
    beta_s,
    gamma_x2
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich_list <- awcls_sandwich(data, models, beta_h_formula, beta_s_formula)
  sandwich <- sandwich_list$sandwich
  pos_theta <- seq(
    length(alpha_s) + length(beta_h) + 1,
    length(vector_estimate)
  )
  var_theta <- sandwich[pos_theta, pos_theta]
  J_theta <- cbind(
    Gamma,
    beta_s[3] * diag(d_r)
  )
  var_beta_r <- J_theta %*% var_theta %*% t(J_theta)
  se_beta_r <- sqrt(diag(var_beta_r))
  beta_r_error <- beta_r - beta_r_true
  beta_r_z_scores <- beta_r_error / se_beta_r
  beta_r_chi2 <- beta_r_error %*% solve(var_beta_r, beta_r_error)
  
  results <- list(
    beta_r=beta_r,
    se_beta_r=se_beta_r,
    var_beta_r=var_beta_r,
    beta_r_chi2=beta_r_chi2,
    beta_r_z_scores=beta_r_z_scores,
    se_beta_s=sqrt(diag(var_theta[1:3, 1:3])),
    sandwich=sandwich,
    bread=sandwich_list$bread,
    meat=sandwich_list$meat
  )
  results
}
