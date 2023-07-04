# The only part of this function I'm not super confident is 
walters_sandwich <- function(data, models, beta_h_formula, beta_s_formula) {
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
  X_beta_r_internal <- model.matrix(formula(models$r), data=data[data$is_internal,])
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_s <- ncol(X_alpha_s)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_r <- ncol(X_beta_r_internal)
  d <- d_alpha_s + d_h + d_s + d_r
  
  # Extract coefficients
  alpha_s <- coef(models$p_s)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  beta_r <- coef(models$r)
  
  # Construct position vectors
  pos_alpha_s <- seq(d_alpha_s)
  pos_beta_h <- d_alpha_s + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_beta_r <- max(pos_beta_s) + seq(d_r)
  
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
  X_beta_s_no_a <- X_beta_s / a_centered
  wcls_s_causal_effects <- c(X_beta_s_no_a %*% beta_s)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_s_fitted_values)
  wcls_weighted_resids <- w * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_s] <- wcls_weighted_resids * X_beta_s
  
  X_beta_hs_scaled <- sqrt(w) * X_beta_hs
  hessian[pos_beta_hs, pos_beta_hs] <- crossprod(X_beta_hs_scaled)
  
  p_s_a_deriv <- (a*(1-p_s) + (1-a)*p_s) / p_s_a
  p_s_deriv <- 1-p_s
  p_s_X_beta_s <- p_s * X_beta_s_raw
  hessian[pos_beta_hs, pos_alpha_s] <- 
    t(X_beta_hs * wcls_weighted_resids) %*% p_s_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_s) * wcls_weighted_resids) %*% p_s_deriv +
    t(X_beta_hs * (p_s * wcls_s_causal_effects * w)) %*% p_s_deriv
  
  # beta_r score
  scores[is_internal, pos_beta_r] <- (wcls_s_causal_effects - c(X_beta_r_internal %*% beta_r)) * X_beta_r_internal
  hessian[pos_beta_r, pos_beta_r] <- crossprod(X_beta_r_internal)
  hessian[pos_beta_r, pos_beta_s] <- -(t(X_beta_r_internal) %*% X_beta_s_no_a[data$is_internal,])
  
  meat <- crossprod(scores)
  half_sandwich <- solve(hessian, t(chol(meat)))
  sandwich <- tcrossprod(half_sandwich)
  sandwich
}

# Need models$r for projection
walters_method <- function(data) {
  beta_s_true <- c(1, 2, -3)
  beta_r_true <- c(-5, -1, 0.9, 0.3)
  theta_true <- c(beta_s_true, beta_r_true)
  
  # Get point estimates
  # Propensity score model
  p_s_mod <- glm(a ~ 1, data=data, family=binomial())
  alpha_s <- coef(p_s_mod)
  data$p_s_hat <- predict(p_s_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_s_hat
  data$p_s_hat_a <- data$a * data$p_s_hat + (1 - data$a) * (1 - data$p_s_hat)
  data$w <- data$p_s_hat_a / data$p_h_a
  data$wcls_s_causal_effects <- 0
  data_internal <- data[data$is_internal,]
  
  # WCLS
  beta_h_formula <- y ~ x1 + x2 + x3
  beta_s_formula <- y ~ 0 + a_centered + I(a_centered * x1) + I(a_centered * x2)
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  d_s <- length(beta_s)
  data_internal$wcls_s_causal_effects <- (model.matrix(beta_s_formula, data_internal) %*% beta_s) / data_internal$a_centered
  
  # beta_r
  r_mod <- lm(wcls_s_causal_effects ~ x1 + I(x1^2) + I(x1^3), data=data_internal)
  beta_r <- coef(r_mod)
  d_r <- length(beta_r)
  
  # Models list
  models <- list(
    p_s=p_s_mod,
    wcls=wcls_mod,
    r=r_mod
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_s,
    beta_h,
    beta_s,
    beta_r
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich <- walters_sandwich(data, models, beta_h_formula, beta_s_formula)
  pos_beta_r <- seq(
    length(alpha_s) + length(beta_h) + length(beta_s) + 1,
    length(vector_estimate)
  )
  var_beta_r <- sandwich[pos_beta_r, pos_beta_r]
  se_beta_r <- sqrt(diag(var_beta_r))
  beta_r_error <- beta_r - beta_r_true
  beta_r_z_scores <- beta_r_error / se_beta_r
  beta_r_chi2 <- beta_r_error %*% solve(var_beta_r, beta_r_error)
  
  results <- list(
    beta_r=beta_r,
    se_beta_r=se_beta_r,
    var_beta_r=var_beta_r,
    beta_r_chi2=beta_r_chi2,
    sandwich=sandwich
  )
  results
}
