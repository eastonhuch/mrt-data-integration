walters_sandwich <- function(data_pooled, data_internal, models, beta_h_formula, beta_s_formula, observational=FALSE) {
  # Extract some columns
  y <- data_pooled$y
  p_h_a <- data_pooled$p_h_a
  a <- data_pooled$a
  a_centered <- data_pooled$a_centered
  is_internal <- data_pooled$is_internal
  x2_internal <- data_pooled$x2[is_internal]
  
  # Construct design matrices
  if (observational) X_alpha_h <- model.matrix(formula(models$p_h), data=data_pooled)
  X_alpha_s <- model.matrix(formula(models$p_s), data=data_pooled)
  X_beta_h <- model.matrix(beta_h_formula, data=data_pooled)
  X_beta_s <- model.matrix(beta_s_formula, data=data_pooled)
  X_beta_s_raw <- X_beta_s / a_centered
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  X_beta_r_internal <- model.matrix(formula(models$r), data=data_internal)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  if (observational) d_alpha_h <- ncol(X_alpha_h)
  d_alpha_s <- ncol(X_alpha_s)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_hs <- d_h + d_s
  d_r <- ncol(X_beta_r_internal)
  d <- d_alpha_s + d_h + d_s + d_r
  if (observational) d <- d_alpha_h + d
  
  # Extract coefficients
  if (observational) alpha_h <- coef(models$p_h)
  alpha_s <- coef(models$p_s)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  beta_r <- coef(models$r)
  
  # Construct position vectors
  pos_alpha_s <- seq(d_alpha_s)
  if (observational) {
    pos_alpha_h <- seq(d_alpha_h)
    pos_alpha_s <- d_alpha_h + pos_alpha_s
  }
  pos_beta_h <- max(pos_alpha_s) + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_beta_r <- max(pos_beta_s) + seq(d_r)
  
  # scores/Hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  
  # p_h_hat score/hessian
  if (observational) {
    p_h_hat <- 1 / (1 + exp(-c(X_alpha_h %*% alpha_h)))  # Assume logit link
    scores[, pos_alpha_h] <- (a - p_h_hat) * X_alpha_h
    sd_p_h_hat <- sqrt(p_h_hat * (1-p_h_hat))
    X_alpha_h_scaled <- sd_p_h_hat * X_alpha_h
    hessian[pos_alpha_h, pos_alpha_h] <- -crossprod(X_alpha_h_scaled)
  }
  
  # p_s_hat score/hessian
  p_s_hat <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))  # Assume logit link
  scores[, pos_alpha_s] <- (a - p_s_hat) * X_alpha_s
  sd_p_s_hat <- sqrt(p_s_hat * (1-p_s_hat))
  X_alpha_s_scaled <- sd_p_s_hat * X_alpha_s
  hessian[pos_alpha_s, pos_alpha_s] <- -crossprod(X_alpha_s_scaled)
  
  # WCLS scores and Hessian
  p_s_hat_a <- a*p_s_hat + (1-a)*(1-p_s_hat)
  w <- p_s_hat_a / p_h_a
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_s_fitted_values <- c(X_beta_s %*% beta_s)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_s_fitted_values)
  wcls_weighted_resids <- w * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_s] <- wcls_weighted_resids * X_beta_s
  
  sqrt_w <- sqrt(w)
  X_beta_hs_scaled <- sqrt_w * X_beta_hs
  GtWG <- crossprod(X_beta_hs_scaled)
  hessian[pos_beta_hs, pos_beta_hs] <- GtWG

  # I think these have negative signs because I'm using (-A^-1) B (-A^-1)
  p_s_hat_a_deriv <- -(2*a-1) * p_s_hat * (1 - p_s_hat)
  log_p_s_hat_a_deriv <- p_s_hat_a_deriv / p_s_hat_a
  log_p_s_deriv <- -(1-p_s_hat)
  p_s_X_beta_s <- p_s_hat * X_beta_s_raw
  hessian[pos_beta_hs, pos_alpha_s] <- 
    t(X_beta_hs * wcls_weighted_resids) %*% log_p_s_hat_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_s) * wcls_weighted_resids) %*% log_p_s_deriv +
    t(X_beta_hs * (p_s_hat * wcls_s_fitted_values / a_centered * w)) %*% log_p_s_deriv
  
  if (observational) {
    p_h_hat_a_deriv <- -(2*a-1) * p_h_hat * (1 - p_h_hat)
    log_p_h_hat_a_deriv <- p_h_hat_a_deriv / p_h_a  # NOTE: This is really p_h_hat_a
    hessian[pos_beta_hs, pos_alpha_h] <- 
      t(X_beta_hs * wcls_weighted_resids) %*% (log_p_h_hat_a_deriv * X_alpha_h)
  }

  # beta_r score
  scores[is_internal, pos_beta_r] <- (data_pooled$wcls_s_causal_effects[data_pooled$is_internal] - c(X_beta_r_internal %*% beta_r)) * X_beta_r_internal
  hessian[pos_beta_r, pos_beta_r] <- crossprod(X_beta_r_internal)
  hessian[pos_beta_r, pos_beta_s] <- -t(X_beta_r_internal) %*% X_beta_s_raw[data_pooled$is_internal,]
  
  # Assemble sandwich
  n_users <- max(data_pooled$user_id)
  t_max <- floor(n / n_users)
  scores_agg <- apply(
    aperm(
      array(scores, dim = c(t_max, n_users, d)),
      c(2,1,3)),
    MARGIN=c(1,3), FUN=sum)
  meat <- crossprod(scores_agg)
  half_sandwich <- solve(hessian, t(chol(meat)))
  sandwich <- tcrossprod(half_sandwich)
  list(
    sandwich=sandwich,
    bread=hessian,
    meat=meat
  )
}

walters_method <- function(data, internal_only=FALSE, observational=FALSE) {
  # Create data_pooled, data_internal
  if (internal_only) {
    data_pooled <- data[data$is_internal,]
  } else {
    data_pooled <- data
  }
  
  # True coefficient vectors
  beta_s_true <- c(1, 2, -3)
  beta_r_true <- c(-5, -1, 0.9, 0.3)
  
  # Get point estimates
  # p_h
  if (observational) {
    if (internal_only) stop("observational==TRUE and internal_only==TRUE not implemented")
    p_h_mod <- glm(a ~ 1 + as.numeric(is_internal) + x1 + x2 + x3, data=data_pooled, family=binomial())
    alpha_h <- coef(p_h_mod)
    data_pooled$p_h <- predict(p_h_mod, newdata=data_pooled, type="response")
    data_pooled$p_h_a <- data_pooled$a * data_pooled$p_h + (1 - data_pooled$a) * (1 - data_pooled$p_h)
  }
  
  # p_s
  p_s_mod <- glm(a ~ 1, data=data_pooled, family=binomial())
  alpha_s <- coef(p_s_mod)
  data_pooled$p_s_hat <- predict(p_s_mod, newdata=data_pooled, type="response")
  data_pooled$a_centered <- data_pooled$a - data_pooled$p_s_hat
  data_pooled$p_s_hat_a <- data_pooled$a * data_pooled$p_s_hat + (1 - data_pooled$a) * (1 - data_pooled$p_s_hat)
  data_pooled$w <- data_pooled$p_s_hat_a / data_pooled$p_h_a
  
  # WCLS
  beta_h_formula <- y ~ x1 + x2 + x3
  beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data_pooled, weights=w)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  data_pooled$wcls_s_causal_effects <- c(model.matrix(beta_s_formula, data=data_pooled) %*% beta_s) / data_pooled$a_centered
  d_s <- length(beta_s)
  
  # beta_r
  data_internal <- data_pooled[data_pooled$is_internal,]
  r_formula <- wcls_s_causal_effects ~ x1 + I(x1^2) + I(x1^3)
  r_mod <- glm(r_formula, data=data_internal)
  beta_r <- coef(r_mod)

  # Models list
  models <- list(
    p_s=p_s_mod,
    wcls=wcls_mod,
    r=r_mod
  )
  if (observational) models$p_h <- p_h_mod
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_s,
    beta_h,
    beta_s,
    beta_r
  )
  if (observational) vector_estimate <- c(alpha_h, vector_estimate)
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich_list <- walters_sandwich(data_pooled, data_internal, models, beta_h_formula, beta_s_formula, observational=observational)
  sandwich <- sandwich_list$sandwich
  pos_beta_r <- length(alpha_s) + length(beta_h) + length(beta_s) + seq_along(beta_r)
  if (observational) pos_beta_r <- length(alpha_h) + pos_beta_r
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
    beta_r_z_scores=beta_r_z_scores,
    sandwich=sandwich,
    bread=sandwich_list$bread,
    meat=sandwich_list$meat
  )
  results
}
