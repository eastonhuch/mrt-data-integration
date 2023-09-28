wcls_sandwich <- function(data, models, beta_h_formula, beta_r_formula, tilt=FALSE) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered

  # Construct design matrices
  X_alpha_r <- model.matrix(formula(models$p_r), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_r <- model.matrix(beta_r_formula, data=data)
  X_beta_r_raw <- X_beta_r / a_centered
  X_beta_hr <- cbind(X_beta_h, X_beta_r)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_r <- ncol(X_alpha_r)
  d_h <- ncol(X_beta_h)
  d_r <- ncol(X_beta_r)
  d_hr <- d_h + d_r
  d <- d_alpha_r + d_h + d_r
  
  # Extract coefficients
  alpha_r <- coef(models$p_r)
  beta_hr <- coef(models$wcls)
  beta_h <- beta_hr[seq(d_h)]
  beta_r <- beta_hr[seq(d_h+1, d_h+d_r)]
  
  # Construct position vectors
  pos_alpha_r <- seq(d_alpha_r)
  pos_beta_h <- d_alpha_r + seq(d_h)
  pos_beta_r <- max(pos_beta_h) + seq(d_r)
  pos_beta_hr <- c(pos_beta_h, pos_beta_r)
  
  # p_r_hat score/hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  p_r_hat <- 1 / (1 + exp(-c(X_alpha_r %*% alpha_r)))  # Assume logit link
  scores[, pos_alpha_r] <- (a - p_r_hat) * X_alpha_r
  sd_p_r_hat <- sqrt(p_r_hat * (1-p_r_hat))
  X_alpha_r_scaled <- sd_p_r_hat * X_alpha_r
  hessian[pos_alpha_r, pos_alpha_r] <- -crossprod(X_alpha_r_scaled)
  
  # WCLS scores and Hessian
  p_r_hat_a <- data$p_r_hat_a
  w <- data$w
  tilt_ratios <- data$tilt_ratios
  w_alt <- data$w_alt
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_r_fitted_values <- c(X_beta_r %*% beta_r)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_r_fitted_values)
  wcls_weighted_resids <- w * tilt_ratios * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_r] <- wcls_weighted_resids * X_beta_r
  
  X_beta_hr_scaled <- sqrt(w_alt) * X_beta_hr
  GtWG <- crossprod(X_beta_hr_scaled)
  hessian[pos_beta_hr, pos_beta_hr] <- GtWG

  p_r_hat_a_deriv <- -(2*a-1) * p_r_hat * (1 - p_r_hat)
  log_p_r_hat_a_deriv <- p_r_hat_a_deriv / p_r_hat_a
  p_r_deriv <- -(1-p_r_hat)
  p_r_X_beta_r <- p_r_hat * X_beta_r_raw
  hessian[pos_beta_hr, pos_alpha_r] <- 
    t(X_beta_hr * wcls_weighted_resids) %*% log_p_r_hat_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_r_X_beta_r) * wcls_weighted_resids) %*% p_r_deriv +
    t(X_beta_hr * (p_r_hat * wcls_r_fitted_values / a_centered * w * tilt_ratios)) %*% p_r_deriv

  # Assemble sandwich
  n_users <- max(data$user_id)
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

wcls <- function(data, tilt=FALSE) {
  beta_r_true <- c(-2, 5)
  
  # Get point estimates
  # Propensity score model
  p_r_mod <- glm(a ~ 1, data=data, family=binomial())
  alpha_r <- coef(p_r_mod)
  data$p_r_hat <- predict(p_r_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_r_hat
  data$p_r_hat_a <- data$a * data$p_r_hat + (1 - data$a) * (1 - data$p_r_hat)
  data$w <- data$p_r_hat_a / data$p_h_a
  
  if (tilt) {
    tilt_mod <- glm(
      is_internal ~ ns(x1, df=5) + ns(x2, df=5) + I(x1 * x2) + I(x1 * x2^2) + I(x1^2 * x2) + I(x1^2 * x2^2),
      family=binomial(), data=data)
    gamma <- coef(tilt_mod)
    internal_prop <- mean(data$is_internal)
    gamma[1] <- gamma[1] - log(internal_prop / (1-internal_prop))
    X_gamma <- model.matrix(tilt_mod)
    raw_tilt_ratios <- c(exp(X_gamma %*% gamma))
    data$tilt_ratios <- data$is_internal + data$is_external * raw_tilt_ratios
  } else {
    data$tilt_ratios <- 1
  }
  data$w_and_tilt <- data$w * data$tilt_ratios
  avg_tilt_ratio_internal <- mean(raw_tilt_ratios[dat$is_internal])
  avg_tilt_ratio_external <- mean(raw_tilt_ratios[dat$is_external])
  if (is.na(avg_tilt_ratio_internal)) avg_tilt_ratio_internal <- 1
  if (is.na(avg_tilt_ratio_external)) avg_tilt_ratio_external <- 1
  data$w_alt <- data$w_and_tilt * 
    (data$is_internal * avg_tilt_ratio_internal +
      data$is_external * avg_tilt_ratio_external)

  # WCLS
  beta_h_formula <- y ~ x1 + x2 + x3
  beta_r_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1)
  beta_r_formula_character <- as.character(update(beta_r_formula, . ~ . + 1))[3]
  beta_r_formula_symbol <- rlang::parse_expr(beta_r_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_r_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w_alt)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_r <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  data$wcls_r_causal_effects <- c(model.matrix(beta_r_formula, data=data) %*% beta_r) / data$a_centered
  d_r <- length(beta_r)
  
  # Models list
  models <- list(
    p_r=p_r_mod,
    wcls=wcls_mod
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_r,
    beta_h,
    beta_r
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich_list <- wcls_sandwich(data, models, beta_h_formula, beta_r_formula, tilt=tilt)
  sandwich <- sandwich_list$sandwich
  pos_beta_r <- length(alpha_r) + length(beta_h) + seq_along(beta_r)
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
