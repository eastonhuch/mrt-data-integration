wcls_sandwich <- function(data, models, beta_h_formula, beta_r_formula, is_balanced=TRUE, estimate_p_r=FALSE, tilt=FALSE) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  is_internal <- data$is_internal
  is_external <- data$is_external

  # Construct design matrices
  if (estimate_p_r) X_alpha_r <- model.matrix(formula(models$p_r), data=data)
  if (tilt) X_delta <- model.matrix(formula(models$tilt), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_r <- model.matrix(beta_r_formula, data=data)
  X_beta_hr <- cbind(X_beta_h, X_beta_r)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  if (estimate_p_r) d_alpha_r <- ncol(X_alpha_r)
  if (tilt) d_delta <- ncol(X_delta)
  d_h <- ncol(X_beta_h)
  d_r <- ncol(X_beta_r)
  d_hr <- d_h + d_r
  d <- d_h + d_r
  if (estimate_p_r) d <- d + d_alpha_r
  if (tilt) d <- d + d_delta
  
  # Extract coefficients
  if (estimate_p_r) alpha_r <- coef(models$p_r)
  if (tilt) delta <- coef(models$tilt)
  beta_hr <- coef(models$wcls)
  beta_h <- beta_hr[seq(d_h)]
  beta_r <- beta_hr[seq(d_h+1, d_h+d_r)]
  
  # Construct position vectors
  curr_idx <- 0
  if (estimate_p_r) {
    pos_alpha_r <- seq(d_alpha_r)
    curr_idx <- curr_idx + d_alpha_r
  }
  if (tilt) {
    pos_delta <- curr_idx + seq(d_delta)
    curr_idx <- curr_idx + d_delta
  }
  pos_beta_h <- curr_idx + seq(d_h)
  pos_beta_r <- max(pos_beta_h) + seq(d_r)
  pos_beta_hr <- c(pos_beta_h, pos_beta_r)
  
  # p_r_hat score/hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  p_r_hat <- data$p_r_hat
  if (estimate_p_r) {
    scores[, pos_alpha_r] <- (a - p_r_hat) * X_alpha_r
    sd_p_r_hat <- sqrt(p_r_hat * (1-p_r_hat))
    X_alpha_r_scaled <- sd_p_r_hat * X_alpha_r
    hessian[pos_alpha_r, pos_alpha_r] <- crossprod(X_alpha_r_scaled)
  }
  
  # Tilt scores and Hessian
  if (tilt) {
    prop_internal <- mean(is_internal)
    rho <- (1 - prop_internal) / prop_internal
    p_delta_num <- rho * c(exp(X_delta %*% delta))
    p_delta <- p_delta_num / (1 + p_delta_num)
    scores[, pos_delta] <- (is_internal - p_delta) * X_delta
    X_delta_weighted <- X_delta * sqrt(p_delta * (1 - p_delta))
    hessian[pos_delta, pos_delta] <- crossprod(X_delta_weighted)
  }

  # WCLS scores and Hessian
  p_r_hat_a <- data$p_r_hat_a
  w <- data$w
  tilt_ratios <- data$tilt_ratios
  w_and_tilt <- data$w_and_tilt
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_r_fitted_values <- c(X_beta_r %*% beta_r)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_r_fitted_values)
  wcls_weighted_resids <- w * tilt_ratios * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_r] <- wcls_weighted_resids * X_beta_r
  
  X_beta_hr_scaled <- sqrt(w_and_tilt) * X_beta_hr
  GtWG <- crossprod(X_beta_hr_scaled)
  hessian[pos_beta_hr, pos_beta_hr] <- GtWG

  if (estimate_p_r) {
    # Note: This doesn't handle the setting of 3+ treatments
    X_beta_r_raw <- X_beta_r / a_centered
    p_r_hat_a_deriv <- -(2*a-1) * p_r_hat * (1 - p_r_hat)
    log_p_r_hat_a_deriv <- p_r_hat_a_deriv / p_r_hat_a
    p_r_deriv <- -(1-p_r_hat)
    p_r_X_beta_r <- p_r_hat * X_beta_r_raw
    hessian[pos_beta_hr, pos_alpha_r] <- 
      t(X_beta_hr * wcls_weighted_resids) %*% log_p_r_hat_a_deriv +
      t(cbind(matrix(0, nrow=n, ncol=d_h), -p_r_X_beta_r) * wcls_weighted_resids) %*% p_r_deriv +
      t(X_beta_hr * (p_r_hat * wcls_r_fitted_values / a_centered * w * tilt_ratios)) %*% p_r_deriv
  }

  if (tilt) {
    hessian[pos_beta_hr, pos_delta] <- -t(is_external * wcls_weighted_resids * X_beta_hr) %*% X_delta
  }
  
  # Assemble sandwich
  n_users <- length(unique(data$user_id))
  if (is_balanced) {
    t_max <- round(n / n_users)
    sandwich <- construct_sandwich_balanced(scores, hessian, n_users, t_max, d)
  } else {
    sandwich <- construct_sandwich_unbalanced(scores, hessian, data$user_id, n_users, d)
  }

  sandwich
}

wcls <- function(data, beta_r_true, beta_h_formula, beta_r_formula, is_balanced=TRUE, p_r_formula=NULL, tilt_formula=NULL) {
  
  estimate_p_r <- !is.null(p_r_formula)
  tilt <- !is.null(tilt_formula)
  # Get point estimates
  # Propensity score model
  if (estimate_p_r) {
    p_r_mod <- glm(p_r_formula, data=data, family=binomial())
    alpha_r <- coef(p_r_mod)
    data$p_r_hat <- predict(p_r_mod, newdata=data, type="response")
  } else {
    data$p_r_hat <- data$p_r
  }
  data$a_centered <- data$a - data$p_r_hat
  data$p_r_hat_a <- data$a * data$p_r_hat + (1 - data$a) * (1 - data$p_r_hat)
  data$w <- data$p_r_hat_a / data$p_h_a
  
  if (tilt) {
    tilt_mod <- glm(
      tilt_formula,
      family=binomial(), data=data)
    delta <- coef(tilt_mod)
    internal_prop <- mean(data$is_internal)
    delta[1] <- delta[1] - log(internal_prop / (1-internal_prop))
    X_delta <- model.matrix(tilt_mod)
    raw_tilt_ratios <- c(exp(X_delta %*% delta))
    data$tilt_ratios <- data$is_internal + data$is_external * raw_tilt_ratios
  } else {
    data$tilt_ratios <- 1
  }
  data$w_and_tilt <- data$w * data$tilt_ratios

  # WCLS
  beta_r_formula_character <- as.character(update(beta_r_formula, . ~ . + 1))[3]
  beta_r_formula_symbol <- rlang::parse_expr(beta_r_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_r_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w_and_tilt)
  last_beta_h_idx <- ncol(model.matrix(beta_h_formula, data=data))
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_r <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  d_r <- length(beta_r)
  
  # Models list
  models <- list(wcls=wcls_mod)
  if (estimate_p_r) models[["p_r"]] <- p_r_mod
  if (tilt) models[["tilt"]] <- tilt_mod
  
  # Total number of parameters
  vector_estimate <- c(beta_h, beta_r)
  if (tilt) vector_estimate <- c(delta, vector_estimate)
  if (estimate_p_r) vector_estimate <- c(alpha_r, vector_estimate)
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich <- wcls_sandwich(data, models, beta_h_formula, beta_r_formula, is_balanced=is_balanced, estimate_p_r=estimate_p_r, tilt=tilt)
  pos_beta_r <- seq(n_params - d_r + 1, n_params)
  var_beta_r <- sandwich[pos_beta_r, pos_beta_r]
  se_beta_r <- sqrt(diag(var_beta_r))
  beta_r_error <- beta_r - beta_r_true
  beta_r_z_scores <- beta_r_error / se_beta_r
  beta_r_chi2 <- beta_r_error %*% solve(var_beta_r, beta_r_error, tol=1e-30)
  
  results <- list(
    beta_r=beta_r,
    se_beta_r=se_beta_r,
    var_beta_r=var_beta_r,
    beta_r_chi2=beta_r_chi2,
    beta_r_z_scores=beta_r_z_scores,
    sandwich=sandwich,
    n=nrow(data),
    p=nrow(sandwich),
    tilt_warning=FALSE
  )
  results
}
