dr_sandwich <- function(data, models, beta_h_formula, beta_s_formula, is_balanced=TRUE) {
  # Extract some columns
  data_internal <- data[data$is_internal,]
  data_external <- data[data$is_external,]
  is_internal <- data$is_internal
  is_external <- data$is_external
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  
  # Construct design matrices
  X_alpha_s <- model.matrix(formula(models$p_s), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_s <- model.matrix(beta_s_formula, data=data)
  X_beta_s_raw <- X_beta_s / a_centered
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  X_omega <- model.matrix(formula(models$tilt), data=data)
  X_beta_r_internal <- model.matrix(formula(models$r), data=data_internal)
  X_beta_r_external <- model.matrix(formula(models$r), data=data_external)
  X_beta_r <- model.matrix(formula(models$r), data=data)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_s <- ncol(X_alpha_s)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_hs <- d_h + d_s
  d_omega <- ncol(X_omega)
  d_pi_internal <- 1
  d_r <- ncol(X_beta_r_internal)
  d <- d_alpha_s + d_h + d_s + d_omega + d_pi_internal + 2 * d_r
  
  # Extract coefficients
  alpha_s <- coef(models$p_s)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  omega <- coef(models$tilt)
  pi_internal <- models$pi_internal
  beta_r <- coef(models$r)
  beta_r_et <- models$r_et
  
  # Construct position vectors
  pos_alpha_s <- seq(d_alpha_s)
  pos_beta_h <- max(pos_alpha_s) + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_omega <- max(pos_beta_hs) + seq(d_omega)
  pos_pi_internal <- max(pos_omega) + 1
  pos_beta_r <- max(pos_pi_internal) + seq(d_r)
  pos_beta_r_et <- max(pos_beta_r) + seq(d_r)
  
  # scores/Hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  
  # p_s_hat score/hessian
  p_s_hat <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))  # Assume logit link
  scores[, pos_alpha_s] <- (a - p_s_hat) * X_alpha_s
  sd_p_s_hat <- sqrt(p_s_hat * (1-p_s_hat))
  X_alpha_s_scaled <- sd_p_s_hat * X_alpha_s
  hessian[pos_alpha_s, pos_alpha_s] <- crossprod(X_alpha_s_scaled) # Should this be negative?
  
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
  
  # These have negative signs because I'm using (-A^-1) B (-A^-1) in my sandwich
  p_s_hat_a_deriv <- -(2*a-1) * p_s_hat * (1 - p_s_hat)
  log_p_s_hat_a_deriv <- p_s_hat_a_deriv / p_s_hat_a
  log_p_s_deriv <- -(1-p_s_hat)
  p_s_X_beta_s <- p_s_hat * X_beta_s_raw
  hessian[pos_beta_hs, pos_alpha_s] <- 
    t(X_beta_hs * wcls_weighted_resids) %*% log_p_s_hat_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_s) * wcls_weighted_resids) %*% log_p_s_deriv +
    t(X_beta_hs * (p_s_hat * wcls_s_fitted_values / a_centered * w)) %*% log_p_s_deriv
  
  # Tilt scores and Hessian
  prop_internal <- mean(is_internal)
  rho <- prop_internal / (1 - prop_internal)
  p_omega_num <- rho * data$raw_tilt_ratios
  p_omega <- p_omega_num / (1 + p_omega_num)
  scores[, pos_omega] <- (is_internal - p_omega) * X_omega
  X_omega_weighted <- X_omega * sqrt(p_omega * (1 - p_omega))
  hessian[pos_omega, pos_omega] <- crossprod(X_omega_weighted) # Should this be negative?
  
  # pi_internal scores and Hessian
  scores[, pos_pi_internal] <- data$is_internal - pi_internal
  hessian[pos_pi_internal, pos_pi_internal] <- nrow(data)
  
  # beta_r score
  scores[is_internal, pos_beta_r] <- (data_internal$y_tilde - c(X_beta_r_internal %*% beta_r)) * X_beta_r_internal
  hessian[pos_beta_r, pos_beta_r] <- crossprod(X_beta_r_internal)
  hessian[pos_beta_r, pos_beta_h] <- t(X_beta_r_internal) %*% (X_beta_h[is_internal,] / data_internal$y_tilde_denom)
  hessian[pos_beta_r, pos_beta_s] <- (t(X_beta_r_internal) %*%
    (X_beta_s[is_internal]/data_internal$y_tilde_denom - X_beta_s_raw[is_internal,]))
  
  # DRET-WCLS
  scores[, pos_beta_r_et] <- (
    data$is_external * data$tilt_ratios * (data$y - data$f_h_a) / (data$y_tilde_denom * (1 - pi_internal)) +
      data$is_internal * (data$f_h_1 - data$f_h_0 - c(X_beta_r %*% beta_r_et)) / pi_internal
  ) * X_beta_r
  hessian[pos_beta_r_et, pos_beta_r_et] <- crossprod(X_beta_r_internal) / pi_internal
  hessian[pos_beta_r_et, pos_omega] <- -t(X_beta_r_external) %*% (
    data_external$tilt_ratios * data_external$y_tilde_frac * X_omega[is_external,]) / (1 - pi_internal)
  hessian[pos_beta_r_et, pos_beta_h] <- t(X_beta_r_external) %*% (
    (data_external$tilt_ratios / data_external$y_tilde_denom) * X_beta_h[is_external,]) / (1 - pi_internal)
  hessian[pos_beta_r_et, pos_beta_s] <- (
    t(X_beta_r_external) %*% (data_external$tilt_ratios / data_external$y_tilde_denom * X_beta_s[is_external,]) / (1 - pi_internal) -
    t(X_beta_r_internal) %*% X_beta_s_raw[is_internal,] / pi_internal)
  hessian[pos_beta_r_et, pos_pi_internal] <- colSums(
    (
      -(1-pi_internal)^(-2) * data$is_external * data$tilt_ratios * (data$y - data$f_h_a) / data$y_tilde_denom
      + pi_internal^(-2) * data$is_internal * (data$f_h_1 - data$f_h_0 - c(X_beta_r %*% beta_r_et))
    ) * X_beta_r
  )
  
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

drwcls <- function(data, beta_r_true, beta_h_formula, beta_s_formula, r_formula, p_s_formula, tilt_formula=NULL, is_balanced=TRUE) {
  # Get point estimates
  # p_s
  p_s_mod <- glm(p_s_formula, data=data, family=binomial())
  alpha_s <- coef(p_s_mod)
  data$p_s_hat <- predict(p_s_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_s_hat
  data$p_s_hat_a <- data$a * data$p_s_hat + (1 - data$a) * (1 - data$p_s_hat)
  data$w <- data$p_s_hat_a / data$p_h_a
  
  # S-moderated model
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_hs <- coef(wcls_mod)
  beta_h <- beta_hs[ seq(last_beta_h_idx)]
  beta_s <- beta_hs[-seq(last_beta_h_idx)]
  data$wcls_s_causal_effects <- c(model.matrix(beta_s_formula, data=data) %*% beta_s) / data$a_centered
  d_s <- length(beta_s)
  data_1 <- data
  data_1$a <- 1
  data_1$a_logical <- TRUE
  data_1$a_centered <- 1 - data_1$p_s_hat
  data_0 <- data
  data_0$a <- 0
  data_0$a_logical <- FALSE
  data_0$a_centered <- - data_1$p_s_hat
  data$f_h_1 <- predict(wcls_mod, newdata = data_1)
  data$f_h_0 <- predict(wcls_mod, newdata = data_0)
  data$f_h_a <- data$a * data$f_h_1 + (1 - data$a) * data$f_h_0
  data$y_tilde_denom <- data$a - (1 - data$p_h)
  data$y_tilde_frac <- (data$y - data$f_h_a) / data$y_tilde_denom
  data$y_tilde <- data$y_tilde_frac + data$wcls_s_causal_effects
  # NOTE: Not currently using weights (\Tilde{\sigma}) for this one
  
  # Tilting Model
  # Try simpler model if there's a warning
  if (is.null(tilt_formula)) {
    tilt_mod <- tryCatch(
      glm(
        is_internal ~ bs(x1, df=3, degree=2)*I(bs(x2, df=3, degree=2)),
        family=binomial(), data=data),
      warning=function(w) tryCatch(
        glm(
          is_internal ~ bs(x1, df=2, degree=2)*I(bs(x2, df=2, degree=2)),
          family=binomial(), data=data),
        warning=function(w) glm(
          is_internal ~ bs(x1, df=1, degree=1)*I(bs(x2, df=1, degree=1)),
          family=binomial(), data=data)
      )
    )
    tilt_warning <- length(coef(tilt_mod)) <= 10
  } else {
    tilt_mod <- glm(tilt_formula, family=binomial(), data=data)
    tilt_warning <- FALSE
  }
  omega <- coef(tilt_mod)
  tilt_warning <- length(omega) <= 10
  pi_internal <- mean(data$is_internal)
  omega[1] <- omega[1] - log(pi_internal / (1-pi_internal))
  X_omega <- model.matrix(tilt_mod)
  data$raw_tilt_ratios <- c(exp(X_omega %*% omega))
  data$tilt_ratios <- data$is_internal + data$is_external * data$raw_tilt_ratios
  data$w_and_tilt <- data$w * data$tilt_ratios
  
  # DRP-WCLS
  data_internal <- data[data$is_internal,]
  r_mod <- glm(r_formula, data=data_internal)
  beta_r <- coef(r_mod)
  
  # DRET-WCLS
  data_external <- data[data$is_external,]
  r_formula_updated <- update(r_formula, y_tilde ~ .)
  X_beta_r_internal <- model.matrix(r_formula_updated, data=data_internal)
  XtX_beta_r_internal <- crossprod(X_beta_r_internal)
  X_beta_r_external <- model.matrix(r_formula_updated, data=data_external)
  beta_r_et <- c(solve(XtX_beta_r_internal / pi_internal, (
    (t(X_beta_r_internal) %*% data_internal$wcls_s_causal_effects) / pi_internal +
      (t(X_beta_r_external) %*% (data_external$tilt_ratios * data_external$y_tilde_frac)) / (1 - pi_internal)
  ), tol=1e-50))

  # Models list
  models <- list(
    p_s=p_s_mod,
    wcls=wcls_mod,
    tilt=tilt_mod,
    pi_internal=pi_internal,
    r=r_mod,
    r_et=beta_r_et
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_s,
    beta_h,
    beta_s,
    omega,
    beta_r,
    beta_r_et
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich <- dr_sandwich(data, models, beta_h_formula, beta_s_formula, is_balanced=is_balanced)
  pos_beta_r <- length(alpha_s) + length(beta_hs) + length(omega) + 1 + seq(2*length(beta_r))
  var_beta_r <- sandwich[pos_beta_r, pos_beta_r]
  Lambda <- chol2inv(chol(var_beta_r))
  half_d_r <- length(beta_r)
  first_half <- 1:half_d_r
  second_half <- (half_d_r+1):(2*half_d_r)
  Lambda_sum <- (
    Lambda[first_half, first_half] + Lambda[first_half, second_half] +
      Lambda[second_half, first_half] + Lambda[second_half, second_half])
  Lambda_sum_inv <- chol2inv(chol(Lambda_sum))
  z <- Lambda %*% c(beta_r, beta_r_et)
  z_sum <- z[first_half] + z[second_half]
  beta_r_pooled <- c(Lambda_sum_inv %*% z_sum)
  Lambda_horiz_sum <- Lambda[first_half,] + Lambda[second_half,]
  var_beta_r_pooled <- Lambda_sum_inv %*% Lambda_horiz_sum %*% var_beta_r %*% t(Lambda_horiz_sum) %*% Lambda_sum_inv
  se_beta_r_pooled <- sqrt(diag(var_beta_r_pooled))
  beta_r_pooled_error <- beta_r_pooled - beta_r_true
  beta_r_pooled_z_scores <- beta_r_pooled_error / se_beta_r_pooled
  beta_r_pooled_chi2 <- beta_r_pooled_error %*% solve(var_beta_r_pooled, beta_r_pooled_error)
  
  results <- list(
    beta_r=beta_r_pooled,
    se_beta_r=se_beta_r_pooled,
    var_beta_r=var_beta_r_pooled,
    beta_r_chi2=beta_r_pooled_chi2,
    beta_r_z_scores=beta_r_pooled_z_scores,
    sandwich=sandwich,
    n=nrow(data),
    p=nrow(sandwich),
    tilt_warning=tilt_warning
  )
  results
}

