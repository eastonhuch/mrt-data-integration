petwcls_sandwich <- function(data, models, beta_h_formula, beta_s_formula, beta_r_formula) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  is_internal <- data$is_internal
  is_external <- data$is_external
  data_internal <- data[is_internal,]
  data_external <- data[is_external,]
  
  # Construct design matrices
  X_alpha_s <- model.matrix(formula(models$p_s), data=data)
  X_delta <- model.matrix(formula(models$tilt), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_s <- model.matrix(beta_s_formula, data=data)
  X_beta_s_raw <- X_beta_s / a_centered
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  X_beta_r_internal <- model.matrix(formula(models$r), data=data_internal)
  X_beta_r <- model.matrix(beta_r_formula, data=data)
  X_beta_hr <- cbind(X_beta_h, X_beta_r)
  X_beta_r_raw <- X_beta_r / a_centered
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_s <- ncol(X_alpha_s)
  d_delta <- ncol(X_delta)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_hs <- d_h + d_s
  d_r <- ncol(X_beta_r_internal)
  d <- d_alpha_s + d_delta + 2*d_h + d_s + 3*d_r # 3*d_r b/c et_wcls has coefs for both internal and external
  
  # Extract coefficients
  alpha_s <- coef(models$p_s)
  delta <- coef(models$tilt)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  beta_hr <- coef(models$r_wcls)
  beta_h_r_wcls <- beta_hr[seq(d_h)]
  beta_r_wcls <- beta_hr[-seq(d_h)]
  beta_r <- coef(models$r)
  
  # Construct position vectors
  pos_alpha_s <- seq(d_alpha_s)
  pos_delta <- d_alpha_s + seq(d_delta)
  pos_beta_h <- max(pos_delta) + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_beta_h_r_wcls <- max(pos_beta_s) + seq(d_h)
  pos_beta_r_wcls <- max(pos_beta_h_r_wcls) + seq(2*d_r)
  pos_beta_hr <- c(pos_beta_h_r_wcls, pos_beta_r_wcls)
  pos_beta_r <- max(pos_beta_hr) + seq(d_r)
  
  # scores/Hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  
  # p_s_hat score/hessian
  p_s_hat <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))  # Assume logit link
  scores[, pos_alpha_s] <- (a - p_s_hat) * X_alpha_s
  sd_p_s_hat <- sqrt(p_s_hat * (1-p_s_hat))
  X_alpha_s_scaled <- sd_p_s_hat * X_alpha_s
  hessian[pos_alpha_s, pos_alpha_s] <- crossprod(X_alpha_s_scaled) # Should this be negative?

  # Tilt scores and Hessian
  prop_internal <- mean(is_internal)
  rho <- prop_internal / (1 - prop_internal)
  p_delta_num <- rho * data$raw_tilt_ratios
  p_delta <- p_delta_num / (1 + p_delta_num)
  scores[, pos_delta] <- (is_internal - p_delta) * X_delta
  X_delta_weighted <- X_delta * sqrt(p_delta * (1 - p_delta))
  hessian[pos_delta, pos_delta] <- crossprod(X_delta_weighted) # Should this be negative?
  
  # beta_hs
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
  
  # beta_hr
  tilt_ratios <- data$tilt_ratios
  w_and_tilt <- data$w_and_tilt
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h_r_wcls)
  wcls_r_fitted_values <- c(X_beta_r %*% beta_r_wcls)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_r_fitted_values)
  wcls_weighted_resids <- w * tilt_ratios * wcls_resids
  scores[, pos_beta_h_r_wcls] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_r_wcls] <- wcls_weighted_resids * X_beta_r
  
  X_beta_hr_scaled <- sqrt(w_and_tilt) * X_beta_hr
  GtWG <- crossprod(X_beta_hr_scaled)
  hessian[pos_beta_hr, pos_beta_hr] <- GtWG
  
  p_s_X_beta_r <- p_s_hat * X_beta_r_raw
  hessian[pos_beta_hr, pos_alpha_s] <-
    hessian[pos_beta_hr, pos_alpha_s] +
    t(X_beta_hr * wcls_weighted_resids) %*% log_p_s_hat_a_deriv +
    t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_r) * wcls_weighted_resids) %*% log_p_s_deriv +
    t(X_beta_hr * (p_s_hat * wcls_r_fitted_values / a_centered * w * tilt_ratios)) %*% log_p_s_deriv
  
  hessian[pos_beta_hr, pos_delta] <- -t(is_external * wcls_weighted_resids * X_beta_hr) %*% X_delta
  
  # beta_r
  scores[is_internal, pos_beta_r] <- (data$wcls_s_causal_effects[data$is_internal] - c(X_beta_r_internal %*% beta_r)) * X_beta_r_internal
  hessian[pos_beta_r, pos_beta_r] <- crossprod(X_beta_r_internal)
  hessian[pos_beta_r, pos_beta_s] <- -t(X_beta_r_internal) %*% X_beta_s_raw[data$is_internal,]
  
  # Assemble sandwich
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

petwcls <- function(data) {
  # True coefficient vectors
  beta_s_true <- c(1, 2, -3)
  beta_r_true <- c(-2, 5)
  
  # p_s
  p_s_mod <- glm(a ~ 1, data=data, family=binomial())
  alpha_s <- coef(p_s_mod)
  data$p_s_hat <- predict(p_s_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_s_hat
  data$p_s_hat_a <- data$a * data$p_s_hat + (1 - data$a) * (1 - data$p_s_hat)
  data$w <- data$p_s_hat_a / data$p_h_a

  # Tilting
  # Try simpler model if there's a warning
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
  delta <- coef(tilt_mod)
  tilt_warning <- length(delta) <= 10
  internal_prop <- mean(data$is_internal)
  delta[1] <- delta[1] - log(internal_prop / (1-internal_prop))
  X_delta <- model.matrix(tilt_mod)
  data$raw_tilt_ratios <- c(exp(X_delta %*% delta))
  data$tilt_ratios <- data$is_internal + data$is_external * data$raw_tilt_ratios
  data$w_and_tilt <- data$w * data$tilt_ratios
  
  # beta_hs
  beta_h_formula <- y ~ x1 + x2 + x3
  beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  data$wcls_s_causal_effects <- c(model.matrix(beta_s_formula, data=data) %*% beta_s) / data$a_centered
  d_s <- length(beta_s)

  # beta_hr
  beta_r_formula <- y ~ 0 + I(is_internal * a_centered) + I(is_internal * a_centered * x1) + I(is_external * a_centered) + I(is_external * a_centered * x1)
  beta_r_formula_character <- as.character(update(beta_r_formula, . ~ . + 1))[3]
  beta_r_formula_symbol <- rlang::parse_expr(beta_r_formula_character)
  r_wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_r_formula_symbol)))
  r_wcls_mod <- lm(r_wcls_formula, data=data, weights=w_and_tilt)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
  beta_h_r_wcls <- coef(r_wcls_mod)[ seq(last_beta_h_idx)]
  beta_r_wcls <- coef(r_wcls_mod)[-seq(last_beta_h_idx)]
  
  # beta_r
  data_internal <- data[data$is_internal,]
  r_formula <- wcls_s_causal_effects ~ x1
  r_mod <- glm(r_formula, data=data_internal)
  beta_r <- coef(r_mod)
  d_r <- length(beta_r)
  
  # Models list
  models <- list(
    p_s=p_s_mod,
    tilt=tilt_mod,
    wcls=wcls_mod,
    r_wcls=r_wcls_mod,
    r=r_mod
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_s,
    delta,
    beta_h,
    beta_s,
    beta_h_r_wcls,
    beta_r_wcls,
    beta_r
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich_list <- petwcls_sandwich(data, models, beta_h_formula, beta_s_formula, beta_r_formula)
  sandwich <- sandwich_list$sandwich
  pos_beta_r <- length(alpha_s) + length(delta) + length(beta_h) + length(beta_s) + length(beta_h) + seq(3*d_r)
  var_beta_r <- sandwich[pos_beta_r, pos_beta_r]
  Lambda <- chol2inv(chol(var_beta_r))
  Lambda_sum <- matrix(0, nrow=2, ncol=2)
  for (j in 1:3) {
    for (k in 1:3) {
      Lambda_sum <- Lambda_sum + Lambda[2*(j-1)+seq(2), 2*(k-1)+seq(2)]
    }
  }
  Lambda_sum_inv <- chol2inv(chol(Lambda_sum))
  z <- Lambda %*% c(beta_r_wcls, beta_r)
  z_sum <- z[1:2] + z[3:4] + z[5:6]
  beta_r_pooled <- c(Lambda_sum_inv %*% z_sum)
  Lambda_horiz_sum <- Lambda[1:2,] + Lambda[3:4,] + Lambda[5:6,]
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
    bread=sandwich_list$bread,
    meat=sandwich_list$meat,
    n=nrow(data),
    p=nrow(sandwich),
    tilt_warning=tilt_warning
  )
  results
}
