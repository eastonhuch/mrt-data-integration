# Need to update below function to account for uncertainty in tilt
etwcls_sandwich <- function(data, models, beta_h_formula, beta_r_formula) {
  # Extract some columns
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  a_centered <- data$a_centered
  is_internal <- data$is_internal
  is_external <- data$is_external
  
  # Construct design matrices
  X_alpha_r <- model.matrix(formula(models$p_r), data=data)
  X_delta <- model.matrix(formula(models$tilt), data=data)
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_r <- model.matrix(beta_r_formula, data=data)
  X_beta_r_raw <- X_beta_r / a_centered
  X_beta_hr <- cbind(X_beta_h, X_beta_r)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  d_alpha_r <- ncol(X_alpha_r)
  d_delta <- ncol(X_delta)
  d_h <- ncol(X_beta_h)
  d_r <- ncol(X_beta_r)
  d_hr <- d_h + d_r
  d <- d_alpha_r + d_delta + d_h + d_r
  
  # Extract coefficients
  alpha_r <- coef(models$p_r)
  delta <- coef(models$tilt)
  beta_hr <- coef(models$wcls)
  beta_h <- beta_hr[seq(d_h)]
  beta_r <- beta_hr[seq(d_h+1, d_h+d_r)]
  
  # Construct position vectors
  pos_alpha_r <- seq(d_alpha_r)
  pos_delta <- d_alpha_r + seq(d_delta)
  pos_beta_h <- max(pos_delta) + seq(d_h)
  pos_beta_r <- max(pos_beta_h) + seq(d_r)
  pos_beta_hr <- c(pos_beta_h, pos_beta_r)
  
  # p_r_hat score/hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  p_r_hat <- 1 / (1 + exp(-c(X_alpha_r %*% alpha_r)))
  scores[, pos_alpha_r] <- (a - p_r_hat) * X_alpha_r
  sd_p_r_hat <- sqrt(p_r_hat * (1-p_r_hat))
  X_alpha_r_scaled <- sd_p_r_hat * X_alpha_r
  hessian[pos_alpha_r, pos_alpha_r] <- crossprod(X_alpha_r_scaled)
  
  # Tilt scores and Hessian
  prop_internal <- mean(is_internal)
  rho <- prop_internal / (1 - prop_internal)
  p_delta_num <- rho * data$raw_tilt_ratios
  p_delta <- p_delta_num / (1 + p_delta_num)
  scores[, pos_delta] <- (is_internal - p_delta) * X_delta
  X_delta_weighted <- X_delta * sqrt(p_delta * (1 - p_delta))
  hessian[pos_delta, pos_delta] <- crossprod(X_delta_weighted)
  
  # WCLS scores and Hessian
  p_r_hat_a <- data$p_r_hat_a
  w <- data$w
  tilt_ratios <- data$tilt_ratios
  w_and_tilt <- data$w_and_tilt
  wcls_h_fitted_values <- c(X_beta_h %*% beta_h)
  wcls_r_fitted_values <- c(X_beta_r %*% beta_r)
  wcls_resids <- c(y - wcls_h_fitted_values - wcls_r_fitted_values)
  wcls_weighted_resids <- w_and_tilt * wcls_resids
  scores[, pos_beta_h] <- wcls_weighted_resids * X_beta_h
  scores[, pos_beta_r] <- wcls_weighted_resids * X_beta_r
  
  X_beta_hr_scaled <- sqrt(w_and_tilt) * X_beta_hr
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
  
  # There must be a big problem with this line:
  hessian[pos_beta_hr, pos_delta] <- -t(is_external * wcls_weighted_resids * X_beta_hr) %*% X_delta
  
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

etwcls <- function(data, pooling_method="full") {
  # pooling_method can be "full", "kronecker", or "equal"
  beta_r_true <- c(-2, 5)
  
  # Propensity score model
  p_r_mod <- glm(a ~ 1, data=data, family=binomial())
  alpha_r <- coef(p_r_mod)
  data$p_r_hat <- predict(p_r_mod, newdata=data, type="response")
  data$a_centered <- data$a - data$p_r_hat
  data$p_r_hat_a <- data$a * data$p_r_hat + (1 - data$a) * (1 - data$p_r_hat)
  data$w <- data$p_r_hat_a / data$p_h_a
  
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
  pi_internal <- mean(data$is_internal)
  delta[1] <- delta[1] - log(pi_internal / (1-pi_internal))
  X_delta <- model.matrix(tilt_mod)
  data$raw_tilt_ratios <- c(exp(X_delta %*% delta))
  data$tilt_ratios <- data$is_internal + data$is_external * data$raw_tilt_ratios
  data$w_and_tilt <- data$w * data$tilt_ratios

  # WCLS
  beta_h_formula <- y ~ 0 + I(as.numeric(is_internal)) + I(is_internal*x1) + I(is_internal*x2) + I(is_internal*x3) + I(as.numeric(is_external)) + I(is_external*x1) + I(is_external*x2) + I(is_external*x3)
  beta_r_formula <- y ~ 0 + I(is_internal * a_centered) + I(is_internal * a_centered * x1) + I(is_external * a_centered) + I(is_external * a_centered * x1)
  beta_r_formula_character <- as.character(update(beta_r_formula, . ~ .))[3]
  beta_r_formula_symbol <- rlang::parse_expr(beta_r_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_r_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data, weights=w_and_tilt)
  last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels"))
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_r <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  d_r <- length(beta_r)
  
  # Models list
  models <- list(
    p_r=p_r_mod,
    wcls=wcls_mod,
    tilt=tilt_mod
  )
  
  # Total number of parameters
  vector_estimate <- c(
    alpha_r,
    delta,
    beta_h,
    beta_r
  )
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich_list <- etwcls_sandwich(data, models, beta_h_formula, beta_r_formula)
  sandwich <- sandwich_list$sandwich
  pos_beta_r <- length(alpha_r) + length(delta) + length(beta_h) + seq_along(beta_r)
  var_beta_r <- sandwich[pos_beta_r, pos_beta_r]
  Lambda <- chol2inv(chol(var_beta_r))
  
  if (pooling_method == "kronecker") {
    Lambda_tilde <- matrix(0, nrow=2, ncol=2)
    Lambda_tilde[1, 1] <- Lambda[1, 1]
    Lambda_tilde[2, 1] <- Lambda[3, 1]
    Lambda_tilde[1, 2] <- Lambda[1, 3]
    Lambda_tilde[2, 2] <- Lambda[3, 3]
    w1 <- sum(Lambda_tilde[,1])
    w2 <- sum(Lambda_tilde[,2])
    w_sum <- w1 + w2
    beta_r_pooled <- (w1 * beta_r[1:2] + w2 * beta_r[3:4]) / (w_sum)
    kron_mat <- t(Lambda_tilde[1,]) %x% diag(2) + t(Lambda_tilde[2,]) %x% diag(2)
    var_beta_r_pooled <- (kron_mat %*% var_beta_r %*% t(kron_mat)) / w_sum^2
  } else if (pooling_method == "full") {
    half_d_r <- round(d_r / 2)
    first_half <- 1:half_d_r
    second_half <- (half_d_r+1):(2*half_d_r)
    Lambda_sum <- (
      Lambda[first_half, first_half] + Lambda[first_half, second_half] +
        Lambda[second_half, first_half] + Lambda[second_half, second_half])
    Lambda_sum_inv <- chol2inv(chol(Lambda_sum))
    z <- Lambda %*% beta_r
    z_sum <- z[first_half] + z[second_half]
    beta_r_pooled <- c(Lambda_sum_inv %*% z_sum)
    Lambda_horiz_sum <- Lambda[first_half,] + Lambda[second_half,]
    var_beta_r_pooled <- Lambda_sum_inv %*% Lambda_horiz_sum %*% var_beta_r %*% t(Lambda_horiz_sum) %*% Lambda_sum_inv
  } else if (pooling_method == "equal") {
    Lambda_tilde <- diag(2) * 0.5
    w1 <- sum(Lambda_tilde[,1])
    w2 <- sum(Lambda_tilde[,2])
    w_sum <- w1 + w2
    beta_r_pooled <- (w1 * beta_r[1:2] + w2 * beta_r[3:4]) / (w_sum)
    kron_mat <- t(Lambda_tilde[1,]) %x% diag(2) + t(Lambda_tilde[2,]) %x% diag(2)
    var_beta_r_pooled <- (kron_mat %*% var_beta_r %*% t(kron_mat)) / w_sum^2
  } else {
    stop("pooling_method must be 'full', 'kronecker', or 'equal'")
  }
  
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
