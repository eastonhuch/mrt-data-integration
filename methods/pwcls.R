pwcls_sandwich <- function(data_pooled, data_internal, models, beta_h_formula, beta_s_formula, wcls_s_causal_effects, X_beta_s_raw, X_beta_s_raw_internal_list, X_beta_r, estimate_p_s=TRUE, observational=FALSE, is_balanced=TRUE) {
  # Extract some columns
  y <- data_pooled$y
  p_h_a <- data_pooled$p_h_a
  w <- data_pooled$w
  a <- data_pooled$a
  a_centered <- data_pooled$a_centered
  is_internal <- data_pooled$is_internal
  
  # Construct design matrices
  if (observational) X_alpha_h <- model.matrix(formula(models$p_h), data=data_pooled)
  if (estimate_p_s) X_alpha_s <- model.matrix(formula(models$p_s), data=data_pooled)
  X_beta_h <- model.matrix(beta_h_formula, data=data_pooled)
  X_beta_s <- model.matrix(beta_s_formula, data=data_pooled)
  X_beta_hs <- cbind(X_beta_h, X_beta_s)
  
  # Store dimensions
  n <- nrow(X_beta_h)
  if (observational) d_alpha_h <- ncol(X_alpha_h)
  if (estimate_p_s) d_alpha_s <- ncol(X_alpha_s)
  d_h <- ncol(X_beta_h)
  d_s <- ncol(X_beta_s)
  d_hs <- d_h + d_s
  d_r <- ncol(X_beta_r)
  n_estimated_treatment_levels <- ncol(wcls_s_causal_effects)
  d <- d_h + d_s + d_r * n_estimated_treatment_levels
  if (observational) d <- d + d_alpha_h
  if (estimate_p_s) d <- d + d_alpha_s
  
  # Extract coefficients
  if (observational) alpha_h <- coef(models$p_h)
  if (estimate_p_s) alpha_s <- coef(models$p_s)
  beta_hs <- coef(models$wcls)
  beta_h <- beta_hs[seq(d_h)]
  beta_s <- beta_hs[seq(d_h+1, d_h+d_s)]
  beta_r <- c(coef(models$r))
  
  # Construct position vectors
  curr_idx <- 0
  if (observational) {
    pos_alpha_h <- seq(d_alpha_h)
    curr_idx <- d_alpha_h
  }
  if (estimate_p_s) {
    pos_alpha_s <- curr_idx + seq(d_alpha_s)
    curr_idx <- curr_idx + d_alpha_s
  }
  pos_beta_h <- curr_idx + seq(d_h)
  pos_beta_s <- max(pos_beta_h) + seq(d_s)
  pos_beta_hs <- c(pos_beta_h, pos_beta_s)
  pos_beta_r <- max(pos_beta_s) + seq(d_r * n_estimated_treatment_levels)
  
  # scores/Hessian
  scores <- matrix(0, nrow=n, ncol=d)
  hessian <- matrix(0, nrow=d, ncol=d)
  
  # p_h_hat score/hessian
  if (observational) {
    p_h_hat <- 1 / (1 + exp(-c(X_alpha_h %*% alpha_h)))
    scores[, pos_alpha_h] <- (a - p_h_hat) * X_alpha_h
    sd_p_h_hat <- sqrt(p_h_hat * (1-p_h_hat))
    X_alpha_h_scaled <- sd_p_h_hat * X_alpha_h
    hessian[pos_alpha_h, pos_alpha_h] <- crossprod(X_alpha_h_scaled)
  }
  
  # p_s_hat score/hessian
  if (estimate_p_s) {
    p_s_hat <- 1 / (1 + exp(-c(X_alpha_s %*% alpha_s)))
    scores[, pos_alpha_s] <- (a - p_s_hat) * X_alpha_s
    sd_p_s_hat <- sqrt(p_s_hat * (1-p_s_hat))
    X_alpha_s_scaled <- sd_p_s_hat * X_alpha_s
    hessian[pos_alpha_s, pos_alpha_s] <- crossprod(X_alpha_s_scaled)
  } else {
    p_s_hat <- data_pooled$p_s_hat
  }
  
  # WCLS scores and Hessian
  p_s_hat_a <- data_pooled$p_s_hat_a
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
  if (estimate_p_s) {
    # NOTE: This portion only works for binary treatments
    p_s_hat_a_deriv <- -(2*a-1) * p_s_hat * (1 - p_s_hat)
    log_p_s_hat_a_deriv <- p_s_hat_a_deriv / p_s_hat_a
    log_p_s_deriv <- -(1-p_s_hat)
    p_s_X_beta_s <- p_s_hat * X_beta_s_raw
    hessian[pos_beta_hs, pos_alpha_s] <- 
      t(X_beta_hs * wcls_weighted_resids) %*% log_p_s_hat_a_deriv +
      t(cbind(matrix(0, nrow=n, ncol=d_h), -p_s_X_beta_s) * wcls_weighted_resids) %*% log_p_s_deriv +
      t(X_beta_hs * (p_s_hat * wcls_s_fitted_values / a_centered * w)) %*% log_p_s_deriv
  }
  
  if (observational) {
    p_h_hat_a_deriv <- -(2*a-1) * p_h_hat * (1 - p_h_hat)
    log_p_h_hat_a_deriv <- p_h_hat_a_deriv / p_h_a  # NOTE: This is really p_h_hat_a
    hessian[pos_beta_hs, pos_alpha_h] <- 
      t(X_beta_hs * wcls_weighted_resids) %*% (log_p_h_hat_a_deriv * X_alpha_h)
  }

  # beta_r score
  idx_beta_r_j <- 0
  idx_beta_s_j <- 0
  d_s_each <- round(d_s / n_estimated_treatment_levels)
  for (j in seq(n_estimated_treatment_levels)) {
    idx_beta_r_j <- max(idx_beta_r_j) + seq(d_r)
    pos_beta_r_j <- pos_beta_r[idx_beta_r_j]
    idx_beta_s_j <- max(idx_beta_s_j) + seq(d_s_each)
    pos_beta_s_j <- pos_beta_s[idx_beta_s_j]
    beta_r_j <- beta_r[idx_beta_r_j]
    scores[is_internal, pos_beta_r_j] <- (wcls_s_causal_effects[,j] - c(X_beta_r %*% beta_r_j)) * X_beta_r
    hessian[pos_beta_r_j, pos_beta_r_j] <- crossprod(X_beta_r)
    hessian[pos_beta_r_j, pos_beta_s_j] <- -t(X_beta_r) %*% X_beta_s_raw_internal_list[[j]]
  }
  
  # Assemble sandwich
  n_users <- length(unique(data_pooled$user_id))
  if (is_balanced) {
    t_max <- round(n / n_users)
    sandwich <- construct_sandwich_balanced(scores, hessian, n_users, t_max, d)
  } else {
    sandwich <- construct_sandwich_unbalanced(scores, hessian, data_pooled$user_id, n_users, d)
  }
  sandwich
}

pwcls <- function(data, beta_r_true, beta_h_formula, beta_s_formula, r_formula, p_s_formula=NULL, internal_only=FALSE, p_h_formula=NULL, is_balanced=TRUE, beta_s_formula_divider_idx=numeric(0)) {
  estimate_p_s <- !is.null(p_s_formula)
  observational <- !is.null(p_h_formula)

  # Create data_pooled, data_internal
  if (internal_only) {
    data_pooled <- data[data$is_internal,]
  } else {
    data_pooled <- data
  }  
  # Get point estimates
  # p_h
  if (observational) {
    if (internal_only) stop("observational==TRUE and internal_only==TRUE not implemented")
    p_h_mod <- glm(p_h_formula, data=data_pooled, family=binomial())
    alpha_h <- coef(p_h_mod)
    # Note: This doens't overwrite the true p_h because data_pooled is a copy of the original data
    data_pooled$p_h <- predict(p_h_mod, newdata=data_pooled, type="response")
    data_pooled$p_h_a <- data_pooled$a * data_pooled$p_h + (1 - data_pooled$a) * (1 - data_pooled$p_h)
  }
  
  # p_s
  if (estimate_p_s) {
    p_s_mod <- glm(p_s_formula, data=data_pooled, family=binomial())
    alpha_s <- coef(p_s_mod)
    data_pooled$p_s_hat <- predict(p_s_mod, newdata=data_pooled, type="response")
    data_pooled$a_centered <- data_pooled$a - data_pooled$p_s_hat
  } else {
    data_pooled$p_s_hat <- data_pooled$p_s
  }
  data_pooled$a_centered <- data_pooled$a - data_pooled$p_s_hat
  
  if (("p_s_a" %in% colnames(data_pooled)) && (!estimate_p_s)) {
    print("Using known p_s_a")
    data_pooled$p_s_hat_a <- data_pooled$p_s_a
  } else {
    data_pooled$p_s_hat_a <- data_pooled$a * data_pooled$p_s_hat + (1 - data_pooled$a) * (1 - data_pooled$p_s_hat) 
  }

  if (("wcls_weight" %in% colnames(data_pooled)) && (!estimate_p_s)) {
    print("Using wcls_weight")
    data_pooled$w <- data_pooled$wcls_weight
  } else {
    data_pooled$w <- data_pooled$p_s_hat_a / data_pooled$p_h_a
  }

  # WCLS
  beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
  beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
  wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
  wcls_mod <- lm(wcls_formula, data=data_pooled, weights=w)
  last_beta_h_idx <- ncol(model.matrix(beta_h_formula, data=data_pooled))
  beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
  beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
  d_s <- length(beta_s)
  wcls_s_causal_effects <- numeric(0)
  X_beta_s_raw <- model.matrix(beta_s_formula, data=data_pooled)
  X_beta_s_raw_internal <- X_beta_s_raw[data_pooled$is_internal,]
  X_beta_s_raw_internal_list <- list()
  list_counter <- 1
  r_formula_start_idx <- 1
  for (divider_idx in c(beta_s_formula_divider_idx, d_s+1)) {
    j_idx <- seq(r_formula_start_idx, divider_idx-1)
    a_j_centered <- X_beta_s_raw[, j_idx[1]]  # Assume first column is intercept
    X_beta_s_raw[, j_idx] <- X_beta_s_raw[, j_idx] / a_j_centered
    X_beta_s_raw_internal[, j_idx] <- X_beta_s_raw[data_pooled$is_internal, j_idx]
    X_beta_s_raw_internal_list[[list_counter]] <- X_beta_s_raw_internal[, j_idx]
    beta_s_j <- beta_s[j_idx]
    wcls_s_causal_effects_j <- X_beta_s_raw_internal[, j_idx] %*% beta_s_j
    wcls_s_causal_effects <- cbind(wcls_s_causal_effects, wcls_s_causal_effects_j)
    r_formula_start_idx <- divider_idx
    list_counter <- list_counter + 1
  }
  
  # beta_r
  data_internal <- data_pooled[data_pooled$is_internal,]
  X_beta_r <- model.matrix(update(r_formula, NULL ~ .), data=data_internal)
  r_mod <- lm(wcls_s_causal_effects ~ 0 + X_beta_r)
  beta_r <- c(coef(r_mod))  # NOTE: We're flattening a matrix of coefficients here
  d_r <- length(beta_r)

  # Models list
  models <- list(
    wcls=wcls_mod,
    r=r_mod
  )
  if (estimate_p_s) models$p_s <- p_s_mod
  if (observational) models$p_h <- p_h_mod
  
  # Total number of parameters
  vector_estimate <- c()
  if (observational) vector_estimate <- c(alpha_h)
  if (estimate_p_s) vector_estimate <- c(vector_estimate, alpha_s)
  vector_estimate <- c(vector_estimate, beta_h, beta_s, beta_r)
  n_params <- length(vector_estimate)
  
  # Standard errors
  sandwich <- pwcls_sandwich(data_pooled, data_internal, models, beta_h_formula, beta_s_formula, wcls_s_causal_effects, X_beta_s_raw, X_beta_s_raw_internal_list, X_beta_r, estimate_p_s=estimate_p_s, observational=observational, is_balanced=is_balanced)
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
