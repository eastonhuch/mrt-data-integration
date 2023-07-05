require(geex)
walters_est_fun <- function(data, models, beta_h_formula, beta_s_formula) {
  # Extract quantities of interest
  y <- data$y
  p_h_a <- data$p_h_a
  a <- data$a
  is_internal <- data$is_internal
  
  # Construct design matrices
  X_alpha_s <- grab_design_matrix(
    data = data,
    rhs_formula = grab_fixed_formula(models$p_s))
  X_beta_h <- model.matrix(beta_h_formula, data=data)
  X_beta_s_raw <- model.matrix(beta_s_formula, data=data) / data$a_centered
  X_beta_r <- grab_design_matrix(
    data = data,
    rhs_formula = grab_fixed_formula(models$r))
  
  # Construct position vectors
  pos_alpha_s <- seq(ncol(X_alpha_s))
  pos_beta_h <- max(pos_alpha_s) + seq(ncol(X_beta_h))
  pos_beta_s <- max(pos_beta_h) + seq(ncol(X_beta_s_raw))
  pos_beta_r <- max(pos_beta_s) + seq(ncol(X_beta_r))
  
  # Grab individual score functions
  p_s_score <- grab_psiFUN(models$p_s, data)
  r_score_raw <- grab_psiFUN(models$r, data)
  r_score <- function(beta_r) r_score_raw(beta_r) * is_internal * summary(models$r)$dispersion
  
  # Create full score function
  score_fun <- function(theta) {
    # Extract individual coefficient vectors
    alpha_s <- theta[pos_alpha_s]
    beta_h <- theta[pos_beta_h]
    beta_s <- theta[pos_beta_s]
    beta_r <- theta[pos_beta_r]
    
    # p_s score
    score <- numeric(length(theta))
    score[pos_alpha_s] <- p_s_score(alpha_s)
    
    # WCLS score
    p_s <- c(1 / (1 + exp(-(X_alpha_s %*% alpha_s))))  # Assume logit link
    X_beta_s <- (a - p_s) * X_beta_s_raw
    p_s_a <- a*p_s + (1-a)*(1-p_s)
    w <- p_s_a / p_h_a
    wcls_resid <- c(y - X_beta_h %*% beta_h - X_beta_s %*% beta_s)
    wcls_weighted_resid <- w * wcls_resid
    score[pos_beta_h] <- wcls_weighted_resid * c(X_beta_h)
    score[pos_beta_s] <- wcls_weighted_resid * c(X_beta_s)
    
    # beta_r score
    score[pos_beta_r] <- r_score(beta_r)
    
    score
  }
}