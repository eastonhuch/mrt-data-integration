eastons_est_fun <- function(data, models, beta_h_formula, beta_s_formula) {
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
  X_gamma_x2 <- grab_design_matrix(
    data = data,
    rhs_formula = grab_fixed_formula(models$s))
  
  # Construct position vectors
  pos_alpha_s <- seq(ncol(X_alpha_s))
  pos_beta_h <- max(pos_alpha_s) + seq(ncol(X_beta_h))
  pos_beta_s <- max(pos_beta_h) + seq(ncol(X_beta_s_raw))
  pos_gamma_x2 <- max(pos_beta_s) + seq(ncol(X_gamma_x2))
  
  # Grab individual score functions
  p_s_score <- grab_psiFUN(models$p_s, data)
  s_score <- grab_psiFUN(models$s, data)
  
  # Create full score function
  score_fun <- function(theta) {
    # Extract individual coefficient vectors
    alpha_s <- theta[pos_alpha_s]
    beta_h <- theta[pos_beta_h]
    beta_s <- theta[pos_beta_s]
    gamma_x2 <- theta[pos_gamma_x2]
    Gamma <- make_Gamma(gamma_x2)
    
    # p_s score
    score <- numeric(length(theta))
    score[pos_alpha_s] <- p_s_score(alpha_s)
    
    # WCLS score
    p_s <- c(1 / (1 + exp(-X_alpha_s %*% alpha_s)))  # Assume logit link
    p_s_a <- a*p_s + (1-a)*(1-p_s)
    X_beta_s <- (a - p_s_a) * X_beta_s_raw
    w <- p_s_a / p_h_a
    wcls_resid <- c(y - X_beta_h %*% beta_h - X_beta_s %*% beta_s)
    wcls_weighted_resid <- w * wcls_resid
    score[pos_beta_h] <- wcls_weighted_resid * c(X_beta_h)
    score[pos_beta_s] <- wcls_weighted_resid * c(X_beta_s)
    
    # Gamma score
    score[pos_gamma_x2] <- s_score(gamma_x2) * is_internal
    
    score
  }
}