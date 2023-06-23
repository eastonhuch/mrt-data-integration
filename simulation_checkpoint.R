# Simulation study for
# "Data integration methods for micro-randomized trials"
set.seed(1)
require(geex)

# Function for simulating data
gen_data <- function(
  t_max=20, n=100, dof=6, n_internal=50, n_external=50,
  ar_param=0.5, plot_simulated_data=FALSE) {
  
  # Internal indicator
  n_obs <- n*t_max
  n <- n_internal + n_external
  is_internal_user <- seq(n) <= n_internal
  is_external_user <- !is_internal_user
  is_internal <- matrix(rep(is_internal_user, t_max), nrow=n, ncol=t_max)
  is_external <- !is_internal
  
  # Creates xs
  # Create one big AR(1) time series for efficiency; then reshape, discard, and split
  gen_ar1_matrix <- function() {
    t(matrix(
      arima.sim(list(order=c(1,0,0), ar=ar_param), n=2*n_obs),
      nrow=t_max*2, ncol=n
    )[1:t_max,])
  }
  r <- x1 <- gen_ar1_matrix()
  x2 <- x1 + is_internal * (2 - 0.3*x1^2 -0.1*x1^3) + rt(n_obs, dof)
  x3 <- -1 + 0.5*x1 - 0.3*x2 + rt(n_obs, dof)
  
  # Plots to check that relationships look right
  par(mfrow=c(1, 3))
  if (plot_simulated_data) {
    blue <- rgb(78, 121, 167, max = 255, alpha = 255, names = "tab_blue")
    green <- rgb(89, 161, 79, max = 255, alpha = 255, names = "tab_green")
    blue50 <- rgb(78, 121, 167, max = 255, alpha = 127, names = "tab_blue_50")
    green50 <- rgb(89, 161, 79, max = 255, alpha = 127, names = "tab_green_50")
    plot(x1[is_internal], x2[is_internal], col=blue50)
    points(x1[is_external], x2[is_external], col=green50)
    plot(x1[is_internal], x3[is_internal], col=blue50)
    points(x1[is_external], x3[is_external], col=green50)
    plot(x2[is_internal], x3[is_internal], col=blue50)
    points(x2[is_external], x3[is_external], col=green50)
    legend("topleft", legend=c("Internal", "External"), col=c(blue, green), pch=1)
    readline(prompt="Press [enter] to continue")
  }
  
  # Treatments
  p_h <- 1 / (1 + exp(
    0.8 - 0.7*is_internal + 0.1*x1 - 0.3*x2 + 0.2*x3))
  a_logical <- runif(n_obs) < p_h
  a <- as.numeric(a_logical)
  
  # Histogram of treatment probabilities
  if (plot_simulated_data) {
    dev.off()
    hist(p_h)
    readline(prompt="Press [enter] to continue")
  }
  
  # Outcomes
  epsilon <- gen_ar1_matrix()
  treatment_effects <- 1 + 2*x1 - 3*x2
  y <- 4 + 2*x1- 1.5*x1*x2 + 0.4*x3^3 + a*treatment_effects + epsilon
  
  # Marginal treatment effects
  marginal_treatment_effects <- -5 - x1 + 0.9*x1^2 + 0.3*x1^3
  if (plot_simulated_data) {
    plot(r[is_external], treatment_effects[is_external], col=green50, cex=0.5,
      xlab="R", ylab="Treatment Effect", xlim=range(r), ylim=range(treatment_effects))
    points(r[is_internal], treatment_effects[is_internal], cex=0.5, col=blue50)
    r_order <- order(r)
    lines(r[r_order], marginal_treatment_effects[r_order], col=blue)
    legend("topleft", legend=c("Internal", "External"), col=c(blue, green), pch=1)
  }
  
  # Make dataframe
  dat <- data.frame(
    "is_internal"=c(is_internal),
    x1=c(x1),
    x2=c(x2),
    x3=c(x3),
    p_h=c(p_h),
    a_logical=c(a_logical),
    a=c(a),
    epsilon=c(epsilon),
    treatment_effect=c(treatment_effects),
    y=c(y),
    user_id=rep(seq(n), t_max)
  )
  dat$ones <- 1
  dat
}

# Generate data
dat <- gen_data(plot=TRUE)

# Get point estimates
p_s_mod <- glm(a ~ 1, data=dat, family=binomial())
alpha_s <- coef(p_s_mod)
dat$p_s_hat <- predict(p_s_mod, newdata=dat, type="response")
dat$a_centered <- dat$a - dat$p_s_hat
dat$w <- dat$p_s_hat / dat$p_h
wcls_mod <- lm(
  y ~ x1 + x2 + x3 + a_centered + I(a_centered * x1) + I(a_centered * x2),
  data=dat, weights=w)
beta_s <- coef(wcls_mod)[5:7]

# Don't need to estimate these
#s_mod_intercept <- lm(ones+0 ~ x1 + I(x1^2) + I(x1^3), data=dat)
#s_mod_x1 <- lm(x1+0 ~ x1 + I(x1^2) + I(x1^3), data=dat)
s_mod_x2 <- lm(x2+0 ~ x1 + I(x1^2) + I(x1^3), data=dat)
Lambda <- cbind(
  c(1, 0, 0, 0),
  c(0, 1, 0, 0),
  coef(s_mod_x2))
Lambda <- zapsmall(Lambda)
row.names(Lambda) <- names(coef(s_mod_x2))
beta_r <- c(Lambda %*% beta_s)
names(beta_r) <- row.names(Lambda)

# Get standard errors
eastons_est_fun <- function(data, models) {
  y <- data$y
  X_p_s <- grab_design_matrix(
    data = data,
    rhs_formula = grab_fixed_formula(models$p_s))
  X_wcls <- grab_design_matrix(
    data = data,
    rhs_formula = grab_fixed_formula(models$wcls))
  s_models <- models$s_models
  X_s_list <- list()
  for (s in names(s_models)) {
    X_s_list[[s]] <- grab_design_matrix(
      data = data,
      rhs_formula = grab_fixed_formula(s_models[[s]]))
  }
  
  p_s_pos <- 1:ncol(X_p_s)
  wcls_pos <- (max(p_s_pos) + 1):(max(p_s_pos) + ncol(X_wcls))
  last_pos <- max(wcls_pos)
  s_pos_list <- list()
  for (s in names(s_models)) {
    s_pos <- (last_pos + 1):(last_pos + ncol(X_s_list[[s]]))
    last_pos <- max(s_pos)
    s_pos_list[[s]] <- s_pos
  }
  
  p_s_scores <- grab_psiFUN(models$p_s, data)
  wcls_scores <- grab_psiFUN(models$wcls, data) # Will need to modify this because of dependence on p_s
  s_scores_list <- list()
  for (s in names(s_models)) {
    s_scores_list[[s]] <- grab_psiFUN(s_models[[s]], data)
  }
  
  score_fun <- function(theta) {
    # Need to do some magic here to get the right scores for WCLS
    
    scores <- c(
      p_s_scores(theta[p_s_pos]),
      wcls_scores(theta[wcls_pos])
    )
    beta_r <- numeric(length(names(s_models)))
    wcls_rev_coef <- rev(coef(models$wcls))
    i <- 0
    for (s in rev(names(s_models))) {
      i <- i + 1
      s_pos <- s_pos_list[[s]]
      s_score <- s_scores_list[[s]]
      Lambda_s <- theta[s_pos]
      scores <- c(scores, s_score(Lambda_s))
      beta_r <- beta_r + theta_s * wcls_coef[i]
    }
    last_pos <- max(s_pos)
    scores <- c(scores, theta[-(1:last_pos)] - beta_r)
    scores
  }
}

models <- list(
  p_s=p_s_mod,
  wcls=wcls_mod,
  s=s_models
)

n_params <- length(alpha_s) +
  length(beta_s) +
  length(Lambda) +
  length(beta_r)

m_estimate(
  estFUN = eastons_est_fun,
  data = dat,
  root_control = setup_root_control(start = rep(0, n_params)),
  outer_args = list(models = models))