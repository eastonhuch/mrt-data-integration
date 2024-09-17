# Function for simulating data
generate_data <- function(
  t_max=20, dof=10, n_internal=100, n_external=100,
  ar_param=0.5, plot_simulated_data=FALSE, x2_coef=-3, x21sq_coef=0) {
  
  # Internal indicator
  n <- n_internal + n_external
  n_obs <- n*t_max
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
  x2 <- is_internal * (1 - x1 + x21sq_coef * x1^2 + 3 * rt(n_obs, dof)) +
        is_external * (2.7 * rt(n_obs, dof))
  # These values make the distributions similar enough that ET-WCLS works but different
  # enough that you see some differences in efficiency based on the pooling method
  x3 <- -1 + 0.5*x1 - 0.8*x2 + rt(n_obs, dof)
  
  # Plots to check that relationships look right
  if (plot_simulated_data) {
    par(mfrow=c(1, 3))
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
    0.2 + 0.3*is_internal + 0.05*x1 - 0.03*x2 + 0.06*x3))
  a_logical <- runif(n_obs) < p_h
  a <- as.numeric(a_logical)
  p_h_a <- a*p_h + (1-a)*(1-p_h)
  
  # Histogram of treatment probabilities
  if (plot_simulated_data) {
    dev.off()
    hist(p_h)
    readline(prompt="Press [enter] to continue")
  }
  
  # Outcomes
  epsilon <- gen_ar1_matrix()
  treatment_effects <- 1 + 2*x1 + x2_coef*x2
  y <- 4 + 2*x1- 1.5*x1*x2 + 0.4*x3^3 + a*treatment_effects + epsilon
  
  # Marginal treatment effects
  marginal_treatment_effects <- -2 + 5* x1
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
    "is_external"=c(is_external),
    x1=c(x1),
    x2=c(x2),
    x3=c(x3),
    p_h=c(p_h),
    p_h_a=c(p_h_a),
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

# set.seed(20)
# dat <- generate_data(
#     t_max=20, dof=10, n_internal=400, n_external=400,
#     ar_param=0.5, plot_simulated_data=FALSE)
# plot(dat$x1[dat$is_internal], dat$x2[dat$is_internal])
# points(dat$x1[dat$is_external], dat$x2[dat$is_external], col=2)
# hist(dat$p_h[dat$is_internal], breaks=seq(0, 1, 0.02))
# hist(dat$p_h[!dat$is_internal], breaks=seq(0, 1, 0.02))
