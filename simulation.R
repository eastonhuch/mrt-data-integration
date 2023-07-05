# Simulation study for
# "Data integration methods for micro-randomized trials"

set.seed(1)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/eastons-method.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
require(geepack)

# Generate data
n <- 1000
run_simulation <- function() {
  dat <- generate_data(dof=1e10, ar_param=1e-10, t_max=2, n_internal=n, n_external=n)
  c(eastons_method(dat)$beta_r_chi2, walters_method(dat)$beta_r_chi2)
}

n_replications <- 1000
chi2_values <- replicate(n_replications, run_simulation())

max_chi2 <- max(chi2_values)
chi2_breaks <- seq(0, max_chi2+1)
chi2_breaks <- seq(0, max(chi2_values[1,])+1)
hist(chi2_values[1,], breaks=chi2_breaks)
hist(chi2_values[2,], breaks=chi2_breaks) # This must be wrong
hist(rchisq(n_replications, 4), breaks=chi2_breaks)

mean(chi2_values[1,])
mean(chi2_values[2,])

median(chi2_values[1,])
median(chi2_values[2,])

dat <- generate_data(dof=1e10, ar_param=1e-10, t_max=2, n_internal=100000, n_external=100000)

dat$p_s <- mean(dat$a)
dat$a_centered <- dat$a - dat$p_s
dat$p_s_a <- dat$a * dat$p_s + (1 - dat$a) * (1 - dat$p_s)
dat$w <- dat$p_s_a / dat$p_h_a
#wcls_mod <- geeglm(
#  y ~ x1 + x2 + x3 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2),
# data=dat,
# id=user_id
#)

eastons_results <- eastons_method(dat)
walters_results <- walters_method(dat)
#sqrt(diag(vcov(wcls_mod)))[5:7]
eastons_results$se_beta_s
sqrt(diag(walters_results$sandwich))[6:8]

eastons_results$beta_r
walters_results$beta_r
eastons_results$se_beta_r
walters_results$se_beta_r
eastons_results$beta_r_z_scores
walters_results$beta_r_z_scores
