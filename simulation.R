# Simulation study for
# "Data integration methods for micro-randomized trials"

#set.seed(1)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/eastons-method.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
require(geepack)

# Generate data
run_simulation <- function() {
  dat <- generate_data(dof=1e10, ar_param=1e-10, t_max=2, n_internal=1000, n_external=1000)
  c(eastons_method(dat)$beta_r_chi2, walters_method(dat)$beta_r_chi2)
}

chi2_values <- replicate(1000, run_simulation())

hist(chi2_values[1,], breaks=seq(0, 300))
hist(chi2_values[2,], breaks=seq(0, 300))
hist(rchisq(1000, 4), breaks=seq(0, 300))

mean(chi2_values[1,])
mean(chi2_values[2,])

median(chi2_values[1,])
median(chi2_values[2,])

dat <- generate_data(dof=1e10, ar_param=1e-10, t_max=2, n_internal=1000, n_external=1000)
eastons_method(dat)$beta_r
walters_method(dat)$beta_r

dat$p_s <- mean(dat$a)
dat$a_centered <- dat$a - dat$p_s
dat$p_s_a <- dat$a * dat$p_s + (1 - dat$a) * (1 - dat$p_s)
dat$w <- dat$p_s_a / dat$p_h_a
wcls_mod <- geeglm(
  y ~ x1 + x2 + x3 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2),
 data=dat,
 id=user_id
)
sqrt(diag(vcov(wcls_mod)))
eastons_method(dat)$se_beta_s
sqrt(diag(walters_method(dat)$sandwich))[6:8]

eastons_method(dat)$se_beta_r
walters_method(dat)$se_beta_r
