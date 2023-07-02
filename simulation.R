# Simulation study for
# "Data integration methods for micro-randomized trials"

#set.seed(1)
source("~/Documents/research/mrt-data-integration/old-r-files/generate_data.R")
source("~/Documents/research/mrt-data-integration/old-r-files/sandwich.R")
require(geepack)

# Generate data
run_simulation <- function() {
dat <- generate_data(dof=1e10, ar_param=1e-200, t_max=2, n_internal=10000, n_external=10000)

make_Gamma <- function(gamma_x2) {
  cbind(
    c(1, 0, 0, 0),
    c(0, 1, 0, 0),
    gamma_x2)
}
beta_s_true <- c(1, 2, -3)
gamma_x2_true <- c(2, 1, -0.3, -0.1)
beta_r_true <- c(-5, -1, 0.9, 0.3)
theta_true <- c(beta_s_true, gamma_x2_true)

# Get point estimates
# Propensity score model
p_s_mod <- glm(a ~ 1, data=dat, family=binomial())
alpha_s <- coef(p_s_mod)
dat$p_s_hat <- predict(p_s_mod, newdata=dat, type="response")
dat$a_centered <- dat$a - dat$p_s_hat
dat$p_s_hat_a <- dat$a * dat$p_s_hat + (1 - dat$a) * (1 - dat$p_s_hat)
dat$w <- dat$p_s_hat_a / dat$p_h_a

# WCLS
beta_h_formula <- y ~ x1 + x2 + x3
beta_s_formula <- y ~ 0 + a_centered + I(a_centered * x1) + I(a_centered * x2)
beta_s_formula_character <- as.character(update(beta_s_formula, . ~ . + 1))[3]
beta_s_formula_symbol <- rlang::parse_expr(beta_s_formula_character)
wcls_formula <- update(beta_h_formula, bquote(. ~ . + .(beta_s_formula_symbol)))
wcls_mod <- lm(wcls_formula, data=dat, weights=w)
last_beta_h_idx <- length(attr(terms(beta_h_formula), "term.labels")) + 1
beta_h <- coef(wcls_mod)[ seq(last_beta_h_idx)]
beta_s <- coef(wcls_mod)[-seq(last_beta_h_idx)]
d_s <- length(beta_s)

# Gamma
s_formula <- x2 ~ x1 + I(x1^2) + I(x1^3)
s_mod <- glm(s_formula, data=dat, subset=dat$is_internal)
gamma_x2 <- coef(s_mod)
Gamma <- make_Gamma(gamma_x2)
row.names(Gamma) <- names(gamma_x2)

# beta_r
beta_r <- c(Gamma %*% beta_s)
d_r <- length(beta_r)
names(beta_r) <- row.names(Gamma)

# Models list
models <- list(
  p_s=p_s_mod,
  wcls=wcls_mod,
  s=s_mod
)

# Total number of parameters
vector_estimate <- c(
  alpha_s,
  beta_h,
  beta_s,
  gamma_x2
)
n_params <- length(vector_estimate)

# Standard errors
sandwich <- eastons_sandwich(dat, models, beta_h_formula, beta_s_formula)
pos_theta <- seq(
  length(alpha_s) + length(beta_h) + 1,
  length(vector_estimate)
)
theta <- c(beta_s, gamma_x2)
theta_error <- theta - theta_true
var_theta <- sandwich[pos_theta, pos_theta]
theta_chi2 <- theta_error %*% solve(var_theta, theta_error) # Checked

beta_s_error <- beta_s - beta_s_true
wcls_gee_mod <- geeglm(wcls_formula, data=dat, weights=w, id=user_id)
beta_s_chi2 <- beta_s_error %*% solve(var_theta[1:3, 1:3], beta_s_error) # Checked

gamma_x2_error <- gamma_x2 - gamma_x2_true
gamma_x2_chi2 <- gamma_x2_error %*% solve(var_theta[4:7, 4:7], gamma_x2_error) # Checked

J_theta <- cbind(
  Gamma,
  beta_s[3] * diag(d_r)
)
var_beta_r <- J_theta %*% var_theta %*% t(J_theta)
se_beta_r <- sqrt(diag(var_beta_r))
beta_r_error <- beta_r - beta_r_true
beta_r_z_scores <- beta_r_error / se_beta_r
beta_r_chi2 <- beta_r_error %*% solve(var_beta_r, beta_r_error)

beta_r_chi2
}

run_simulation()
chi2_values <- replicate(100, run_simulation())
mean(chi2_values)
median(chi2_values)
hist(chi2_values, breaks=seq(0, 70)) # This is close but not great
hist(rchisq(100, 4), breaks=seq(0, 70))