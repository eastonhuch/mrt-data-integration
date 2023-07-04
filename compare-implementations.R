# Simulation study for
# "Data integration methods for micro-randomized trials"

#set.seed(1)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/eastons-method.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
require(geepack)

dat <- generate_data(dof=1e10, ar_param=1e-10, t_max=2, n_internal=20, n_external=20)
eastons_method(dat)$beta_r
walters_method(dat)$beta_r

dat$p_s_hat <- mean(dat$a)
dat$a_centered <- dat$a - dat$p_s
dat$p_s_hat_a <- dat$a * dat$p_s_hat + (1 - dat$a) * (1 - dat$p_s_hat)
dat$w <- dat$p_s_hat_a / dat$p_h_a
wcls_mod <- geeglm(
  y ~ x1 + x2 + x3 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2),
 data=dat,
 id=user_id
)
sqrt(diag(vcov(wcls_mod)))
eastons_method(dat)$se_beta_s
sqrt(diag(walters_method(dat)$sandwich))[6:8]

eastons_results <- eastons_method(dat)
walters_results <- walters_method(dat)

##############
# OLD WAY
source("~/Documents/research/mrt-data-integration/geex-implementation.R")

# Generate data
make_Gamma <- function(gamma_x2) {
  cbind(
    c(1, 0, 0, 0),
    c(0, 1, 0, 0),
    gamma_x2)
}

# Get point estimates
# Propensity score model
p_s_mod <- glm(a ~ 1, data=dat, family=binomial())
alpha_s <- coef(p_s_mod)
dat$p_s_hat <- predict(p_s_mod, newdata=dat, type="response")
dat$p_s_hat_a <- dat$a * dat$p_s_hat + (1 - dat$a) * (1 - dat$p_s_hat)
dat$w <- dat$p_s_hat_a / dat$p_h_a

# WCLS
beta_h_formula <- y ~ x1 + x2 + x3
beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
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
beta_r

# Models list
models_list <- list(
  p_s=p_s_mod,
  wcls=wcls_mod,
  s=s_mod
)
models$p_s

# Total number of parameters
vector_estimate <- c(
  alpha_s,
  beta_h,
  beta_s,
  gamma_x2
)
n_params <- length(vector_estimate)

# Standard errors
fitted_model <- m_estimate(
  estFUN = eastons_est_fun, 
  data   = dat,
  root_control = setup_root_control(start = rep(0, n_params)),
  roots=vector_estimate,
  compute_roots = FALSE,
  outer_args = list(models=models_list, beta_h_formula=beta_h_formula, beta_s_formula=beta_s_formula)
)
# These differ slightly in the 5, 1 element
geex_bread <- grab_bread(fitted_model@sandwich_components)
max(abs(geex_bread - eastons_results$bread))

# These are essentially identical
geex_meat <- grab_meat(fitted_model@sandwich_components)
max(abs(geex_meat - eastons_results$meat))

# geex sandwich is formed as expected
geex_reconstructed_sandwich <- solve(geex_bread) %*% geex_meat %*% t(solve(geex_bread))
vcov(fitted_model) / geex_reconstructed_sandwich