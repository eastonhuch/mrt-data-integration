# Simulation study for
# "Data integration methods for micro-randomized trials"

set.seed(1)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/eastons-method.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
source("~/Documents/research/mrt-data-integration/wcls.R")
require(geepack)
require(abind)

# True value of beta_r
beta_r_true <- c(-5, -1, 0.9, 0.3)
names(beta_r_true) <- c("A - p(A=1 | s)", "[A - p(A=1 | s)] * x1", "[A - p(A=1 | s)] * x1^2", "[A - p(A=1 | s)] * x1^3")

# Function for applying WCLS to data
apply_wcls <- function(data) {
  p_s_mod <- glm(a ~ 1, data=data, family=binomial())
  data$p_s_hat <- predict(p_s_mod, newdata=data, type="response")
  data$p_s_hat_a <- data$a * data$p_s_hat + (1 - data$a) * (1 - data$p_s_hat)
  data$w <- data$p_s_hat_a / data$p_h_a
  data$a_centered <- data$a - data$p_s_hat
  
  # TODO: Update this method to account for uncertainty in p_mod
  model_wcls <- geeglm(
    y ~ x1 + x2 + x3 + I(a_centered) + I(a_centered * x1) + I(a_centered * x1^2) + + I(a_centered * x1^3),
    data=data, weights=w, id=user_id
  )
  model_summary_wcls <- summary(model_wcls)
  vcov_wcls <- vcov(model_wcls)
  wcls_beta_r_idx <- seq(5, 8)
  beta_r_wcls <- model_wcls$coefficients[wcls_beta_r_idx]
  se_wcls <- sqrt(diag(vcov_wcls))[wcls_beta_r_idx]
  covered_wcls <- (
    (beta_r_true >= beta_r_wcls - 1.96 * se_wcls) & 
      (beta_r_true <= beta_r_wcls + 1.96 * se_wcls))
  
  results <- cbind(
    estimate=beta_r_wcls,
    se=se_wcls,
    covered=covered_wcls
  )
  results
}

# Function for simulating one data set and applying methods
simulate_one <- function(n_internal, n_external) {
  dat <- generate_data(t_max=20, dof=10, n_internal=n_internal, n_external=n_external, ar_param=0.5)
  
  # Walter method
  model_walter <- walters_method(dat)
  covered_walter <- (
    (beta_r_true >= model_walter$beta_r - 1.96 * model_walter$se_beta_r) & 
    (beta_r_true <= model_walter$beta_r + 1.96 * model_walter$se_beta_r))
  results_walter <- cbind(
    estimate=model_walter$beta_r,
    se=model_walter$se_beta_r,
    covered=covered_walter
  )
  
  # Internal only
  dat_internal <- dat[dat$is_internal, ]
  results_internal <- apply_wcls(dat_internal)
  
  # Naive pooling
  results_naive <- apply_wcls(dat)
  
  # Bind results together
  results <- abind(
    results_walter,
    results_internal,
    results_naive,
    along=3
  )
  dimnames(results)[[3]] <- c("P-WCLS", "WCLS, Internal Only", "WCLS, Naive Pooling")
  results
}

# Run simulation
simulate_all <- function(n_internal, n_external, n_replications) {
  results <- replicate(n_replications, simulate_one(n_internal=n_internal, n_external=n_external))
  
  # Process results
  coverage <- apply(results[,"covered",,], MARGIN=c(1, 2), FUN=mean)
  avg_estimate <- apply(results[,"estimate",,], MARGIN=c(1, 2), FUN=mean)
  bias <- avg_estimate - beta_r_true
  empirical_se <- sqrt(apply(
    (results[,"estimate",,] - (avg_estimate %x% array(1, dim=c(1, 1, n_replications))))^2,
    MARGIN=c(1, 2),
    FUN=mean))
  empirical_relative_efficiency <- empirical_se / empirical_se[, "P-WCLS"]
  avg_analytical_se <- apply(results[,"se",,], MARGIN=c(1, 2), FUN=mean)
  analytical_relative_efficiency <- avg_analytical_se / avg_analytical_se[, "P-WCLS"]
  mse <- apply((results[,"estimate",,] - beta_r_true)^2, MARGIN=c(1, 2), FUN=mean)
  rmse <- sqrt(mse)
  
  # Create results list
  list(
    results=results,
    coverage=coverage,
    avg_estimate=avg_estimate,
    bias=bias,
    empirical_se=empirical_se,
    empirical_relative_efficiency=empirical_relative_efficiency,
    avg_analytical_se=avg_analytical_se,
    analytical_relative_efficiency=analytical_relative_efficiency,
    mse=mse,
    rmse=rmse
  )
}

# Make table
# Let's do this one parameter at a time
create_pretty_table <- function(result_list) {
  n_methods <- ncol(result_list$coverage)
  n_coefs <- nrow(result_list$coverage)
  coef_table <- data.frame()
  
  coef_table <- cbind(
    rep(beta_r_true, each=n_methods),
    c(t(result_list$avg_estimate)),
    c(t(result_list$empirical_se)),
    c(t(result_list$empirical_relative_efficiency)),
    c(t(result_list$avg_analytical_se)),
    c(t(result_list$analytical_relative_efficiency)),
    c(t(result_list$rmse)),
    c(t(result_list$coverage))
  )
  coef_table <- round(coef_table, digits=3) %>% as.data.frame()
  coef_table <- cbind(
    rep(names(beta_r_true), each=n_methods),
    rep(colnames(result_list$coverage), times=n_coefs),
    coef_table)
  row.names(coef_table) <- NULL
  colnames(coef_table) <- c(
    "Coefficient Name",
    "Method",
    "True Value",
    "Avg Estimate",
    "Empirical Standard Error",
    "Empirical Relative Efficiency",
    "Analytical Standard Error",
    "Analytical Relative Efficiency",
    "rMSE",
    "Coverage (95% Nominal)"
  )
  coef_table  
}

n_replications <- 400

# 25 external data points
results_25 <- simulate_all(100, 25, n_replications)
View(create_pretty_table(results_25))

# 50 external data points
results_50 <- simulate_all(100, 50, n_replications)
View(create_pretty_table(results_50))

# 100 external data points
results_100 <- simulate_all(100, 100, n_replications)
View(create_pretty_table(results_100))

# 200 external data points
results_200 <- simulate_all(100, 200, n_replications)
View(create_pretty_table(results_200))
