# Simulation study for
# "Data integration methods for micro-randomized trials"

set.seed(1)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/eastons-method.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
source("~/Documents/research/mrt-data-integration/wcls.R")
require(geepack)
require(abind)
require(tidyverse)

# True value of beta_r
beta_r_true <- c(-5, -1, 0.9, 0.3)
names(beta_r_true) <- c("A - p(A=1 | s)", "[A - p(A=1 | s)] * x1", "[A - p(A=1 | s)] * x1^2", "[A - p(A=1 | s)] * x1^3")
method_names <- c("P-WCLS, Pooled", "P-WCLS, Internal Only", "WCLS, Pooled", "WCLS, Internal Only")

process_results <- function(model) {
  covered <- (
    (beta_r_true >= model$beta_r - 1.96 * model$se_beta_r) & 
      (beta_r_true <= model$beta_r + 1.96 * model$se_beta_r))
  results <- cbind(
    estimate=model$beta_r,
    se=model$se_beta_r,
    covered=covered
  )
  results
}

# Function for simulating one data set and applying methods
simulate_one <- function(n_internal, n_external) {
  dat <- generate_data(t_max=20, dof=10, n_internal=n_internal, n_external=n_external, ar_param=0.5)
  
  # P-WCLS, Pooled
  model_pwcls_pooled <- walters_method(dat)
  results_pwcls_pooled <- process_results(model_pwcls_pooled)
  
  # P-WCLS, Internal Only
  model_pwcls_internal <- walters_method(dat, internal_only=TRUE)
  results_pwcls_internal <- process_results(model_pwcls_internal)
  
  # WCLS, Pooled
  model_wcls_pooled <- wcls(dat)
  results_wcls_pooled <- process_results(model_wcls_pooled)
  
  # WCLS, Internal Only
  dat_internal <- dat[dat$is_internal, ]
  model_wcls_internal <- wcls(dat_internal)
  results_wcls_internal <- process_results(model_wcls_internal)
  
  # Bind results together
  results <- abind(
    results_pwcls_pooled,
    results_pwcls_internal,
    results_wcls_pooled,
    results_wcls_internal,
    along=3
  )
  dimnames(results)[[3]] <- method_names
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
  empirical_relative_efficiency <- empirical_se[, "P-WCLS, Pooled"] / empirical_se
  avg_analytical_se <- apply(results[,"se",,], MARGIN=c(1, 2), FUN=mean)
  analytical_relative_efficiency <- avg_analytical_se[, "P-WCLS, Pooled"] / avg_analytical_se
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
    rmse=rmse,
    n_internal=n_internal,
    n_external=n_external
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
    rep(result_list$n_internal, each=n_methods*n_coefs),
    rep(result_list$n_external, each=n_methods*n_coefs),
    rep(names(beta_r_true), each=n_methods),
    rep(colnames(result_list$coverage), times=n_coefs),
    coef_table)
  row.names(coef_table) <- NULL
  colnames(coef_table) <- c(
    "Internal Sample Size",
    "External Sample Size",
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

# Run simulation across many sample sizes
n_replications <- 100
sample_sizes <- c(25, 100, 400, 1600)
result_df <- NULL
for (n_internal in sample_sizes) {
  for (n_external in sample_sizes) {
    results_i <- simulate_all(n_internal, n_external, n_replications)
    result_df_i <- create_pretty_table(results_i)
    if (is.null(result_df)) {
      result_df <- result_df_i
    } else {
      result_df <- rbind(result_df, result_df_i)
    }
  }
}
colnames(result_df) <- colnames(result_df_i)

# Checkpoint result dataframe
# result_df_file <- "~/Documents/research/mrt-data-integration/simulation_results.csv"
# write.csv(result_df, file=result_df_file, row.names=FALSE)
# result_df <- read.csv(result_df_file)

# Plot effect of increasing external sample size
unbiased_method_names <- method_names[method_names != "WCLS, Pooled"]
result_df_100_internal <- result_df %>% filter(`Internal Sample Size` == 100)
pdf(file="Documents/research/mrt-data-integration/internal_sample_size_efficiency.pdf",
    width=8, height=2.7)
par(mfrow=c(1, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1.1, cex.lab=1.2)
subplot_names <- c("(a) Intercept", "(b) Linear Term", "(c) Quadratic Term", "(d) Cubic Term")
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- names(beta_r_true)[coef_counter]
  plot(NULL, type="n", xlim=c(20, 1700), ylim=c(0.05, 1.1), log="xy", xaxt="n",
       xlab="External Sample Size", ylab="", main=subplot_names[coef_counter])
  if (coef_counter == 1) {
    title(ylab="Relative Efficiency", line=2.7)
  }
  axis(side=1, at=sample_sizes, labels=sample_sizes)
  if (coef_counter == 1) {
    legend("bottomleft", legend=unbiased_method_names, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names))
  }
  method_counter <- 0
  for (method in unbiased_method_names) {
    method_counter <- method_counter + 1
    result_df_100_internal_i <- result_df_100_internal %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`External Sample Size`)  
    lines(
      result_df_100_internal_i$`External Sample Size`,
      result_df_100_internal_i$`Empirical Relative Efficiency`,
      type="b", pch=method_counter, col=method_counter
    )
  }
}
dev.off()

# Plot effect of increasing internal sample size
result_df_100_external <- result_df %>% filter(`External Sample Size` == 100)
pdf(file="Documents/research/mrt-data-integration/external_sample_size_efficiency.pdf",
    width=8, height=2.7)
par(mfrow=c(1, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1.1, cex.lab=1.2)
subplot_names <- c("(a) Intercept", "(b) Linear Term", "(c) Quadratic Term", "(d) Cubic Term")
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- names(beta_r_true)[coef_counter]
  plot(NULL, type="n", xlim=c(20, 1700), ylim=c(0.05, 1.1), log="xy", xaxt="n",
       xlab="Internal Sample Size", ylab="", main=subplot_names[coef_counter])
  if (coef_counter == 1) {
    title(ylab="Relative Efficiency", line=2.7)
  }
  axis(side=1, at=sample_sizes, labels=sample_sizes)
  if (coef_counter == 1) {
    legend("bottomleft", legend=unbiased_method_names, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names))
  }
  method_counter <- 0
  for (method in unbiased_method_names) {
    method_counter <- method_counter + 1
    result_df_100_external_i <- result_df_100_external %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`Internal Sample Size`)  
    lines(
      result_df_100_internal_i$`Internal Sample Size`,
      result_df_100_internal_i$`Empirical Relative Efficiency`,
      type="b", pch=method_counter, col=method_counter
    )
  }
}
dev.off()
