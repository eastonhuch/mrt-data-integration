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
require(ggplot2)
require(forcats)
require(xtable)
require(scales)

# True value of beta_r
beta_r_true <- c(-5, -1, 0.9, 0.3)
names(beta_r_true) <- c("Intercept", "Linear", "Quadratic", "Cubic")
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
# n_replications <- 400
# sample_sizes <- c(25, 100, 400, 1600, 6400)
# result_df <- NULL
# results_25_25 <- NULL
# for (n_internal in sample_sizes) {
#   for (n_external in sample_sizes) {
#     results_i <- simulate_all(n_internal, n_external, n_replications)
#     result_df_i <- create_pretty_table(results_i)
#     if (is.null(result_df)) {
#       result_df <- result_df_i
#     } else {
#       result_df <- rbind(result_df, result_df_i)
#     }
#     
#     if ((n_internal == 25) && (n_external == 25)) {
#       results_25_25 <- results_i
#     }
#   }
# }
# colnames(result_df) <- colnames(result_df_i)

# Checkpoint result dataframe
result_df_file <- "~/Documents/research/mrt-data-integration/simulation_results.csv"
write.csv(result_df, file=result_df_file, row.names=FALSE)
result_df <- read.csv(result_df_file)
colnames(result_df) <- colnames(result_df_i)

# Checkpoint results_25_25
results_25_25_file <- "~/Documents/research/mrt-data-integration/results_25_25.RData"
save(results_25_25, file=results_25_25_file)
load(results_25_25_file)

# Plot effect of increasing external sample size
unbiased_method_names <- method_names[method_names != "WCLS, Pooled"]
result_df_25_internal <- result_df %>% filter(`Internal Sample Size` == 25)
pdf(file="~/Documents/research/mrt-data-integration/sample_size_se.pdf",
    width=8, height=5.4)
par(mfrow=c(2, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1, cex.lab=1.2)
subplot_names <- names(beta_r_true)
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- names(beta_r_true)[coef_counter]
  max_se <- result_df_25_internal %>% filter(
    Method == method,
    `Coefficient Name` == coef_name) %>%
    pull(`Empirical Standard Error`) %>%
    max()
  plot(NULL, type="n", xlim=c(20, 6800), ylim=c(-(0.3 * max_se), 1.2*max_se), log="x", xaxt="n",
       xlab="External Sample Size", ylab="", main=subplot_names[coef_counter])
  abline(h=0, col="gray", lty=2)
  if (coef_counter == 1) {
    title(ylab="Standard Error", line=2.7)
  }
  axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=sample_sizes[c(1, 3, 5)])
  axis(side=1, at=sample_sizes[c(2, 4)], labels=sample_sizes[c(2, 4)])
  if (coef_counter == 1) {
    legend("bottomleft", legend=unbiased_method_names_r, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names), bg="white")
  }
  method_counter <- 0
  for (method in unbiased_method_names) {
    method_counter <- method_counter + 1
    result_df_25_internal_i <- result_df_25_internal %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`External Sample Size`)  
    lines(
      result_df_25_internal_i$`External Sample Size`,
      result_df_25_internal_i$`Empirical Standard Error`,
      type="b", pch=method_counter, col=method_counter
    )
  }
}
#dev.off()

# Plot effect of increasing internal sample size
result_df_25_external <- result_df %>% filter(`External Sample Size` == 25)
#pdf(file="Documents/research/mrt-data-integration/internal_sample_size_se.pdf", width=8, height=2.7)
#par(mfrow=c(1, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1, cex.lab=1.2)
subplot_names <- c("(a) Intercept", "(b) Linear Term", "(c) Quadratic Term", "(d) Cubic Term")
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- names(beta_r_true)[coef_counter]
  max_se <- result_df_25_external %>% filter(
    Method == method,
    `Coefficient Name` == coef_name) %>%
    pull(`Empirical Standard Error`) %>%
    max()
  plot(NULL, type="n", xlim=c(20, 6800), ylim=c(-(0.3 * max_se), 1.2*max_se), log="x", xaxt="n",
       xlab="Internal Sample Size", ylab="")#main=subplot_names[coef_counter])
  abline(h=0, col="gray", lty=2)
  if (coef_counter == 1) {
    title(ylab="Standard Error", line=2.7)
  }
  axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=sample_sizes[c(1, 3, 5)])
  axis(side=1, at=sample_sizes[c(2, 4)], labels=sample_sizes[c(2, 4)])
  # if (coef_counter == 1) {
  #   legend("bottomleft", legend=unbiased_method_names, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names), bg="white")
  # }
  method_counter <- 0
  for (method in unbiased_method_names) {
    method_counter <- method_counter + 1
    result_df_25_external_i <- result_df_25_external %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`Internal Sample Size`)  
    lines(
      result_df_25_external_i$`Internal Sample Size`,
      result_df_25_external_i$`Empirical Standard Error`,
      type="b", pch=method_counter, col=method_counter
    )
  }
}
dev.off()

# Plot effect of increasing both sample sizes
result_df_balanced_samples <- result_df %>% filter(`External Sample Size` == `Internal Sample Size`)
pdf(file="Documents/research/mrt-data-integration/sample_size_efficiency.pdf",
    width=8, height=2.7)
par(mfrow=c(1, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1, cex.lab=1.2)
subplot_names <- c("(a) Intercept", "(b) Linear Term", "(c) Quadratic Term", "(d) Cubic Term")
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- names(beta_r_true)[coef_counter]
  plot(NULL, type="n", xlim=c(20, 6800), ylim=c(-0.3, 1.35), log="x", xaxt="n",
       xlab="Sample Sizes", ylab="", main=subplot_names[coef_counter])
  abline(h=0, col="gray", lty=2)
  if (coef_counter == 1) {
    title(ylab="Relative Efficiency", line=2.7)
  }
  axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=sample_sizes[c(1, 3, 5)])
  axis(side=1, at=sample_sizes[c(2, 4)], labels=sample_sizes[c(2, 4)])
  if (coef_counter == 1) {
    legend("bottomleft", legend=unbiased_method_names, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names), bg="white")
  }
  method_counter <- 0
  for (method in unbiased_method_names) {
    method_counter <- method_counter + 1
    result_df_balanced_samples_i <- result_df_balanced_samples %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`Internal Sample Size`)  
    lines(
      result_df_balanced_samples_i$`Internal Sample Size`,
      result_df_balanced_samples_i$`Empirical Relative Efficiency`,
      type="b", pch=method_counter, col=method_counter
    )
  }
}
dev.off()

# Add grouped side-by-side boxplots
estimates_25_25 <- results_25_25$results[,"estimate",,]
readable_coef_names <- c("Intercept", "Linear", "Quadratic", "Cubic")
adjust_method_names <- function(x) {
  method_names <- c("P-WCLS, Pooled", "P-WCLS, Internal", "WCLS, Pooled", "WCLS, Internal")
  names(method_names) <- c("P-WCLS, Pooled", "P-WCLS, Internal Only", "WCLS, Pooled", "WCLS, Internal Only")
  method_names[x]
}
dimnames(estimates_25_25)[[1]] <- readable_coef_names
estimates_25_25_df <- data.frame()
for (coef_name in readable_coef_names) {
  for (method in method_names) {
    estimates_25_25_df_i <- data.frame(Estimate=estimates_25_25[coef_name, method,])
    estimates_25_25_df_i["Method"] <- method
    estimates_25_25_df_i["MethodNumber"] <- which.max(method_names == method)
    estimates_25_25_df_i["Coefficient"] <- coef_name
    estimates_25_25_df_i["CoefficientNumber"] <- which.max(readable_coef_names == coef_name)
    estimates_25_25_df <- rbind(
      estimates_25_25_df,
      estimates_25_25_df_i
    )
  }
}
boxplots_25_25_df <- estimates_25_25_df %>% mutate(
  Method = adjust_method_names(Method)
) %>% mutate(
  Method = fct_reorder(Method, MethodNumber),
  Coefficient = fct_reorder(Coefficient, CoefficientNumber)
)
pdf("~/Documents/research/mrt-data-integration/estimates_25_25.pdf", width=10, height=3)
ggplot(boxplots_25_25_df, aes(x=Coefficient, y=Estimate, fill=Method)) + 
  geom_boxplot()
dev.off()

# Tables
print_exact_number_nicely <- function(x, digits=1) {
  x_rounded <- round(x)
  if (abs(x - x_rounded) < 1e-6) {
    result <- as.character(x_rounded)
  } else {
    result <- as.character(round(x, 1))
  }
  result
}


result_table <- result_df %>% filter(
    `Internal Sample Size` == 100,
    `External Sample Size` == 100,
  ) %>% mutate(
    `Coefficient Name` = sapply(`Coefficient Name`, function(x_i) paste0("$", x_i, "$")),
    Method = adjust_method_names(Method),
    `True Value` = sapply(`True Value`, print_exact_number_nicely),
    `Relative Efficiency` = label_percent(accuracy=0.1)(`Empirical Relative Efficiency`),
    `Coverage` = label_percent(accuracy=0.1)(`Coverage (95% Nominal)`),
    `Avg Estimate` = round(`Avg Estimate`, 2),
    `rMSE` = round(`rMSE`, 2),
  ) %>% select(
    `True Value`,
    Method,
    `Avg Estimate`,
    `Relative Efficiency`,
    `rMSE`,
    `Coverage`
  ) %>% xtable(
    label="integration:tab:simulation_results",
    caption="Results from the simulation with 100 individuals in both the internal and external studies. As expected, all methods except WCLS, Pooled generate approximately unbiased estimates. P-WCLS, Pooled is the most efficient estimator for all coefficients except the intercept. All methods except WCLS, Pooled achieve near-nominal coverage, with the two versions of P-WCLS outperforming WCLS in most comparisons."
  )
print(
  result_table,
  sanitize.text.function = function(x) gsub("\\%", "\\\\\\%", x),
  include.rownames = FALSE
)  
