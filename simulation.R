# Simulation study for
# "Data integration methods for micro-randomized trials"

random_seed <- 0 # Other seeds are set within function simulate_all()
set.seed(random_seed)
source("./generate_data.R")
source("./extras/helpers.R")
source("./methods/wcls.R")
source("./methods/pwcls.R")
source("./methods/etwcls.R")
source("./methods/drwcls.R")
source("./methods/petwcls.R")
require(geepack)
require(abind)
require(tidyverse)
require(ggplot2)
require(forcats)
require(xtable)
require(scales)
require(stringr)
require(splines)

# True value of beta_r
beta_r_true <- c(-2, 5)
coef_names <- c("Intercept", "Slope")
names(beta_r_true) <- coef_names
method_names <- c("WCLS-Internal", "WCLS-Pooled", "P-WCLS-Internal", "P-WCLS-Pooled", "P-WCLS-Pooled-Obs", "ET-WCLS-Equal", "ET-WCLS-Kron", "ET-WCLS", "DR-WCLS", "PET-WCLS")
beta_h_formula <- y ~ x1 + x2 + x3
beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
beta_r_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1)
et_beta_h_formula <- y ~ 0 + I(as.numeric(is_internal)) + I(is_internal*x1) + I(is_internal*x2) + I(is_internal*x3) + I(as.numeric(is_external)) + I(is_external*x1) + I(is_external*x2) + I(is_external*x3)
et_beta_r_formula <- y ~ 0 + I(is_internal * a_centered) + I(is_internal * a_centered * x1) + I(is_external * a_centered) + I(is_external * a_centered * x1)
pwcls_r_formula <- wcls_s_causal_effects ~ x1
a_intercept_formula <- a ~ 1
p_h_formula <- a ~ 1 + as.numeric(is_internal) + x1 + x2 + x3

process_results <- function(model) {
  dof <- model$n - model$p
  t_quantile <- qt(0.975, dof)
  covered <- (
    (beta_r_true >= model$beta_r - t_quantile * model$se_beta_r) & 
      (beta_r_true <= model$beta_r + t_quantile * model$se_beta_r))
  results <- cbind(
    estimate=model$beta_r,
    se=model$se_beta_r,
    covered=covered,
    tilt_warning=model$tilt_warning
  )
  results
}

# Function for simulating one data set and applying methods
simulate_one <- function(n_internal, n_external) {
  dat <- generate_data(t_max=20, dof=10, n_internal=n_internal, n_external=n_external, ar_param=0.5)
  
  # WCLS-Internal
  dat_internal <- dat[dat$is_internal, ]
  model_wcls_internal <- wcls(dat_internal, beta_r_true, beta_h_formula, beta_r_formula, p_r_formula=a_intercept_formula)
  results_wcls_internal <- process_results(model_wcls_internal)
  
  # WCLS-Pooled
  model_wcls_pooled <- wcls(dat, beta_r_true, beta_h_formula, beta_r_formula, p_r_formula=a_intercept_formula)
  results_wcls_pooled <- process_results(model_wcls_pooled)

  # P-WCLS-Internal
  model_pwcls_internal <- pwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula, internal_only=TRUE)
  results_pwcls_internal <- process_results(model_pwcls_internal)
    
  # P-WCLS-Pooled
  model_pwcls_pooled <- pwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula)
  results_pwcls_pooled <- process_results(model_pwcls_pooled)
  
  # P-WCLS-Pooled-OBS
  model_pwcls_pooled_obs <- pwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula, p_h_formula=p_h_formula)
  results_pwcls_pooled_obs <- process_results(model_pwcls_pooled_obs)

  # ET-WCLS-Equal
  model_et_wcls_equal <- etwcls(dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="equal")
  results_et_wcls_equal <- process_results(model_et_wcls_equal)
  
  # ET-WCLS-Kron
  model_et_wcls_kron <- etwcls(dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="kronecker")
  results_et_wcls_kron <- process_results(model_et_wcls_kron)

  # ET-WCLS
  model_et_wcls <- etwcls(dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="full")
  results_et_wcls <- process_results(model_et_wcls)
  
  # DR-WCLS
  model_dr_wcls <- drwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, a_intercept_formula)
  results_dr_wcls <- process_results(model_dr_wcls)
  
  # PET-WCLS
  model_pet_wcls <- petwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, et_beta_r_formula, a_intercept_formula)
  results_pet_wcls <- process_results(model_pet_wcls)
  
  # Bind results together
  results <- abind(
    results_wcls_internal,
    results_wcls_pooled,
    results_pwcls_internal,
    results_pwcls_pooled,
    results_pwcls_pooled_obs,
    results_et_wcls_equal,
    results_et_wcls_kron,
    results_et_wcls,
    results_dr_wcls,
    results_pet_wcls,
    along=3
  )
  dimnames(results)[[3]] <- method_names
  results
}

# Run simulation
simulate_all <- function(n_internal, n_external, n_replications) {
  results <- sapply(random_seed + seq(n_replications), simplify="array", FUN=function(seed) {
    print(paste("Random seed:", seed))
    set.seed(seed)
    simulate_one(n_internal=n_internal, n_external=n_external)
  })

  # Process results
  coverage <- apply(results[,"covered",,], MARGIN=c(1, 2), FUN=mean)
  tilt_warnings <- apply(results[,"tilt_warning",,], MARGIN=c(1, 2), FUN=sum)
  avg_estimate <- apply(results[,"estimate",,], MARGIN=c(1, 2), FUN=mean)
  bias <- avg_estimate - beta_r_true
  empirical_se <- sqrt(apply(
    (results[,"estimate",,] - (avg_estimate %x% array(1, dim=c(1, 1, n_replications))))^2,
    MARGIN=c(1, 2),
    FUN=mean))
  empirical_relative_efficiency <- empirical_se[, "WCLS-Internal"] / empirical_se
  avg_analytical_se <- apply(results[,"se",,], MARGIN=c(1, 2), FUN=mean)
  analytical_relative_efficiency <- avg_analytical_se[, "WCLS-Internal"] / avg_analytical_se
  mse <- apply((results[,"estimate",,] - beta_r_true)^2, MARGIN=c(1, 2), FUN=mean)
  rmse <- sqrt(mse)
  
  # Create results list
  list(
    results=results,
    tilt_warnings=tilt_warnings,
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
    c(t(result_list$coverage)),
    c(t(result_list$tilt_warnings))
  )
  coef_table <- round(coef_table, digits=3) %>% as.data.frame()
  coef_table <- cbind(
    rep(result_list$n_internal, each=n_methods*n_coefs),
    rep(result_list$n_external, each=n_methods*n_coefs),
    rep(coef_names, each=n_methods),
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
    "Coverage (95% Nominal)",
    "# Tilt Warnings"
  )
  coef_table  
}

# Run simulation across many sample sizes
loop_start_time <- Sys.time()
n_replications <- 2
result_df <- NULL
results_25_25 <- NULL
sample_sizes <- c(25, 100, 400, 1600, 6400)
sample_size_pairs <- list(
  c(25, 25), c(100, 100), c(400, 400), c(1600, 1600), c(6400, 6400),
  c(100, 25), c(100, 400), c(100, 1600), c(100, 6400),
  c(25, 100), c(400, 100), c(1600, 100), c(6400, 100)
)
# For checking results fast
# sample_size_pairs <- list(c(100, 100), c(20, 200))
# sample_size_pairs <- list(c(100, 100))
# n_replications <- 20

for (sample_size_pair in sample_size_pairs) {
  print(sample_size_pair)
  n_internal <- sample_size_pair[1]
  n_external <- sample_size_pair[2]
  results_i <- simulate_all(n_internal, n_external, n_replications)
  result_df_i <- create_pretty_table(results_i)
  if (is.null(result_df)) {
    result_df <- result_df_i
  } else {
    result_df <- rbind(result_df, result_df_i)
  }

  if ((n_internal == 25) && (n_external == 25)) {
    results_25_25 <- results_i
  }
  
  if ((n_internal == 100) && (n_external == 100)) {
    results_100_100 <- results_i
  }
  
  if ((n_internal == 400) && (n_external == 400)) {
    results_400_400 <- results_i
  }
  
  if ((n_internal == 1600) && (n_external == 1600)) {
    results_1600_1600 <- results_i
  }
  
  if ((n_internal == 6400) && (n_external == 6400)) {
    results_6400_6400 <- results_i
  }
  
}
colnames(result_df) <- colnames(result_df_i)
loop_end_time <- Sys.time()
loop_time <- loop_end_time - loop_start_time
estimated_time_full <- loop_time * 400 / n_replications
print(paste("Estimated time for full run:", estimated_time_full))

# For recreating just results_25_25
# results_25_25 <- simulate_all(25, 25, n_replications)

# Checkpoint result dataframe
result_df_file <- "./results/simulation_results.csv"
write.csv(result_df, file=result_df_file, row.names=FALSE)
result_df <- read.csv(result_df_file, check.names=FALSE)
colnames(result_df) <- colnames(result_df_i)

# Checkpoint results for specific sample sizes
results_25_25_file <- "./results/results_25_25.RData"
save(results_25_25, file=results_25_25_file)
load(results_25_25_file)

results_100_100_file <- "./results/results_100_100.RData"
save(results_100_100, file=results_100_100_file)
load(results_100_100_file)

results_400_400_file <- "./results/results_400_400.RData"
save(results_400_400, file=results_400_400_file)
load(results_400_400_file)

results_1600_1600_file <- "./results/results_1600_1600.RData"
save(results_1600_1600, file=results_1600_1600_file)
load(results_1600_1600_file)

results_6400_6400_file <- "./results/results_6400_6400.RData"
save(results_6400_6400, file=results_6400_6400_file)
load(results_6400_6400_file)

# Plot effect of increasing external sample size
unbiased_method_names <- method_names[method_names != "WCLS-Pooled"]
methods_for_se_plot <- c(
  "WCLS-Internal",
  # "WCLS-Pooled"
  # "P-WCLS-Internal",
  "P-WCLS-Pooled",
  # "P-WCLS-Pooled-Obs",
  # "ET-WCLS-Equal",
  # "ET-WCLS-Kron",
  "ET-WCLS",
  "DR-WCLS",
  "PET-WCLS"
)
method_colors <- c(
  "#5778a4",
  "#e49444",
  "#d1615d",
  "#85b6b2",
  "#6a9f58",
  "#e7ca60",
  "#a87c9f",
  "#f1a2a9",
  "#967662"#,
  #"#b8b0ac"
)
names(method_colors) <- c(
  "WCLS-Internal",
  "WCLS-Pooled",
  "P-WCLS-Internal",
  "P-WCLS-Pooled",
  "P-WCLS-Pooled-Obs",
  # "ET-WCLS-Equal",
  "ET-WCLS-Kron",
  "ET-WCLS",
  "DR-WCLS",
  "PET-WCLS"
)
  
result_df_100_internal <- result_df %>% filter(`Internal Sample Size` == 100)
pdf(file="./figures/sample_size_se.pdf", width=12, height=2.5)
k1 <- 10
k2 <- 1
k3 <- 10
layout(matrix(c(rep(1, k1), rep(2, k2), rep(seq(3,6), each=k3)), 1, k1+k2+4*k3, byrow = TRUE))
par(mai=c(0, 0, 0.2, 0), cex.main=1.5, cex.axis=1, cex.lab=1.2)
plot.new()
legend(x=0.03, y=0.86, legend=methods_for_se_plot, col=method_colors[methods_for_se_plot], lty=1, lwd=2, cex=1.5)
plot.new()
text(x=0.76, y=0.57, "Standard Error", adj=c(0.5, 0.5), srt=90, cex=1.5)

par(mai=c(0.6, 0.4, 0.5, 0.1), cex.main=1.5, cex.axis=1.2, cex.lab=1.5)
subplot_names <- coef_names
y_axticks <- c(0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
y_axtick_labels <- c(".25", ".5", "1", "2", "4", "8")
mains <- list(
  expression(paste("(a) ", "Intercept, ", n[internal], "=", 100)),
  expression(paste("(b) ", "Slope, ", n[internal], "=", 100)),
  expression(paste("(c) ", "Intercept, ", n[external], "=", 100)),
  expression(paste("(d) ", "Slope, ", n[external], "=", 100))
)
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- coef_names[coef_counter]
  all_ses <- result_df_100_internal %>%
    filter(`Coefficient Name` == coef_name) %>%
    pull(`Empirical Standard Error`)
    # pull(`Analytical Standard Error`)
  min_se <- min(all_ses)
  max_se <- max(all_ses)
  plot(NULL, type="n", xlim=c(20, 6800), 
       # ylim=c(0.5*min_se, 1.2*max_se),
       ylim=c(0.2, 8),
       log="xy", xaxt="n", yaxt="n",
       xlab=expression(n[external]), ylab="", main=mains[[coef_counter]])
  axis(side=2, at=y_axticks[c(1, 3, 5)], labels=y_axtick_labels[c(1, 3, 5)])
  axis(side=2, at=y_axticks[c(2, 4, 6)], labels=y_axtick_labels[c(2, 4, 6)])
  axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=comma(sample_sizes[c(1, 3, 5)]))
  axis(side=1, at=sample_sizes[c(2, 4)], labels=comma(sample_sizes[c(2, 4)]))
  # if (coef_counter == 1) {
  #   legend("bottomleft", legend=methods_for_se_plot, col=seq_along(methods_for_se_plot), lty=1, bg="white")
  # }
  method_counter <- 0
  for (method in methods_for_se_plot) {
    method_counter <- method_counter + 1
    result_df_100_internal_i <- result_df_100_internal %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`External Sample Size`)
    lines(
      result_df_100_internal_i$`External Sample Size`,
      result_df_100_internal_i$`Empirical Standard Error`,
      # result_df_100_internal_i$`Analytical Standard Error`,
      type="b",
      col=method_colors[method],
      lwd=2
    )
  }
}

# Plot effect of increasing internal sample size
result_df_100_external <- result_df %>% filter(`External Sample Size` == 100)
for (coef_counter in seq_along(beta_r_true)) {
  coef_name <- coef_names[coef_counter]
  all_ses <- result_df_100_external %>%
    filter(`Coefficient Name` == coef_name) %>%
    pull(`Empirical Standard Error`)
    # pull(`Analytical Standard Error`)
  min_se <- min(all_ses)
  max_se <- max(all_ses)
  plot(NULL, type="n", xlim=c(20, 6800),
       # ylim=c(0.5*min_se, 1.2*max_se), 
       ylim=c(0.2, 8),
       log="xy", xaxt="n", yaxt="n",
       xlab=expression(n[internal]), ylab="", main=mains[2+coef_counter])
  axis(side=2, at=y_axticks[c(1, 3, 5)], labels=y_axtick_labels[c(1, 3, 5)])
  axis(side=2, at=y_axticks[c(2, 4, 6)], labels=y_axtick_labels[c(2, 4, 6)])
  axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=comma(sample_sizes[c(1, 3, 5)]))
  axis(side=1, at=sample_sizes[c(2, 4)], labels=comma(sample_sizes[c(2, 4)]))
  method_counter <- 0
  for (method in methods_for_se_plot) {
    method_counter <- method_counter + 1
    result_df_100_external_i <- result_df_100_external %>% filter(
      Method == method,
      `Coefficient Name` == coef_name) %>% arrange(`Internal Sample Size`)
    lines(
      result_df_100_external_i$`Internal Sample Size`,
      result_df_100_external_i$`Empirical Standard Error`,
      # result_df_100_external_i$`Analytical Standard Error`,
      type="b",
      col=method_colors[method],
      lwd=2
    )
  }
}
dev.off()

# Add grouped side-by-side boxplots
estimates_400_400 <- results_400_400$results[,"estimate",,]
dimnames(estimates_400_400)[[1]] <- coef_names
estimates_400_400_df <- data.frame()
for (coef_name in coef_names) {
  for (method in method_names) {
    estimates_400_400_df_i <- data.frame(Estimate=estimates_400_400[coef_name, method,])
    estimates_400_400_df_i["Estimation Error"] <- estimates_400_400_df_i["Estimate"] - beta_r_true[coef_name]
    estimates_400_400_df_i["Method"] <- method
    estimates_400_400_df_i["MethodNumber"] <- which.max(method_names == method)
    estimates_400_400_df_i["Coefficient"] <- coef_name
    estimates_400_400_df_i["CoefficientNumber"] <- which.max(coef_names == coef_name)
    estimates_400_400_df <- rbind(
      estimates_400_400_df,
      estimates_400_400_df_i
    )
  }
}
boxplots_400_400_df <- estimates_400_400_df %>% mutate(
  Method = fct_reorder(Method, MethodNumber),
  Coefficient = fct_reorder(Coefficient, CoefficientNumber)
) %>% filter(Method != "ET-WCLS-Equal")
# NOTE: I filtered out ET-WCLS-Equal because it makes the y-axis too large

pdf("./figures/estimates_400_400.pdf", width=10, height=3)
boxplot_linewidth <- 0.3
ggplot() +
  geom_boxplot(data=boxplots_400_400_df, aes(x=Coefficient, y=`Estimation Error`, fill=Method)) +
  geom_line(
    data=data.frame(x=c(0.4, 2.6), y=rep(0, 2)),
    aes(x=x, y=y, group=1),
    linewidth=boxplot_linewidth,
    color="black") +
  scale_fill_manual(values=method_colors) +
  ylab(expression(paste("Estimation Error:  ", hat(beta)[r] - beta[r]))) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=13),
  )
dev.off()

# Plot confidence intervals for first two methods
# estimates_100_100 <- results_100_100$results[,"estimate",,]
# 
# pdf("./figures/x1_se_plot.pdf",
#     width=6.5, height=2.5)
# par(mfrow=c(1, 3), mai=c(0.6, 0.6, 0.5, 0.07), cex.main=1.5, cex.axis=0.9, cex.lab=1.2)
# 
# dat <- generate_data(n_internal=100000, n_external=100000)
# max_abs_x1 <- round(max(abs(dat$x1)), 1) + 0.2
# x1_values <- seq(-max_abs_x1, max_abs_x1, 0.2)
# hist(dat$x1, breaks=x1_values, probability=TRUE, xlab=expression(X[1]), main="(a) Histogram")
# 
# x1_design_matrix <- cbind(1, x1_values)
# plot(NULL, type="n", xlim=c(-max_abs_x1, max_abs_x1), ylim=c(-50, 1100),
#      xlab=expression(X[1]), ylab="Causal Effect",
#      main="(b) 95% CIs")
# dash_lty <- 2
# x1_plot_method_numbers <- c(1, 3)
# for (method_number in x1_plot_method_numbers) {
#   method_name <- method_names[method_number]
#   method_estimates <- estimates_100_100[,method_name,]
#   method_fitted_values <- x1_design_matrix %*% method_estimates
#   mean_method_fitted_values <- rowMeans(method_fitted_values)
#   lines(x1_values, mean_method_fitted_values, col=method_number)
#   lines(x1_values, apply(method_fitted_values, MARGIN=1, function(x) quantile(x, probs=0.0100)), col=method_number, lty=dash_lty)
#   lines(x1_values, apply(method_fitted_values, MARGIN=1, function(x) quantile(x, probs=0.975)), col=method_number, lty=dash_lty)
# }
# legend("topleft", legend=method_names[x1_plot_method_numbers], lty=1, col=x1_plot_method_numbers, cex=0.85)
# 
# plot(NULL, type="n", xlim=c(-max_abs_x1, max_abs_x1), ylim=c(0, 29),
#      xlab=expression(X[1]), ylab="Standard Error",
#      main="(c) SE Comparison")
# for (method_number in x1_plot_method_numbers) {
#   method_name <- method_names[method_number]
#   method_estimates <- estimates_100_100[,method_name,]
#   method_fitted_values <- x1_design_matrix %*% method_estimates
#   method_ses <- apply(method_fitted_values, MARGIN=1, sd)
#   lines(x1_values, method_ses, col=method_number)
# }
# legend("topleft", legend=method_names[x1_plot_method_numbers], lty=1, col=x1_plot_method_numbers)
# 
# dev.off()

# Number of warnings by method/sample size
result_df %>%
  group_by(Method, `Internal Sample Size`, `External Sample Size`) %>%
  summarise(`# Tilt Warnings`=sum(`# Tilt Warnings`))
  

# Table
print_exact_number_nicely <- function(x, digits=1) {
  x_rounded <- round(x)
  if (abs(x - x_rounded) < 1e-6) {
    result <- as.character(x_rounded)
  } else {
    result <- as.character(round(x, 1))
  }
  result
}

make_table <- function(table_sample_size, method_vector=method_names) {
  result_table <- result_df %>% filter(
    `Internal Sample Size` == table_sample_size,
    `External Sample Size` == table_sample_size,
    Method %in% method_vector
  ) %>% mutate(
    `True Value (Numeric)` = `True Value`,
    `Relative Efficiency (Numeric)` = `Empirical Relative Efficiency`,
    `Coverage (Numeric)` = `Coverage (95% Nominal)`,
    `Avg Estimate (Numeric)` = `Avg Estimate`,
    `rMSE (Numeric)` = `rMSE`,
    `True Value` = sapply(`True Value`, print_exact_number_nicely),
    `Relative Efficiency` = label_percent(accuracy=0.1)(`Empirical Relative Efficiency`),
    `Coverage` = label_percent(accuracy=0.1)(`Coverage (95% Nominal)`),
    `Avg Estimate` = sprintf("%.2f", `Avg Estimate`),
    `rMSE` = sprintf("%.2f", `rMSE`),
  ) %>% mutate(
    `Relative Efficiency` = ifelse(`Method` == "WCLS-Pooled", "N/A", `Relative Efficiency`),
    `Relative Efficiency (Numeric)` = ifelse(`Method` == "WCLS-Pooled", 0, `Relative Efficiency (Numeric)`)
  )
  result_table$idx <- seq(nrow(result_table))
  
  make_bold <- function(x) paste0("\\textbf{", x, "}")
  
  # Coefficient-specific formatting
  for (coef_name in coef_names) {
    # Best performance
    coef_table <- filter(result_table, `Coefficient Name` == coef_name)
    # Not filtering properly
    
    # Relative efficiency
    best_eff_idx <- result_table %>%
      filter(
        `Coefficient Name` == coef_name,
        `Relative Efficiency (Numeric)` == max(coef_table$`Relative Efficiency (Numeric)`)) %>%
      pull(idx)
    result_table[best_eff_idx, "Relative Efficiency"] <- make_bold(result_table[best_eff_idx, "Relative Efficiency"])
    
    # rMSE
    best_rmse_idx <- result_table %>%
      filter(
        `Coefficient Name` == coef_name,
        `rMSE (Numeric)` == min(coef_table$`rMSE (Numeric)`)) %>%
      pull(idx)
    result_table[best_rmse_idx, "rMSE"] <- make_bold(result_table[best_rmse_idx, "rMSE"])
    
    # Multi-row coefficient names
    first_coef_row <- TRUE
    for (i in seq(nrow(result_table))) {
      coef_name_i <- result_table[i, "Coefficient Name"]
      if (coef_name_i == coef_name) {
        if (first_coef_row) {
          result_table[i, "Coefficient Name"] <- paste0("\\multirow{10}{*}{", coef_name_i, "}") 
        } else {
          result_table[i, "Coefficient Name"] <- ""
        }
        first_coef_row <- FALSE
      }
    }
  }
  
  coverage_mc_error <- 3 * sqrt(0.05 * 0.95 / n_replications)
  for (i in seq(nrow(result_table))) {
    # Bold columns within MC error
    coverage_i <- result_table[i, "Coverage (Numeric)"]
    if (abs(coverage_i - 0.95) < coverage_mc_error) {
      result_table[i, "Coverage"] <- make_bold(result_table[i, "Coverage"])
    }
    
    avg_estimate_i <- result_table[i, "Avg Estimate (Numeric)"]
    true_value_i <- result_table[i, "True Value (Numeric)"]
    se_i <- result_table[i, "Empirical Standard Error"] / sqrt(n_replications)
    if (abs((avg_estimate_i - true_value_i) / se_i) < 3) {
      result_table[i, "Avg Estimate"] <- make_bold(result_table[i, "Avg Estimate"])
    }
  }
  
  result_table <- result_table %>% select(
    `Coefficient Name`,
    `True Value`,
    Method,
    `Avg Estimate`,
    `Relative Efficiency`,
    `rMSE`,
    `Coverage`
  )
  
  colnames(result_table) <- c(
    "\\multirow{2}{*}{\\parbox{1pt}{Coefficient Name}}",
    "\\multirow{2}{*}{\\parbox{25pt}{True Value}}",
    "\\multirow{2}{*}{\\parbox{1pt}{Method}}",
    "\\multirow{2}{*}{\\parbox{42pt}{Avg\\\\Estimate}}",
    "\\multirow{2}{*}{\\parbox{48pt}{Relative\\\\Efficiency}}",
    "\\multirow{2}{*}{\\parbox{28pt}{rMSE}}",
    "\\multirow{2}{*}{\\parbox{40pt}{Coverage}}"
  )
  result_table <- rbind(rep("", ncol(result_table)), result_table)
  
  xtable_results <- xtable(
    result_table,
    label="integration:tab:simulation_results",
    caption=paste(
      "Results from the simulation with",
      table_sample_size,
      "individuals in both the internal and external studies.
For the Avg Estimate and Coverage columns, the boldface indicates values within Monte Carlo error ($3\\sigma$) of the truth.
For the Relative Efficiency and rMSE columns, the boldface indicates the best performance for each coefficient (PET-WCLS in both cases)."
      )
    ) %>% print(
      sanitize.text.function = function(x) gsub("\\%", "\\\\\\%", x),
      include.rownames = FALSE) %>%
    str_replace("&  &  \\\\\\\\ \\n", "&  &  \\\\\\\\\n\\\\hline\n") %>% 
    str_replace("\\\\hline\\n &  &  ", " &  &  ") %>%
    str_replace("\\n  \\\\multirow\\{10\\}\\{\\*\\}\\{Slope\\}", "\n \\\\hline \n  \\\\multirow{10}{*}{Slope}")
  cat(xtable_results)
}

make_table(25)
make_table(100)
table_method_names <- c("WCLS-Internal", "WCLS-Pooled", "P-WCLS-Internal", "P-WCLS-Pooled", "ET-WCLS", "DR-WCLS", "PET-WCLS")
make_table(400, method_vector=table_method_names)  # This is the one shown in the paper
make_table(400)
make_table(1600)
make_table(6400)
