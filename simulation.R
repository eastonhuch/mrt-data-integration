# Simulation study for
# "Data integration methods for micro-randomized trials"

random_seed <- 0 # Other seeds are set within function simulate_all()
set.seed(random_seed)
source("~/Documents/research/mrt-data-integration/generate_data.R")
source("~/Documents/research/mrt-data-integration/walters-method.R")
source("~/Documents/research/mrt-data-integration/wcls.R")
source("~/Documents/research/mrt-data-integration/et-wcls.R")
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
method_names <- c("P-WCLS-Pooled", "P-WCLS-Pooled-OBS", "P-WCLS-Internal", "WCLS-Pooled", "WCLS-Internal", "ET-WCLS-Equal", "ET-WCLS")

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
  
  # P-WCLS-Pooled
  model_pwcls_pooled <- walters_method(dat)
  results_pwcls_pooled <- process_results(model_pwcls_pooled)
  
  # P-WCLS-Pooled-OBS
  model_pwcls_pooled_obs <- walters_method(dat, observational=TRUE)
  results_pwcls_pooled_obs <- process_results(model_pwcls_pooled_obs)

  # P-WCLS-Internal
  model_pwcls_internal <- walters_method(dat, internal_only=TRUE)
  results_pwcls_internal <- process_results(model_pwcls_internal)
  
  # WCLS-Pooled
  model_wcls_pooled <- wcls(dat)
  results_wcls_pooled <- process_results(model_wcls_pooled)
  
  # WCLS-Internal
  dat_internal <- dat[dat$is_internal, ]
  model_wcls_internal <- wcls(dat_internal)
  results_wcls_internal <- process_results(model_wcls_internal)
  
  # ET-WCLS-Equal
  model_et_wcls_equal <- wcls(dat, tilt=TRUE)
  results_et_wcls_equal <- process_results(model_et_wcls_equal)

  # ET-WCLS
  model_et_wcls <- etwcls(dat)
  results_et_wcls <- process_results(model_et_wcls)
  
  # Bind results together
  results <- abind(
    results_pwcls_pooled,
    results_pwcls_pooled_obs,
    results_pwcls_internal,
    results_wcls_pooled,
    results_wcls_internal,
    results_et_wcls_equal,
    results_et_wcls,
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
  avg_estimate <- apply(results[,"estimate",,], MARGIN=c(1, 2), FUN=mean)
  bias <- avg_estimate - beta_r_true
  empirical_se <- sqrt(apply(
    (results[,"estimate",,] - (avg_estimate %x% array(1, dim=c(1, 1, n_replications))))^2,
    MARGIN=c(1, 2),
    FUN=mean))
  empirical_relative_efficiency <- empirical_se[, "P-WCLS-Pooled"] / empirical_se
  avg_analytical_se <- apply(results[,"se",,], MARGIN=c(1, 2), FUN=mean)
  analytical_relative_efficiency <- avg_analytical_se[, "P-WCLS-Pooled"] / avg_analytical_se
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
    "Coverage (95% Nominal)"
  )
  coef_table  
}

# Run simulation across many sample sizes
n_replications <- 400
# sample_sizes <- c(25, 100, 400, 1600, 6400)
sample_sizes <- c(400)
result_df <- NULL
results_25_25 <- NULL
for (n_internal in sample_sizes) {
  for (n_external in sample_sizes) {
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
  }
}
colnames(result_df) <- colnames(result_df_i)

# For recreating just results_25_25
# results_25_25 <- simulate_all(25, 25, n_replications)

# Checkpoint result dataframe
result_df_file <- "~/Documents/research/mrt-data-integration/simulation_results.csv"
write.csv(result_df, file=result_df_file, row.names=FALSE)
result_df <- read.csv(result_df_file, check.names=FALSE)
colnames(result_df) <- colnames(result_df_i)

# Checkpoint results_25_25
# results_25_25_file <- "~/Documents/research/mrt-data-integration/results_25_25.RData"
# save(results_25_25, file=results_25_25_file)
# load(results_25_25_file)

# Checkpoint results_100_100
# results_100_100_file <- "~/Documents/research/mrt-data-integration/results_100_100.RData"
# save(results_100_100, file=results_100_100_file)
# load(results_100_100_file)

results_400_400_file <- "~/Documents/research/mrt-data-integration/results_400_400.RData"
save(results_400_400, file=results_400_400_file)
load(results_400_400_file)

# Plot effect of increasing external sample size
unbiased_method_names <- method_names[method_names != "WCLS-Pooled"]
# result_df_25_internal <- result_df %>% filter(`Internal Sample Size` == 25)
# pdf(file="~/Documents/research/mrt-data-integration/sample_size_se.pdf",
#     width=8, height=5.4)
# par(mfrow=c(2, 4), mai=c(0.6, 0.5, 0.5, 0.07), cex.main=1.5, cex.axis=1, cex.lab=1.2)
# subplot_names <- coef_names
# for (coef_counter in seq_along(beta_r_true)) {
#   coef_name <- coef_names[coef_counter]
#   max_se <- result_df_25_internal %>% 
#     filter(`Coefficient Name` == coef_name) %>%
#     pull(`Empirical Standard Error`) %>%
#     max()
#   plot(NULL, type="n", xlim=c(20, 6800), ylim=c(-(0.3 * max_se), 1.2*max_se), log="x", xaxt="n",
#        xlab="External Sample Size", ylab="", main=subplot_names[coef_counter])
#   abline(h=0, col="gray", lty=2)
#   if (coef_counter == 1) {
#     title(ylab="Standard Error", line=2.7)
#   }
#   axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=sample_sizes[c(1, 3, 5)])
#   axis(side=1, at=sample_sizes[c(2, 4)], labels=sample_sizes[c(2, 4)])
#   if (coef_counter == 1) {
#     legend("bottomleft", legend=unbiased_method_names, pch=seq_along(unbiased_method_names), col=seq_along(unbiased_method_names), bg="white")
#   }
#   method_counter <- 0
#   for (method in unbiased_method_names) {
#     method_counter <- method_counter + 1
#     result_df_25_internal_i <- result_df_25_internal %>% filter(
#       Method == method,
#       `Coefficient Name` == coef_name) %>% arrange(`External Sample Size`)  
#     lines(
#       result_df_25_internal_i$`External Sample Size`,
#       result_df_25_internal_i$`Empirical Standard Error`,
#       type="b", pch=method_counter, col=method_counter
#     )
#   }
# }

# Plot effect of increasing internal sample size
# result_df_25_external <- result_df %>% filter(`External Sample Size` == 25)
# for (coef_counter in seq_along(beta_r_true)) {
#   coef_name <- coef_names[coef_counter]
#   max_se <- result_df_25_external %>% 
#     filter(`Coefficient Name` == coef_name) %>%
#     pull(`Empirical Standard Error`) %>%
#     max()
#   plot(NULL, type="n", xlim=c(20, 6800), ylim=c(-(0.3 * max_se), 1.2*max_se), log="x", xaxt="n",
#        xlab="Internal Sample Size", ylab="")
#   abline(h=0, col="gray", lty=2)
#   if (coef_counter == 1) {
#     title(ylab="Standard Error", line=2.7)
#   }
#   axis(side=1, at=sample_sizes[c(1, 3, 5)], labels=sample_sizes[c(1, 3, 5)])
#   axis(side=1, at=sample_sizes[c(2, 4)], labels=sample_sizes[c(2, 4)])
#   method_counter <- 0
#   for (method in unbiased_method_names) {
#     method_counter <- method_counter + 1
#     result_df_25_external_i <- result_df_25_external %>% filter(
#       Method == method,
#       `Coefficient Name` == coef_name) %>% arrange(`Internal Sample Size`)  
#     lines(
#       result_df_25_external_i$`Internal Sample Size`,
#       result_df_25_external_i$`Empirical Standard Error`,
#       type="b", pch=method_counter, col=method_counter
#     )
#   }
# }
# dev.off()

# Add grouped side-by-side boxplots
# estimates_25_25 <- results_25_25$results[,"estimate",,]
# dimnames(estimates_25_25)[[1]] <- coef_names
# estimates_25_25_df <- data.frame()
# for (coef_name in coef_names) {
#   for (method in method_names) {
#     estimates_25_25_df_i <- data.frame(Estimate=estimates_25_25[coef_name, method,])
#     estimates_25_25_df_i["Method"] <- method
#     estimates_25_25_df_i["MethodNumber"] <- which.max(method_names == method)
#     estimates_25_25_df_i["Coefficient"] <- coef_name
#     estimates_25_25_df_i["CoefficientNumber"] <- which.max(coef_names == coef_name)
#     estimates_25_25_df <- rbind(
#       estimates_25_25_df,
#       estimates_25_25_df_i
#     )
#   }
# }
# boxplots_25_25_df <- estimates_25_25_df %>% mutate(
#   Method = fct_reorder(Method, MethodNumber),
#   Coefficient = fct_reorder(Coefficient, CoefficientNumber)
# )
# pdf("~/Documents/research/mrt-data-integration/estimates_25_25.pdf", width=10, height=3)
# boxplot_linewidth <- 0.3
# ggplot() + 
#   geom_boxplot(data=boxplots_25_25_df, aes(x=Coefficient, y=Estimate, fill=Method)) +
#   geom_line(
#     data=data.frame(x=c(0.5, 1.5), y=rep(-2, 2)),
#     aes(x=x, y=y, group=1),
#     linewidth=boxplot_linewidth) +
#   geom_line(
#     data=data.frame(x=c(1.5, 2.5), y=rep(-5, 2)),
#     aes(x=x, y=y, group=1),
#     linewidth=boxplot_linewidth)
# dev.off()

# Plot confidence intervals for first two methods
estimates_400_400 <- results_400_400$results[,"estimate",,]

pdf("~/Documents/research/mrt-data-integration/x1_se_plot.pdf",
    width=6.5, height=2.5)
par(mfrow=c(1, 3), mai=c(0.6, 0.6, 0.5, 0.07), cex.main=1.5, cex.axis=0.9, cex.lab=1.2)

# dat <- generate_data(n_internal=100000, n_external=100000)
# max_abs_x1 <- round(max(abs(dat$x1)), 1) + 0.2
# x1_values <- seq(-max_abs_x1, max_abs_x1, 0.2)
# hist(dat$x1, breaks=x1_values, probability=TRUE, xlab=expression(X[1]), main="(a) Histogram")
# 
# x1_design_matrix <- cbind(1, x1_values)
# plot(NULL, type="n", xlim=c(-max_abs_x1, max_abs_x1), ylim=c(-50, 125),
#      xlab=expression(X[1]), ylab="Causal Effect",
#      main="(b) 95% CIs")
# dash_lty <- 2
# x1_plot_method_numbers <- c(1, 3)
# for (method_number in x1_plot_method_numbers) {
#   method_name <- method_names[method_number]
#   method_estimates <- estimates_400_400[,method_name,]
#   method_fitted_values <- x1_design_matrix %*% method_estimates
#   mean_method_fitted_values <- rowMeans(method_fitted_values)
#   lines(x1_values, mean_method_fitted_values, col=method_number)
#   lines(x1_values, apply(method_fitted_values, MARGIN=1, function(x) quantile(x, probs=0.025)), col=method_number, lty=dash_lty)
#   lines(x1_values, apply(method_fitted_values, MARGIN=1, function(x) quantile(x, probs=0.975)), col=method_number, lty=dash_lty)
# }
# legend("topleft", legend=method_names[x1_plot_method_numbers], lty=1, col=x1_plot_method_numbers, cex=0.85)
# 
# plot(NULL, type="n", xlim=c(-max_abs_x1, max_abs_x1), ylim=c(0, 29),
#      xlab=expression(X[1]), ylab="Standard Error",
#      main="(c) SE Comparison")
# for (method_number in x1_plot_method_numbers) {
#   method_name <- method_names[method_number]
#   method_estimates <- estimates_400_400[,method_name,]
#   method_fitted_values <- x1_design_matrix %*% method_estimates
#   method_ses <- apply(method_fitted_values, MARGIN=1, sd)
#   lines(x1_values, method_ses, col=method_number)
# }
# legend("topleft", legend=method_names[x1_plot_method_numbers], lty=1, col=x1_plot_method_numbers, cex=0.85)
# 
# dev.off()

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

result_table <- result_df %>% filter(
    `Internal Sample Size` == 400,
    `External Sample Size` == 400,
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
        result_table[i, "Coefficient Name"] <- paste0("\\multirow{4}{*}{", coef_name_i, "}") 
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
    caption="Results from the simulation with 400 individuals in both the internal and external studies. As expected, all methods except WCLS-Pooled generate approximately unbiased estimates. P-WCLS-Pooled is the most efficient estimator for all coefficients except the intercept. All methods except WCLS-Pooled achieve near-nominal coverage, with the two versions of P-WCLS outperforming WCLS in most comparisons."
  ) %>%
  print(
    sanitize.text.function = function(x) gsub("\\%", "\\\\\\%", x),
    include.rownames = FALSE) %>%
  str_replace("&  &  \\\\\\\\ \\n", "&  &  \\\\\\\\\n\\\\hline\n") %>% 
  str_replace("\\\\hline\\n &  &  ", " &  &  ") %>%
  str_replace("\\n  \\\\multirow\\{4\\}\\{\\*\\}\\{Slope\\}", "\n \\\\hline \n  \\\\multirow{4}{*}{Slope}")
cat(xtable_results)
