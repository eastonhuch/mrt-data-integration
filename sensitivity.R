# Sensitivity analysis for
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

sensitivity_analysis <- function(beta_r_true, sens_label, x2_coef=-3, x21sq_coef=0) {
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
    dat <- generate_data(t_max=20, dof=10, n_internal=n_internal, n_external=n_external, ar_param=0.5, x2_coef=x2_coef, x21sq_coef=x21sq_coef)
    
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
    model_pet_wcls <- petwcls(dat, beta_r_true, beta_h_formula, beta_s_formula, et_beta_r_formula, pwcls_r_formula, a_intercept_formula)
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
      # print(paste("Random seed:", seed))
      cat("\r", seed)
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
  n_replications <- 400
  result_df <- NULL
  sample_sizes <- c(400)
  sample_size_pairs <- list(c(400, 400))
  
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
    results_400_400 <- results_i
  }
  colnames(result_df) <- colnames(result_df_i)
  loop_end_time <- Sys.time()
  loop_time <- loop_end_time - loop_start_time
  estimated_time_full <- loop_time * 400 / n_replications
  print(paste("Full run time:", estimated_time_full))
  
  # Checkpoint result dataframe
  result_df_file <- "./results/simulation_results.csv"
  write.csv(result_df, file=result_df_file, row.names=FALSE)
  result_df <- read.csv(result_df_file, check.names=FALSE)
  colnames(result_df) <- colnames(result_df_i)
  
  results_400_400_file <- "./results/results_400_400.RData"
  save(results_400_400, file=results_400_400_file)
  load(results_400_400_file)
  
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
      "\\multirow{2}{*}{\\parbox{1pt}{Coefficient name}}",
      "\\multirow{2}{*}{\\parbox{25pt}{True value}}",
      "\\multirow{2}{*}{\\parbox{1pt}{Method}}",
      "\\multirow{2}{*}{\\parbox{42pt}{Avg\\\\estimate}}",
      "\\multirow{2}{*}{\\parbox{48pt}{Relative\\\\efficiency}}",
      "\\multirow{2}{*}{\\parbox{28pt}{rMSE}}",
      "\\multirow{2}{*}{\\parbox{40pt}{Coverage}}"
    )
    result_table <- rbind(rep("", ncol(result_table)), result_table)
    
    xtable_results <- xtable(
      result_table,
      label=sens_label,
      caption=paste(
        "Results from the simulation with",
        table_sample_size,
        "individuals in both the internal and external studies.
  For the ``Avg estimate'' and ``Coverage'' columns, the boldface indicates values within Monte Carlo error ($3\\sigma$) of the truth.
  For the ``Relative efficiency'' and ``rMSE'' columns, the boldface indicates the best performance for each coefficient (PET-WCLS in both cases)."
        )
      ) %>% print(
        sanitize.text.function = function(x) gsub("\\%", "\\\\\\%", x),
        include.rownames = FALSE,
        table.placement=NULL,
        floating.environment = "table*",
        print.results=TRUE) %>%
      str_replace("&  &  \\\\\\\\ \\n", "&  &  \\\\\\\\\n\\\\hline\n") %>% 
      str_replace("\\\\hline\\n &  &  ", " &  &  ") %>%
      str_replace("\\n  \\\\multirow\\{10\\}\\{\\*\\}\\{Slope\\}", "\n \\\\hline \n  \\\\multirow{10}{*}{Slope}")
    cat(xtable_results)
  }
  
  table_method_names <- c("WCLS-Internal", "WCLS-Pooled", "P-WCLS-Internal", "P-WCLS-Pooled", "ET-WCLS", "DR-WCLS", "PET-WCLS")
  make_table(400, method_vector=table_method_names)  # This is the one shown in the paper
  # make_table(400)  # All 10 methods
}

sensitivity_analysis(c(1, 2), "simulation-x2-0", x2_coef=0, x21sq_coef=0)
sensitivity_analysis(c(3, 0), "simulation-x2-2", x2_coef=2, x21sq_coef=0)
sensitivity_analysis(c(-2, 5), "simulation-x12sq", x2_coef=-3, x21sq_coef=0.3)