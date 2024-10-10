# Example analysis using simulated data for
# "Data integration methods for micro-randomized trials"

# Set random seed
random_seed <- 1
set.seed(random_seed)

# Import helper functions & helpful packages
source("./generate_data.R")
source("./extras/helpers.R")
source("./methods/wcls.R")
source("./methods/pwcls.R")
source("./methods/etwcls.R")
source("./methods/drwcls.R")
source("./methods/petwcls.R")
require(splines)

# Copy the setup from the simulation and simulate data
beta_r_true <- c(-2, 5)
coef_names <- c("Intercept", "Slope")
names(beta_r_true) <- coef_names
num_coefs <- length(coef_names)
method_names <- c("WCLS-Internal", "WCLS-Pooled", "P-WCLS-Internal", "P-WCLS-Pooled", "P-WCLS-Pooled-Obs", "ET-WCLS-Equal", "ET-WCLS-Kron", "ET-WCLS", "DR-WCLS", "PET-WCLS")
num_methods <- length(method_names)
beta_h_formula <- y ~ x1 + x2 + x3
beta_s_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1) + I(a_centered * x2)
beta_r_formula <- y ~ 0 + I(a_centered) + I(a_centered * x1)
et_beta_h_formula <- y ~ 0 + I(as.numeric(is_internal)) + I(is_internal*x1) + I(is_internal*x2) + I(is_internal*x3) + I(as.numeric(is_external)) + I(is_external*x1) + I(is_external*x2) + I(is_external*x3)
et_beta_r_formula <- y ~ 0 + I(is_internal * a_centered) + I(is_internal * a_centered * x1) + I(is_external * a_centered) + I(is_external * a_centered * x1)
pwcls_r_formula <- wcls_s_causal_effects ~ x1
a_intercept_formula <- a ~ 1
p_h_formula <- a ~ 1 + as.numeric(is_internal) + x1 + x2 + x3
example_dat <- generate_data(t_max=20, dof=10, n_internal=400, n_external=400, ar_param=0.5)
write.csv(example_dat, "example_dat.csv")
example_dat_internal <- example_dat[example_dat$is_internal, ]

# Apply methods
all_models <- list()
all_models[[1]] <- wcls(example_dat_internal, beta_r_true, beta_h_formula, beta_r_formula, p_r_formula=a_intercept_formula)
all_models[[2]] <- wcls(example_dat, beta_r_true, beta_h_formula, beta_r_formula, p_r_formula=a_intercept_formula)
all_models[[3]] <- pwcls(example_dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula, internal_only=TRUE)
all_models[[4]] <- pwcls(example_dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula)
all_models[[5]] <- pwcls(example_dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, p_s_formula=a_intercept_formula, p_h_formula=p_h_formula)
all_models[[6]] <- etwcls(example_dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="equal")
all_models[[7]] <- etwcls(example_dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="kronecker")
all_models[[8]] <- etwcls(example_dat, beta_r_true, et_beta_h_formula, et_beta_r_formula, a_intercept_formula, pooling_method="full")
all_models[[9]] <- drwcls(example_dat, beta_r_true, beta_h_formula, beta_s_formula, pwcls_r_formula, a_intercept_formula)
all_models[[10]] <- petwcls(example_dat, beta_r_true, beta_h_formula, beta_s_formula, et_beta_r_formula, pwcls_r_formula, a_intercept_formula)

# Inspect results
all_coefs <- sapply(all_models, function(d) d$beta_r)
all_ses <- sapply(all_models, function(d) d$se_beta_r)
coefs_ses_flat <- paste0(round(all_coefs, 3), " (", round(all_ses, 3), ")")
coefs_ses_mat <- matrix(coefs_ses_flat, nrow=num_coefs, ncol=num_methods)
coefs_ses_df <- as.data.frame(coefs_ses_mat)
row.names(coefs_ses_df) <- coef_names
colnames(coefs_ses_df) <- method_names
write.csv(coefs_ses_df, "example_estimates.csv")
