library(dplyr)
library(tidyr)

# set working directory
setwd("~/Documents/03_Paper/03_in_preparation/2019_rare_events_ma/06_Data_Prep")

# define functions and variables
# assess convergence for lme4 objects
conv_glmer <- function(model_summary) {
  # returns TRUE when model converged
  conv <- (model_summary$optinfo$conv$opt == 0)
  return(conv)
}

# assess convergence for glmmTMB objects
conv_glmmTMB <- function(model_summary) {
  # returns TRUE when model converged
  conv <- !anyNA(model_summary$coefficients$cond[, "Std. Error"]) &
    !anyNA(model_summary$coefficients$zi[, "Std. Error"])
  return(conv)
}

# assess convergence for Cai model
conv_cai <- function(model) {
  # returns TRUE when model converged
  conv <- (model$convergence == 0)
  return(conv)
}

# mean squared error
compute_mse <- function(estimates, true_param) {
  n <- length(estimates[!is.na(estimates)])
  mse <- sum((estimates - true_param)^2, na.rm = TRUE) / n
  return(mse)
}

# root mean squared error
compute_rmse <- function(estimates, true_param) {
  n <- length(estimates[!is.na(estimates)])
  mse <- sum((estimates - true_param)^2, na.rm = TRUE) / n
  return(sqrt(mse))
}

# mean absolute error
compute_mae <- function(estimates, true_param) {
  n <- length(estimates[!is.na(estimates)])
  mae <- sum(abs(estimates - true_param), na.rm = TRUE) / n
  return(mae)
}

# max error
compute_me <- function(estimates, true_param) {
  if (any(!is.na(estimates))) {
    me <- max(abs(estimates - true_param), na.rm = TRUE)
  } else {
    me <- NA
  }
  return(me)
}

# bias
compute_bias <- function(estimates, true_param) {
  bias <- mean(estimates, na.rm = TRUE) - true_param
  return(bias)
}

# median bias
compute_median_bias <- function(estimates, true_param) {
  bias <- median(estimates - true_param, na.rm = TRUE)
  return(bias)
}

# coverage
compute_coverage <- function(estimates, ses, true_param) {
  n <- length(estimates[!is.na(estimates)])
  lower_CI <- estimates - 1.96 * ses
  upper_CI <- estimates + 1.96 * ses
  covered <- ((true_param > lower_CI) & (true_param < upper_CI))
  coverage <- sum(covered, na.rm = TRUE) / n
  return(coverage)
}

# coverage for log transformed CI
compute_log_coverage <- function(estimates, ses, true_param) {
  n <- length(estimates[!is.na(estimates)])
  lower_CI_prep <- estimates - 1.96 * ses
  lower_CI_prep <- ifelse(
    lower_CI_prep < 0,
    0,
    lower_CI_prep
  )
  lower_CI <- log(lower_CI_prep)
  upper_CI <- log(estimates + 1.96 * ses)
  covered <- ((true_param > lower_CI) & (true_param < upper_CI))
  coverage <- sum(covered, na.rm = TRUE) / n
  return(coverage)
}

# compute Monte Carlo SE
compute_mc_se <- function(estimates) {
  n <- length(estimates[!is.na(estimates)])
  mc_se <- sqrt(var(estimates, na.rm = TRUE) / n)
  return(mc_se)
}

results_table <- expand.grid(
  true_RR = c(0.5, 1, 2),
  baseline_p = c(0.05, 0.1), 
  n_studies = c(5, 30, 100), 
  n_control = 50, 
  R = c(0.5, 1, 2),
  tau = c(0, 0.6, 1)
)

results_table$condition <- paste0(
  "RR", results_table$true_RR,
  "_p", results_table$baseline_p,
  "_nstudies", results_table$n_studies,
  "_ncontrol", results_table$n_control,
  "_R", results_table$R,
  "_tau", results_table$tau
)

n_conditions <- length(results_table$condition)

n_simulations <- 1000

model_names <- c(
  "m_REML", "m_REML_tcc", "m_DL", "m_DL_tcc", "m_SJ", "m_SJ_tcc",
  "m_poiss", "m_zip_rifs", "m_zip_fifs", "m_zip_ri", "m_zip_fi",
  "m_cond_binom",  "m_beta_binom", "m_kuss_binom", "m_cai_binom"
)

n_models <- length(model_names)

# load param_summary 
load("extracted_parameters.rds")

# make parameter table
param_table <- param_list[[1]]
param_table$condition <- rep(names(param_list)[1], length(param_table$sim_trial))
for (i in 2:length(param_list)) {
  param_list[[i]]$condition <- rep(names(param_list)[i], length(param_list[[i]]$sim_trial))
  param_table <- rbind(param_table, param_list[[i]])
}

# save 
save(param_table, file = "param_table.rds")

# make parameter summary
param_summary <- param_table %>%
  gather(c(log_RR_m_REML:se_log_br_m_poiss),
         key = "measure_model", value = "value") %>%
  extract(col = "measure_model", into = c("measure", "model"), 
          regex = "([[:graph:]]+)_m_([[:graph:]]+)")

# save
save(param_summary, file = "param_summary.rds")

# compute datasets for plots


#### model summaries ####
# Normal models
sum_normal <- param_summary %>%
  filter(grepl("REML", model) | grepl("DL", model) | grepl("SJ", model)) %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(log_RR)),
    mean_log_RR = mean(log_RR, na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log_RR),
    mean_RR = exp(mean_log_RR),
    median_log_RR = median(log_RR, na.rm = TRUE),
    median_RR = exp(median_log_RR),
    var_log_RR = var(log_RR, na.rm = TRUE),
    sd_log_RR = sd(log_RR, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log_RR, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log_RR, true_param = log(true_RR[1])),
    coverage_log_RR = compute_coverage(estimates = log_RR, 
                                       ses = se_log_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log_RR, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log_RR, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log_RR, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log_RR, true_param = log(true_RR[1])),
    mean_tau2 = mean(tau2, na.rm = TRUE),
    median_tau2 = median(tau2, na.rm = TRUE),
    var_tau2 = var(tau2, na.rm = TRUE),
    sd_tau2 = sd(tau2, na.rm = TRUE),
    bias_tau2 = compute_bias(estimates = tau2, true_param = true_tau[1]^2),
    median_bias_tau2 = compute_median_bias(estimates = tau2, true_param = true_tau[1]^2),
    mae_tau2 = compute_mae(estimates = tau2, true_param = true_tau[1]^2),
    me_tau2 = compute_me(estimates = tau2, true_param = true_tau[1]^2),
    rmse_tau2 = compute_rmse(estimates = tau2, true_param = true_tau[1]^2)
  )
save(sum_normal, file = "sum_normal.rds")

# Poisson model
sum_poiss <- param_summary %>%
  filter(model == "poiss") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(log_br)),
    mean_log_br = mean(log_br, na.rm = TRUE),
    mean_br = mean(exp(log_br), na.rm = TRUE),
    median_log_br = median(log_br, na.rm = TRUE),
    median_br = median(exp(log_br), na.rm = TRUE),
    var_log_br = var(log_br, na.rm = TRUE),
    sd_log_br = sd(log_br, na.rm = TRUE),
    bias_log_br = compute_bias(estimates = log_br, true_param = log(true_baseline_p[1])),
    median_bias_log_br = compute_median_bias(estimates = log_br, true_param = log(true_baseline_p[1])),
    coverage_log_br = compute_coverage(estimates = log_br, 
                                       ses = se_log_br, true_param = log(true_baseline_p[1])),
    mae_log_br = compute_mae(estimates = log_br, true_param = log(true_baseline_p[1])),
    me_log_br = compute_me(estimates = log_br, true_param = log(true_baseline_p[1])),
    mse_log_br = compute_mse(estimates = log_br, true_param = log(true_baseline_p[1])),
    rmse_log_br = compute_rmse(estimates = log_br, true_param = log(true_baseline_p[1])),
    mean_log_RR = mean(log_RR, na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log_RR),
    mean_RR = mean(exp(log_RR), na.rm = TRUE),
    median_log_RR = median(log_RR, na.rm = TRUE),
    median_RR = median(exp(log_RR), na.rm = TRUE),
    var_log_RR = var(log_RR, na.rm = TRUE),
    sd_log_RR = sd(log_RR, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log_RR, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log_RR, true_param = log(true_RR[1])),
    coverage_log_RR = compute_coverage(estimates = log_RR, 
                                       ses = se_log_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log_RR, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log_RR, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log_RR, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log_RR, true_param = log(true_RR[1])),
    mean_tau2 = mean(tau2, na.rm = TRUE),
    median_tau2 = median(tau2, na.rm = TRUE),
    var_tau2 = var(tau2, na.rm = TRUE),
    sd_tau2 = sd(tau2, na.rm = TRUE),
    bias_tau2 = compute_bias(estimates = tau2, true_param = true_tau[1]^2),
    median_bias_tau2 = compute_median_bias(estimates = tau2, true_param = true_tau[1]^2),
    mae_tau2 = compute_mae(estimates = tau2, true_param = true_tau[1]^2),
    me_tau2 = compute_me(estimates = tau2, true_param = true_tau[1]^2),
    rmse_tau2 = compute_rmse(estimates = tau2, true_param = true_tau[1]^2)
  )
save(sum_poiss, file = "sum_poiss.rds")

# ZIP model with slope in zero-inflation arm 
sum_zip_slope <- param_summary %>%
  filter(model == "zip_fifs" | model == "zip_rifs") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(mod_RR)),
    mean_log_RR = mean(log(mod_RR), na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log(mod_RR)),
    median_log_RR = median(log(mod_RR), na.rm = TRUE),
    var_log_RR = var(log(mod_RR), na.rm = TRUE),
    sd_log_RR = sd(log(mod_RR), na.rm = TRUE),
    bias_log_RR = compute_bias(log(mod_RR), true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(log(mod_RR), true_param = log(true_RR[1])),
    coverage_log_RR = compute_log_coverage(estimates = mod_RR, 
                                       ses = se_mod_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(log(mod_RR), true_param = log(true_RR[1])),
    me_log_RR = compute_me(log(mod_RR), true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(log(mod_RR), true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(log(mod_RR), true_param = log(true_RR[1])),
    mean_RR = mean(mod_RR, na.rm = TRUE),
    mc_SE_RR = compute_mc_se(mod_RR),
    median_RR = median(mod_RR, na.rm = TRUE),
    var_RR = var(mod_RR, na.rm = TRUE),
    sd_RR = sd(mod_RR, na.rm = TRUE),
    bias_RR = compute_bias(mod_RR, true_param = true_RR[1]),
    median_bias_RR = compute_median_bias(mod_RR, true_param = true_RR[1]),
    coverage_RR = compute_coverage(estimates = mod_RR, ses = se_mod_RR, 
                                   true_param = true_RR[1]),
    mae_RR = compute_mae(mod_RR, true_param = true_RR[1]),
    me_RR = compute_me(mod_RR, true_param = true_RR[1]),
    mse_RR = compute_mse(mod_RR, true_param = true_RR[1]),
    rmse_RR = compute_rmse(mod_RR, true_param = true_RR[1])
  )
save(sum_zip_slope, file = "sum_zip_slope.rds")

# ZIP model with only intercept in zero-inflation arm 
sum_zip <- param_summary %>%
  filter(model == "zip_fi" | model == "zip_ri") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(b_cond)),
    mean_log_RR = mean(b_cond, na.rm = TRUE),
    # for no slope in zi part, the (1-p)/(1-q) = 1 and modified relative risk
    # is just the b weight from the conditional model
    mc_SE_log_RR = compute_mc_se(b_cond),
    mean_RR = mean(exp(b_cond), na.rm = TRUE),
    median_log_RR = median(b_cond, na.rm = TRUE),
    median_RR = median(exp(b_cond), na.rm = TRUE),
    var_log_RR = var(b_cond, na.rm = TRUE),
    sd_log_RR = sd(b_cond, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = b_cond, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = b_cond, true_param = log(true_RR[1])),
    coverage_log_RR = compute_coverage(estimates = b_cond, 
                                       ses = se_b_cond, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = b_cond, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = b_cond, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = b_cond, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = b_cond, true_param = log(true_RR[1])),
    mean_tau2 = mean(var_b_cond, na.rm = TRUE),
    median_tau2 = median(var_b_cond, na.rm = TRUE),
    var_tau2 = var(var_b_cond, na.rm = TRUE),
    sd_tau2 = sd(var_b_cond, na.rm = TRUE),
    # see above (comment mean_log_RR) -> why we now have estimate for tau
    bias_tau2 = compute_bias(estimates = var_b_cond, true_param = true_tau[1]^2),
    median_bias_tau2 = compute_median_bias(estimates = var_b_cond, true_param = true_tau[1]^2),
    mae_tau2 = compute_mae(estimates = var_b_cond, true_param = true_tau[1]^2),
    me_tau2 = compute_me(estimates = var_b_cond, true_param = true_tau[1]^2),
    rmse_tau2 = compute_rmse(estimates = var_b_cond, true_param = true_tau[1]^2)
  )
save(sum_zip, file = "sum_zip.rds")

# Conditional binomial model
sum_binom <- param_summary %>%
  filter(model == "cond_binom") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(log_RR)),
    mean_log_RR = mean(log_RR, na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log_RR),
    mean_RR = mean(exp(log_RR), na.rm = TRUE),
    median_log_RR = median(log_RR, na.rm = TRUE),
    median_RR = median(exp(log_RR), na.rm = TRUE),
    var_log_RR = var(log_RR, na.rm = TRUE),
    sd_log_RR = sd(log_RR, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log_RR, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log_RR, true_param = log(true_RR[1])),
    coverage_log_RR = compute_coverage(estimates = log_RR, 
                                       ses = se_log_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log_RR, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log_RR, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log_RR, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log_RR, true_param = log(true_RR[1])),
    mean_tau2 = mean(tau2, na.rm = TRUE),
    median_tau2 = median(tau2, na.rm = TRUE),
    var_tau2 = var(tau2, na.rm = TRUE),
    sd_tau2 = sd(tau2, na.rm = TRUE),
    bias_tau2 = compute_bias(estimates = tau2, true_param = true_tau[1]^2),
    median_bias_tau2 = compute_median_bias(estimates = tau2, true_param = true_tau[1]^2),
    mae_tau2 = compute_mae(estimates = tau2, true_param = true_tau[1]^2),
    me_tau2 = compute_me(estimates = tau2, true_param = true_tau[1]^2),
    rmse_tau2 = compute_rmse(estimates = tau2, true_param = true_tau[1]^2)
  )
save(sum_binom, file = "sum_binom.rds")

# Beta binomial model
sum_betabinom <- param_summary %>%
  filter(model == "beta_binom") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(log_RR)),
    mean_log_RR = mean(log_RR, na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log_RR),
    mean_RR = mean(exp(log_RR), na.rm = TRUE),
    median_log_RR = median(log_RR, na.rm = TRUE),
    median_RR = median(exp(log_RR), na.rm = TRUE),
    var_log_RR = var(log_RR, na.rm = TRUE),
    sd_log_RR = sd(log_RR, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log_RR, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log_RR, true_param = log(true_RR[1])),
    coverage_log_RR = compute_coverage(estimates = log_RR, 
                                       ses = se_log_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log_RR, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log_RR, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log_RR, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log_RR, true_param = log(true_RR[1]))
  )
save(sum_betabinom, file = "sum_betabinom.rds")

# Kuss beta-binomial model
sum_kuss <- param_summary %>%
  filter(model == "kuss_binom") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(log_RR)),
    mean_log_RR = mean(log_RR, na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log_RR),
    mean_RR = mean(exp(log_RR), na.rm = TRUE),
    median_log_RR = median(log_RR, na.rm = TRUE),
    median_RR = median(exp(log_RR), na.rm = TRUE),
    coverage_log_RR = compute_coverage(estimates = log_RR, ses = se_log_RR,
                                       true_param = log(true_RR[1])),
    var_log_RR = var(log_RR, na.rm = TRUE),
    sd_log_RR = sd(log_RR, na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log_RR, true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log_RR, true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log_RR, true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log_RR, true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log_RR, true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log_RR, true_param = log(true_RR[1]))
  )
save(sum_kuss, file = "sum_kuss.rds")

# Cai conditional binomial model
sum_cai <- param_summary %>%
  filter(model == "cai_binom") %>%
  spread(key = "measure", value = "value") %>%
  group_by(condition, model) %>%
  summarize(
    converged = sum(!is.na(pooled_RR_cai)),
    mean_log_RR = mean(log(pooled_RR_cai), na.rm = TRUE),
    mc_SE_log_RR = compute_mc_se(log(pooled_RR_cai)),
    median_log_RR = median(log(pooled_RR_cai), na.rm = TRUE),
    coverage_log_RR = compute_log_coverage(estimates = pooled_RR_cai,
                                           ses = se_pooled_RR_cai, 
                                           true_param = log(true_RR[1])),
    var_log_RR = var(log(pooled_RR_cai), na.rm = TRUE),
    bias_log_RR = compute_bias(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    median_bias_log_RR = compute_median_bias(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    mae_log_RR = compute_mae(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    me_log_RR = compute_me(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    mse_log_RR = compute_mse(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    rmse_log_RR = compute_rmse(estimates = log(pooled_RR_cai), true_param = log(true_RR[1])),
    mean_RR = mean(pooled_RR_cai, na.rm = TRUE),
    mc_SE_RR = compute_mc_se(pooled_RR_cai),
    median_RR = median(pooled_RR_cai, na.rm = TRUE),
    var_RR = var(pooled_RR_cai, na.rm = TRUE),
    bias_RR = compute_bias(estimates = pooled_RR_cai, true_param = true_RR[1]),
    median_bias_RR = compute_median_bias(estimates = pooled_RR_cai, true_param = true_RR[1]),
    coverage_RR = compute_coverage(estimates = pooled_RR_cai, ses = se_pooled_RR_cai,
                                   true_param = true_RR[1]),
    mae_RR = compute_mae(estimates = pooled_RR_cai, true_param = true_RR[1]),
    me_RR = compute_me(estimates = pooled_RR_cai, true_param = true_RR[1]),
    mse_RR = compute_mse(estimates = pooled_RR_cai, true_param = true_RR[1]),
    rmse_RR = compute_rmse(estimates = pooled_RR_cai, true_param = true_RR[1])
  )
save(sum_cai, file = "sum_cai.rds")

#### data for violin plots ####
violin_plot_data <- param_summary %>%
  filter(measure == "log_RR" | measure == "pooled_RR_cai" | measure == "mod_RR" |
           (measure == "b_cond" & (model == "zip_ri" | model == "zip_fi"))) %>%
  mutate(value_RR_log_scale = ifelse((measure == "pooled_RR_cai" | measure == "mod_RR"), 
                                     log(value), value),
         factor_model = as.factor(model)) %>%
  select(-c(measure, value)) %>%
  filter(value_RR_log_scale != -Inf)

save(violin_plot_data, file = "violin_plot_data.rds")

#### data for convergence plots ####
conv_plot_data <- param_summary %>%
  filter(measure == "log_RR" | measure == "pooled_RR_cai" | measure == "mod_RR" |
           (measure == "b_cond" & (model == "zip_ri" | model == "zip_fi"))) %>%
  group_by(model, condition) %>%
  summarize(not_converged = sum(is.na(value))) %>%
  spread(key = "model", value = "not_converged") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau, group_ratio, 
             n_studies, n_control, condition),
    by = "condition"
  ) %>%
  gather(beta_binom:zip_rifs, key = "model", value = "not_converged")

save(conv_plot_data, file = "conv_plot_data.rds")

#### data for zero plots ####

zeros <- param_table %>%
  group_by(condition) %>%
  summarize(
    avg_n_zeros_treat = mean(n_zeros_treat),
    per_n_zeros_treat = mean(n_zeros_treat / n_studies[1]) ,
    sd_n_zeros_treat = sd(n_zeros_treat),
    avg_n_zeros_control = mean(n_zeros_control),
    per_n_zeros_control = mean(n_zeros_control / n_studies[1]),
    sd_n_zeros_control = sd(n_zeros_control),
    avg_n_zeros = mean(n_zeros_treat + n_zeros_control),
    sd_n_zeros = sd(n_zeros_treat + n_zeros_control)
  )

zero_plots <- param_table %>%
  filter(sim_trial == 1) %>%
  select(true_RR, true_baseline_p, true_tau,
         group_ratio, n_studies, n_control, condition) %>%
  mutate(
    true_RR_f = as.factor(true_RR), 
    true_tau_f = as.factor(true_tau)
  ) %>%
  inner_join(zeros, by = "condition")

# save
save(zero_plots, file = "zero_plots_data.rds")

#### data for coverage plots ####

coverage_plot <- sum_betabinom %>%
  select(condition, model, coverage_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_cai %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, coverage_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, coverage_log_RR)) %>%
  spread(key = "model", value = "coverage_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = c("condition" = "condition")
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "coverage_log_RR")

# save
save(coverage_plot, file = "coverage_plot_data.rds")

# for the models which don't give pooled RR estimates on log scale
coverage_plot_not_log <- sum_cai %>%
  select(condition, model, coverage_RR) %>%
  rbind(sum_zip_slope %>% select(condition, model, coverage_RR)) %>%
  spread(key = "model", value = "coverage_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = c("condition" = "condition")
  ) %>%
  gather(cai_binom:zip_rifs, key = "model", value = "coverage_RR")

# save
save(coverage_plot_not_log, file = "coverage_plot_not_log_data.rds")

#### data for bias plots ####

bias_plot <- sum_betabinom %>%
  select(condition, model, bias_log_RR, median_bias_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, bias_RR, median_bias_RR) %>%
          mutate(
            bias_log_RR = bias_RR,
            median_bias_log_RR = median_bias_RR
          ) %>%
          select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_cai %>% select(condition, model, bias_RR, median_bias_RR) %>%
          mutate(
            bias_log_RR = bias_RR,
            median_bias_log_RR = median_bias_RR
          ) %>%
          select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  # note that we don't compute bias for log RR but the RR for Cai and ZIP with slope model,
  # we just renamed it so that it works with this code
  select(condition, model, bias_log_RR) %>%
  spread(key = "model", value = "bias_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "bias_log_RR")

all_biases_plot <- sum_betabinom %>%
  select(condition, model, bias_log_RR, median_bias_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, bias_RR, median_bias_RR) %>%
          mutate(
            bias_log_RR = bias_RR,
            median_bias_log_RR = median_bias_RR
          ) %>%
          select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  rbind(sum_cai %>% select(condition, model, bias_RR, median_bias_RR) %>%
          mutate(
            bias_log_RR = bias_RR,
            median_bias_log_RR = median_bias_RR
          ) %>%
          select(condition, model, bias_log_RR, median_bias_log_RR)) %>%
  # note that we don't compute bias for log RR but the RR for Cai model
  select(condition, model, median_bias_log_RR) %>%
  spread(key = "model", value = "median_bias_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "median_bias_log_RR") %>%
  inner_join(bias_plot %>% select(condition, model, bias_log_RR),
             by = c("condition" = "condition", "model" = "model"))

save(all_biases_plot, file = "all_biases_plot.rds")

#### data for monte carlo se plots ####
mc_se_plot <- sum_betabinom %>%
  select(condition, model, mc_SE_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, mc_SE_RR) %>%
          mutate(mc_SE_log_RR = mc_SE_RR) %>%
          select(condition, model, mc_SE_log_RR)) %>%
  rbind(sum_cai %>% select(condition, model, mc_SE_RR) %>%
          mutate(mc_SE_log_RR = mc_SE_RR) %>%
          select(condition, model, mc_SE_log_RR)) %>%
  # note that we selected actually the MC for the RR not the log RR
  # we just renamed the variable to work better with the code
  spread(key = "model", value = "mc_SE_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "mc_SE_log_RR")

save(mc_se_plot, file = "mc_se_plot.rds")

#### data for tau2 bias plots ####

bias_plot_tau2 <- sum_binom %>%
  select(condition, model, bias_tau2, median_bias_tau2) %>%
  rbind(sum_normal %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  rbind(sum_poiss %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  rbind(sum_zip %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  select(condition, model, bias_tau2) %>%
  spread(key = "model", value = "bias_tau2") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = "condition"
  ) %>%
  gather(DL:zip_ri, key = "model", value = "bias_tau2")

all_biases_plot_tau2 <- sum_binom %>%
  select(condition, model, bias_tau2, median_bias_tau2) %>%
  rbind(sum_normal %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  rbind(sum_poiss %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  rbind(sum_zip %>% select(condition, model, bias_tau2, median_bias_tau2)) %>%
  select(condition, model, median_bias_tau2) %>%
  spread(key = "model", value = "median_bias_tau2") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR)),
    by = "condition"
  ) %>%
  gather(DL:zip_ri, key = "model", value = "median_bias_tau2") %>%
  inner_join(bias_plot_tau2 %>% select(condition, model, bias_tau2),
             by = c("condition" = "condition", "model" = "model"))

# save
save(all_biases_plot_tau2, file = "all_biases_plot_tau2_realdata.rds")

#### data for RMSE (log RR) plots ####
plot_rmse_data <- sum_betabinom %>%
  select(condition, model, rmse_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, rmse_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, rmse_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, rmse_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, rmse_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, rmse_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, rmse_log_RR)) %>%
  spread(key = "model", value = "rmse_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "rmse_log_RR")

save(plot_rmse_data, file = "plot_rmse_data.rds")

#### data for MAE (log RR) plots ####
plot_mae_data <- sum_betabinom %>%
  select(condition, model, mae_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, mae_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, mae_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, mae_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, mae_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, mae_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, mae_log_RR)) %>%
  spread(key = "model", value = "mae_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "mae_log_RR")

save(plot_mae_data, file = "plot_mae_data.rds")

#### data for ME (log RR) plots ####
plot_me_data <- sum_betabinom %>%
  select(condition, model, me_log_RR) %>%
  rbind(sum_binom %>% select(condition, model, me_log_RR)) %>%
  rbind(sum_normal %>% select(condition, model, me_log_RR)) %>%
  rbind(sum_poiss %>% select(condition, model, me_log_RR)) %>%
  rbind(sum_kuss %>% select(condition, model, me_log_RR)) %>%
  rbind(sum_zip %>% select(condition, model, me_log_RR)) %>%
  rbind(sum_zip_slope %>% select(condition, model, me_log_RR)) %>%
  spread(key = "model", value = "me_log_RR") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_rifs, key = "model", value = "me_log_RR")

save(plot_me_data, file = "plot_me_data.rds")

#### data for RMSE, MAE, ME for Cai model ####
plot_measures_data_cai <- sum_cai %>%
  select(condition, model, rmse_RR, mae_RR, me_RR) %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  ungroup() %>%
  gather(rmse_RR:me_RR, key = "measure", value = "value") %>%
  mutate(measure = recode(measure, rmse_RR = "RMSE", mae_RR = "MAE", me_RR = "ME"))

# save  
save(plot_measures_data_cai, file = "plot_measures_data_cai.rds")

### data for tau2 RMSE ####

plot_rmse_data_tau2 <- sum_binom %>%
  select(condition, model, rmse_tau2) %>%
  rbind(sum_normal %>% select(condition, model, rmse_tau2)) %>%
  rbind(sum_poiss %>% select(condition, model, rmse_tau2)) %>%
  rbind(sum_zip %>% select(condition, model, rmse_tau2)) %>%
  spread(key = "model", value = "rmse_tau2") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_ri, key = "model", value = "rmse_tau2")

# save
save(plot_rmse_data_tau2, file = "plot_rmse_data_tau2.rds")


### data for tau2 MAE ####

plot_mae_data_tau2 <- sum_binom %>%
  select(condition, model, mae_tau2) %>%
  rbind(sum_normal %>% select(condition, model, mae_tau2)) %>%
  rbind(sum_poiss %>% select(condition, model, mae_tau2)) %>%
  rbind(sum_zip %>% select(condition, model, mae_tau2)) %>%
  spread(key = "model", value = "mae_tau2") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_ri, key = "model", value = "mae_tau2")

# save
save(plot_mae_data_tau2, file = "plot_mae_data_tau2.rds")

#### data for tau2 ME ####

plot_me_data_tau2 <- sum_binom %>%
  select(condition, model, me_tau2) %>%
  rbind(sum_normal %>% select(condition, model, me_tau2)) %>%
  rbind(sum_poiss %>% select(condition, model, me_tau2)) %>%
  rbind(sum_zip %>% select(condition, model, me_tau2)) %>%
  spread(key = "model", value = "me_tau2") %>%
  inner_join(
    param_table %>%
      filter(sim_trial == 1) %>%
      select(true_RR, true_baseline_p, true_tau,
             group_ratio, n_studies, n_control, condition) %>%
      mutate(true_RR_f = as.factor(true_RR),
             true_tau_f = as.factor(true_tau)),
    by = "condition"
  ) %>%
  gather(DL:zip_ri, key = "model", value = "me_tau2")

# save
save(plot_me_data_tau2, file = "plot_me_data_tau2.rds")












