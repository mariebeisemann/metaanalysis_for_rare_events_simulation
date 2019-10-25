# this script was run on the computing cluster PALMA II, University of Muenster


# packages
library(dplyr)
library(tidyr)
library(glmmTMB)

# set working directory
#setwd("/scratch/tmp/m_beis01")

# functions
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

# compute Monte Carlo SE
compute_mc_se <- function(estimates) {
  n <- length(estimates[!is.na(estimates)])
  mc_se <- sqrt(var(estimates, na.rm = TRUE) / n)
  return(mc_se)
}

# compute modified relative risk
compute_mod_RR <- function(int_zi, b_zi, b_pois) {
  p <- exp(int_zi + b_zi) / (1 + exp(int_zi + b_zi))
  q <- exp(int_zi) / (1 + exp(int_zi))
  MRR <- ((1 - p)/(1 - q))*exp(b_pois)
  return(MRR)
}

# compute se for modified rr
compute_se_mod_rr <- function(int_zi, b_zi, b_pois, cov_matrix) {
  jac_int_zi <- ((1 - exp(b_zi)) * exp(b_pois + int_zi)) / (1 + exp(int_zi + b_zi))^2
  jac_b_zi <- (-1) * ((1 + exp(int_zi)) * exp(int_zi + b_zi + b_pois)) /
    (1 + exp(int_zi + b_zi))^2
  jac_b_pois <- (exp(b_pois) + exp(int_zi + b_pois)) / (1 + exp(int_zi + b_zi))
  jac <- matrix(c(jac_int_zi, jac_b_zi, jac_b_pois), nrow = 1, ncol = 3, byrow = TRUE)
  variance <- jac %*% cov_matrix %*% t(jac)
  se <- sqrt(variance)
  return(as.numeric(se))
}

# compute pooled RR for Cai model
pooled_rr_cai <- function(fit, W) {
  # define function for interval
  integral_function <- function(tau, gamma, psi, Wi) {
    out <- exp(tau) * (((exp(tau) / (Wi + exp(tau)))^(psi*gamma) * 
                          (Wi / (Wi + exp(tau)))^(psi*Wi)) / 
                         (beta((psi * gamma), (psi * Wi))))
    return(out)
  }
  integral_sum <- 0
  error <- FALSE
  for (i in seq_along(W)) {
    int <- try(
      integrate(integral_function, -100, 100,
                gamma = fit$par["gamma"], psi = fit$par["psi"], Wi = W[i],
                subdivisions=2000),
      silent = TRUE
    )
    if (class(int) != "try-error") {
      integral_sum <- integral_sum + int$value
    } else {
      error <- TRUE
    }
  }
  if (error | integral_sum == Inf) {
    RR <- NA
  } else {
    RR <- integral_sum / length(W)
  }
  return(RR)
}

# compute SE for pooled RR for Cai model
se_pooled_rr_cai <- function(fit, W) {
  gamma <- fit$par["gamma"]
  psi <- fit$par["psi"]
  cov_matrix <- solve(fit$hessian)
  
  # first partial derivative after psi
  jac_psi_int_function <- function(theta, gamma, psi, W) {
    out <- (exp(theta) * (W/(W + exp(theta)))^(W * psi) * 
              (exp(theta)/(W + exp(theta)))^(psi * gamma) *
              (gamma * log(exp(theta)/(W + exp(theta))) + W * log(W/(W + exp(theta))) + 
                 W * digamma(psi * (W + gamma)) + gamma * digamma(psi * (W + gamma)) -
                 W * digamma(W * psi) - gamma * digamma(psi * gamma)))/beta(psi * gamma, W * psi)
    return(out)
  }
  integral_sum_psi <- 0
  error <- FALSE
  for (i in seq_along(W)) {
    int_psi <- try(
      integrate(jac_psi_int_function, -100, 100,
                gamma = gamma, psi = psi, W = W[i],
                subdivisions=3000),
      silent = TRUE
    )
    if (class(int_psi) != "try-error") {
      integral_sum_psi <- integral_sum_psi + int_psi$value
    } else {
      error <- TRUE
    }
  }
  if (error | integral_sum_psi == Inf) {
    jac_psi <- NA
  } else {
    jac_psi <- integral_sum_psi / length(W)
  }
  
  # first partial derivative after gamma
  jac_gamma_int_function <- function(theta, gamma, psi, W) {
    out <- (exp(theta) * psi * (W/(W + exp(theta)))^(W * psi) *
              (exp(theta)/(W + exp(theta)))^(gamma * psi) *
              (log(exp(theta)/(W + exp(theta))) * beta(gamma * psi, W * psi) -
                 beta(gamma * psi, W * psi) * (digamma(gamma * psi) -
                                                 digamma(gamma * psi + W * psi)))) /
      beta(gamma * psi, W * psi)^2
     return(out)
  }
  integral_sum_gamma <- 0
  error <- FALSE
  for (i in seq_along(W)) {
    int_gamma <- try(
      integrate(jac_gamma_int_function, -100, 100,
                gamma = gamma, psi = psi, W = W[i],
                subdivisions=3000),
      silent = TRUE
    )
    if (class(int_gamma) != "try-error") {
      integral_sum_gamma <- integral_sum_gamma + int_gamma$value
    } else {
      error <- TRUE
    }
  }
  if (error | integral_sum_gamma == Inf) {
    jac_gamma <- NA
  } else {
    jac_gamma <- integral_sum_gamma / length(W)
  }
  
  jac <- matrix(c(jac_gamma, jac_psi), nrow = 1, ncol = 2, byrow = TRUE)
  variance <- jac %*% cov_matrix %*% t(jac) 
  se <- sqrt(variance)
  return(as.numeric(se))  
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

### prep data ###

param_list <- list()

for (i in 1:n_conditions) {
  # read in current sim results
  file_name <- paste0("sim_res_", i, ".rds")
  results <- readRDS(file_name)
  
  # define constants and declare variables
  condition <- results_table$condition[i]
  true_RR <- results_table$true_RR[i]
  true_baseline_p <- results_table$baseline_p[i]
  true_tau <- results_table$tau[i]
  R <- results_table$R[i]
  n_studies <- results_table$n_studies[i]
  n_control <- results_table$n_control[i]
  params <- data.frame(
    sim_trial = 1:n_simulations,
    true_RR = rep(true_RR, n_simulations),
    true_baseline_p = rep(true_baseline_p, n_simulations),
    true_tau = rep(true_tau, n_simulations),
    group_ratio = rep(R, n_simulations),
    n_studies = rep(n_studies, n_simulations),
    n_control = rep(n_control, n_simulations),
    n_zeros_treat = rep(NA, n_simulations),
    n_zeros_control = rep(NA, n_simulations),
    stringsAsFactors = FALSE
  )
  
  # add columns for parameters & different models to dataframe
  for (m in 1:n_models) {
    # columns for metafor models
    if (grepl("REML", model_names[m]) | grepl("DL", model_names[m]) |
        grepl("SJ", model_names[m])) {
      col_name <- paste0("log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("tau2_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_tau2_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if metafor model
    # columns for Poisson model
    if (model_names[m] == "m_poiss") {
      col_name <- paste0("log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("tau2_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if Poisson model
    if (model_names[m] == "m_cond_binom") {
      col_name <- paste0("log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("tau2_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if conditional binomial model
    if (model_names[m] == "m_zip_rifs" | model_names[m] == "m_zip_fifs") {
      col_name <- paste0("int_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_int_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("b_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_b_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("int_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_int_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("var_b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("mod_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_mod_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if slope zip model
    if (model_names[m] == "m_zip_ri" | model_names[m] == "m_zip_fi") {
      col_name <- paste0("int_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_int_zi_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("int_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_int_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("var_b_cond_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if intercept only zip model
    if (model_names[m] == "m_cai_binom") {
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("pooled_RR_cai_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_pooled_RR_cai_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if Cai model
    if (model_names[m] == "m_beta_binom") {
      col_name <- paste0("log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if beta-binomial model
    if (model_names[m] == "m_kuss_binom") {
      col_name <- paste0("log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
      col_name <- paste0("se_log_RR_", model_names[m])
      params[[col_name]] <- rep(NA, times = n_simulations)
    } # end if Kuss model
  }
  
  # loop through simulation trials and extract parameters
  for (j in 1:length(results)) {
    # extract how many zeros there are per trial
    data_set <- results[[j]]$data
    params[j, c("n_zeros_treat")] <- sum(data_set$treat_event == 0)
    params[j, c("n_zeros_control")] <- sum(data_set$control_event == 0)
    
    # loop through models
    for (k in 1:n_models) {
      model_name <- names(results[[j]]$model_list)[k]
      current_model <- results[[j]]$model_list[[model_name]]$model
      
      ### conventional normal-normal models ###
      if (class(current_model)[1] == "rma.uni") {
        params[j, paste0("log_RR_", model_name)] <- current_model$b
        params[j, paste0("se_log_RR_", model_name)] <- current_model$se
        params[j, paste0("tau2_", model_name)] <- current_model$tau2
        params[j, paste0("se_tau2_", model_name)] <- current_model$se.tau2
      }
      
      ### other models ###
      if (class(current_model)[1] == "list") {
        # when list, then either glmer or glmmTMB or Cai model
        
        ### GLmer Models (Conditional Binomial and Poisson) ###
        if (length(current_model) == 17) {
          # model is glmer model
          # check if model has converged and only extract parameters if 
          # that is the case
          if (conv_glmer(current_model)) {
            # model has converged, parameters should be extracted
            if (current_model$family == "poisson") {
              params[j, paste0("log_RR_", model_name)] <- 
                current_model$coefficients["grouptreat", "Estimate"]
              params[j, paste0("se_log_RR_", model_name)] <- 
                current_model$coefficients["grouptreat", "Std. Error"]
              params[j, paste0("tau2_", model_name)] <- 
                current_model$varcor$study["grouptreat", "grouptreat"]
              params[j, paste0("log_br_", model_name)] <- 
                current_model$coefficients["(Intercept)", "Estimate"]
              params[j, paste0("se_log_br_", model_name)] <- 
                current_model$coefficients["(Intercept)", "Std. Error"]
            } 
            if (current_model$family == "binomial") {
              params[j, paste0("log_RR_", model_name)] <- 
                current_model$coefficients["(Intercept)", "Estimate"]
              params[j, paste0("se_log_RR_", model_name)] <- 
                current_model$coefficients["(Intercept)", "Std. Error"]
              params[j, paste0("tau2_", model_name)] <- 
                current_model$varcor$study["(Intercept)", "(Intercept)"]
            }
          } 
        } # end if glmer
        
        ### glmmTMB models (ZIP models and beta binomial) ###
        if (length(current_model) == 10) {
          # model is glmmTMB
          if (conv_glmmTMB(current_model)) {
            # check whether ZIP or beta binomial model 
            if (current_model$family == "poisson") {
              # for all ZIP models, we want to extract the intercept and b-weight
              # from the zi arm of the model and the b-weight from the conditional
              # model to be able to compute the modified RR
              params[j, paste0("int_zi_", model_name)] <- 
                current_model$coefficients$zi["(Intercept)", "Estimate"]
              params[j, paste0("se_int_zi_", model_name)] <- 
                current_model$coefficients$zi["(Intercept)", "Std. Error"]
              if (model_name == "m_zip_rifs" | model_name == "m_zip_fifs") {
                # check if we have a slope in the zero-inflation arm,
                # extract that and compute the modified RR
                params[j, paste0("b_zi_", model_name)] <- 
                  current_model$coefficients$zi["grouptreat", "Estimate"]
                params[j, paste0("se_b_zi_", model_name)] <- 
                  current_model$coefficients$zi["grouptreat", "Std. Error"]
                params[j, paste0("mod_RR_", model_name)] <- 
                  compute_mod_RR(int_zi = current_model$coefficients$zi["(Intercept)", "Estimate"],
                                 b_zi = current_model$coefficients$zi["grouptreat", "Estimate"], 
                                 b_pois = current_model$coefficients$cond["grouptreat", "Estimate"])
                cov_matrix <- results[[j]]$model_list[[model_name]]$cov_matrix
                params[j, paste0("se_mod_RR_", model_name)] <- 
                  compute_se_mod_rr(int_zi = current_model$coefficients$zi["(Intercept)", "Estimate"],
                                    b_zi = current_model$coefficients$zi["grouptreat", "Estimate"],
                                    b_pois = current_model$coefficients$cond["grouptreat", "Estimate"],
                                    cov_matrix = cov_matrix)
                
              }
              params[j, paste0("int_cond_", model_name)] <- 
                current_model$coefficients$cond["(Intercept)", "Estimate"]
              params[j, paste0("se_int_cond_", model_name)] <- 
                current_model$coefficients$cond["(Intercept)", "Std. Error"]
              params[j, paste0("b_cond_", model_name)] <- 
                current_model$coefficients$cond["grouptreat", "Estimate"]
              params[j, paste0("se_b_cond_", model_name)] <- 
                current_model$coefficients$cond["grouptreat", "Std. Error"]
              params[j, paste0("var_b_cond_", model_name)] <- 
                current_model$varcor$cond$study["grouptreat", "grouptreat"]
            }
            if (current_model$family == "betabinomial") {
              # beta binomial model
              params[j, paste0("log_RR_", model_name)] <- 
                current_model$coefficients$cond["(Intercept)", "Estimate"]
              params[j, paste0("se_log_RR_", model_name)] <- 
                current_model$coefficients$cond["(Intercept)", "Std. Error"]
            }
          } 
        } # end if glmmTMB
        
        ### Cai & Kuss model ###
        if (length(current_model) == 6) {
          # either Cai or Kuss model
          if (conv_cai(current_model)) {
            # check if model (either Cai or Kuss) has converged
            if (names(current_model$par)[1] == "gamma") {
              # Cai model
              # model has converged and we extract parameters
              params[j, paste0("pooled_RR_cai_", model_name)] <-
                pooled_rr_cai(fit = current_model, W = rep(1/R, n_studies))
              params[j, paste0("se_pooled_RR_cai_", model_name)] <-
                se_pooled_rr_cai(fit = current_model, W = rep(1/R, n_studies))
            }
            if (names(current_model$par)[1] == "b0") {
              # Kuss model
              params[j, paste0("log_RR_", model_name)] <- current_model$par["b1"]
              params[j, paste0("se_log_RR_", model_name)] <-
                sqrt(solve(current_model$hessian)["b1", "b1"])
            }
            
          } 
        } # end if Cai model
      } # end if class == list
    } # end loop through models for the jth simulation trial
  } # end loop through simulation trials for the ith condition
  param_list[[condition]] <- params
} # end loop through conditions

# save
save(param_list, file = "extracted_parameters.rds")
     
     
     
     