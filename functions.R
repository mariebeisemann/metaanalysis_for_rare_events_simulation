library(XML)
library(tidyverse)
library(lme4)
library(metafor)
library(glmmTMB)
library(extraDistr)

#### functions ####

extract_data <- function(data, no_outcome, subgroup_from, subgroup_to) {
  # extract datasets
  subgroups <- subgroup_from:subgroup_to
  datasets <- list()
  count <- 1
  for (group in subgroups) {
    studies <- which(names(data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]]) == "DICH_DATA")
    n_studies <- length(studies)
    dataset <- data.frame(
      study_id = vector(mode = "character", length = n_studies),
      events_1 = vector(mode = "numeric", length = n_studies),
      events_2 = vector(mode = "numeric", length = n_studies),
      total_1 = vector(mode = "numeric", length = n_studies),
      total_2 = vector(mode = "numeric", length = n_studies),
      stringsAsFactors = FALSE
    )
    if (length(dataset$study_id) > 0) {
      count2 <- 1
      for (study in studies) {
        if (is.character(data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]])) {
          dataset$study_id[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]["STUDY_ID"]
          dataset$events_1[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]["EVENTS_1"]
          dataset$events_2[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]["EVENTS_2"]
          dataset$total_1[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]["TOTAL_1"]
          dataset$total_2[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]["TOTAL_2"]
        }
        if (is.list(data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]])) {
          dataset$study_id[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]$.attrs["STUDY_ID"]
          dataset$events_1[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]$.attrs["EVENTS_1"]
          dataset$events_2[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]$.attrs["EVENTS_2"]
          dataset$total_1[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]$.attrs["TOTAL_1"]
          dataset$total_2[count2] <- 
            data$ANALYSES_AND_DATA$COMPARISON[[no_outcome]][[group]][[study]]$.attrs["TOTAL_2"]
        }
        count2 <- count2 + 1
      }
      datasets[[count]] <- dataset
      count <- count + 1
    } else {
      datasets[[count]] <- NA
      count <- count + 1
    }
  } 
  
  # rearrange data into one data frame
  result <- datasets[[1]]
  for (i in 2:length(datasets)) {
    if (!is.na(datasets[i])) {
      result <- rbind(result, datasets[[i]])
    }
  }
  
  # output
  return(as.data.frame(result))
}

# RR with treatment arm correction
treat_arm_correction_RR <- function(treat_event, treat_no_event,
                                    control_event, control_no_event, k) {
  treat_which <- which(treat_event == 0)
  control_which <- which(control_event == 0)
  if (!is.null(treat_which) | !is.null(control_which)) {
    to_correct <- unique(c(treat_which, control_which))
    for (i in to_correct) {
      corr_treat <- k / (treat_event[i] + treat_no_event[i])
      corr_control <- k / (control_event[i] + control_no_event[i])
      treat_event[i] <- treat_event[i] + corr_treat
      treat_no_event[i] <- treat_no_event[i] + corr_treat
      control_event[i] <- control_event[i] + corr_control
      control_no_event[i] <- control_no_event[i] + corr_control
    }
  }
  risk_treat <- treat_event / (treat_event + treat_no_event)
  risk_control <- control_event / (control_event + control_no_event)
  RR <- risk_treat / risk_control
  return(RR)
}

# var for log RR with treatment arm correction
treat_arm_correction_logRR_var <- function(treat_event, treat_no_event,
                                           control_event, control_no_event, k) {
  treat_which <- which(treat_event == 0)
  control_which <- which(control_event == 0)
  if (!is.null(treat_which) | !is.null(control_which)) {
    to_correct <- unique(c(treat_which, control_which))
    for (i in to_correct) {
      corr_treat <- k / (treat_event[i] + treat_no_event[i])
      corr_control <- k / (control_event[i] + control_no_event[i])
      treat_event[i] <- treat_event[i] + corr_treat
      treat_no_event[i] <- treat_no_event[i] + corr_treat
      control_event[i] <- control_event[i] + corr_control
      control_no_event[i] <- control_no_event[i] + corr_control
    }
  }
  log_RR_var <- (1 / treat_event) - (1 / (treat_event + treat_no_event)) +
    (1 / control_event) - (1 / (control_event + control_no_event))
  return(log_RR_var)
}

# likelihood function for conditional model proposed by Cai et al.
likelihood_cai <- function(theta, Y1, Y2, W) {
  # W is the sample size ratio: N_{i2}/N_{i1}
  gamma <- theta[1]
  psi <- theta[2]
  num <- beta((psi * gamma + Y1), (psi * W + Y2))
  denom <- beta((psi * gamma), (psi * W))
  lik <- num / denom
  return(-1*sum(log(lik)))
}

# likelihood for Kuss beta-binomial model
likelihood_kuss <- function(theta,
                            treat_event, treat_no_event,
                            control_event, control_no_event) {
  b0 <- theta[1]
  b1 <- theta[2]
  prec <- theta[3]
  treat_total <- treat_event + treat_no_event
  control_total <- control_event + control_no_event
  mu_treat <- exp(b0 + b1)
  mu_control <- exp(b0)
  
  # mean precision parameterization of beta-binomial distribution
  alpha_treat <- mu_treat * prec
  beta_treat <- (1 - mu_treat) * prec
  alpha_control <- mu_control * prec
  beta_control <- (1 - mu_control) * prec
  
  to_sum_treat <- dbbinom(
    treat_event, size = treat_total, alpha = alpha_treat,
    beta = beta_treat, log = TRUE
  )
  to_sum_control <- dbbinom(
    control_event, size = control_total, alpha = alpha_control,
    beta = beta_control, log = TRUE
  )
  
  sum_total <- sum(to_sum_treat) + sum(to_sum_control)
  
  return(-1 * sum_total)
}

# prep data for standard models

prep_data_standard <- function(data, 
                               treat_event, treat_no_event,
                               control_event, control_no_event, 
                               k) {
  data2 <- escalc("RR", ai = data[[treat_event]], bi = data[[treat_no_event]],
                  ci = data[[control_event]], di = data[[control_no_event]],
                  data = data, 
                  add = k, to = "only0", drop00 = TRUE) %>%
    filter(!is.na(yi))
  # double zero studies are dropped from data2 to be consistent with
  # binomial models
  
  # treatment arm continuity correction
  data2$RR_tcc <- treat_arm_correction_RR(
    data2[[treat_event]], data2[[treat_no_event]],
    data2[[control_event]], data2[[control_no_event]],
    k = k)
  data2$log_RR_tcc <- log(data2$RR_tcc)
  
  data2$log_RR_tcc_var <- treat_arm_correction_logRR_var(
    data2[[treat_event]], data2[[treat_no_event]],
    data2[[control_event]], data2[[control_no_event]],
    k = k)
  
  data2$group_ratio <- data2$total_1 / data2$total_2
  
  # output
  return(data2)
}

# describe data

describe_data <- function(data, data_long, data2, data2_long) {
  data_name <- deparse(substitute(data))
  n_studies <- length(data[["study_id"]])
  n_single_zeroes <- sum(((data[["events_1"]] == 0) | (data[["events_2"]] == 0))) -
    sum((data[["events_1"]] == 0) & (data[["events_2"]] == 0))
  n_double_zeroes <- sum((data[["events_1"]] == 0) & (data[["events_2"]] == 0))
  min_sample_size <- min(data[["total_1"]] + data[["total_2"]])
  max_sample_size <- max(data[["total_1"]] + data[["total_2"]])
  mean_sample_size <- mean(data[["total_1"]] + data[["total_2"]])
  sd_sample_size <- sd(data[["total_1"]] + data[["total_2"]])
  print(paste0("For dataset ", data_name, ":"))
  print(paste0("No. of studies: ", n_studies))
  print(paste0("Out of those, no. of single-zero studies: ", n_single_zeroes))
  print(paste0("Out of those, no. of double-zero studies: ", n_double_zeroes))  
  print(paste0("Average sample size was ", round(mean_sample_size, 2),
               " with SD ", round(sd_sample_size, 2)))  
  print( paste0("Minimal sample size is ", min_sample_size, 
                ", maximal sample size is ", max_sample_size))
}

# run analyses

run_analyses <- function(data, data_long, data2, data2_long) {
  model_list <- list()
  
  m_REML <- try(
    rma(yi = yi, vi = vi, data = data2, method = "REML")
  )
  m_REML_tcc <- try(
    rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "REML")
  )
  
  # estimation with DL
  m_DL <- try(
    rma(yi = yi, vi = vi, data = data2, method = "DL")
  )
  m_DL_tcc <- try(
    rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "DL")
  )
  
  # estimation with SJ
  m_SJ <- try(
    rma(yi = yi, vi = vi, data = data2, method = "SJ")
  )
  m_SJ_tcc <- try(
    rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "SJ")
  )
  
  # save models in list
  model_list[["m_REML"]] <- m_REML
  model_list[["m_REML_tcc"]] <- m_REML_tcc
  model_list[["m_DL"]] <- m_DL
  model_list[["m_DL_tcc"]] <- m_DL_tcc
  model_list[["m_SJ"]] <- m_SJ
  model_list[["m_SJ_tcc"]] <- m_SJ_tcc
  
  ### Poisson models ###
  # Poisson random-effects model
  m_poiss <- try(
    glmer(count ~ 1 + group + (1 + group | study_id) + offset(log(n)),
          data = data_long, family = "poisson")
  )
  
  # zero-inflated Poisson, random intercept, fixed slope in ZI
  m_zip_rifs <- try(
    glmmTMB(count ~ 1 + group + (1 + group | study_id) + offset(log(n)),
            data = data_long, family = poisson,
            ziformula = ~ 1 + group + (1 | study_id))
  )
  
  # zero-inflated Poisson, fixed intercept, fixed slope in ZI
  m_zip_fifs <- try(
    glmmTMB(count ~ 1 + group + (1 + group | study_id) + offset(log(n)),
            data = data_long, family = poisson,
            ziformula = ~ 1 + group)
  )
  
  # zero-inflated Poisson, random intercept in ZI
  m_zip_ri <- try(
    glmmTMB(count ~ 1 + group + (1 + group | study_id) + offset(log(n)),
            data = data_long, family = poisson,
            ziformula = ~ 1 + (1 | study_id))
  )
  
  # zero-inflated Poisson, fixed intercept in ZI
  m_zip_fi <- try(
    glmmTMB(count ~ 1 + group + (1 + group | study_id) + offset(log(n)),
            data = data_long, family = poisson,
            ziformula = ~ 1)
  )
  
  # save models to list
  model_list[["m_poiss"]] <- m_poiss
  model_list[["m_zip_rifs"]] <- m_zip_rifs
  model_list[["m_zip_fifs"]] <- m_zip_fifs
  model_list[["m_zip_ri"]] <- m_zip_ri
  model_list[["m_zip_fi"]] <- m_zip_fi
  
  ### Binomial models ###
  # conditional binomial model (Boehning et al., 2015; Cai et al., 2010; Stijnen et al., 2010)
  m_cond_binom <- try(
    glmer(cbind(events_1, events_2) ~ 1 + (1 | study_id) + offset(log(group_ratio)),
          data = data2, family = "binomial")
  )
  
  # beta binomial (with logit link)
  m_beta_binom <- try(
    glmmTMB(cbind(events_1, events_2) ~ 1 + offset(log(group_ratio)),
            data = data2, family = betabinomial(link = "logit"))
  )
  
  # actual beta binomial model in Kuss, 2015 (with log link)
  # can include double-zero studies
  theta_init_beta <- c(-1, 0, 1)
  names(theta_init_beta) <- c("b0", "b1", "prec")
  m_kuss_binom <- try(
    optim(
      par = theta_init_beta, 
      fn = likelihood_kuss, 
      treat_event = data$events_1,
      treat_no_event = data$no_events_1,
      control_event = data$events_2,
      control_no_event = data$no_events_2,
      hessian = TRUE
    )
  )
  
  # actual model by Cai et al., 2010
  theta_init <- c(0.1, 0.1)
  names(theta_init) <- c("gamma", "psi")
  m_cai_binom <- try(
    optim(
      par = theta_init,
      fn = likelihood_cai,
      Y1 = data2$events_1,
      Y2 = data2$events_2,
      W = (1 / data2$group_ratio),
      hessian = TRUE
    )
  )
  # pooled RR can then be computed using respective function based on fit
  # object and sample size ratio which is included in design table (R)
  
  
  # save models to list
  model_list[["m_cond_binom"]] <- m_cond_binom
  model_list[["m_beta_binom"]] <- m_beta_binom
  model_list[["m_kuss_binom"]] <- m_kuss_binom
  model_list[["m_cai_binom"]] <- m_cai_binom
  
  # return model list
  return(model_list)
}

# show RR for standard models (with CI)
RR_for_stand_models <- function(results, model_name) {
  if (class(results[[model_name]])[1] == "rma.uni") {
      out <- c(
        exp(summary(results[[model_name]])$b),
        exp(summary(results[[model_name]])$ci.lb),
        exp(summary(results[[model_name]])$ci.ub)
      )
      names(out) <- c("RR", "CI_lower", "CI_upper")
  } else {
    out <- "There was an error in computing the model. No RR available."
  }
  return(out)
}

RR_for_poiss <- function(results) {
  if (class(results[["m_poiss"]])[1] == "glmerMod") {
    out <- c(
      exp(summary(results[["m_poiss"]])$coefficients["groupevents_1", "Estimate"]),
      exp(
        summary(results[["m_poiss"]])$coefficients["groupevents_1", "Estimate"] -
          1.96 * summary(results[["m_poiss"]])$coefficients["groupevents_1", "Std. Error"]
      ),
      exp(
        summary(results[["m_poiss"]])$coefficients["groupevents_1", "Estimate"] +
          1.96 * summary(results[["m_poiss"]])$coefficients["groupevents_1", "Std. Error"]
      )
    )
    names(out) <- c("RR", "CI_lower", "CI_upper")
  } else {
    out <- "There was an error in computing the model. No RR available."
  }
  return(out)
}

RR_for_cond_binom <- function(results) {
  if (class(results[["m_cond_binom"]])[1] == "glmerMod") {
    out <- c(
      exp(summary(results[["m_cond_binom"]])$coefficients["(Intercept)", "Estimate"]),
      exp(
        summary(results[["m_cond_binom"]])$coefficients["(Intercept)", "Estimate"] -
          1.96 * summary(results[["m_cond_binom"]])$coefficients["(Intercept)", "Std. Error"]
      ),
      exp(
        summary(results[["m_cond_binom"]])$coefficients["(Intercept)", "Estimate"] +
          1.96 * summary(results[["m_cond_binom"]])$coefficients["(Intercept)", "Std. Error"]
      )
    )
    names(out) <- c("RR", "CI_lower", "CI_upper")
  } else {
    out <- "There was an error in computing the model. No RR available."
  }
  return(out)
}

RR_for_beta_binom <- function(results) {
  if (class(results[["m_beta_binom"]])[1] == "glmmTMB") {
    out <- c(
      exp(summary(results[["m_beta_binom"]])$coefficients$cond["(Intercept)", "Estimate"]),
      exp(
        summary(results[["m_beta_binom"]])$coefficients$cond["(Intercept)", "Estimate"] -
          1.96 * summary(results[["m_beta_binom"]])$coefficients$cond["(Intercept)", "Std. Error"]
      ),
      exp(
        summary(results[["m_beta_binom"]])$coefficients$cond["(Intercept)", "Estimate"] +
          1.96 * summary(results[["m_beta_binom"]])$coefficients$cond["(Intercept)", "Std. Error"]
      )
    )
    names(out) <- c("RR", "CI_lower", "CI_upper")
  } else {
    out <- "There was an error in computing the model. No RR available."
  }
  return(out)
}

RR_for_kuss_binom <- function(results) {
  if (class(results[["m_kuss_binom"]]) == "list") {
    if (names(results[["m_kuss_binom"]][1]) == "par") {
      se <- try(sqrt(solve(results[["m_kuss_binom"]]$hessian)["b1", "b1"]), silent = TRUE)
      if (class(se) == "numeric") {
        out <- c(
          exp(results[["m_kuss_binom"]]$par["b1"]),
          exp(results[["m_kuss_binom"]]$par["b1"] - 1.96 * se),
          exp(results[["m_kuss_binom"]]$par["b1"] + 1.96 * se)
        )
        names(out) <- c("RR", "CI_lower", "CI_upper")
      } else {
        out <- c(
          exp(results[["m_kuss_binom"]]$par["b1"]),
          NA,
          NA
        )
        names(out) <- c("RR", "CI_lower", "CI_upper")
        out <- list(out, "Hessian matrix included only zeroes.")
      }
    }
  } else {
    out <- "There was an error in computing the model. No RR available."
  }
  return(out)
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
                gamma = fit[["par"]]["gamma"], psi = fit[["par"]]["psi"], Wi = W[i],
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
  gamma <- fit[["par"]]["gamma"]
  psi <- fit[["par"]]["psi"]
  cov_matrix <- try(solve(fit$hessian), silent = TRUE)
  
  if (class(cov_matrix) == "matrix") {
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
    out <- as.numeric(se) 
  } else {
    out <- "Hessian matrix contained only zeroes and coould not be reversed."
  }
  return(out)
}

RR_for_cai_binom <- function(results, data2) {
  if (length(data2[[1]]) == 0) {
    out <- "All studies were double-zero studies and excluded for the computation
            of this model. Thus, this model can't provide a sensible estimate for the RR."
  } else {
    if (class(results[["m_cai_binom"]]) == "list") {
      if (names(results[["m_cai_binom"]][1]) == "par") {
        RR <- pooled_rr_cai(fit = results[["m_cai_binom"]], W = (1 / data2[["group_ratio"]]))
        se <- se_pooled_rr_cai(fit = results[["m_cai_binom"]], W = (1 / data2[["group_ratio"]]))
        if (class(se) == "numeric") {
          out <- c(RR, RR - 1.96 * se, RR + 1.96 * se)
          names(out) <- c("RR", "CI_lower", "CI_upper")
        } else {
          out <- c(RR, NA, NA)
          names(out) <- c("RR", "CI_lower", "CI_upper")
          out <- list(out, se)
        }
      }  
    } else {
      out <- "There was an error in computing the model. No RR available."
    }
  }
  return(out)
}

# show results for all models per data frame

show_results <- function(results, data2) {
  text <- "Results for "
  
  print(paste0(text, "m_REML:"))
  print(summary(results[["m_REML"]]))
  print("RR for m_REML:")
  print(RR_for_stand_models(results, model_name = "m_REML"))
  
  print(paste0(text, "m_REML_tcc:"))
  print(summary(results[["m_REML_tcc"]]))
  print("RR for m_REML_tcc:")
  print(RR_for_stand_models(results, model_name = "m_REML_tcc"))
  
  print(paste0(text, "m_DL:"))
  print(summary(results[["m_DL"]]))
  print("RR for m_DL:")
  print(RR_for_stand_models(results, model_name = "m_DL"))
  
  print(paste0(text, "m_DL_tcc:"))
  print(summary(results[["m_DL_tcc"]]))
  print("RR for m_DL_tcc:")
  print(RR_for_stand_models(results, model_name = "m_DL_tcc"))
  
  print(paste0(text, "m_SJ:"))
  print(summary(results[["m_SJ"]]))
  print("RR for m_SJ:")
  print(RR_for_stand_models(results, model_name = "m_SJ"))
  
  print(paste0(text, "m_SJ_tcc:"))
  print(summary(results[["m_SJ_tcc"]]))
  print("RR for m_SJ_tcc:")
  print(RR_for_stand_models(results, model_name = "m_SJ_tcc"))
  
  print(paste0(text, "m_poiss:"))
  print(summary(results[["m_poiss"]]))
  # warning: is singular
  print("RR for m_poiss:")
  print(RR_for_poiss(results))
  
  print(paste0(text, "m_zip_rifs:"))
  print(summary(results[["m_zip_rifs"]]))
  if (is.na(summary(results[["m_zip_rifs"]])$coefficients$cond["(Intercept)", "Std. Error"]) |
      is.na(summary(results[["m_zip_rifs"]])$coefficients$zi["(Intercept)", "Std. Error"])) {
    print("m_zip_rifs did not converge properly (no standard errors were computed.")
  }
  # not converged (no standard errors)
  
  print(paste0(text, "m_zip_fifs:"))
  print(summary(results[["m_zip_fifs"]]))
  if (is.na(summary(results[["m_zip_fifs"]])$coefficients$cond["(Intercept)", "Std. Error"]) |
      is.na(summary(results[["m_zip_fifs"]])$coefficients$zi["(Intercept)", "Std. Error"])) {
    print("m_zip_fifs did not converge properly (no standard errors were computed.")
  }
  # not converged (no standard errors)
  
  print(paste0(text, "m_zip_ri:"))
  print(summary(results[["m_zip_ri"]]))
  if (is.na(summary(results[["m_zip_ri"]])$coefficients$cond["(Intercept)", "Std. Error"]) |
      is.na(summary(results[["m_zip_ri"]])$coefficients$zi["(Intercept)", "Std. Error"])) {
    print("m_zip_ri did not converge properly (no standard errors were computed.")
  }
  # not converged (no standard errors)
  
  print(paste0(text, "m_zip_fi:"))
  print(summary(results[["m_zip_fi"]]))
  if (is.na(summary(results[["m_zip_fi"]])$coefficients$cond["(Intercept)", "Std. Error"]) |
      is.na(summary(results[["m_zip_fi"]])$coefficients$zi["(Intercept)", "Std. Error"])) {
    print("m_zip_fi did not converge properly (no standard errors were computed.")
  }
  # not converged (no standard errors)
  
  print(paste0(text, "m_cond_binom:"))
  print(summary(results[["m_cond_binom"]]))
  # warning: is singular
  print("RR for m_cond_binom:")
  print(RR_for_cond_binom(results))
  
  print(paste0(text, "m_beta_binom:"))
  print(summary(results[["m_beta_binom"]]))
  print("RR for m_beta_binom:")
  print(RR_for_beta_binom(results))
  
  print(paste0(text, "m_kuss_binom:"))
  print(results[["m_kuss_binom"]])
  if (class(results[["m_kuss_binom"]]) == "list") {
    if (names(results[["m_kuss_binom"]][1]) == "par") {
      if (results[["m_kuss_binom"]]$convergence == 0) {
        print("RR for m_kuss_binom:")
        print(RR_for_kuss_binom(results))
      } else {
        print("The model m_kuss_binom has not converged.")
      }
    }
  } else {
    print("Error in computing m_kuss_binom.")
  }
  
  print(paste0(text, "m_cai_binom:"))
  print(results[["m_cai_binom"]])
  if (class(results[["m_cai_binom"]]) == "try-error") {
    print("Trying to compute m_cai_binom failed.")
  }
  print("RR for cai_binom:")
  print(RR_for_cai_binom(results, data2))
}



