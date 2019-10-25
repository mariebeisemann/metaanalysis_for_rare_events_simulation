# this script was run on the computing cluster PALMA II, University of Muenster

# packages
library(doParallel)
library(dplyr)
library(tidyr)
library(lme4)
library(metafor)
library(glmmTMB)
library(extraDistr)
library(truncnorm)
library(nleqslv)

# set working directory
# set working directory
setwd("/scratch/tmp/m_beis01")

n_cluster <- 72
# n_cluster <- 1
cluster <- makeCluster(n_cluster)
registerDoParallel(cluster)

# define functions
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

# function to solve for mu and sigma for truncnorm
solve_for_mu_and_sigma <- function(par, tau, log_RR, upper_bound) {
  mu <- par[1]
  sigma <- par[2]
  beta <- (upper_bound - mu) / sigma
  sol_mu <- log_RR + sigma * (dnorm(beta) / pnorm(beta)) - mu
  sol_sigma <- sigma^2 * (1 - (beta * (dnorm(beta)/pnorm(beta))) - (dnorm(beta)/pnorm(beta))^2) - tau^2
  sol <- c(sol_mu, sol_sigma)
  return(sol)
}

# define conditions
design <- expand.grid(
  true_RR = c(0.5, 1, 2),
  baseline_p = c(0.05, 0.1),
  n_studies = c(5, 30, 100), 
  n_control = 50,
  R = c(0.5, 1, 2),
  tau = c(0, 0.6, 1) 
)

n_conditions <- length(design$true_RR)

design$condition <- paste0(
  "RR", design$true_RR,
  "_p", design$baseline_p,
  "_nstudies", design$n_studies,
  "_ncontrol", design$n_control,
  "_R", design$R,
  "_tau", design$tau
)

n_simulations <- 1000

results <- vector(mode = "list", length = n_conditions)
names(results) <- design$condition

# sim_data <- vector(mode = "list", length = n_conditions)
# names(sim_data) <- design$condition

# simulate
I <- 1:n_conditions
for (i in I) {
  # make list to save simulation trials and data in
  # results_simulation_trials <- vector(mode = "list", length = n_simulations)
  
  # print where we are
  message("We are in the ", i, "th condition of ", n_conditions, ".")
  
  # set seed
  set.seed(17)
  
  # save current state
  file_name <- paste0("Okt_2019_random_state_cond_", i, ".RData")
  save(".Random.seed", file = file_name)

  # run different simulation trials
  results_simulation_trials <- foreach (j = 1:n_simulations) %dopar% {
    require(dplyr)
    require(tidyr)
    require(lme4)
    require(metafor)
    require(glmmTMB)
    require(extraDistr)
    require(truncnorm)
    require(nleqslv)

    ### simulate data set for the trial ###
    
    # solve for parameters mu and sigma of truncnorm distribution so that
    # log RR and tau are what the design says they are
    init <- c(1, 1)
    names(init) <- c("mu", "sigma")
    sol <- nleqslv(
      x = init,
      fn = solve_for_mu_and_sigma,
      tau = design$tau[i],
      log_RR = log(design$true_RR[i]),
      upper_bound = log(1/design$baseline_p[i])
    )

    # draw the effects for the studies
    RR <- exp(rtruncnorm(design$n_studies[i], mean = sol$x[["mu"]], sd = sol$x[["sigma"]],
                         a = -Inf, b = log(1/design$baseline_p[i])))

    # compute event probabilities in treatment group
    treat_p <- design$baseline_p[i] * RR

    # compute sizes of treatment groups
    n_treat <- design$n_control[i] * design$R[i]

    # draw observations for treatment and control group for each of
    # the n studies for this meta-analysis
    obs_treat <- vector(mode = "list", length = length(treat_p))
    obs_control <- vector(mode = "list", length = length(design$baseline_p[i]))
    for (k in 1:length(obs_treat)) {
      obs_treat[[k]] <- rbinom(n_treat, 1, treat_p[k])
      obs_control[[k]] <- rbinom(design$n_control[i], 1, design$baseline_p[i])
    }

    # make dataframe for the meta-analysis of this trial
    data <- data.frame(
      study = 1:length(obs_treat),
      treat_event = sapply(obs_treat, sum),
      treat_no_event = sapply(obs_treat, length) - sapply(obs_treat, sum),
      treat_sample_size = sapply(obs_treat, length),
      control_event = sapply(obs_control, sum),
      control_no_event = sapply(obs_control, length) - sapply(obs_control, sum),
      control_sample_size = sapply(obs_control, length),
      group_ratio = rep(design$R[i], length(obs_treat)),
      group_ratio_rev = rep(1 / design$R[i], length(obs_treat))
    )

    # make list of models to save
    model_list <- list()

    # prep data
    data_long <- data %>%
      gather(c(treat_event,control_event), key = "group", value = "event_count") %>%
      mutate(
        group = as.factor(group),
        group = recode_factor(group, control_event = "control", treat_event = "treat"),
        dummy_group = ifelse(group == "treat", 1, 0),
        n = ifelse(group == "treat", treat_sample_size, control_sample_size)
      ) %>%
      select(study, group, event_count, dummy_group, n)

    # calculate RR (with standard continuity correction as default in metafor:
    # add 0.5 to zero-studies)
    data2 <- escalc("RR", ai = treat_event, bi = treat_no_event,
                   ci = control_event, di = control_no_event,
                   data = data, add = 0.5, to = "only0", drop00 = TRUE) %>%
      filter(!is.na(yi))
    # double zero studies are dropped from data2 to be consistent with
    # binomial models

    # treatment arm continuity correction
    data2$RR_tcc <- treat_arm_correction_RR(data2$treat_event, data2$treat_no_event,
                                           data2$control_event, data2$control_no_event,
                                           k = 0.5)
    data2$log_RR_tcc <- log(data2$RR_tcc)

    data2$log_RR_tcc_var <- treat_arm_correction_logRR_var(data2$treat_event, data2$treat_no_event,
                                                          data2$control_event, data2$control_no_event,
                                                          k = 0.5)
    
    # create dataframe in long format without double-zero studies
    # for neg binomial model
    data_long2 <- data2 %>%
      gather(c(treat_event,control_event), key = "group", value = "event_count") %>%
      mutate(
        group = as.factor(group),
        group = recode_factor(group, control_event = "control", treat_event = "treat"),
        dummy_group = ifelse(group == "treat", 1, 0),
        n = ifelse(group == "treat", treat_sample_size, control_sample_size)
      ) %>%
      select(study, group, event_count, dummy_group, n)

    # estimate models

    ### conventional normal-normal models ###
    # default estimation
    m_REML <- try(
      rma(yi = yi, vi = vi, data = data2, method = "REML"),
      silent = TRUE
    )
    m_REML_tcc <- try(
      rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "REML"),
      silent = TRUE
    )

    # estimation with DL
    m_DL <- try(
      rma(yi = yi, vi = vi, data = data2, method = "DL"),
      silent = TRUE
    )
    m_DL_tcc <- try(
      rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "DL"),
      silent = TRUE
    )

    # estimation with SJ
    m_SJ <- try(
      rma(yi = yi, vi = vi, data = data2, method = "SJ"),
      silent = TRUE
    )
    m_SJ_tcc <- try(
      rma(yi = log_RR_tcc, vi = log_RR_tcc_var, data = data2, method = "SJ"),
      silent = TRUE
    )

    # save models in list
    model_list[["m_REML"]] <- list(model = m_REML)
    model_list[["m_REML_tcc"]] <- list(model = m_REML_tcc)
    model_list[["m_DL"]] <- list(model = m_DL)
    model_list[["m_DL_tcc"]] <- list(model = m_DL_tcc)
    model_list[["m_SJ"]] <- list(model = m_SJ)
    model_list[["m_SJ_tcc"]] <- list(model = m_SJ_tcc)

    ### Poisson models ###
    # Poisson random-effects model
    m_poiss <- try(
      glmer(event_count ~ 1 + group + (1 + group | study) + offset(log(n)),
            data = data_long, family = "poisson"),
      silent = TRUE
    )

    # zero-inflated Poisson, random intercept, fixed slope in ZI
    m_zip_rifs <- try(
      glmmTMB(event_count ~ 1 + group + (1 + group | study) + offset(log(n)),
              data = data_long, family = poisson,
              ziformula = ~ 1 + group + (1 | study)),
      silent = TRUE
    )
    
    # zero-inflated Poisson, fixed intercept, fixed slope in ZI
    m_zip_fifs <- try(
      glmmTMB(event_count ~ 1 + group + (1 + group | study) + offset(log(n)),
              data = data_long, family = poisson,
              ziformula = ~ 1 + group),
      silent = TRUE
    )
    
    # zero-inflated Poisson, random intercept in ZI
    m_zip_ri <- try(
      glmmTMB(event_count ~ 1 + group + (1 + group | study) + offset(log(n)),
              data = data_long, family = poisson,
              ziformula = ~ 1 + (1 | study)),
      silent = TRUE
    )
    
    # zero-inflated Poisson, fixed intercept in ZI
    m_zip_fi <- try(
      glmmTMB(event_count ~ 1 + group + (1 + group | study) + offset(log(n)),
              data = data_long, family = poisson,
              ziformula = ~ 1),
      silent = TRUE
    )
    
    # save models to list
    model_list[["m_poiss"]] <- list(model = summary(m_poiss)[-15])
    model_list[["m_zip_rifs"]] <- list(
      model = summary(m_zip_rifs)[1:10],
      cov_matrix = if (!is(m_zip_rifs, "try-error")) {
        vcov(m_zip_rifs, full = TRUE)[
          c("zi~(Intercept)", "zi~grouptreat", "grouptreat"),
          c("zi~(Intercept)", "zi~grouptreat", "grouptreat")
          ]
      }
    )
    model_list[["m_zip_fifs"]] <- list(
      model = summary(m_zip_fifs)[1:10],
      cov_matrix = if (!is(m_zip_fifs, "try-error")) {
        vcov(m_zip_fifs, full = TRUE)[
          c("zi~(Intercept)", "zi~grouptreat", "grouptreat"),
          c("zi~(Intercept)", "zi~grouptreat", "grouptreat")
          ]
      }
    )
    model_list[["m_zip_ri"]] <- list(model = summary(m_zip_ri)[1:10])
    model_list[["m_zip_fi"]] <- list(model = summary(m_zip_fi)[1:10])

    ### Binomial models ###
    # conditional binomial model (Boehning et al., 2015; Cai et al., 2010; Stijnen et al., 2010)
    m_cond_binom <- try(
      glmer(cbind(treat_event, control_event) ~ 1 + (1 | study) + offset(log(group_ratio)),
            data = data2, family = "binomial"),
      silent = TRUE
    )

    # beta binomial (with logit link)
    m_beta_binom <- try(
      glmmTMB(cbind(treat_event, control_event) ~ 1 + offset(log(group_ratio)),
              data = data2, family = betabinomial(link = "logit")),
      silent = TRUE
    )
    
    # actual beta binomial model in Kuss, 2015 (with log link)
    # can include double-zero studies
    theta_init_beta <- c(-1, 0, 1)
    names(theta_init_beta) <- c("b0", "b1", "prec")
    m_kuss_binom <- try(
      optim(
        par = theta_init_beta, 
        fn = likelihood_kuss, 
        treat_event = data$treat_event,
        treat_no_event = data$treat_no_event,
        control_event = data$control_event,
        control_no_event = data$control_no_event,
        hessian = TRUE
      ),
      silent = TRUE
    )
    
    # actual model by Cai et al., 2010
    theta_init <- c(0.1, 0.1)
    names(theta_init) <- c("gamma", "psi")
    m_cai_binom <- try(
      optim(
        par = theta_init,
        fn = likelihood_cai,
        Y1 = data2$treat_event,
        Y2 = data2$control_event,
        W = data2$group_ratio_rev,
        hessian = TRUE
      ),
      silent = TRUE
    )
    # pooled RR can then be computed using respective function based on fit
    # object and sample size ratio which is included in design table (R)
      

    # save models to list
    model_list[["m_cond_binom"]] <- list(model = summary(m_cond_binom)[-15])
    model_list[["m_beta_binom"]] <- list(model = summary(m_beta_binom)[1:10])
    model_list[["m_kuss_binom"]] <- list(model = m_kuss_binom)
    model_list[["m_cai_binom"]] <- list(model = m_cai_binom)

    # # save model list to simulation trial results
    # results_simulation_trials[[j]] <- model_list
    #
    # # save data
    # sim_data_trials[[j]] <- data

    return(list(model_list = model_list, data = data))
  } # end loop simulation trial
  
  # save simulation trial 
  save_to_file <- paste0("Okt_2019_sim_res_", i, ".rds")
  saveRDS(results_simulation_trials, file = save_to_file)
}

session_info <- sessionInfo()
saveRDS(session_info, file = "session_info.rds")




