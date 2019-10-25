library(tidyverse)


# set working directory
setwd("/Users/mariebeisemann/Documents/03_Paper/03_in_preparation/2019_rare_events_ma/03_Plots")

path_to_data <- "/Users/mariebeisemann/Documents/03_Paper/03_in_preparation/2019_rare_events_ma/06_Data_Prep"

#### Violin plots ####

file <- paste0(path_to_data, "/violin_plot_data.rds")
load(file)

### R = 0.5
mean_data_violin <- violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 0.5) %>%
  group_by(model_numbers, n_studies, true_tau, true_baseline_p, true_RR) %>%
  summarize(mean = mean(value_RR_log_scale, na.rm = TRUE))

violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 0.5) %>%
  ggplot(aes(x = model_numbers, y = value_RR_log_scale)) +
  geom_violin(aes(fill = factor_model)) +
  scale_fill_viridis_d(name = "Model",
                       breaks = c("DL", "DL_tcc", "REML", "REML_tcc", "SJ", "SJ_tcc", 
                                  "beta_binom", "cai_binom", "cond_binom", "kuss_binom", 
                                  "poiss", "zip_fi", "zip_fifs", "zip_ri", "zip_rifs"),
                       labels = c("DL (1)", "DL_tcc (2)", "REML (3)", "REML_tcc (4)", 
                                  "SJ (5)", "SJ_tcc (6)", "beta_binom (7)", "cai_binom (8)", 
                                  "cond_binom (9)", "kuss_binom (10)", "poiss (11)", 
                                  "zip_fi (12)", "zip_fifs (13)", "zip_ri (14)", "zip_rifs (15)")) +
  geom_point(data = mean_data_violin, aes(x = model_numbers, y = mean), 
             shape = 17, color = "red", size = 1) +
  geom_hline(aes(yintercept = log(true_RR)), color = "gray60") +
  facet_grid(n_studies + true_tau ~ true_baseline_p + true_RR) +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(0.6)),
        panel.grid = element_blank()) +
  ylim(-2.5, 2.5) +
  ylab("log RR") +
  xlab("Models")

ggsave("violin_R05.pdf", width = 11, height = 8.5)

### R = 1
mean_data_violin <- violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 1) %>%
  group_by(model_numbers, n_studies, true_tau, true_baseline_p, true_RR) %>%
  summarize(mean = mean(value_RR_log_scale, na.rm = TRUE))

violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 1) %>%
  ggplot(aes(x = model_numbers, y = value_RR_log_scale)) +
  geom_violin(aes(fill = factor_model)) +
  scale_fill_viridis_d(name = "Model",
                       breaks = c("DL", "DL_tcc", "REML", "REML_tcc", "SJ", "SJ_tcc", 
                                  "beta_binom", "cai_binom", "cond_binom", "kuss_binom", 
                                  "poiss", "zip_fi", "zip_fifs", "zip_ri", "zip_rifs"),
                       labels = c("DL (1)", "DL_tcc (2)", "REML (3)", "REML_tcc (4)", 
                                  "SJ (5)", "SJ_tcc (6)", "beta_binom (7)", "cai_binom (8)", 
                                  "cond_binom (9)", "kuss_binom (10)", "poiss (11)", 
                                  "zip_fi (12)", "zip_fifs (13)", "zip_ri (14)", "zip_rifs (15)")) +
  geom_point(data = mean_data_violin, aes(x = model_numbers, y = mean), 
             shape = 17, color = "red", size = 1) +
  geom_hline(aes(yintercept = log(true_RR)), color = "gray60") +
  facet_grid(n_studies + true_tau ~ true_baseline_p + true_RR) +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(0.6)),
        panel.grid = element_blank()) +
  ylim(-2.5, 2.5) +
  ylab("log RR") +
  xlab("Models")

ggsave("violin_R1.pdf", width = 11, height = 8.5)

### R = 2
mean_data_violin <- violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 2) %>%
  group_by(model_numbers, n_studies, true_tau, true_baseline_p, true_RR) %>%
  summarize(mean = mean(value_RR_log_scale, na.rm = TRUE))

violin_plot_data %>%
  mutate(model_numbers = as.factor(as.numeric(as.factor(model)))) %>%
  filter(group_ratio == 2) %>%
  ggplot(aes(x = model_numbers, y = value_RR_log_scale)) +
  geom_violin(aes(fill = factor_model)) +
  scale_fill_viridis_d(name = "Model",
                       breaks = c("DL", "DL_tcc", "REML", "REML_tcc", "SJ", "SJ_tcc", 
                                  "beta_binom", "cai_binom", "cond_binom", "kuss_binom", 
                                  "poiss", "zip_fi", "zip_fifs", "zip_ri", "zip_rifs"),
                       labels = c("DL (1)", "DL_tcc (2)", "REML (3)", "REML_tcc (4)", 
                                  "SJ (5)", "SJ_tcc (6)", "beta_binom (7)", "cai_binom (8)", 
                                  "cond_binom (9)", "kuss_binom (10)", "poiss (11)", 
                                  "zip_fi (12)", "zip_fifs (13)", "zip_ri (14)", "zip_rifs (15)")) +
  geom_point(data = mean_data_violin, aes(x = model_numbers, y = mean), 
             shape = 17, color = "red", size = 1) +
  geom_hline(aes(yintercept = log(true_RR)), color = "gray60") +
  facet_grid(n_studies + true_tau ~ true_baseline_p + true_RR) +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(0.6)),
        panel.grid = element_blank()) +
  ylim(-2.5, 2.5) +
  ylab("log RR") +
  xlab("Models")

ggsave("violin_R2.pdf", width = 11, height = 8.5)

#### Convergence plots ####

file <- paste0(path_to_data, "/conv_plot_data.rds")
load(file)

### for R = 0.5
conv_plot_data %>%
  filter(group_ratio == 0.5) %>%
  mutate(model_names = ifelse(model == "beta_binom",
                              "beta",
                              ifelse(model == "cai_binom",
                                     "cai",
                                     ifelse(model == "cond_binom",
                                            "binom",
                                            ifelse(model == "kuss_binom",
                                                   "kuss",
                                                   model
                                                   )
                                            )
                                     )
                              ),
         condition2 = gsub("ncontrol50_", "", condition),
         condition2 = gsub("_R0.5_", "_", condition2),
         condition2 = gsub("_", ", ", condition2),
         condition2 = gsub("nstudies", "no. of studies = ", condition2),
         condition2 = gsub("RR", "RR = ", condition2),
         condition2 = gsub("p", "p = ", condition2),
         condition2 = gsub("tau", "tau = ", condition2)
         
  ) %>%
  ggplot(aes(model_names, condition2)) + 
  geom_tile(aes(fill = not_converged), colour = "white") + 
  scale_fill_gradient(low = "green", high = "red", name = "No. of trials with \nfailed convergence") +
  xlab("Model") +
  ylab("Condition")

ggsave("convergence_R05.pdf")

### for R = 1
conv_plot_data %>%
  filter(group_ratio == 1) %>%
  mutate(model_names = ifelse(model == "beta_binom",
                              "beta",
                              ifelse(model == "cai_binom",
                                     "cai",
                                     ifelse(model == "cond_binom",
                                            "binom",
                                            ifelse(model == "kuss_binom",
                                                   "kuss",
                                                   model
                                            )
                                     )
                              )
  ),
  condition2 = gsub("ncontrol50_", "", condition),
  condition2 = gsub("_R1_", "_", condition2),
  condition2 = gsub("_", ", ", condition2),
  condition2 = gsub("nstudies", "no. of studies = ", condition2),
  condition2 = gsub("RR", "RR = ", condition2),
  condition2 = gsub("p", "p = ", condition2),
  condition2 = gsub("tau", "tau = ", condition2)
  
  ) %>%
  ggplot(aes(model_names, condition2)) + 
  geom_tile(aes(fill = not_converged), colour = "white") + 
  scale_fill_gradient(low = "green", high = "red", name = "No. of trials with \nfailed convergence") +
  xlab("Model") +
  ylab("Condition")

ggsave("convergence_R1.pdf")

### for R = 2
conv_plot_data %>%
  filter(group_ratio == 2) %>%
  mutate(model_names = ifelse(model == "beta_binom",
                              "beta",
                              ifelse(model == "cai_binom",
                                     "cai",
                                     ifelse(model == "cond_binom",
                                            "binom",
                                            ifelse(model == "kuss_binom",
                                                   "kuss",
                                                   model
                                            )
                                     )
                              )
  ),
  condition2 = gsub("ncontrol50_", "", condition),
  condition2 = gsub("_R2_", "_", condition2),
  condition2 = gsub("_", ", ", condition2),
  condition2 = gsub("nstudies", "no. of studies = ", condition2),
  condition2 = gsub("RR", "RR = ", condition2),
  condition2 = gsub("p", "p = ", condition2),
  condition2 = gsub("tau", "tau = ", condition2)
  
  ) %>%
  ggplot(aes(model_names, condition2)) + 
  geom_tile(aes(fill = not_converged), colour = "white") + 
  scale_fill_gradient(low = "green", high = "red", name = "No. of trials with \nfailed convergence") +
  xlab("Model") +
  ylab("Condition")

ggsave("convergence_R2.pdf")

#### Zero plots ####

file <- paste0(path_to_data, "/zero_plots_data.rds")
load(file)

# Treatment group zeroes
zero_plots %>%
  gather(c(per_n_zeros_treat, per_n_zeros_control), key = "group", value = per_zeroes) %>%
  mutate(group = ifelse(group == "per_n_zeros_treat", "treat", "control")) %>%
  ggplot(aes(y = per_zeroes, x = true_RR_f, fill = true_tau_f)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(n_control + n_studies + group ~ group_ratio + true_baseline_p) +
  scale_fill_viridis_d(name = "True heterogenity \n(SD, log scale)") +
  theme_bw() +
  ylim(c(0, 1)) +
  xlab("True RR") +
  ylab("Relative Frequency of Zeros")

ggsave("zeroes.pdf")

#### Coverage plots ####

file <- paste0(path_to_data, "/coverage_plot_data.rds")
load(file)

model_names <- c(
  "m_REML", "m_REML_tcc", "m_DL", "m_DL_tcc", "m_SJ", "m_SJ_tcc",
  "m_poiss", "m_zip_rifs", "m_zip_fifs", "m_zip_ri", "m_zip_fi",
  "m_cond_binom",  "m_beta_binom", "m_kuss_binom", "m_cai_binom"
)

n_models <- length(model_names)

# RR = 0.5
coverage_plot %>%
  mutate(true_tau_f = as.factor(true_tau),
         true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))) %>%
  ggplot(aes(x = true_tau_f, y = coverage_log_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(linetype = model, group = model)) +
  geom_hline(aes(yintercept = 0.95), color = "black") +
  facet_grid(n_studies + true_baseline_p ~ group_ratio + true_RR_f) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") +
  ylab("Coverage of the 95% CI for the pooled log RR")

ggsave("coverage_plot.pdf", width = 11, height = 8.5)

# for models which don't estimate a log RR, but a RR (Cai and ZIP with slope)

file <- paste0(path_to_data, "/coverage_plot_not_log_data.rds")
load(file)

coverage_plot_not_log %>%
#  filter(true_RR == 2) %>%
  mutate(true_tau_f = as.factor(true_tau),
         true_RR_f = as.factor(true_RR)) %>%
  ggplot(aes(x = true_tau_f, y = coverage_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(linetype = model, group = model)) +
  geom_hline(aes(yintercept = 0.95), color = "black") +
  facet_grid(n_control + n_studies + true_RR_f ~ group_ratio + true_baseline_p) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") +
  ylab("Coverage of the 95% CI for the pooled log RR")

ggsave("coverage_plot_not_log.pdf")

#### Bias plots ####

file <- paste0(path_to_data, "/all_biases_plot.rds")
load(file)

# bias plots for models with pooled log RR
col_data <- unique(
  all_biases_plot %>% 
    ungroup() %>%
    gather(c(median_bias_log_RR, bias_log_RR), key = "measure_type", value = "value") %>%
    mutate(measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias")) %>%
    select(n_studies, measure_type2, true_RR, true_baseline_p)
)

# R = 0.5, restricted value range
plot_data <- all_biases_plot %>%
  filter(group_ratio == 0.5,
         model != "cai_binom",
         model != "zip_fifs",
         model != "zip_rifs") %>%
  gather(c(median_bias_log_RR, bias_log_RR), key = "measure_type", value = "value") %>%
  filter(value > -2) %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                                         group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ true_RR + true_baseline_p,
             scales = "free_y") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_R05_restricted_value_range.pdf")

# R = 0.5, full value range
plot_data <- all_biases_plot %>%
  filter(group_ratio == 0.5,
         model != "cai_binom",
         model != "zip_fifs",
         model != "zip_rifs") %>%
  gather(c(median_bias_log_RR, bias_log_RR), key = "measure_type", value = "value") %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                  group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ true_RR + true_baseline_p,
             scales = "free_y") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_R05_full_value_range.pdf")

# R = 1, full value range
plot_data <- all_biases_plot %>%
  filter(group_ratio == 1,
         model != "cai_binom",
         model != "zip_fifs",
         model != "zip_rifs") %>%
  gather(c(median_bias_log_RR, bias_log_RR), key = "measure_type", value = "value") %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                                         group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ true_RR + true_baseline_p,
             scales = "free_y") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_R1_full_value_range.pdf")

# R = 2, full value range
plot_data <- all_biases_plot %>%
  filter(group_ratio == 2,
         model != "cai_binom",
         model != "zip_fifs",
         model != "zip_rifs") %>%
  gather(c(median_bias_log_RR, bias_log_RR), key = "measure_type", value = "value") %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                  group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ true_RR + true_baseline_p) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_R2_full_value_range.pdf")

# bias plots for models with pooled RR (cai and zip with slope)
col_data <- unique(
  all_biases_plot %>% 
    gather(c(bias_log_RR, median_bias_log_RR), key = "measure_type", value = "value") %>%
    mutate(measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias")) %>%
    select(n_studies, measure_type2, group_ratio, true_baseline_p)
)

plot_data <- all_biases_plot %>%
  filter(
         model == "zip_fifs" | model == "zip_rifs"| model == "cai_binom") %>%
  gather(c(bias_log_RR, median_bias_log_RR), key = "measure_type", value = "value") %>%
  mutate(true_tau_f = as.factor(true_tau),
         true_RR_f = as.factor(ifelse(true_RR == 0.5, "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")
                                      )
                               ),
         measure_type2 = ifelse(measure_type == "bias_log_RR", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                  group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ group_ratio + true_baseline_p + true_RR_f) +
  ylim(-0.2, 0.8) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies) ")

ggsave("bias_plot_cai_and_zip_restricted_value_range.pdf", 
       width = 297, height = 210, units = "mm")

#### MC SE plots ####

file <- paste0(path_to_data, "/mc_se_plot.rds")
load(file)

# MC SE plot for models with pooled log RR estimate, restricted value range
mc_se_plot %>%
  filter(model != "cai_binom",
         model != "zip_rifs",
         model != "zip_fifs",
         mc_SE_log_RR < 0.5) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
         ) %>%
  ggplot(aes(x = true_tau_f, y = mc_SE_log_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(color = model, group = model, linetype = model)) +
  facet_grid(true_RR_f + n_studies ~ group_ratio + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("Monte Carlo Standard Error") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mc_se_plot_restricted_value_range.pdf")

# MC SE plot for models with pooled log RR estimate, full value range
mc_se_plot %>%
  filter(model != "cai_binom",
         model != "zip_rifs",
         model != "zip_fifs") %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  ggplot(aes(x = true_tau_f, y = mc_SE_log_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(color = model, group = model, linetype = model)) +
  facet_grid(true_RR_f + n_studies ~ group_ratio + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("Monte Carlo Standard Error") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mc_se_plot_full_value_range.pdf")

# MC SE plot for the models with pooled RR (not log scale), i.e. 
# Cai and zip models with slope

# restricted value range
mc_se_plot %>%
  filter(model == "zip_rifs" | model == "zip_fifs" | model == "cai_binom",
         mc_SE_log_RR < 1) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  ggplot(aes(x = true_tau_f, y = mc_SE_log_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(color = model, group = model, linetype = model)) +
  facet_grid(true_RR_f + n_studies ~ group_ratio + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("Monte Carlo standard error") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mc_se_plot_cai_and_zip_restricted_value_range.pdf")

# full value range
mc_se_plot %>%
  filter(model == "zip_rifs" | model == "zip_fifs" | model == "cai_binom") %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  ggplot(aes(x = true_tau_f, y = mc_SE_log_RR, color = model, shape = model)) +
  geom_point() +
  geom_line(aes(color = model, group = model, linetype = model)) +
  facet_grid(true_RR_f + n_studies ~ group_ratio + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("Monte Carlo standard error") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mc_se_plot_cai_and_zip_full_value_range.pdf")

#### Bias plots for tau2 (online supplement only) ####

file <- paste0(path_to_data, "/all_biases_plot_tau2_realdata.rds")
load(file)

col_data <- unique(
  all_biases_plot_tau2 %>% 
    gather(bias_tau2:median_bias_tau2, key = "measure_type", value = "value") %>%
    mutate(measure_type2 = ifelse(measure_type == "bias_tau2", "Mean Bias", "Median Bias")) %>%
    select(n_studies, measure_type2, group_ratio, true_baseline_p)
)

# for RR = 0.5, restricted value range
plot_data <- all_biases_plot_tau2 %>%
  filter(true_RR == 0.5) %>%
  gather(bias_tau2:median_bias_tau2, key = "measure_type", value = "value") %>%
  filter(value < 5) %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_tau2", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.05) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                                         group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ group_ratio + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_tau2_RR05_restricted_value_range.pdf")
# full value range shown in table

# for RR = 1, restricted value range
plot_data <- all_biases_plot_tau2 %>%
  filter(true_RR == 1) %>%
  gather(bias_tau2:median_bias_tau2, key = "measure_type", value = "value") %>%
  filter(value < 2) %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_tau2", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.05) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                  group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ group_ratio + true_baseline_p) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_tau2_RR1_restricted_value_range.pdf")
# full value range shown in table

# for RR = 2, restricted value range
plot_data <- all_biases_plot_tau2 %>%
  filter(true_RR == 2) %>%
  gather(bias_tau2:median_bias_tau2, key = "measure_type", value = "value") %>%
  filter(value < 2) %>%
  mutate(true_tau_f = as.factor(true_tau),
         measure_type2 = ifelse(measure_type == "bias_tau2", "Mean Bias", "Median Bias"))

ggplot() +
  geom_rect(data = col_data, aes(fill = measure_type2),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.05) +
  scale_fill_manual(values = c("white", "grey"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, aes(x = true_tau_f, y = value, color = model, 
                                  group = model, linetype = model)) +
  facet_grid(n_studies + measure_type2 ~ group_ratio + true_baseline_p) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("(Median) Bias") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("bias_plot_tau2_RR2_restricted_value_range.pdf")
# full value range shown in table

#### RMSE plots ####

# for models with provide an estimate for the pooled log RR

file <- paste0(path_to_data, "/plot_rmse_data.rds")
load(file)

# restricted value range 
plot_data <- plot_rmse_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs", 
         rmse_log_RR < 2)

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = rmse_log_RR, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("rmse_plot_restricted_value_range.pdf")

# larger value range 
plot_data <- plot_rmse_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs",
         rmse_log_RR < 10)

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = rmse_log_RR, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("rmse_plot_larger_value_range.pdf")

# full value range 
plot_data <- plot_rmse_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs")

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = rmse_log_RR, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("rmse_plot_full_value_range.pdf")

# RMSE for models which estimate the pooled RR, i.e., Cai and zip with slope

file <- paste0(path_to_data, "/plot_measures_data_cai.rds")
load(file)

# join data for zip with slopes and data for cai
plot_data <- plot_measures_data_cai %>% 
  filter(measure == "RMSE") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

plot_data <- plot_data[,c("condition", "true_RR", "true_baseline_p", "true_tau",
                          "group_ratio", "n_studies", "n_control", "model", "value")]

plot_data <- plot_data %>%
  rbind(
    plot_rmse_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = rmse_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(rmse_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))),
    true_tau_f = as.factor(true_tau)
  )

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

# restricted value range
plot_data_r <- plot_data %>%
  filter(value < 10)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data_r,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data_r, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("rmse_cai_and_zip_restricted_value_range.pdf")

# full value range
ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("rmse_cai_and_zip_full_value_range.pdf")
# also make a table for full value range

#### MAE plots ####

file <- paste0(path_to_data, "/plot_mae_data.rds")
load(file)

# for models which provide estimate for the pooled log RR

# restricted value range
plot_data <- plot_mae_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs",
         mae_log_RR < 3)

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = mae_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = mae_log_RR, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mae_plot_restricted_value_range.pdf")

# full value range
plot_data <- plot_mae_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs")

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = mae_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = mae_log_RR, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("mae_plot_full_value_range.pdf")

# MAE for models which estimate the pooled RR, i.e., Cai and zip with slope

# join data for zip with slopes and data for cai
plot_data <- plot_measures_data_cai %>% 
  filter(measure == "MAE") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

plot_data <- plot_data[,c("condition", "true_RR", "true_baseline_p", "true_tau",
                          "group_ratio", "n_studies", "n_control", "model", "value")]

plot_data <- plot_data %>%
  rbind(
    plot_mae_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = mae_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(mae_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))),
    true_tau_f = as.factor(true_tau)
  )

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

# restricted value range
plot_data_r <- plot_data %>%
  filter(value < 10)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data_r,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data_r, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("mae_cai_and_zip_restricted_value_range.pdf")

# full value range
ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("mae_cai_and_zip_full_value_range.pdf")
# also make a table for full value range

#### ME plots ####

file <- paste0(path_to_data, "/plot_me_data.rds")
load(file)

# for models which provide estimate for the pooled log RR

# restricted value range
plot_data <- plot_me_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs",
         me_log_RR < 10)

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_log_RR, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("me_plot_restricted_value_range.pdf")

# fuller value range
plot_data <- plot_me_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs",
         me_log_RR < 50)

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_log_RR, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("me_plot_fuller_value_range.pdf")

# full value range
plot_data <- plot_me_data %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                      "RR: 0.5",
                                      ifelse(true_RR == 1, "RR: 1", "RR: 2")))
  ) %>%
  filter(model != "zip_rifs", 
         model != "zip_fifs")

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_log_RR, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_log_RR, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)")

ggsave("me_plot_full_value_range.pdf")

# MAE for models which estimate the pooled RR, i.e., Cai and zip with slope

# join data for zip with slopes and data for cai
plot_data <- plot_measures_data_cai %>% 
  filter(measure == "ME") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

plot_data <- plot_data[,c("condition", "true_RR", "true_baseline_p", "true_tau",
                          "group_ratio", "n_studies", "n_control", "model", "value")]

plot_data <- plot_data %>%
  rbind(
    plot_me_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = me_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(me_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5, 
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))),
    true_tau_f = as.factor(true_tau)
  )

col_data2 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

# restricted value range
plot_data_r <- plot_data %>%
  filter(value < 50)

ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data_r,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data_r, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("me_cai_and_zip_restricted_value_range.pdf")

# full value range
ggplot() +
  geom_rect(data = col_data2, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = value, color = model, shape = model)) +
  geom_line(data = plot_data, 
            aes(x = true_tau_f, y = value, color = model,
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (standard deviation of true log RR for different studies)") 

ggsave("me_cai_and_zip_full_value_range.pdf")
# also make a table for full value range

#### RMSE for tau2 ####

file <- paste0(path_to_data, "/plot_rmse_data_tau2.rds")
load(file)

# restricted value range
plot_data <- plot_rmse_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  filter(rmse_tau2 < 2) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = rmse_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("rmse_plot_tau2_restricted_value_range.pdf")

# fuller value range
plot_data <- plot_rmse_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  filter(rmse_tau2 < 10) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = rmse_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("rmse_plot_tau2_fuller_value_range.pdf")

# full value range
plot_data <- plot_rmse_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = rmse_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = rmse_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("RMSE") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("rmse_plot_tau2_full_value_range.pdf")

#### MAE for tau2 ####

file <- paste0(path_to_data, "/plot_mae_data_tau2.rds")
load(file)

# restricted value range
plot_data <- plot_mae_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  filter(mae_tau2 < 2) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = mae_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = mae_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("mae_plot_tau2_restricted_value_range.pdf")

# full value range
plot_data <- plot_mae_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = mae_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = mae_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("MAE") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("mae_plot_tau2_full_value_range.pdf")

#### ME for tau2 ####

file <- paste0(path_to_data, "/plot_me_data_tau2.rds")
load(file)

# restricted value range
plot_data <- plot_me_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  filter(me_tau2 < 10) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("me_plot_tau2_restricted_value_range.pdf")

# fuller value range
plot_data <- plot_me_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  filter(me_tau2 < 80) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("me_plot_tau2_fuller_value_range.pdf")

# restricted value range
plot_data <- plot_me_data_tau2 %>%
  ungroup() %>%
  select(-true_RR_f) %>%
  mutate(
    true_RR_f = as.factor(ifelse(true_RR == 0.5,
                                 "RR: 0.5",
                                 ifelse(true_RR == 1, "RR: 1", "RR: 2"))))

col_data4 <- unique(
  plot_data %>% 
    select(n_studies, group_ratio, true_RR_f, true_baseline_p) 
)

ggplot() +
  geom_rect(data = col_data4, aes(fill = true_RR_f),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("slategray1","white", "palegoldenrod"), guide = FALSE) +
  geom_point(data = plot_data,
             aes(x = true_tau_f, y = me_tau2, color = model, shape = model)) +
  geom_line(data = plot_data,
            aes(x = true_tau_f, y = me_tau2, color = model, 
                group = model, linetype = model)) +
  facet_grid(n_studies + group_ratio ~ true_RR_f + true_baseline_p,
             scales = "free") +
  theme_bw() +
  scale_color_viridis_d() +
  scale_shape_manual(values = 1:n_models) +
  ylab("ME") +
  xlab("Degree of heterogeneity (tau on log scale)")

ggsave("me_plot_tau2_full_value_range.pdf")




