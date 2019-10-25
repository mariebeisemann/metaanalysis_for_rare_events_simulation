library(tidyverse)

# set working directory
setwd("/Users/mariebeisemann/Documents/03_Paper/03_in_preparation/2019_rare_events_ma/07_Tables")

path_to_data <- "/Users/mariebeisemann/Documents/03_Paper/03_in_preparation/2019_rare_events_ma/06_Data_Prep"

#### table for bias Cai model ####

file <- paste0(path_to_data, "/all_biases_plot.rds")
load(file)

bias_table_cai <- all_biases_plot %>%
  ungroup() %>%
  filter(model == "cai_binom") %>%
  select(-c(condition, model))

write.csv2(bias_table_cai, file = "bias_table_cai.csv")

#### table for MC se Cai model ####

file <- paste0(path_to_data, "/mc_se_plot.rds")
load(file)

mc_se_table_cai <- mc_se_plot %>%
  ungroup() %>%
  filter(model == "cai_binom") %>%
  select(-c(condition, model))

write.csv2(mc_se_table_cai, file = "mc_se_table_cai.csv")  

### table for bias in tau2 ###

file <- paste0(path_to_data, "/all_biases_plot_tau2_realdata.rds")
load(file)
  
median_bias_tau2_table <- all_biases_plot_tau2 %>%
  ungroup() %>%
  select(-c(condition, bias_tau2)) %>%
  spread(key = model, value = median_bias_tau2)

write.csv2(median_bias_tau2_table, file = "median_bias_tau2_table.csv") 
  
mean_bias_tau2_table <- all_biases_plot_tau2 %>%
  ungroup() %>%
  select(-c(condition, median_bias_tau2)) %>%
  spread(key = model, value = bias_tau2)  

write.csv2(mean_bias_tau2_table, file = "mean_bias_tau2_table.csv") 
  
### table for RMSE for RR for Cai and zip with slope ####

file <- paste0(path_to_data, "/plot_rmse_data.rds")
load(file)

file <- paste0(path_to_data, "/plot_measures_data_cai.rds")
load(file)

# join data for zip with slopes and data for cai
table_rmse_cai_zip <- plot_measures_data_cai %>% 
  filter(measure == "RMSE") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

table_rmse_cai_zip <- table_rmse_cai_zip[,c("condition", "true_RR", "true_baseline_p", 
                                            "true_tau","group_ratio", "n_studies", 
                                            "n_control", "model", "value")]

table_rmse_cai_zip <- table_rmse_cai_zip %>%
  rbind(
    plot_rmse_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = rmse_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(rmse_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(RMSE = value) %>%
  select(-value)

write.csv2(table_rmse_cai_zip, file = "table_rmse_cai_zip.csv")

#### table for MAE for RR for Cai and zip with slope ####

file <- paste0(path_to_data, "/plot_mae_data.rds")
load(file)

file <- paste0(path_to_data, "/plot_measures_data_cai.rds")
load(file)

# join data for zip with slopes and data for cai
table_mae_cai_zip <- plot_measures_data_cai %>% 
  filter(measure == "MAE") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

table_mae_cai_zip <- table_mae_cai_zip[,c("condition", "true_RR", "true_baseline_p", 
                                            "true_tau","group_ratio", "n_studies", 
                                            "n_control", "model", "value")]

table_mae_cai_zip <- table_mae_cai_zip %>%
  rbind(
    plot_mae_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = mae_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(mae_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(MAE = value) %>%
  select(-value)

write.csv2(table_mae_cai_zip, file = "table_mae_cai_zip.csv")

### table for ME for RR for Cai and zip with slope ####

file <- paste0(path_to_data, "/plot_me_data.rds")
load(file)

file <- paste0(path_to_data, "/plot_measures_data_cai.rds")
load(file)

# join data for zip with slopes and data for cai
table_me_cai_zip <- plot_measures_data_cai %>% 
  filter(measure == "ME") %>% 
  select(-c(measure, true_RR_f, true_tau_f))

table_me_cai_zip <- table_me_cai_zip[,c("condition", "true_RR", "true_baseline_p", 
                                          "true_tau","group_ratio", "n_studies", 
                                          "n_control", "model", "value")]

table_me_cai_zip <- table_me_cai_zip %>%
  rbind(
    plot_me_data %>%
      ungroup() %>%
      filter(model == "zip_rifs" | model == "zip_fifs") %>%
      mutate(value = me_log_RR,
             condition = as.character(condition),
             model = as.character(model)) %>%
      select(-c(me_log_RR, true_RR_f, true_tau_f))
  ) %>%
  mutate(ME = value) %>%
  select(-value)

write.csv2(table_me_cai_zip, file = "table_me_cai_zip.csv")






  
  
  
  
  
  
  


  