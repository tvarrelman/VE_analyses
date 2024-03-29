library(survival)
library(broom)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggrepel)
library(viridis)  
library(ggbeeswarm)

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

# function that calculates the regression
ve.function <- function(data) {
  main_df <- data %>%
    drop_na() %>%
    filter(e3 %in% c(1,2)) %>%
    mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))
  final_df <- data.frame()
  for (sympt in unique(main_df$label)){
    data_sub <- main_df %>% 
    filter(label == sympt) %>% 
    mutate(strata_id = paste0(iso_3,"_",e3,"_",e4))
    cl_output <- clogit(synd_id ~ vacc_id + strata(strata_id), data_sub, method="approximate")
    tidy_out <- tidy(cl_output, conf.int = TRUE) %>%
      mutate(estimate = estimate) %>%
      mutate(conf.low = conf.low) %>%
      mutate(conf.high = conf.high) %>%
      mutate(model_name = "clogig") %>%
      mutate(label = sympt) %>%
      mutate(sympt1 = ifelse(grepl('&', sympt), str_split(sympt, " & ")[[1]][[1]], sympt)) %>%
      mutate(sympt2 = ifelse(grepl('&', sympt), str_split(sympt, " & ")[[1]][[2]], NA))
    
    final_df <- rbind(final_df, tidy_out)
    
  }
  
  inter_coef <- final_df %>%
    filter(term %in% c("vacc_id:wave")) %>%
    arrange(desc(estimate))
  
  main_coef <- final_df %>%
    filter(term %in% c("vacc_id")) %>%
    mutate(ve = 1 - exp(estimate)) %>%
    mutate(or = exp(estimate)) %>%
    mutate(conf.low.ve = 1-exp(conf.low)) %>%
    mutate(conf.high.ve = 1-exp(conf.high)) 
  return(main_coef)
}

# Below are the various datasets that are used in the analyses
# uncomment the country filters to produce the individual country estimates
main_data <- read.csv('../data/output/GTM_MEX_ZAF_regression_df_2dose.csv') #%>% 
  #filter(iso_3 == 'MEX')
  #filter(iso_3 == 'ZAF')
  #filter(iso_3 == 'GTM')
#main_data <- read.csv('../data/output/GTM_MEX_ZAF_regression_df_2dose_mild.csv')
#main_data <- read.csv('../data/output/GTM_MEX_ZAF_regression_df_2dose_severe.csv')  

# subset the data to only the omicron wave
df_omi <- main_data %>%
  filter(wave==1)

# subset the data to only the delta wave
df_delta <- main_data %>%
  filter(wave==0) 

# determine the counts for the omicron wav
df_omi_count <- main_data %>%
  filter(wave==1) %>%
  drop_na() %>%
  filter(e3 %in% c(1,2)) %>%
  mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))

# determine the counts for the delta wave
df_delta_count <- main_data %>%
  filter(wave==0) %>%
  drop_na() %>%
  filter(e3 %in% c(1,2)) %>%
  mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))

# perform the regression on the omicron wave 
omi_result <- ve.function(df_omi) %>%
  mutate(ref = 'Omicron')

# perform the regression on the delta wave
delta_result <- ve.function(df_delta) %>%
  mutate(ref = 'Delta')

# merge the output of the omicron analysis with the delta analysis for the final df
final_result <- rbind(omi_result, delta_result)

# save the output that will be used to create figures
#saveRDS(final_result, '../data/output/GTM_MEX_ZAF_clogit_output_df_2dose.rds')
#saveRDS(final_result, '../data/output/GTM_MEX_ZAF_clogit_output_df_2dose_mild.rds')
#saveRDS(final_result, '../data/output/GTM_MEX_ZAF_clogit_output_df_2dose_severe.rds')

#save indivudal country results
#saveRDS(final_result, '../data/output/MEX_clogit_output_df_2dose.rds')
#saveRDS(final_result, '../data/output/ZAF_clogit_output_df_2dose.rds')
#saveRDS(final_result, '../data/output/GTM_clogit_output_df_2dose.rds')

