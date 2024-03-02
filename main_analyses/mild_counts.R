library(tidyverse)
setwd("/Users/tannervarrelman/Documents/Comms_med_VE/data/output/")

main_data <- read.csv('GTM_MEX_ZAF_regression_df_2dose_mild_9_24_23.csv')

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

nrow(df_omi_count %>% filter(synd_id == 1))
nrow(df_omi_count %>% filter(synd_id == 0))


# determine the counts for the delta wave
df_delta_count <- main_data %>%
  filter(wave==0) %>%
  drop_na() %>%
  filter(e3 %in% c(1,2)) %>%
  mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))

nrow(df_delta_count %>% filter(synd_id == 1))
nrow(df_delta_count %>% filter(synd_id == 0))
