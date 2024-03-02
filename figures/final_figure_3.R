library(tidyverse)
library(ggplot2)
setwd("/Users/tannervarrelman/Documents/Comms_med_VE/data/output/")

mild_df <- readRDS('GTM_MEX_ZAF_clogit_output_df_2dose_mild_8_27_23.rds') %>%
  mutate(clinical_def = 'mild')
sev_df <- readRDS('GTM_MEX_ZAF_clogit_output_df_2dose_severe_8_27_23.rds') %>%
  mutate(clinical_def = 'severe')

df <- rbind(sev_df, mild_df)

severity_fig <- ggplot(data=df, aes(x=clinical_def, y=ve)) +
  geom_boxplot(aes(fill=ref)) +
  ylim(0,1) +
  scale_shape_manual(values=c(17, 16)) +
  labs(x="", y="Vaccine Effectiveness", fill="Wave Period") +
  theme_classic()

severity_fig

sev_out <- data.frame(ggplot_build(severity_fig)$data)

pdf('/Users/tannervarrelman/Documents/Comms_med_VE/publication_figures/Figure_3.pdf',
    width=3.54331, height=2.6574825)
print(severity_fig)
dev.off()

mild_df <- read.csv('GTM_MEX_ZAF_regression_df_2dose_mild_8_27_23.csv') %>%
  drop_na() %>%
  filter(e3 %in% c(1,2)) %>%
  mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))

severe_df <- read.csv('GTM_MEX_ZAF_regression_df_2dose_8_27_23.csv') %>%
  drop_na() %>%
  filter(e3 %in% c(1,2)) %>%
  mutate(e4 = ifelse(e4 %in% c(1,2,3), 1, ifelse(e4 %in% c(4,5,6,7), 2, NA)))

omi_mild <- mild_df %>% 
  filter(wave==1)

delta_mild <- mild_df %>%
  filter(wave==0)

omi_severe <- severe_df %>%
  filter(wave==1)

delta_severe <- severe_df %>%
  filter(wave==0)

