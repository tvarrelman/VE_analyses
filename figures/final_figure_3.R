library(tidyverse)
library(ggplot2)

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

mild_df <- readRDS('../data/output/GTM_MEX_ZAF_clogit_output_df_2dose_mild.rds') %>%
  mutate(clinical_def = 'mild')
sev_df <- readRDS('../data/output/GTM_MEX_ZAF_clogit_output_df_2dose_severe.rds') %>%
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

pdf('../figure_pdfs/Figure_3.pdf',
    width=3.54331, height=2.6574825)
print(severity_fig)
dev.off()

