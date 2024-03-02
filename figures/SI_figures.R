library(tidyverse)
library(ggplot2)
setwd("/Users/tannervarrelman/Documents/Comms_med_VE/data/output/")

main_df <- readRDS('GTM_MEX_ZAF_clogit_output_df_2dose_8_27_23.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'MEX+ZAF+GTM')

gtm_df <- readRDS('GTM_clogit_output_df_2dose_8_27_23.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'GTM')

mex_df <- readRDS('MEX_clogit_output_df_2dose_8_27_23.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'MEX')

zaf_df <- readRDS('ZAF_clogit_output_df_2dose_8_27_23.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'ZAF')

process_change <- function(df, region, phase) {
  delt_df <- data.frame()
  for (label in unique(df$label)){
    delt_sub <- df[(df$label==label) & (df$ref=='Delta'),]
    omi_sub <- df[(df$label==label) & (df$ref=='Omicron'),]
    diff <- omi_sub$ve - delt_sub$ve
    sympt1 <- delt_sub$sympt1
    sympt2 <- delt_sub$sympt2
    inter_df <- data.frame(sympt1, sympt2, label, diff)
    delt_df <- rbind(delt_df, inter_df)
  }
  sympt_list <- c("Loss of smell or taste", "Fever", "Headache",
                  "Aches or muscle pain","Chest pain","Fatigue",               
                  "Nausea", "Stuffy or runny nose", "Chills",
                  "Difficulty breathing", "Cough", "Sore throat")  
  dist_df <- data.frame()
  for (symptom in sympt_list){
    sympt_sub <- delt_df %>%
      filter(label == symptom | sympt1 == symptom | sympt2==symptom) %>%
      mutate(sympt_ref = symptom) %>%
      mutate(secondary_sympt = ifelse(sympt1==symptom, sympt2, ifelse(sympt2==symptom, sympt1, NA)))
    dist_df <- rbind(dist_df, sympt_sub)
  }
  
  dist_df2 <- dist_df %>%
    mutate(synd_def=ifelse(sympt2 %in% NA, 'Individual Symptom', 'CLI (Pairwise Combination)')) %>%
    filter(synd_def == 'CLI (Pairwise Combination)') %>%
    mutate(country = region)
    
  
  med_df <- data.frame()
  for (ref_sympt in unique(dist_df2$sympt_ref)){
    med_sub <- dist_df2 %>%
      filter(sympt_ref == ref_sympt)
    diff_med <- median(med_sub$diff)
    diff_inter <- data.frame(ref_sympt, diff_med)
    med_df <- rbind(med_df, diff_inter)
  }
  
  sympt_order <- med_df %>%
    arrange(desc(diff_med))
  
  dist_df2$sympt_ref <- factor(dist_df2$sympt_ref, levels=sympt_order$ref_sympt)
  
  delt_df2 <- delt_df %>%
    mutate(country = region)
  
  return(list(dist_df2, delt_df2))
}

main_delt <- process_change(main_df, 'MEX+ZAF+GTM', 'Delta Surge-Omicron Surge')  
main_delt_dist_df <- data.frame(main_delt[[1]])
main_delt_df <- data.frame(main_delt[[2]])

mex_delt <- process_change(mex_df, 'MEX', 'Delta Surge-Omicron Surge') 
mex_delt_dist_df <- data.frame(mex_delt[[1]])
mex_delt_df <- data.frame(mex_delt[[2]])

gtm_delt <- process_change(gtm_df, 'GTM', 'Delta Surge-Omicron Surge')
gtm_delt_dist_df <- data.frame(gtm_delt[[1]])
gtm_delt_df <- data.frame(gtm_delt[[2]])

zaf_delt <- process_change(zaf_df, 'ZAF', 'Delta Surge-Omicron Surge')
zaf_delt_dist_df <- data.frame(zaf_delt[[1]])
zaf_delt_df <- data.frame(zaf_delt[[2]])

fig_2b_df <- rbind(main_delt_dist_df, mex_delt_dist_df, gtm_delt_dist_df, zaf_delt_dist_df)
fig_2b_df$country <- factor(fig_2b_df$country, levels=c('MEX+ZAF+GTM', 'MEX', 'ZAF', 'GTM', 'MEX-Sample'))
fig_2a_df <- rbind(main_delt_df, mex_delt_df, gtm_delt_df, zaf_delt_df)
fig_2a_df$country <- factor(fig_2a_df$country, levels=c('MEX+ZAF+GTM', 'MEX', 'ZAF', 'GTM'))

fig2a <- ggplot(data=fig_2a_df, aes(y=diff, x=country, color=country)) +
  geom_boxplot(width=0.2, show.legend=TRUE) +
  theme_classic(base_size = 10) +
  scale_color_discrete(guide = guide_legend()) +
  theme(
    legend.position="top",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0, 'cm')) +#,
  labs(x='', y=expression(Change~In~VE~(VE['Omicron']-VE['Delta'])), color="", shape="Secondary Symptom") +
  ylim(-.6, 0)
fig2a

png('/Users/tannervarrelman/Documents/Comms_med_VE/publication_figures/SI_figure1_S1_8_27_23.png',
    width=10, height=10, units="cm", res=500)
print(fig2a)
dev.off()

fig2b <- ggplot(data = fig_2b_df, aes(x=sympt_ref, y=diff, color=country)) +
  geom_boxplot(show.legend=TRUE, width=0.5) +
  theme_classic(base_size = 10) +
  scale_color_discrete(guide = guide_legend()) +
  theme(
    axis.text.x = element_text(angle = 45, hjust=0.95),
    legend.position="top",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0, 'cm')) +
  labs(x='Anchor Symptom', y=expression(Change~In~VE~(VE['Omicron']-VE['Delta'])), color="", shape="Secondary Symptom") +
  ylim(-.6, 0)
fig2b

png('/Users/tannervarrelman/Documents/Comms_med_VE/publication_figures/Figure_S2_8_27_23.png',
    width=30, height=10, units="cm", res=500)
print(fig2b)
dev.off()

output_2a <- data.frame(ggplot_build(fig2a)$data)



