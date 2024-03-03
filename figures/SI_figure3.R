library(tidyverse)
library(ggplot2)

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

main_df <- readRDS('../data/output/GTM_MEX_ZAF_clogit_output_df_2dose.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'MEX+ZAF+GTM')
  
gtm_df <- readRDS('../data/output/GTM_clogit_output_df_2dose.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'GTM')

mex_df <- readRDS('../data/output/MEX_clogit_output_df_2dose.rds') %>%
  filter(!sympt2 %in% NA) %>%
  mutate(region = 'MEX')

zaf_df <- readRDS('../data/output/ZAF_clogit_output_df_2dose.rds') %>%
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
    filter(synd_def == 'CLI (Pairwise Combination)') 
  
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
  
  fig2b <- ggplot(data = dist_df2, aes(sympt_ref, y=diff)) +
    geom_boxplot(show.legend=FALSE, color='black', width=0.2) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust=0.95),
          legend.position="top",
          legend.box = "vertical",
          legend.direction = "vertical",
          legend.spacing.x = unit(0, 'cm'))+
    guides(shape = guide_legend(title.position="bottom", title.hjust = 0.5)) +
    labs(x='Anchor Symptom', y=expression(Change~In~VE~(VE['Omicron']-VE['Delta'])), color="", shape="Secondary Symptom") +
    ylim(-1, 0)
  
  fig2a <- ggplot(data=delt_df, aes(y=diff)) +
    geom_boxplot(show.legend=FALSE, color='black', width=0.2) +
    theme_classic(base_size = 10) +
    theme(
      legend.position="top",
      legend.box = "vertical",
      legend.direction = "vertical",
      legend.spacing.x = unit(0, 'cm'),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    guides(shape = guide_legend(title.position="bottom", title.hjust = 0.5)) +
    labs(x='', y=expression(Change~In~VE~(VE['Omicron']-VE['Delta'])), color="", shape="Secondary Symptom") +
    ylim(-.6, 0)
  print(fig2a)
  output_2b <- data.frame(ggplot_build(fig2b)$data)
  output_2b$ref_symptom <- sympt_order$ref_sympt
  output_2b$region <- rep(region, nrow(output_2b))
  output_2b$phase <- rep(phase, nrow(output_2b))
  return(output_2b)
}

main_delt <- process_change(main_df, 'MEX+ZAF+GTM', 'Delta Surge-Omicron Surge')  
mex_delt <- process_change(mex_df, 'MEX', 'Delta Surge-Omicron Surge') 
gtm_delt <- process_change(gtm_df, 'GTM', 'Delta Surge-Omicron Surge')
zaf_delt <- process_change(zaf_df, 'ZAF', 'Delta Surge-Omicron Surge')

output_df <- rbind(main_delt, mex_delt, gtm_delt, zaf_delt)
output_df$region <- factor(output_df$region, levels=c('MEX+ZAF+GTM',
                                                      'MEX',
                                                      'ZAF',
                                                      'GTM'))

final_output <- output_df
final_output$ref_symptom <- factor(final_output$ref_symptom, levels=main_delt$ref_symptom)
pd <- position_dodge(0.1)
si_fig <- ggplot(data=final_output) +
  geom_jitter(aes(x=region, y=middle, fill=ref_symptom), shape=21, color="black", alpha=0.9, position=pd) +
  theme_classic() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x='', y='Median Change in VE', fill='Anchor Symptom') +
  theme(axis.text.x = element_text(angle = 45, hjust=0.95),
        plot.margin = margin(t = 2, r = .25, b = 0.25, l = 1.75, "cm")) +
  ylim(-.6, 0.0) #+
si_fig

png('../figure_pdfs/SI_figure3.png',
   width=12, height=10, units="cm", res=500)
print(si_fig)
dev.off()


