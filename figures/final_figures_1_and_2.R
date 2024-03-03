library(tidyverse)
library(ggplot2)
library(grid)
library(ggbeeswarm)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(circlize)

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

df <- readRDS('../data/output/GTM_MEX_ZAF_clogit_output_df_2dose.rds') %>%
  filter(!sympt2 %in% NA)

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

wave_df <- data.frame()
sympt_list2 <- c("Loss of smell or taste")
for (symptom in sympt_list2){
  wave_sub <- df %>%
    filter(label == symptom | sympt1 == symptom | sympt2==symptom) %>%
    mutate(sympt_ref = symptom) %>%
    mutate(secondary_sympt = ifelse(sympt1==symptom, sympt2, ifelse(sympt2==symptom, sympt1, NA)))
  wave_df <- rbind(wave_df, wave_sub)
}

dist_df2 <- dist_df %>%
  mutate(synd_def=ifelse(sympt2 %in% NA, 'Individual Symptom', 'CLI (Pairwise Combination)')) %>%
  filter(synd_def == 'CLI (Pairwise Combination)') 

wave_df2 <- wave_df %>%
  mutate(synd_def=ifelse(sympt2 %in% NA, 'Individual Symptom', 'CLI (Pairwise Combination)')) %>%
  filter(synd_def == 'CLI (Pairwise Combination)')

wave_df2$ref <- factor(wave_df2$ref, levels = c('Delta', 'Omicron'))
delt_wave <- wave_df2 %>%
  filter(ref=='Delta')
delt_med <- median(delt_wave$ve)
omi_wave <- wave_df2 %>%
  filter(ref=='Omicron')
omi_med <- median(omi_wave$ve)
summary(omi_wave$ve)
summary(delt_wave$ve)

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

dist_df3 <- dist_df2 %>%
  group_by(sympt_ref) %>%
  mutate(med_diff = median(diff) - (-0.3995696)) %>%
  ungroup()

dist_df3$sympt_ref <- factor(dist_df3$sympt_ref, levels=sympt_order$ref_sympt)

fig2b <- ggplot(data = dist_df3, aes(x=sympt_ref, y=diff)) +
  geom_hline(yintercept = -0.3995696, color='gray', linetype = "dashed") + 
  geom_boxplot(show.legend=FALSE, color='black', width=0.2) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.95),
        legend.position="right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'cm'))+
  guides(shape = guide_legend(title.position="bottom", title.hjust = 0.5)) +
  labs(x='Anchor Symptom', y=expression(Change~In~VE~(VE['Omicron']-VE['Delta'])), color="", shape="Secondary Symptom") +
  ylim(-.6, 0)
fig2b

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
fig2a

wave_order <- wave_df2 %>% 
  filter(ref=='Delta') %>%
  arrange(desc(ve))

wave_df2$secondary_sympt <- factor(wave_df2$secondary_sympt, levels=wave_order$secondary_sympt)
fig1a <- ggplot(data = wave_df2, aes(y=ve, x=secondary_sympt)) +
  geom_errorbar(aes(ymax=conf.low.ve, ymin=conf.high.ve), width=0.0) +
  geom_point(aes(shape=ref), show.legend=TRUE, size=2) +
  scale_shape_manual(values=c(17, 16)) +
  scale_color_viridis(discrete=TRUE, guide="legend") +
  labs(x='Symptom Paired with Loss of Smell or Taste', y=expression(Vaccine~Effectiveness~(VE)), title="", shape = "", color="", fill="") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.95),
        legend.position="none") +
  guides(shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  ylim(0, 1)
fig1a

fig1b <- ggplot(data=df, aes(x=ref, y=ve)) +
  geom_boxplot(show.legend=FALSE, color='black', width=0.2) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.95),
        legend.position="top",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'cm'))+
  guides(shape = guide_legend(title.position="bottom", title.hjust = 0.5)) +
  labs(x='Period', y=expression(Vaccine~Effectiveness~(VE))) +
  ylim(0, 1)
fig1b

output_1a <- data.frame(ggplot_build(fig1a)$data)
output_1b <- data.frame(ggplot_build(fig1b)$data)
output_2a <- data.frame(ggplot_build(fig2a)$data)
output_2b <- data.frame(ggplot_build(fig2b)$data[2])

figure1 <- ggarrange(fig1a, fig1b, nrow=1, align='h', labels=c('a', 'b'), widths=c(1,0.3))
figure1

figure2 <- ggarrange(fig2a, fig2b, nrow=1, align='h', labels=c('a', 'b'), widths=c(0.15, 1))
figure2

test_figure1 <- ggarrange(fig1a, fig1b, fig2b, fig2a, nrow=2, ncol=2, labels=c('a', 'b', 'c', 'd'), widths=c(1,0.5, 0.125, 1))

fig2 <- annotate_figure(figure2, top = text_grob("Median Date, No Filters: South Africa", 
                                      color = "red", face = "bold", size = 14))

pdf('../figure_pdfs/Figure1.pdf',
    width=7.08661, height=5.3149575)
print(figure1)
dev.off()

pdf('../figure_pdfs/Figure2.pdf',
    width=7.08661, height=5.3149575)
print(figure2)
dev.off()
