### plotting correlation plots for hepatocystis transcriptome
### intensities between spleen and liver

### loading the libraries
library(ggplot2)
library(GGally)
library(dplyr)
library(MASS)
library(tidyverse)
library(ggeffects)
library(ggrain)

redoMetadat <- FALSE
redoDE <- FALSE
redoAnnotation <- FALSE

if(redoMetadat){ 
    source("R/2_metadata.R")
}else{
    metadata <- read.csv("intermediateData/metadata_expanded.csv")
}

## Figure 1: distribution of parasitemia (blood smears)
ggplot(metadata, aes(x = Age_2category, y = Parasitemia_in_percent,
                     fill = Age_2category)) +
        geom_rain(rain.side="l", alpha=.5)

ggplot(metadata, aes(x = Age_2category, y = Parasitemia_in_percent+1,
                     fill = Age_2category)) +
    geom_rain(rain.side="l", alpha=.5) +
    scale_y_log10("parasitemia (log10 +1)")


if(redoAnnotation){
    source("R/2_annotation.R")
} else{
   allLocusGO <- readRDS("intermediateData/GOtermAnnot.RDS")
}

if(redoDE){
    source("R/4_DE_analysis.R")
} else {
  res_ALL  <- readRDS("intermediateData/DETs_ALL.RDS")
}

### P-VALUE THRESHOLDS!!!
## a 5% ADJUSTED p-value threshold
adjpval_thresh <- 0.05
## no fold-change threshold
fc_thresh <- 0


## adding a column for significance TRUE/FALSE
res_ALL <- lapply(res_ALL, function(x){
    x$significance <- abs(x$log2FoldChange) > fc_thresh &
        x$padj < adjpval_thresh
    x
})


## ## the enrichment alnalysis is currently only in tables, directly
## ## produced in 5_Enrichment.R

## redoEnrichment <- FALSE
## if(redoEnrichment){ source("R/5_Enrichment.R") } else { DETs_ALL <-
## readRDS("intermediateData/GOtermAnnot.RDS") }

## set up plotting
source("R/plot_setup.R")
## special fonts!
extrafont::loadfonts(device = "all") ## run every time

### The first basic DE tests using blood stage intenisty.

Liver_Pas_Plot <- res_ALL[["liverPas:Parasitemia_in_percent"]] %>%
    as.data.frame() %>%
    ggplot(aes(log2FoldChange, -log10(padj),
               fill = significance)) +
    scale_y_continuous("Negative logarithm of adjusted p-value (-log10(adj-pval))",
                       limits=c(0, 5.2)) +
    scale_x_continuous("Log2 Fold-Change (per % blood parasitemia change)",
                       limits=c(-4.2, 6.1)) +
    geom_point(shape = 21, color = "white", size = 5, stroke = 0.7) +
    scale_fill_manual(values = colors_significance, labels = c("> 0.05", "<= 0.05")) +
    ggtitle("Differential expression in the liver") +
    theme(legend.position = c(1, 0.84),  #  position within the plot
          legend.background = element_rect(color = NA, fill = "white"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_legend(title = "adjusted\np-value")) 



Spleen_Pas_Plot <- res_ALL[["spleenPas:Parasitemia_in_percent"]] %>%
    as_tibble() %>%
    dplyr::filter(!is.na(significance))%>%
    ggplot(aes(log2FoldChange, -log10(padj),
               fill = significance)) +
    scale_y_continuous("", limits=c(0, 5.2)) +
    scale_x_continuous("Log2 Fold-Change (per % blood parasitemia change)",
                       limits=c(-4.2, 6.1)) +
    scale_fill_manual(values = colors_significance, labels = c("> 0.05", "<= 0.05")) +
    geom_point(shape = 21, color = "white", size = 5, stroke = 0.7) +
    ggtitle("Differential expression in the spleen") +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5))    



Liver_Trans_Plot <- res_ALL[["liver:rpmh_scaled"]] %>%
    as_tibble() %>%
    dplyr::filter(!is.na(significance))%>%
    ggplot(aes(log2FoldChange, -log10(padj),
               fill = significance)) +
    scale_y_continuous("Negative logarithm of adjusted p-value (-log10(adj-pval))",
                       limits=c(0, 5.4)) +
    scale_x_continuous("Log2 Fold-Change (per scaled rpmh)",
                       limits=c(-2.5, 2.5)) +
    scale_fill_manual(values = colors_significance, labels = c("> 0.05", "<= 0.05")) +
    geom_point(shape = 21, color = "white", size = 5, stroke = 0.7) +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5))    


Spleen_Trans_Plot <- res_ALL[["spleen:rpmh_scaled"]] %>%
    as_tibble() %>%
    dplyr::filter(!is.na(significance))%>%
    ggplot(aes(log2FoldChange, -log10(padj),
               fill = significance)) +
    scale_y_continuous("", limits=c(0, 5.4)) +
    scale_x_continuous("Log2 Fold-Change (per scaled rpmh)",
                       limits=c(-2.5, 2.5)) +
    scale_fill_manual(values = colors_significance, labels = c("> 0.05", "<= 0.05")) +
    geom_point(shape = 21, color = "white", size = 5, stroke = 0.7) +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5))    


Fig3 <- ggpubr::ggarrange(Liver_Pas_Plot, Spleen_Pas_Plot,
                          Liver_Trans_Plot, Spleen_Trans_Plot,
                          labels = c("a", "b", "c", "d"),
                          ncol = 2, nrow = 2)


# Scatter Plot of transcriptome infection estimates coloured by organ
ggsave("figures/Fig3.png", Fig3, width = 16, height = 16, units = "in",
       bg = 'white')

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
                  y = hepatocystis_transcriptome_parasitemia, 
                  color = Organ)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm") +
  labs(title = "Scatter Plot of transcriptome infection estimates coloured by organ", 
       x = "blood parasitemia in % 
       (percent of infected erythrocytes)", 
       y = "# reads of Hepatocystis transcripts)") +
  theme_bw()

dev.off()

# Boxplot of transcriptome infection estimates coloured by organ
pdf("plots/Boxplot_of_transcriptome_infection_estimates_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           fill = Organ)) +
  geom_boxplot() +
  labs(title = "Boxplot of transcriptome infection estimates coloured by organ", 
       y = "transcriptome infection estimate (#reads)") +
  theme_bw()

dev.off()

#parasitemia model
model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ 
                              Parasitemia_in_percent * Organ + 
                            offset(log(Sequencing_depth)),
       metadata)

summary(model_parasitemia)


pdf("plots/model_parasitemia.pdf")
ggpredict(model_parasitemia, 
              terms = c("Parasitemia_in_percent", "Organ"),
          condition = c(Sequencing_depth = mean(metadata$Sequencing_depth))) %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
              alpha = .1) +
  #ylim(c(0, 2e+07)) +
  geom_point(data = metadata, aes(x = Parasitemia_in_percent,
             y = hepatocystis_transcriptome_parasitemia*
             mean_correction_factor,
             color = Organ)) +
  scale_y_log10() +
  labs(x = "blood parasitemia in %
       (percentage of infected erythrocytes)",
       y = "# of reads of Hepatocystis transcripts") +
  theme_bw()

dev.off()
 

# boxplot of transcripts correlations by infection status coloured by organ
pdf("plots/boxplot_of_parasitemia_infected_erythrocytes_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Infection_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "infection status
       (*by microscopy of blood smears and PCR)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations coloured by organ") +
  theme_bw()

dev.off()

# scatter plot of transcripts correlation by organ coloured by infection status
pdf("plots/scatter_plot_of_organ_paired_parasitemia.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Infection_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations by organ coloured by infection status") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold"))
  
dev.off()

# scatter plot of transcripts correlation by organ coloured by infection status and using the rpmh values
pdf("plots/scatter_plot_of_organ_paired_parasitemia_using_rpmh_values.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = rpmh, 
           color = Infection_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations by organ coloured 
       by infection status (using rpmh values)") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold"))

dev.off()


# scatter plot of blood parasitemia count coloured by organ
pdf("plots/scatter_plot_of_blood_parasitemia_count.pdf")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "blood parasitemia in % 
       (percent of infected erythrocytes)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Blood parasitemia count coloured by organ") +
  theme_bw()

dev.off()


# box plot of infected vs undetected coloured by organ
pdf("plots/boxplot_of_transcriptome_intensity_in_liver_vs_spleen_samples.pdf")

ggplot(metadata, 
       aes(x = Infection_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "infection status
       (* by microscopy of blood smears and PCR)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcriptome intensity in liver vs spleen samples based on infection status") +
  theme_bw()

dev.off()


# scatter plot of liver vs spleen parasitemia coloured by infected and uninfected 
#for samples with both spleen and liver tissues
pdf("plots/scatter_plot_for_spleen_and_liver_infection_status_in_samples_with_both_tissues.pdf")
as_tibble(metadata) %>% 
  dplyr::select(hepatocystis_transcriptome_parasitemia,
                Organ, SampleID, Infection_status_blood) %>%
  pivot_wider(values_from = c(hepatocystis_transcriptome_parasitemia),
              names_from = Organ,
              values_fill = NA) %>% 
  filter(!is.na(Liver) & !is.na(Spleen)) %>%
  ggplot(aes(Liver, Spleen, color = Infection_status_blood)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray") +  # Reference line
  scale_y_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  scale_x_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  labs(title = "Liver vs spleen parasitemia in samples with both tissues") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold")) #+
# coord_fixed(ratio = 1) +
# xlim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE))) +
# ylim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE)))

dev.off()


# scatter plot of liver vs spleen parasitemia coloured by infected and uninfected using rpmh values 
#for samples with both spleen and liver tissues
pdf("plots/scatter_plot_for_spleen_and_liver_infection_status_in_samples_with_both_tissues_using_rpmh_values.pdf")
as_tibble(metadata) %>% 
  dplyr::select(rpmh,
                Organ, SampleID, Infection_status_blood) %>%
  pivot_wider(values_from = c(rpmh),
              names_from = Organ,
              values_fill = NA) %>% 
  filter(!is.na(Liver) & !is.na(Spleen)) %>%
  ggplot(aes(Liver, Spleen, color = Infection_status_blood)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray") +  # Reference line
  scale_y_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  scale_x_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  labs(title = "Liver vs spleen parasitemia in samples with both tissues using rpmh values") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold")) #+
# coord_fixed(ratio = 1) +
# xlim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE))) +
# ylim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE)))

dev.off()



