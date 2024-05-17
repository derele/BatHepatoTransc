### All plotting and other output for the analysis is produced in this
### script

### loading the libraries
library(ggplot2)
library(GGally)
library(dplyr)
library(MASS)
library(tidyverse)
library(ggeffects)
library(ggrain)
library(interpretCI)
library(coin)

redoMetadat <- FALSE
redoDE <- FALSE
redoAnnotation <- FALSE

if(redoMetadat){ 
    source("R/2_metadata.R")
}else{
    metadata <- read.csv("intermediateData/metadata_expanded.csv")
}

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


## ## the enrichment alnalysis is currently only in tables, directly
## ## produced in 5_Enrichment.R

## redoEnrichment <- FALSE
## if(redoEnrichment){ source("R/5_Enrichment.R") } else { DETs_ALL <-
## readRDS("intermediateData/GOtermAnnot.RDS") }

## set up plotting
source("R/plot_setup.R")
## special fonts!
extrafont::loadfonts(device = "all") ## run every time

## Figure 1: distribution of parasitemia (blood smears)
DistPara <- ggplot(metadata, aes(x = Age_2category, y = Parasitemia_in_percent,
                                 fill = Age_2category)) +
    geom_rain(rain.side="l", alpha=.5) +
    facet_wrap(~Season, nrow=2, scales="free_y",
               labeller =  labeller(Season = ~ paste("Season:", .x)))+
    scale_y_continuous("Parasitemia in % blood cells") +
    scale_x_discrete("Age category") +
    guides(fill="none") +
    theme_minimal() 
ggsave("figures/Fig1.png", DistPara, bg="white")

### Man Whitney U tests for the statements about parasitemia in the
### results section ###############################################

 ## one of liver or spleen
MetaUniq <- metadata[!duplicated(metadata$SampleID),]

## Intensity of infection > 0 
MetaUniqInf <- MetaUniq[MetaUniq$Parasitemia_in_percentO>0 &
                        !is.na(MetaUniq$Parasitemia_in_percentO),] 
                                       
## 84 were positive on blood smears, but only for 
## 78 it was possible to record the intensity ... for later on

propCI(n = nrow(MetaUniq),
       p = 84/nrow(MetaUniq),
       alpha = 0.05)

84/nrow(MetaUniq)
84/nrow(MetaUniq) - 0.8176755
## 73.68% +/- 8.08%

propCI(n = nrow(MetaUniq),
       p =  99/nrow(MetaUniq),
       alpha = 0.05)

99/nrow(MetaUniq)
99/nrow(MetaUniq) - 0.8063693

median(MetaUniqInf$Parasitemia_in_percentO)
max(MetaUniqInf$Parasitemia_in_percentO)

##  Adult vs. Young

tapply(MetaUniqInf$Parasitemia_in_percentO,
       MetaUniqInf$Age_2category,
       median, na.rm=TRUE)

table(MetaUniqInf$Age_2category)

coin::wilcox_test(Parasitemia_in_percentO ~ as.factor(Age_2category),
                  data=MetaUniqInf)

## Male vs. Female
table(MetaUniqInf$Sex)

tapply(MetaUniqInf$Parasitemia_in_percent, MetaUniqInf$Sex, median)

coin::wilcox_test(Parasitemia_in_percentO ~ as.factor(Sex),
                  data=MetaUniqInf)

coin::wilcox_test(Parasitemia_in_percentO ~ as.factor(Sex),
                  data=MetaUniqInf)



## Reproductive activity
table(MetaUniqInf$Reproductive_status)

tapply(MetaUniqInf$Parasitemia_in_percent, MetaUniqInf$Reproductive_status, median)

coin::wilcox_test(Parasitemia_in_percentO ~ as.factor(Reproductive_status),
                  data=MetaUniqInf)

## Season
table(MetaUniqInf$Season)

tapply(MetaUniqInf$Parasitemia_in_percent, MetaUniqInf$Season, median)

coin::wilcox_test(Parasitemia_in_percentO ~ as.factor(Season),
                  data=MetaUniqInf)

wilcox.test(Parasitemia_in_percentO ~ as.factor(Season),
                  data=MetaUniqInf)




### Figure 2 the relation of blood parasitemia and transcriptom intensity


Fig2a <- ggplot(metadata, 
       aes(x = Infection_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.25))+
    scale_y_log10()+
    labs(x = "infection status
       (*by microscopy of blood smears and PCR)", 
       y = "# of reads of Hepatocystis transcripts") +
    theme_bw()

#parasitemia model
model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ 
                              Parasitemia_in_percent * Organ + 
                            offset(log(Sequencing_depth)),
       metadata)

summary(model_parasitemia)
 

Fig2b <- ggpredict(model_parasitemia, 
              terms = c("Parasitemia_in_percent", "Organ"),
          condition = c(Sequencing_depth = mean(metadata$Sequencing_depth))) %>%
  ggplot(aes(x = x, y = predicted)) +
##   geom_line(aes(color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
              alpha = .1) +
  #ylim(c(0, 2e+07)) +
  geom_point(data = metadata, aes(x = Parasitemia_in_percent,
             y = hepatocystis_transcriptome_parasitemia*
             mean_correction_factor,
             color = Organ)) +
    scale_y_log10() +
    labs(x = "Parasitemia in % blood cells", 
       y = "# reads for Hepatocystis transcripts") +
    theme_bw()


Fig2c <- ggplot(metadata, 
                aes(x = Organ, 
                    y = rpmh, 
                    color = Infection_status_blood)) +
    geom_point()+
    geom_line(aes(group = SampleID))+
    scale_y_log10()+
    labs(x = "organ", 
         y = "# of reads for Hepatocystis transcripts")+
    theme_bw() +
    scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold"))



Fig2d <- as_tibble(metadata) %>% 
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

Fig2 <- ggpubr::ggarrange(Fig2a, Fig2b,
                          Fig2c, Fig2d,
                          labels = c("a", "b", "c", "d"),
                          ncol = 2, nrow = 2)


ggsave("figures/Fig2.png", Fig2)


### P-VALUE THRESHOLDS!!!
### SEE Script 5_Enrichment.R!!

### Figure 3 The first basic DE tests using blood stage intenisty.

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


### Linking the organs!!!!
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




