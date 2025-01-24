setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/Shared drives/Crop soil microbiome/Manuscripts/Manuscript 1 - crop species/R analyses/Nhu's_analyses")
library(tidyverse)
library(vegan)
library(phyloseq)
library(qiime2R)
library(multcompView)
library(patchwork)
library(pairwiseAdonis)
library(RColorBrewer)
library(agricolae)
library(nlme)

display.brewer.all()
brewer.pal(n=9,"YlGnBu")
brewer.pal(n=9,"PuBuGn")

##===================================================================================================
## Read in 16S data
##===================================================================================================
#Read in diversity metrics. Make sure that each of the TSV files have a "SampleID" header.
Observed_otus_16S <- read.table(file="16S/observed-otus.tsv", header=TRUE)
Shannon_div_16S <- read.table(file="16S/shannon-diversity.tsv", header=TRUE)
Faith_PD_16S <- read.table(file="16S/faith-PD.tsv", header=TRUE)
PCoA_16S <- read.table(file="16S/pcoa-aitchison.tsv", header=TRUE)

#Read in metadata file and merge other metrics with the metadata file
metadata_table_16S <- read.table(file="16S/metadata-16S-LKM.tsv", header=TRUE) %>% 
  left_join(Observed_otus_16S, by="SampleID") %>% 
  left_join(Shannon_div_16S, by="SampleID") %>% 
  left_join(Faith_PD_16S, by="SampleID") %>% 
  left_join(PCoA_16S, by="SampleID") %>% 
  drop_na() %>% 
  filter(observed_features < 500) %>% #Removing one sample set of extreme outliers at C1
  as_tibble()

write.table(metadata_table_16S, "16S/metadata_table_16S.tsv", sep="\t")

#Keeping all samples
metadata_table2_16S <- read.table(file="16S/metadata-16S-LKM.tsv", header=TRUE) %>% 
  left_join(Observed_otus_16S, by="SampleID") %>% 
  left_join(Shannon_div_16S, by="SampleID") %>% 
  left_join(Faith_PD_16S, by="SampleID") %>% 
  left_join(PCoA_16S, by="SampleID") %>% 
  drop_na() %>% 
  as_tibble()

#Read in rarefied OTU table; import from tsv file
otu_table_16S <- read.table(file="16S/table-rarefied-with-taxonomy.tsv", header=TRUE, sep='\t', check.names=FALSE, row.names=1) %>%
  select(-taxonomy) %>%
  t() %>% 
  as_tibble()

#Format data tables for phyloseq
otu_16S <- otu_table(otu_table_16S, taxa_are_rows=FALSE)
metadata_16S <- sample_data(metadata_table_16S)
#metadata2_16S <- sample_data(metadata_table2_16S)
tree <- read_tree("16S/rooted-tree.nwk")
 
# #tree can also be imported directly from QIIME
# tree2 <- read_qza("16S/rooted_tree.qza")
# tree2 <- tree2[["data"]]

#Merge data tables into a phyloseq object
merged_tables_16S <- phyloseq(otu_16S, metadata_16S, tree)#, taxa)

##===================================================================================================
## Read in ITS data
##===================================================================================================
Observed_otus_ITS <- read.table(file="ITS/observed-otus.tsv", header=TRUE)
Shannon_div_ITS <- read.table(file="ITS/shannon-diversity.tsv", header=TRUE)
PCoA_ITS <- read.table(file="ITS/pcoa_aitchison.tsv", header=TRUE)

#Read in metadata file and merge Obseverd OTUs and Shannon Diversity with the metadata
metadata_table_ITS <- read.table(file="ITS/metadata-ITS-LKM.tsv", header=TRUE) %>% 
  left_join(Observed_otus_ITS, by="SampleID") %>% 
  left_join(Shannon_div_ITS, by="SampleID") %>% 
  left_join(PCoA_ITS, by="SampleID") %>% 
  drop_na() %>% 
  as_tibble()

#Read in rarefied OTU table
otu_table_ITS <- read_qza("ITS/rarefied_table.qza")
otu_table_ITS <- otu_table_ITS[["data"]]
otu_table_ITS <- as_tibble(t(otu_table_ITS))

#Format data tables for phyloseq
otu_ITS <- otu_table(otu_table_ITS, taxa_are_rows=FALSE)
metadata_ITS <- sample_data(metadata_table_ITS)

#Merge data tables into a phyloseq object
merged_tables_ITS <- phyloseq(otu_ITS, metadata_ITS)#, taxa)

##===================================================================================================
## FIGURE 1 PCoA ORDIATION PLOTS: community differentiation between the three crop species
##===================================================================================================
#Perform the ordination
set.seed(2002)
pcoa_unifrac <- ordinate(merged_tables_16S, method="PCoA", distance="unifrac", weighted=TRUE)
set.seed(2002)
pcoa_bray <- ordinate(merged_tables_ITS, method="PCoA", distance="bray")

text_theme <-  theme(axis.text.x = element_text(size = 11), 
                     axis.title.x = element_text(size = 14), 
                     axis.text.y = element_text(size = 11), 
                     axis.title.y = element_text(size = 14),
                     text=element_text(size=16), 
                     plot.tag=element_text(size=18), 
                     legend.title = element_text(size = 14), 
                     legend.position="right")


#Plot PCoA -- (A) Bacteria/Archaea colored by crop species
(pcoa_plot_16S_crop <- plot_ordination(merged_tables_16S, pcoa_unifrac, type="samples", color="Crop") +
    theme_classic() +
    geom_point(size=2.5, show.legend=TRUE) +
    stat_ellipse(level=0.75) +
    scale_color_brewer(palette="Dark2", aesthetics ="color", name="Crop Species", labels=c("Lettuce", "Maize", "Mustard cabbage")) +
    labs(tag="A") + text_theme + theme(legend.position="none")
)

#Plot PCoA -- (B) Fungi colored by crop species
(pcoa_plot_ITS_crop <- plot_ordination(merged_tables_ITS, pcoa_bray, type="samples", color="Crop") +
    theme_classic() +
    geom_point(size=2.5, show.legend=FALSE) +
    stat_ellipse(level=0.75) +
    scale_color_brewer(palette="Dark2", aesthetics ="color", name="Crop Species", labels=c("Lettuce", "Maize", "Mustard cabbage")) +
    labs(tag="B") + text_theme
)

#Plot PCoA -- (C) Bacteria/Archaea colored by CropCycle (time)
(pcoa_plot_16S_cycle <- plot_ordination(merged_tables_16S, pcoa_unifrac, type="samples", color="CropCycle") +
    theme_classic() +
    geom_point(size=2.5, show.legend=FALSE) +
    stat_ellipse(level=0.75) +
    scale_color_manual(values=c("#A6BDDB", "#3690C0", "#225EA8"), name="Crop Cycle", labels=c("Cycle 1", "Cycle 2", "Cycle 3")) +
    labs(tag="C") + text_theme + theme(legend.position="none")
)

#Plot PCoA -- (D) Fungi colored by crop species
(pcoa_plot_ITS_cycle <- plot_ordination(merged_tables_ITS, pcoa_bray, type="samples", color="CropCycle") +
    theme_classic() +
    geom_point(size=2.5, show.legend=FALSE) +
    stat_ellipse(level=0.75) +
    scale_color_manual(values=c("#A6BDDB", "#3690C0", "#225EA8"), name="Crop Cycle", labels=c("Cycle 1", "Cycle 2", "Cycle 3")) +
    labs(tag="D") + text_theme
)


#Plot PCoA -- (E) Bacteria/Archaea colored by Soil
(pcoa_plot_16S_soil <- plot_ordination(merged_tables_16S, pcoa_unifrac, type="samples", color="Soil") +
    theme_classic() +
    geom_point(size=2.5, show.legend=TRUE) +
    stat_ellipse(level=0.75) +
    scale_color_manual(values=c("#1b9e77", "#db9511"), name="Soil", labels=c("Cultivated", "Uncultivated")) +
    labs(tag="E") + text_theme + theme(legend.position="none")
)

#Plot PCoA -- (F) Fungi colored by Soil
(pcoa_plot_ITS_soil <- plot_ordination(merged_tables_ITS, pcoa_bray, type="samples", color="Soil") +
    theme_classic() +
    geom_point(size=2.5, show.legend=FALSE)+
    stat_ellipse(level=0.75) +
    scale_color_manual(values=c("#1b9e77", "#db9511"), name="Soil", labels=c("Cultivated", "Uncultivated")) +
    labs(tag="F") + text_theme)

#Plot together
pcoa_plot_16S_crop + pcoa_plot_ITS_crop + pcoa_plot_16S_cycle + pcoa_plot_ITS_cycle + pcoa_plot_16S_soil + pcoa_plot_ITS_soil + plot_layout(ncol = 2)
#(pcoa_plot_16S_crop + pcoa_plot_ITS_crop) / (pcoa_plot_16S_cycle + pcoa_plot_ITS_cycle) / (pcoa_plot_16S_soil + pcoa_plot_ITS_soil)

#Save the plot
ggsave("Figure1-ordination.pdf", device="pdf", width=8.5, height=8.5)

##===================================================================================================
## FIGURE 1 STATISTICS: Community composition shifts with crop species, crop cycle, and soils
##===================================================================================================
#PERMANOVA comparisons among community members
#Make distance matrices
unifracdist_w <- phyloseq::distance(merged_tables_16S, method="unifrac", weighted=TRUE)
braydist <- phyloseq::distance(merged_tables_ITS, method="bray")

#Run PERMANOVA for 16S data
set.seed(2023)
adonis2(unifracdist_w ~ Soil, data=metadata_table2_16S, strata=metadata_table2_16S$CropCycle, permutations=999) #S, P=0.001, R2=0.035 | sig. between cultivated and uncultivated
set.seed(2023)
adonis2(unifracdist_w ~ Crop, data=metadata_table2_16S, strata=metadata_table2_16S$CropCycle, permutations=999) #S, P=0.001, R2=0.069 | sig. between crop species, across different cycles
set.seed(2023)
adonis2(unifracdist_w ~ CropCycle, data=metadata_table2_16S, permutations=999) #S, P=0.002, R2=0.050 | sig. between crop cycles taking into account cultivation

#Run pair-wise PERMANOVA across different crop species
set.seed(2023)
pairwise.adonis2(unifracdist_w ~ Crop, data=metadata_table2_16S, permutations=999)

#Interactions effect
set.seed(2023)
adonis2(unifracdist_w ~ Soil*Crop*CropCycle, data=metadata_table2_16S, permutations=999)
#Soil:Crop  P=0.001; significant interactions between crop cycles and cultivated vs. uncultivated soils
#Crop:Cycle P= 0.242; no significant effect between crop species and crop cycles

#Run PERMANOVA for ITS data
set.seed(2023)
adonis2(braydist ~ Soil, data=metadata_table_ITS, strata=metadata_table_ITS$Cycle, permutations=999) #S, P=0.001, R2=0.062 | sig. between cultivated and uncultivated
set.seed(2023)
adonis2(braydist ~ Crop, data=metadata_table_ITS, strata=metadata_table_ITS$Cycle, permutations=999) #NS, P=0.234, R2=0.034 | no sig. between crop species, across different cycles
set.seed(2023)
adonis2(braydist ~ CropCycle, data=metadata_table_ITS, permutations=999) #S, P=0.001, R2=0.372

#Interactions effect
set.seed(2023)
adonis2(braydist ~ Soil*Crop*Cycle, data=metadata_table_ITS, permutations=999)
#Soil:Crop  P=0.948; no significant interactions between crop cycles and cultivated vs. uncultivated soils
#Crop:Cycle P= 0.496; no significant effect between crop species and crop cycles

#Run pair-wise PERMANOVA across different crop species
set.seed(2023)
pairwise.adonis2(unifracdist_w ~ Crop, data=metadata_table2_16S, permutations=999)

##===================================================================================================
## FIGURE 2 BOXPLOTS: legacies of soils on richness and community differentiation over time
##===================================================================================================
#Prokaryotic richness based on Soil
(observed_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=observed_features, color=Soil)) +
  geom_boxplot(outlier.color="gray", linewidth=1) +
  theme_minimal() +
  scale_color_manual(values=c("#1b9e77", "#db9511")) +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  labs(x="", y = "Observed ASVs (16S)", tag="A") +
  theme(plot.tag = element_text(size=24)) +
  theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) + # aes(group = Day, colour =Day)) +
  guides(alpha = "none", color = "none")
)

#Prokaryotic PCoA based on Soil
(pcoa_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=Axis1, color=Soil)) +
  geom_boxplot(outlier.color="gray", linewidth=1) +
  theme_minimal() +
  scale_color_manual(values=c("#1b9e77", "#db9511")) +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  labs(x="", y = "Principal Coordinate \nAxis 1 (16S)", tag="B") +
  theme(plot.tag = element_text(size=24)) +
  theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) +# aes(group = Day, colour =Day)) +
  guides(alpha = "none", color = "none")
)

#Fungal richness based on Soil
(observed_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=observed_features, color=Soil)) +
  geom_boxplot(outlier.color="gray", linewidth=1) +
  theme_minimal() +
  scale_color_manual(values=c("#1b9e77", "#db9511")) +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  labs(x="", y = "Observed OTUs (ITS)", tag="C") +
  theme(plot.tag = element_text(size=24)) +
  theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6)+ # aes(group = Day, colour =Day)) +
  guides(alpha = "none", color = "none")
)

#Fungal PCoA based on Soil
(pcoa_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=Axis1, color=Soil)) +
  geom_boxplot(outlier.color="gray", linewidth=1) +
  theme_minimal() +
  scale_color_manual(values=c("#1b9e77", "#db9511")) +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  labs(x="", y = "Principal Coordinate \nAxis 1 (ITS)", tag="D") +
  theme(plot.tag = element_text(size=24)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=15)) +
  theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) # aes(group = Day, colour =Day))
)

#Plot together
observed_Soil_16S + pcoa_Soil_16S + observed_Soil_ITS + pcoa_Soil_ITS

#Save the plot
ggsave("Figure1-Soil-treatment-boxplot.pdf", device="pdf", width=14, height=9.5)

##===================================================================================================
## FIGURE 2 STATISTICS: comparisons between community parameters over time
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal == OK for ANOVA
shapiro.test(metadata_table_16S$observed_features)#NS
shapiro.test(metadata_table_16S$Axis1)#S
shapiro.test(metadata_table_ITS$observed_features)#NS
shapiro.test(metadata_table_ITS$Axis1)#NS

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(metadata_table_16S$observed_features ~ metadata_table_16S$CropCycle2)#NS
bartlett.test(metadata_table_16S$Axis1 ~ metadata_table_16S$CropCycle2)#NS
bartlett.test(metadata_table_ITS$observed_features ~ metadata_table_ITS$CropCycle2)#NS
bartlett.test(metadata_table_ITS$Axis1 ~ metadata_table_ITS$CropCycle2)#NS

# (A) Observed ASVs vs Cycle
# ANOVA assumptions OK
anova_observed_asv <- aov(metadata_table_16S$observed_features ~ metadata_table_16S$CropCycle2)
summary(anova_observed_asv)#S

# Post-hoc test
tukey_observed_asv <- TukeyHSD(anova_observed_asv, ordered="TRUE")
tukey_observed_asv

# Extract letters associated with significance with the "multcompView" package
observed_letters_asv <- multcompLetters4(anova_observed_asv, tukey_observed_asv)
as.data.frame.list(observed_letters_asv$`metadata_table_16S$CropCycle2`)

# (B) PCoA Axis 1 vs Cycle
# ANOVA assumptions failed; using non-parametric test
library(agricolae)
kruskal(metadata_table_16S$Axis1, metadata_table_16S$CropCycle2, group=TRUE, p.adj="bonferroni")$groups

# (C) Observed OTUs vs Cycle
# ANOVA assumptions OK
anova_observed_otu <- aov(metadata_table_ITS$observed_features ~ metadata_table_ITS$CropCycle2)
summary(anova_observed_otu)#S
tukey_observed_otu <- TukeyHSD(anova_observed_otu, ordered="TRUE")
tukey_observed_otu

# Extract letters associated with significance with the "multcompView" package
observed_letters_otu <- multcompLetters4(anova_observed_otu, tukey_observed_otu)
as.data.frame.list(observed_letters_otu$`metadata_table_ITS$CropCycle2`)

# (D) PCoA Axis 1 vs Cycle
# ANOVA assumptions OK; but using non-parametric test to be comparable with (B)
kruskal(metadata_table_ITS$Axis1, metadata_table_ITS$CropCycle2, group=TRUE, p.adj="bonferroni")$groups

##===================================================================================================
## SUPPLEMENTAL FIGURE 1 NUTRIENT CONTENT: Observed features and principal coordinates over time
##===================================================================================================

## Read in soil nutrient data
nutrient_data <- read.table(file="soil_nutrients/soil_nutrient_data.txt", header=TRUE, sep="\t")

## Draw box & whisker plots
(pH_plot <- ggplot(nutrient_data, aes(x=CropCycle2, y=pH, color=Soil)) +
    geom_boxplot(outlier.color="gray", linewidth=1) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "pH", tag="A") +
    theme(plot.tag = element_text(size=24)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) +# aes(group = Day, colour =Day)) +
    guides(alpha = "none", color = "none")
)

(TotalC_plot <- ggplot(nutrient_data, aes(x=CropCycle2, y=PercentC, color=Soil)) +
    geom_boxplot(outlier.color="gray", linewidth=1) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Total C (%)", tag="B") +
    theme(plot.tag = element_text(size=24)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) +# aes(group = Day, colour =Day)) +
    guides(alpha = "none", color = "none")
)

(TotalN_plot <- ggplot(nutrient_data, aes(x=CropCycle2, y=PercentN, color=Soil)) +
    geom_boxplot(outlier.color="gray", linewidth=1) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Total N (%)", tag="C") +
    theme(plot.tag = element_text(size=24)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) +# aes(group = Day, colour =Day)) +
    guides(alpha = "none", color = "none")
)

#Plot together
pH_plot/TotalC_plot/TotalN_plot

#Save the plot
ggsave("FigureS1-nutrients-over-time.pdf", device="pdf", width=6.5, height=8)

##===================================================================================================
## SUPPLEMENTAL FIGURE 1 STATISTICS: comparisons of soil nutrient content
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal == OK for ANOVA
shapiro.test(nutrient_data$pH)#NS
shapiro.test(nutrient_data$PercentC)#NS
shapiro.test(nutrient_data$PercentN)#S

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(nutrient_data$pH ~ nutrient_data$CropCycle2)#NS
bartlett.test(nutrient_data$PercentC ~ nutrient_data$CropCycle2)#NS
bartlett.test(nutrient_data$PercentN ~ nutrient_data$CropCycle2)#NS

# (A) pH vs Cycle
# ANOVA assumptions OK
anova_pH <- aov(nutrient_data$pH ~ nutrient_data$CropCycle2)
summary(anova_pH)#NS

# (B) Percent C vs Cycle
# ANOVA assumptions OK
anova_C <- aov(nutrient_data$PercentC ~ nutrient_data$CropCycle2)
summary(anova_C)#NS

# (C) Percent N vs Cycle
# ANOVA assumptions failed; using non-parametric test
library(agricolae)
kruskal(nutrient_data$PercentN, nutrient_data$CropCycle2, group=TRUE, p.adj="bonferroni")$groups


##===================================================================================================
## SUPPLEMENTAL FIGURE 2 DOTPLOT REGRESSION: Observed features and principal coordinates across CROP SPECIES
##===================================================================================================
#Prokaryotic Richness based on Crop species
(observed_Crop_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=observed_features, color=Crop, group=Crop)) +
   geom_point(alpha=0.8, size=4) +
   geom_smooth(method=lm, alpha=0.2) +
   theme_minimal() +
   scale_color_brewer(palette="Dark2") +
   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
   labs(x="", y = "Observed ASVs (16S)", tag="A") +
   theme(plot.tag = element_text(size=20)) +
   theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
         axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
   guides(alpha = "none", color = "none") 
)

#Prokaryotic PCoA based on Crop species
(pcoa_Crop_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=Axis1, color=Crop, group=Crop)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_brewer(palette="Dark2") +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 1 (16S)", tag="B") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    guides(alpha = "none", color = "none")
)

#Fungal Richness based on Crop species
(observed_Crop_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=observed_features, color=Crop, group=Crop)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_brewer(palette="Dark2") +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Observed OTUs (ITS)", tag="C") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    guides(alpha = "none", color = "none")
)

#Fungal PCoA based on Crop species
(pcoa_Crop_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=Axis1, color=Crop, group=Crop)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_brewer(palette="Dark2") +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 1 (ITS)", tag="D") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    theme(legend.position="bottom", text=element_text(size=15)) +
    theme(legend.title = element_blank())
)

#Plot together
observed_Crop_16S + pcoa_Crop_16S + observed_Crop_ITS + pcoa_Crop_ITS

#Save the plot
ggsave("FigureS2-Crop-regresssion-2.pdf", device="pdf", width=14, height=9)




##===================================================================================================
## SUPPLEMENTAL FIGURE 2 STATISTICS: Linear effects of community parameters and CROP SPECIES
##===================================================================================================
#Linear mixed-effect models; accounting for random effects
#Testing for effects between different plant species over time
#(A)
lme(observed_features ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #NS; the model predicts no significance among the three crop species
#(B)
lme(Axis1 ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #NS; the model predicts no significance among the three crop species
#(C)
lme(observed_features ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance among the three crop species
#(D)
lme(Axis1 ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance among the three crop species


#Testing Axis 2
lme(Axis2 ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #NS; the model predicts no significance among the three crop species
lme(Axis2 ~ Crop, random = ~1 | Soil/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance among the three crop species


##===================================================================================================
## SUPPLEMENTAL FIGURE 3 DOTPLOT REGRESSION: Observed features and principal coordinates across SOILS
## Slope and intercepts for both soils together is found in the analyses in the next section
##===================================================================================================
#(A) Prokaryotic richness based on Soil
(observed_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=observed_features, color=Soil, group=Soil)) +
   geom_point(alpha=0.8, size=4) +
   geom_smooth(method=lm, alpha=0.2) +
   theme_minimal() +
   scale_color_manual(values=c("#1b9e77", "#db9511")) +
   geom_abline(intercept = 332.432, slope = -8.5814, color="#807DBA", linetype="dashed", linewidth=1.5) +
   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
   labs(x="", y = "Observed ASVs (16S)", tag="A") +
   theme(plot.tag = element_text(size=20)) +
   theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
         axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
   guides(alpha = "none", color = "none")
)
  
#(B) Prokaryotic PCoA based on Soil
(pcoa_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=Axis1, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    geom_abline(intercept=-8.925672, slope = 2.358429, color="#807DBA", linetype="dashed", linewidth=1.5) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 1 (16S)", tag="B") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  guides(alpha = "none", color = "none")
)

#(C) Fungal richness based on Soil
(observed_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=observed_features, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    geom_abline(intercept = 211.55000, slope = -11.23333, color="#807DBA", linetype="dashed", linewidth=1.5) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Observed OTUs (ITS)", tag="C") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    guides(alpha = "none", color = "none")
)

#(D) Fungal PCoA based on Soil
(pcoa_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=Axis1, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    geom_abline(intercept = -10.451757, slope = 2.850479, color="#807DBA", linetype="dashed", linewidth=1.5) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 1 (ITS)", tag="D") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  theme(legend.position="bottom") +
    theme(legend.title = element_blank())
)

#Plot together
observed_Soil_16S + pcoa_Soil_16S + observed_Soil_ITS + pcoa_Soil_ITS

#Save the plot
ggsave("Figure2-Soil-treatment-dotplot.pdf", device="pdf", width=14, height=7)


##===================================================================================================
## SUPPLEMENTAL FIGURE 3 STATISTICS: Linear effects of community parameters and SOILS
##===================================================================================================
#Split the dataset to be able to analyze them separately
metadata_table_16S_c <- filter(metadata_table_16S, Soil == "Cultivated")
metadata_table_16S_uc <- filter(metadata_table_16S, Soil == "Uncultivated")

metadata_table_ITS_c <- filter(metadata_table_ITS, Soil == "Cultivated")
metadata_table_ITS_uc <- filter(metadata_table_ITS, Soil == "Uncultivated")

#Linear mixed-effect models; accounting for random effects
#Testing for effects of cultivated over time; using numeric time value
#(A)
lme(observed_features ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_c) %>% 
  summary() #S; the model predicts significance across the cycle
lme(observed_features ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_uc) %>% 
  summary() #NS; the model predicts no significance across the cycle; this checks out

#(B)
lme(Axis1 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_c) %>% 
  summary() #S; the model predicts significance across the cycle
lme(Axis1 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_uc) %>% 
  summary() #S; the model predicts significance across the cycle

#(C)
lme(observed_features ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_c) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
lme(observed_features ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_uc) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils

#(D)
lme(Axis1 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_c) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
lme(Axis1 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_uc) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils


#Testing Axis 2
lme(Axis2 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_c) %>% 
  summary() #NS; the model predicts no significance across the cycle
lme(Axis2 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_16S_uc) %>%
  summary() #S; the model predicts no significance across the cycle

lme(Axis2 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_c) %>%
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
lme(Axis2 ~ Cycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS_uc) %>%
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils


#Linear mixed-effect models; accounting for random effects
#Testing for significance between cultivated and uncultivated soils => these results reflect the results from the models above.
#(A)
lme(observed_features ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; the model predicts significance between cultivated vs uncultivated soils
#(B)
lme(Axis1 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
lme(Axis2 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_16S) %>%
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils

#(C)
lme(observed_features ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
#(D)
lme(Axis1 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
lme(Axis2 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_ITS) %>%
  summary() #S; the model predicts significance between cultivated vs uncultivated soils

#Linear mixed-effect models; accounting for random effects
#Testing for the effects of ALL SOILS together; regression slope and intercept used in building Supplemental Figure 2.
#(A)
lme(observed_features ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; P = 0.0354; the model predicts significance across the cycle
  #intercept = 332.432, slope = -8.5814

#(B)
lme(Axis1 ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; P < 0.001; the model predicts strong significance across the cycle
  #intercept = -8.925672, slope = 2.358429
# lme(Axis2 ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>%
#   summary() #S; P=0.1; the model predicts no significance across the cycle

#(C)
lme(observed_features ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #S; P < 0.001; the model predicts strong significance across the cycle
  #intercept = 211.55000, slope = -11.23333

#(D)
lme(Axis1 ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #S; P < 0.001; the model predicts strong significance across the cycle
  #intercept = -10.451757, slope = 2.850479
# lme(Axis2 ~ TimePoint, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>%
#   summary() #S; P = 0.0378; the model predicts weak significance across the cycle



#Linear mixed-effect models; accounting for random effects
#Testing for significance different CYCLES (1,2,3)
#(A)
lme(observed_features ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; the model predicts significance among the 3 cycles, (P<0.034)

#(B)
lme(Axis1 ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; the model predicts significance among the 3 cycles, (P<0.001)
lme(Axis2 ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_16S) %>%
  summary() #NS; the model predicts no significance among the 3 cycles, (P>0.072)

#(C)
lme(observed_features ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #S; the model predicts significance among the 3 cycles, (P<0.014)

#(D)
lme(Axis1 ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #S; the model predicts significance among the 3 cycles, (P<0.001)
lme(Axis2 ~ CropCycle, random = ~1 | Replication/StudyID, data=metadata_table_ITS) %>%
  summary() #S; the model predicts significance with cycles 1 & 2 cycles (P<0.0475) but not significant with cycle 3 (P=0.217)


##===================================================================================================
## SUPPLEMENTAL FIGURE 4 BOXPLOT (using PCoA Axis 2): legacies of soils on richness and community differentiation over time
##===================================================================================================
#Prokaryotic PCoA based on Soil (Axis2)
(pcoa_Soil_16S2 <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=Axis2, color=Soil)) +
   geom_boxplot(outlier.color="gray", linewidth=1) +
   theme_minimal() +
   scale_color_manual(values=c("#1b9e77", "#db9511")) +
   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
   labs(x="", y = "Principal Coordinate \nAxis 2 (16S)", tag="A") +
   theme(plot.tag = element_text(size=24)) +
   theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
         axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
   stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) +# aes(group = Day, colour =Day)) +
   guides(alpha = "none", color = "none")
)

#Fungal PCoA based on Soil
(pcoa_Soil_ITS2 <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=Axis2, color=Soil)) +
    geom_boxplot(outlier.color="gray", linewidth=1) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 2 (ITS)", tag="B") +
    theme(plot.tag = element_text(size=24)) +
    theme(legend.position="bottom") +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=15)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    stat_summary(fun=median, geom='line', aes(group=1), color="gray40", linetype=2, linewidth=0.6) # aes(group = Day, colour =Day))
)

pcoa_Soil_16S2 + pcoa_Soil_ITS2

#Save the plot
ggsave("FigureS3-Soil-treatment-boxplot-Axis2.pdf", device="pdf", width=14, height=4.25)

##===================================================================================================
## SUPPLEMENTAL FIGURE 4 STATISTICS: comparisons between community parameters over time for Axis 2
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal == OK for ANOVA
shapiro.test(metadata_table_16S$Axis2)#NS
shapiro.test(metadata_table_ITS$Axis2)#NS

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(metadata_table_16S$Axis2 ~ metadata_table_16S$CropCycle2)#NS
bartlett.test(metadata_table_ITS$Axis2 ~ metadata_table_ITS$CropCycle2)#NS

# (A) PCoA Axis 2 vs Cycle for prokaryotes
# ANOVA assumptions OK
anova_pcoa_16S2 <- aov(metadata_table_16S$Axis2 ~ metadata_table_16S$CropCycle2)
summary(anova_pcoa_16S2)#NS

# (B) PCoA Axis 2 vs Cycle for fungi
# ANOVA assumptions OK
anova_pcoa_ITS2 <- aov(metadata_table_ITS$Axis2 ~ metadata_table_ITS$CropCycle2)
summary(anova_pcoa_ITS2)#S
tukey_pcoa_ITS2 <- TukeyHSD(anova_pcoa_ITS2, ordered="TRUE")
tukey_pcoa_ITS2

# Extract letters associated with significance with the "multcompView" package
observed_letters_pcoa <- multcompLetters4(anova_pcoa_ITS2, tukey_pcoa_ITS2)
as.data.frame.list(observed_letters_pcoa$`metadata_table_ITS$CropCycle2`)








# ##===================================================================================================
# ## ADDITIONAL STATISTICS NOT USED IN MANUSCRIPT
# ##===================================================================================================
# ##===================================================================================================
# ## SUPPLEMENTAL FIGURE 1 STATISTICS: Linear effects of community parameters and SOILS (not used in manuscript)
# ##===================================================================================================
# #Split the dataset to be able to analyze them separately
# nutrient_data_c <- filter(nutrient_data, Soil == "Cultivated")
# nutrient_data_uc <- filter(nutrient_data, Soil == "Uncultivated")
# 
# #Linear mixed-effect models; accounting for random effects
# #Testing for effects of cultivated over time; using numeric time value
# library(nlme)
# #(A)
# lme(pH ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_c) %>% 
#   summary() #NS; the model predicts no significance across the cycle
# lme(pH ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_uc) %>% 
#   summary() #NS; the model predicts no significance across the cycle; this checks out
# 
# #(B)
# lme(PercentC ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_c) %>% 
#   summary() #NS; the model predicts no significance across the cycle
# lme(PercentC ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_uc) %>% 
#   summary() #NS; the model predicts no significance across the cycle
# 
# #(C)
# lme(PercentN ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_c) %>% 
#   summary() #NS; the model predicts no significance across the cycle
# lme(PercentN ~ Cycle, random = ~1 | Replication/Sample_name, data=nutrient_data_uc) %>% 
#   summary() #NS; the model predicts no significance across the cycle
# 
# 
# #Testing for the effects of all soils together
# #(A) pH
# lme(pH ~ TimePoint, random = ~1 | Replication/Sample_name, data=nutrient_data) %>% 
#   summary() #NS; P = 0.501; the model predicts no significance across the cycle
# 
# #(B) Total C
# lme(PercentC ~ TimePoint, random = ~1 | Replication/Sample_name, data=nutrient_data) %>% 
#   summary() #NS; P = 0.577 ; the model predicts no significance across the cycle
# 
# #(C) Total N
# lme(PercentN ~ TimePoint, random = ~1 | Replication/Sample_name, data=nutrient_data) %>% 
#   summary() #NS; P = 0.220; the model predicts no significance across the cycle


##===================================================================================================
#Linear fixed-effect models (not used in manuscript)
##===================================================================================================
#Testing for effects between cultivated and uncultivated soils
#(A)
# lm(observed_features ~ Soil, data=metadata_table_16S) %>% 
#   summary() #S; the model predicts significance between cultivated vs uncultivated soils
# #(B)
# lm(Axis1 ~ Soil, data=metadata_table_16S) %>% 
#   summary() #S; the model predicts significance between cultivated vs uncultivated soils
# #(C)
# lm(observed_features ~ Soil, data=metadata_table_ITS) %>% 
#   summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
# #(D)
# lm(Axis1 ~ Soil, data=metadata_table_ITS) %>% 
#   summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
# 
# #Testing for effects between different crop species
# #(A)
# lm(observed_features ~ Crop, data=metadata_table_16S) %>% 
#   summary() #NS; the model predicts no significance among the three crop species
# #(B)
# lm(Axis1 ~ Crop, data=metadata_table_16S) %>% 
#   summary() #NS; the model predicts no significance among the three crop species
# #(C)
# lm(observed_features ~ Crop, data=metadata_table_ITS) %>% 
#   summary() #NS; the model predicts no significance among the three crop species
# #(D)
# lm(Axis1 ~ Crop, data=metadata_table_ITS) %>% 
#   summary() #NS; the model predicts no significance among the three crop species

##===================================================================================================
## Dotplot grouped by crop cycle (not used in manuscript)
##===================================================================================================
# #Richness based on Crop species
# ggplot(metadata_table, aes(x=CropCycle, y=observed_features, color=Crop, group=Crop)) +
#   geom_point(alpha=0.8, size=4) +
#   geom_smooth(method=lm, alpha=0.2) +
#   theme_minimal() +
#   scale_color_brewer(palette="Dark2") +
#   scale_x_discrete(limit=c("C1", "C2", "C3"),
#                    labels = c("Cycle 1", "Cycle 2", "Cycle 3")) +
#   labs(x="", y = "Observed ASVs (16S)", tag="A") +
#   theme(plot.tag = element_text(size=20)) +
#   theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))

# #PC based on Crop species (divergence of community)
# ggplot(metadata_table, aes(x=CropCycle, y=Axis1, color=Crop, group=Crop)) +
#   geom_point(alpha=0.8, size=4) +
#   geom_smooth(method=lm, alpha=0.2) +
#   theme_minimal() +
#   scale_color_brewer(palette="Dark2") +
#   scale_x_discrete(limit=c("C1", "C2", "C3"), 
#                    labels = c("Cycle 1", "Cycle 2", "Cycle 3")) +
#   labs(x="", y = "Principal Coordinate Axis 1 (16S)", tag="B") +
#   theme(plot.tag = element_text(size=20)) +
#   theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))

#PC based on Crop species (divergence of community)
# ggplot(metadata_table, aes(x=CropCycle, y=Axis1, color=Soil, group=Soil)) +
#   geom_point(alpha=0.8, size=4) +
#   geom_smooth(method=lm, alpha=0.2) +
#   theme_minimal() +
#   scale_color_manual(values=c("#66a61e", "#a6761d")) +
#   scale_x_discrete(limit=c("C1", "C2", "C3"), 
#                    labels = c("Cycle 1", "Cycle 2", "Cycle 3")) +
#   labs(x="", y = "Principal Coordinate Axis 1 (16S)", tag="D") +
#   theme(plot.tag = element_text(size=20)) +
#   theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), 
#         axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
