setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/Shared drives/Crop soil microbiome/Manuscripts/Manuscript 1 - crop species/Env Microbiology submission/EMI 2nd submission/R commands/Nhu's_analyses")
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

#display.brewer.all()

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

#write.table(metadata_table_16S, "16S/metadata_table_16S.tsv", sep="\t")

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
metadata2_16S <- sample_data(metadata_table2_16S)
tree <- read_tree("16S/rooted-tree.nwk")
 
# #tree can also be imported directly from QIIME
# tree2 <- read_qza("16S/rooted_tree.qza")
# tree2 <- tree2[["data"]]

#Merge data tables into a phyloseq object
merged_tables_16S <- phyloseq(otu_16S, metadata_16S, tree)#, taxa)
merged_tables2_16S <- phyloseq(otu_16S, metadata2_16S, tree)#, taxa)

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
## Build Figure 1 BOXPLOT - Observed features and principal coordinates across Cycles
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

#Prokaryotic PCoA based on Soil
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
## Build Figure 2 DOTPLOT REGRESSION - Observed features and principal coordinates across cycles
##===================================================================================================
#Prokaryotic PCoA based on Soil
(observed_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=observed_features, color=Soil, group=Soil)) +
   geom_point(alpha=0.8, size=4) +
   geom_smooth(method=lm, alpha=0.2) +
   theme_minimal() +
   scale_color_manual(values=c("#1b9e77", "#db9511")) +
   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
   labs(x="", y = "Observed ASVs (16S)", tag="A") +
   theme(plot.tag = element_text(size=20)) +
   theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
         axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
   guides(alpha = "none", color = "none")
)

#Prokaryotic PCoA based on Soil
(pcoa_Soil_16S <- ggplot(metadata_table_16S, aes(x=CropCycle2, y=Axis1, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Principal Coordinate \nAxis 1 (16S)", tag="B") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
  guides(alpha = "none", color = "none")
)

#Fungal PCoA based on Soil
(observed_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=observed_features, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
    scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c","C3uc","C3c"), 
                     labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
    labs(x="", y = "Observed OTUs (ITS)", tag="C") +
    theme(plot.tag = element_text(size=20)) +
    theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14)) +
    guides(alpha = "none", color = "none")
)

#Prokaryotic PCoA based on Soil
(pcoa_Soil_ITS <- ggplot(metadata_table_ITS, aes(x=CropCycle2, y=Axis1, color=Soil, group=Soil)) +
    geom_point(alpha=0.8, size=4) +
    geom_smooth(method=lm, alpha=0.2) +
    theme_minimal() +
    scale_color_manual(values=c("#1b9e77", "#db9511")) +
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
## Build Figure 3 PCoA PLOT between the three crop species
##===================================================================================================
#Perform the ordination
set.seed(2002)
pcoa_unifrac <- ordinate(merged_tables2_16S, method="PCoA", distance="unifrac", weighted=TRUE)
set.seed(2002)
pcoa_bray <- ordinate(merged_tables_ITS, method="PCoA", distance="bray")

#Plot PCoA -- (A) Bacteria/Archaea
(pcoa_plot_16S <- plot_ordination(merged_tables2_16S, pcoa_unifrac, type="samples", color="Crop") +
  theme_classic() +
  geom_point(size=3, show.legend=TRUE) +
  stat_ellipse(level=0.75) +
  #scale_color_manual(values = c("lawn"="darkorange", "soil" = "#7a81ff"), name="Soil", labels=c("Lawn", "Soil")) +
  scale_color_brewer(palette="Dark2", aesthetics ="color", name="", labels=c("Lettuce", "Maize", "Mustard cabbage")) +
    labs(tag="A") +
    theme(axis.text.x = element_text(size = 11), 
          axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), 
          axis.title.y = element_text(size = 14),
          text=element_text(size=17), 
          plot.tag=element_text(size=24), 
          legend.position="bottom")
)

(pcoa_plot_ITS <- plot_ordination(merged_tables_ITS, pcoa_bray, type="samples", color="Crop") +
    theme_classic() +
    geom_point(size=3, show.legend=FALSE) +
    stat_ellipse(level=0.75) +
    #scale_color_manual(values = c("lawn"="darkorange", "soil" = "#7a81ff"), name="Soil", labels=c("Lawn", "Soil")) +
    scale_color_brewer(palette="Dark2", aesthetics ="color", name="", labels=c("Lettuce", "Maize", "Mustard cabbage")) +
    labs(tag="B") +
    theme(axis.text.x = element_text(size = 11), 
          axis.title.x = element_text(size = 14), 
          axis.text.y = element_text(size = 11), 
          axis.title.y = element_text(size = 14),
          text=element_text(size=17), 
          plot.tag=element_text(size=24), 
          legend.position="bottom") +
  guides(alpha = "none", color = "none")
)


#Plot together
pcoa_plot_16S + pcoa_plot_ITS

#Save the plot
ggsave("Figure3-CropSpecies-ordination.pdf", device="pdf", width=10, height=5.5)


##===================================================================================================
## Build Supplemental Figure 1 DOTPLOT REGRESSION- Observed features and principal coordinates across CROP SPECIES
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
ggsave("FigureS1-Crop-regresssion.pdf", device="pdf", width=14, height=7)


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

  
##===================================================================================================
##  FIG. 1 STATISTICS: COMPARISONS BETWEEN COMMUNITY PARAMETERS OVER TIME
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal = OK for ANOVA
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
##  FIG. 2 STATISTICS: LINEAR EFFECTS OF COMMUNITY PARAMETERS OVER TIME
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

#----------
#Linear mixed-effect models; accounting for random effects
#Testing for significance between cultivated and uncultivated soils
#(A)
lme(observed_features ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; the model predicts significance between cultivated vs uncultivated soils
#(B)
lme(Axis1 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_16S) %>% 
  summary() #S; the model predicts significance between cultivated vs uncultivated soils
#(C)
lme(observed_features ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils
#(D)
lme(Axis1 ~ Soil, random = ~1 | Crop/Replication/StudyID, data=metadata_table_ITS) %>% 
  summary() #NS; the model predicts no significance between cultivated vs uncultivated soils

##===================================================================================================
##  FIG. 3 STATISTICS: COMPARISONS BETWEEN COMMUNITY PARAMETERS AND PLANT SPECIES
##===================================================================================================
#PERMANOVA comparisons among community members
#Make distance matrices
unifracdist_w <- phyloseq::distance(merged_tables2_16S, method="unifrac", weighted=TRUE)
braydist <- phyloseq::distance(merged_tables_ITS, method="bray")


#Run PERMANOVA for 16S data
set.seed(2023)
summary(adonis2(unifracdist_w ~ Soil, data=metadata_table2_16S, permutations=999)) #S, P=0.003, R2=0.035
set.seed(2023)
adonis2(unifracdist_w ~ Crop, data=metadata_table2_16S, strata=metadata_table2_16S$Cycle, permutations=999) #NS, P=0.003, R2=0.069
set.seed(2023)
adonis2(unifracdist_w ~ CropCycle2, data=metadata_table2_16S, permutations=999) #S, P=0.001, R2=0.112

#Run pair-wise PERMANOVA
set.seed(2002)
pairwise.adonis2(unifracdist_w ~ Crop, data=metadata_table2_16S, permutations=999)

# adonis2(unifracdist_w ~ CropSoil, data=metadata_table_ITS, permutations=999) #S, R2=0.121
# adonis2(unifracdist_w ~ CycleSoil, data=metadata_table_ITS, permutations=999) #S, R2=0.221
# adonis2(unifracdist_w ~ CropSoil, strata=metadata_table_ITS$Cycle, data=metadata_table_ITS, permutations=999) #S, R2=0.121
# adonis2(unifracdist_w ~ Soil*Crop*CropCycle2, data=metadata_table_ITS, permutations=999)

#Run PERMANOVA for ITS data
set.seed(2023)
adonis2(braydist ~ Crop, data=metadata_table_ITS, strata=metadata_table_ITS$Cycle, permutations=999) #NS, P=0.234, R2=0.034


# set.seed(2023)
# adonis2(braydist ~ Soil, data=metadata_table_ITS, permutations=999) #S, P=0.001, R2=0.062
# set.seed(2023)
# adonis2(braydist ~ CropCycle2, data=metadata_table_ITS, permutations=999) #S, P=0.001, R2=0.372
# adonis2(braydist ~ CropSoil, data=metadata_table_ITS, permutations=999) #S, R2=0.107
# adonis2(braydist ~ CycleSoil, data=metadata_table_ITS, permutations=999) #S, R2=0.372
# adonis2(braydist ~ Crop, strata=metadata_table_ITS$Cycle, data=metadata_table_ITS, permutations=999) #S, R2=0.107
# adonis2(braydist ~ Soil*Crop*Cycle, data=metadata_table_ITS, permutations=999)


##===================================================================================================
##  SUPPLEMENTAL FIGURE 1 STATISTICS: COMPARISONS BETWEEN COMMUNITY PARAMETERS AND PLANT SPECIES
##===================================================================================================
#Linear mixed-effect models; accounting for random effects
#Testing for effects between different plant species
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

#Linear fixed-effect models (not used in MS)
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




