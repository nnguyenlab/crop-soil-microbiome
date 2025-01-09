setwd("/Volumes/GoogleDrive/Shared drives/Crop soil microbiome/Manuscripts/Manuscript 1 - crop species/Env Microbiology submission/EMI 2nd submission/R commands/Nhu's_analyses")
library(metafor)

#geneColors=c("16S"="#d95f02", "ITS"="#1c6d77")
geneColors=c("16S"="#7570b3", "ITS"="#d95f02")

##########################################
#Comparing observed features
##########################################
observed_features_stats <- read.csv("response_ratio/observed_features_stats.csv", header=TRUE)

LRR_effect_sizes_richness <- escalc(data=observed_features_stats, "ROM",
                          m1i=mean_baseline, 
                          sd1i=SD_baseline, 
                          n1i=sample_size_baseline,
                          m2i=mean_treatment,
                          sd2i=SD_treatment,
                          n2i=sample_size_treatment)

LRR_richness <- ggplot(LRR_effect_sizes_richness, aes(x=CropCycle2, y=yi, color=gene, group=gene)) +
  geom_point(size=5) +
  geom_line(linetype=1, linewidth=0.3, alpha=0.6) +
  geom_errorbar(aes(ymin=yi-vi, ymax=yi+vi), width=.1, color="white") +
  theme_minimal() +
  theme(legend.title= element_blank()) +
  theme(legend.justification=c(-0.25,1.2), legend.position=c(0,1)) +
  #theme(legend.justification=c(1.15,-0.25), legend.position=c(1,0)) +
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="gray")) +
  theme(plot.tag = element_text(size=20)) +
  scale_color_manual(values=geneColors, name="Organism", labels=c("Bacteria & Archaea","Fungi")) +
  labs(x="", y="Log Response Ratio (richness)", tag="A") +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c", "C3uc", "C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  guides(color=guide_legend(override.aes = list(size = 3)))
  

##########################################
#Comparing Aitchison Distances
##########################################
aitchison_stats <- read.csv("response_ratio/aitchison_stats.csv", header=TRUE)

LRR_effect_sizes_aitchison <- escalc(data=aitchison_stats, "ROM",
                                   m1i=mean_baseline, 
                                   sd1i=SD_baseline, 
                                   n1i=sample_size_baseline,
                                   m2i=mean_treatment,
                                   sd2i=SD_treatment,
                                   n2i=sample_size_treatment)

LRR_aitchison <- ggplot(LRR_effect_sizes_aitchison, aes(x=CropCycle2, y=yi, color=gene, group=gene)) +
  geom_point(size=5) +
  geom_line(linetype=1, linewidth=0.3, alpha=0.6) +
  geom_errorbar(aes(ymin=yi-vi, ymax=yi+vi), width=.1, color="white") +
  theme_minimal() +
  theme(plot.tag = element_text(size=20)) +
#  theme(legend.justification=c(1.15,-0.25), legend.position=c(1,0)) +
#  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="gray")) +
  scale_color_manual(values=geneColors, name="Organism", labels=c("Bacteria & Archaea","Fungi")) +
  labs(x="", y="Log Response Ratio (community distance)", tag="B") +
  scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c", "C3uc", "C3c"), 
                   labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
  guides(alpha="none", color="none")

#Plot together
LRR_richness + LRR_aitchison

#Save the plot
ggsave("Figure3-response.pdf", device="pdf", width=10, height=4)






#Tried with aitchison axis 1 stats but due to negative values, LLR cannot be calculated
# aitchison_axis1_stats <- read.csv("aitchison_axis1_stats.csv", header=TRUE)
# LRR_effect_sizes_composition <- escalc(data=aitchison_axis1_stats, "ROM",
#                            m1i=mean_baseline, 
#                            sd1i=SD_baseline, 
#                            n1i=sample_size_baseline,
#                            m2i=mean_treatment,
#                            sd2i=SD_treatment,
#                            n2i=sample_size_treatment)
# 
# LRR_composition <- ggplot(LRR_effect_sizes_composition, aes(x=CropCycle2, y=yi, color=gene, group=gene)) +
#   geom_point(size=5) +
#   geom_line(linetype=1, linewidth=0.3, alpha=0.6) +
#   #geom_smooth(method=lm, alpha=0.2, se=FALSE) +
#   geom_errorbar(aes(ymin=yi-vi, ymax=yi+vi), width=.1, color="white") +
#   theme_minimal() +
#   scale_color_manual(values=geneColors, name="Organism", labels=c("Bacteria & Archaea","Fungi")) +
#   labs(x="", y="Log Likelihood Ratio (community composition)", tag="B") +
#   #theme(legend.justification=c(1,0), legend.position=c(1,0)) +
#   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c", "C3uc", "C3c"), 
#                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
#   guides(alpha="none", color="none")


##########################################
#Comparing Shannon Diversity -- not used in manuscript
##########################################
# shannon_stats <- read.csv("response_ratio/shannon_stats.csv", header=TRUE)
# 
# LRR_effect_sizes_shannon <- escalc(data=shannon_stats, "ROM",
#                                    m1i=mean_baseline, 
#                                    sd1i=SD_baseline, 
#                                    n1i=sample_size_baseline,
#                                    m2i=mean_treatment,
#                                    sd2i=SD_treatment,
#                                    n2i=sample_size_treatment)
# 
# LRR_shannon <- ggplot(LRR_effect_sizes_shannon, aes(x=CropCycle2, y=yi, color=gene, group=gene)) +
#   geom_point(size=5) +
#   geom_line(linetype=1, linewidth=0.3, alpha=0.6) +
#   #geom_smooth(method=lm, alpha=0.2, se=FALSE) +
#   geom_errorbar(aes(ymin=yi-vi, ymax=yi+vi), width=.1, color="gray") +
#   theme_minimal() +
#   theme(plot.tag = element_text(size=20)) +
#   theme(legend.justification=c(1.15,-0.25), legend.position=c(1,0)) +
#   theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="gray")) +
#   scale_color_manual(values=geneColors, name="Organism", labels=c("Bacteria & Archaea","Fungi")) +
#   labs(x="", y="Log Response Ratio (Shannon diversity)", tag="B") +
#   scale_x_discrete(limit=c("C1uc", "C1c", "C2uc", "C2c", "C3uc", "C3c"), 
#                    labels = c("Cycle 1 \nuncultivated", "Cycle 1 \ncultivated","Cycle 2 \nuncultivated","Cycle 2 \ncultivated", "Cycle 3 \nuncultivated", "Cycle 3 \ncultivated")) +
#   guides(color=guide_legend(override.aes = list(size = 3)))  
# 
