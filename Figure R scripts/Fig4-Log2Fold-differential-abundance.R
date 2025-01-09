#load libraries
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
theme_set(theme_bw())

setwd("/Volumes/GoogleDrive/Shared drives/Crop soil microbiome/Manuscripts/Manuscript 1 - crop species/Env Microbiology submission/EMI 2nd submission/R commands/Nhu's_analyses")

#Color Brewer Pallette
#RColorBrewer::display.brewer.all()
#Manual color scales
PhylaColors <- c("Acidobacteriota"="#1c9e77", "Actinobacteriota"="#756fb4", "Bacteroidota"="#f45c85", 
                 "Chloroflexi"="#1c67dd", "Crenarchaeota"="gold", "Cyanobacteria"="seagreen", "Desulfobacterota"="gold1", 
                 "Firmicutes"="#d671aa", "Gemmatimonadota"="#ed4131", "Myxococcota"="#c2e871", "Methylomirabilota"="#1c6d77", "Nitrospirota"="deepskyblue",
                 "Patescibacteria"="#3f5d82", "Planctomycetota"="#7a81ff", "Proteobacteria"="orange", "Verrucomicrobiota"="#66a621")

FungalColors = c("Ascomycota"="#1c6d77", "Basidiomycota"="#FD8D3C", Chytridiomycota="deepskyblue", Glomeromycota="#66a621", Kickxellomycota="#7a81ff", Mortierellomycota="#ed4131", Mucoromycota="#1c67dd")

####################################################################################
### IMPORT AND PREPARE 16S DATASET
####################################################################################

#Import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
Crop_microbe <- import_biom("16S/table-LKM-with-taxonomy.biom")#, treefilename="tree-Crop-bact.nwk", refseqfilename="seqs-Crop-bact.fasta")

#Rename columns
colnames(tax_table(Crop_microbe)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Save otu table (optional) -- this will save a simple OTU table without taxonomy
#write.csv(as.data.frame(otu_table(Crop_microbe)),"otu_Crop_microbe_initial.csv")

#Read in metadata; standard QIIME2 metadata table can be used here
metadata <- read.csv("16S/metadata-16S-LKM.csv")

#Create another row for "SampleID"; create a dataframe for the metadata; add metadata to main phyloseq object
row.names(metadata) <- metadata$SampleID
metadata$SampleID <- metadata$SampleID
sample_data(Crop_microbe) <- metadata
#View(data.frame(sample_data(Crop_microbe)))

#De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_Crop_microbe_0 <- as.data.frame(tax_table(Crop_microbe))

#OPTIONAL export taxa table
# write.csv(tax_table_Crop_microbe_0, "tax_table_Crop_microbe_0.csv")

#Make a copy of the table
tax_table_Crop_microbe <- tax_table_Crop_microbe_0

#Renaming the taxonomy if not standard (this may not be necessary depending on the taxonomic database used)
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Acidobacteriota", "Acidobacteria", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Actinobacteriota", "Actinobacteria", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Armatimonadota", "Armatimonadetes", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Bacteroidota", "Bacteroidetes", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Gemmatimonadota", "Gemmatimonadetes", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Halobacterota", "Halobacteria", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Planctomycetota", "Planctomycetes", x)}))
#tax_table_Crop_microbe <- data.frame(lapply(tax_table_Crop_microbe, function(x) {gsub("Verrucomicrobiota", "Verrucomicrobia", x)}))

#Remove taxonomic notations
tax_table_Crop_microbe$Kingdom<- gsub("d__", "", tax_table_Crop_microbe$Kingdom)#sometimes "k" is replaced by "d" so make sure that is is properly removed.
tax_table_Crop_microbe$Phylum <- gsub("p__", "", tax_table_Crop_microbe$Phylum)
tax_table_Crop_microbe$Class <- gsub("c__", "", tax_table_Crop_microbe$Class)
tax_table_Crop_microbe$Order <- gsub("o__", "", tax_table_Crop_microbe$Order)
tax_table_Crop_microbe$Family <- gsub("f__", "", tax_table_Crop_microbe$Family)
tax_table_Crop_microbe$Genus <- gsub("g__", "", tax_table_Crop_microbe$Genus)
tax_table_Crop_microbe$Species <- gsub("s__", "", tax_table_Crop_microbe$Species)
#tax_table_Crop_microbe$Gen_Fam <- paste(tax_table_Crop_microbe$Genus, " (", tax_table_Crop_microbe$Family,")",sep="")

#View(tax_table_Crop_microbe_0)
#View(tax_table_Crop_microbe)

row.names(tax_table_Crop_microbe) <- row.names(tax_table_Crop_microbe_0)
tax_table(Crop_microbe) <- as.matrix(tax_table_Crop_microbe)
# View(data.frame(tax_table(Crop_microbe)))

#subsetting your datasets (often it will requires a slow narrowing down of each category until you get the samples you want)
#can also use for subsetting taxa: Crop_microbe_sub0 = subset_taxa(Crop_microbe, Kingdom=="Bacteria")
#Crop_microbe_cultivated = subset_samples(Crop_microbe, Soil == "Cultivated")#keep only cultivated samples
#Crop_microbe_uncultivated = subset_samples(Crop_microbe, Soil == "Uncultivated")#keep only uncultivated samples
#Crop_microbe = subset_samples(Crop_microbe, SampleID != "02.2.5.1.16S.a")#exclude the following sample based on ID

##Crop_microbe <- Crop_microbe_treatment

# View(data.frame(otu_table(Crop_microbe)))

# Check rarefaction of the data
# rarecurve(t(otu_table(Crop_microbe)), step=50, cex=0.5)

#removing samples that didn't work -- low read counts (this would best be done in QIIME)
#Crop_microbe_0 <- prune_samples(sample_sums(Crop_microbe) >= 100, Crop_microbe)#if done here, pass this object onto downstream workflow instead of Crop_microbe
# #View(data.frame(sample_data(Crop_microbe_0)))
# #View(data.frame(otu_table(Crop_microbe_0)))
# rarecurve(t(otu_table(Crop_microbe_0)), step=100, cex=0.5)
# Crop_microbe_0_otu <- data.frame(otu_table(Crop_microbe_0))
# write.csv(Crop_microbe_0_otu,"Crop_microbe_0_otu.csv")

###Note: This approach requires you to make a new phyloseq object for each comparison, which is shown below (using "subset_samples")

#Subsetting your dataset to make various comparisons
Crop_microbe_cultivated <- subset_samples(Crop_microbe, Soil=="Cultivated")#keep only cultivated samples
Crop_microbe_cultivated <- subset_samples(Crop_microbe_cultivated, CropCycle !="C2")#keep only C1 and C3 samples

Crop_microbe_uncultivated <- subset_samples(Crop_microbe, Soil=="Uncultivated")#keep only uncultivated samples
Crop_microbe_uncultivated <- subset_samples(Crop_microbe_uncultivated, CropCycle !="C2")#keep only C1 and C3 samples

Crop_microbe_first_last <-subset_samples(Crop_microbe, CycleSoil %in% c("C1_uncultivated", "C3_cultivated"))#keep the first and last sampling points


####################################################################################
###DESeq Comparison [C1, uncultivated (left) vs C3, cultivated (right)] -- USED IN MANUSCRIPT3f
###This comparison is to show how the very first sample compares to the very last sample, 
###without the noise of everything in between since it does seem that the community shifts in progression
####################################################################################
#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
Crop_microbe_first_last <- filter_taxa(Crop_microbe_first_last, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the Crop_microbe phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
Crop_microbe_first_last <- phyloseq((otu_table(Crop_microbe_first_last)+1), sample_data(Crop_microbe_first_last), tax_table(Crop_microbe_first_last))

#Convert phyloseq to DESeq Data Set object (dds)
Crop_microbe_dds <- phyloseq_to_deseq2(Crop_microbe_first_last, ~CropCycle)

#Determine which level should the dataset be set as the REFERENCE sample class
Crop_microbe_dds$CropCycle <- relevel(Crop_microbe_dds$CropCycle, "C1")

#Perform the contrast using standard and recognized parameters and tests
Crop_microbe_dds = DESeq(Crop_microbe_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(Crop_microbe_dds)

#Performing the final calculations and extracting the points
Crop_microbe_dds_results = results(Crop_microbe_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#Crop_microbe_dds_results = lfcShrink(dds=Crop_microbe_dds, contrast = c("SampleType","Soil","Mat"), res=Crop_microbe_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = Crop_microbe_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table = Crop_microbe_dds_results[which(Crop_microbe_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table = cbind(as(sig_table, "data.frame"), as(tax_table(Crop_microbe)[rownames(sig_table), ], "matrix"))

#View(sig_table_C1)
write.csv(sig_table, "16S/DESeq_sig_table_C1_C3.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_C1 = Crop_microbe_C1dds_results[which(Crop_microbe_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_C1 = cbind(as(notsig_table_C1, "data.frame"), as(tax_table(Crop_microbe_C1)[rownames(notsig_table_C1), ], "matrix"))
#head(notsig_table_C1)
# write.csv(notsig_table_C1, "notsig_table_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_Phyla <- unique(sig_table$Phylum)
sig_table_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_sub <- subset(sig_table, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("_", " ", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_sub$Genus)

#Plot the logfold changes
sig_table_subp <- ggplot(sig_table_sub, aes(x=log2FoldChange, y=Genus, color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Cycle 1            Cycle 3") + 
  theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  labs(x="Log2 Fold Change", y ="", tag="A") + 
  theme(plot.tag = element_text(size=20))

#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot_16S <- sig_table_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = PhylaColors)

Diff_Abund_Plot_16S

#Save a PDF of your file. Final edits can be made in Inkscape.
#ggsave("Figure3-DESeq-cultivated.pdf", device="pdf", width=8, height=9.5)



####################################################################################
###DESeq Comparison [C1 (left) vs C3 (right), cultivated soil] -- did not use in manuscript
####################################################################################

#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
Crop_microbe_filtered <- filter_taxa(Crop_microbe_cultivated, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the Crop_microbe phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
Crop_microbe_filtered <- phyloseq((otu_table(Crop_microbe_filtered)+1), sample_data(Crop_microbe_filtered), tax_table(Crop_microbe_filtered))

#Convert phyloseq to DESeq Data Set object (dds)
Crop_microbe_dds <- phyloseq_to_deseq2(Crop_microbe_filtered, ~CropCycle)

#Determine which level should the dataset be set as the REFERENCE sample class
Crop_microbe_dds$CropCycle <- relevel(Crop_microbe_dds$CropCycle, "C1")

#Perform the contrast using standard and recognized parameters and tests
Crop_microbe_dds = DESeq(Crop_microbe_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(Crop_microbe_dds)

#Performing the final calculations and extracting the points
Crop_microbe_dds_results = results(Crop_microbe_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#Crop_microbe_dds_results = lfcShrink(dds=Crop_microbe_dds, contrast = c("SampleType","Soil","Mat"), res=Crop_microbe_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = Crop_microbe_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table_cultivated = Crop_microbe_dds_results[which(Crop_microbe_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table_cultivated = cbind(as(sig_table_cultivated, "data.frame"), as(tax_table(Crop_microbe)[rownames(sig_table_cultivated), ], "matrix"))

#View(sig_table_cultivated_C1)
write.csv(sig_table_cultivated, "16S/ sig_table_cultivated.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_cultivated_C1 = Crop_microbe_C1dds_results[which(Crop_microbe_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_cultivated_C1 = cbind(as(notsig_table_cultivated_C1, "data.frame"), as(tax_table(Crop_microbe_C1)[rownames(notsig_table_cultivated_C1), ], "matrix"))
#head(notsig_table_cultivated_C1)
# write.csv(notsig_table_cultivated_C1, "notsig_table_cultivated_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_cultivated_Phyla <- unique(sig_table_cultivated$Phylum)
sig_table_cultivated_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_cultivated_sub <- subset(sig_table_cultivated, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_cultivated_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_cultivated_sub$Genus)
sig_table_cultivated_sub$Genus <- gsub("_", " ", sig_table_cultivated_sub$Genus)
sig_table_cultivated_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_cultivated_sub$Genus)
sig_table_cultivated_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_cultivated_sub$Genus)
sig_table_cultivated_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_cultivated_sub$Genus)

#Plot the logfold changes
sig_table_cultivated_subp <- ggplot(sig_table_cultivated_sub, aes(x=log2FoldChange, y=Genus, color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Cycle 1                                   Cycle 3") + 
  theme(plot.title = element_text(hjust = 0.5, size=11))

#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot <- sig_table_cultivated_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = PhylaColors) #+scale_color_brewer(palette="Dark2")
Diff_Abund_Plot

#Save a PDF of your file. Final edits can be made in Inkscape.
ggsave("Figure3-DESeq-cultivated.pdf", device="pdf", width=6, height=8)



####################################################################################
###DESeq Comparison [C1 (left) vs C3 (right), uncultivated soil] -- did not use in manuscript
####################################################################################

#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
Crop_microbe_unfiltered <- filter_taxa(Crop_microbe_uncultivated, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the Crop_microbe phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
Crop_microbe_unfiltered <- phyloseq((otu_table(Crop_microbe_unfiltered)+1), sample_data(Crop_microbe_unfiltered), tax_table(Crop_microbe_unfiltered))

#Convert phyloseq to DESeq Data Set object (dds)
Crop_microbe_dds <- phyloseq_to_deseq2(Crop_microbe_unfiltered, ~CropCycle)

#Determine which level should the dataset be set as the REFERENCE sample class
Crop_microbe_dds$CropCycle <- relevel(Crop_microbe_dds$CropCycle, "C1")

#Perform the contrast using standard and recognized parameters and tests
Crop_microbe_dds = DESeq(Crop_microbe_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(Crop_microbe_dds)

#Performing the final calculations and extracting the points
Crop_microbe_dds_results = results(Crop_microbe_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#Crop_microbe_dds_results = lfcShrink(dds=Crop_microbe_dds, contrast = c("SampleType","Soil","Mat"), res=Crop_microbe_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = Crop_microbe_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table_uncultivated = Crop_microbe_dds_results[which(Crop_microbe_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table_uncultivated = cbind(as(sig_table_uncultivated, "data.frame"), as(tax_table(Crop_microbe)[rownames(sig_table_uncultivated), ], "matrix"))

#View(sig_table_uncultivated_C1)
write.csv(sig_table_uncultivated, "16S/ sig_table_uncultivated.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_uncultivated_C1 = Crop_microbe_C1dds_results[which(Crop_microbe_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_uncultivated_C1 = cbind(as(notsig_table_uncultivated_C1, "data.frame"), as(tax_table(Crop_microbe_C1)[rownames(notsig_table_uncultivated_C1), ], "matrix"))
#head(notsig_table_uncultivated_C1)
# write.csv(notsig_table_uncultivated_C1, "notsig_table_uncultivated_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_uncultivated_Phyla <- unique(sig_table_uncultivated$Phylum)
sig_table_uncultivated_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_uncultivated_sub <- subset(sig_table_uncultivated, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_uncultivated_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_uncultivated_sub$Genus)
sig_table_uncultivated_sub$Genus <- gsub("_", " ", sig_table_uncultivated_sub$Genus)
sig_table_uncultivated_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_uncultivated_sub$Genus)
sig_table_uncultivated_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_uncultivated_sub$Genus)
sig_table_uncultivated_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_uncultivated_sub$Genus)

#Plot the logfold changes
sig_table_uncultivated_subp <- ggplot(sig_table_uncultivated_sub, aes(x=log2FoldChange, y=Genus, color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Cycle 1                      Cycle 3") + 
  theme(plot.title = element_text(hjust = 0.5, size=11)) +
  labs(x="", y ="", tag="B") +
  theme(plot.tag = element_text(size=20))
  
#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot_uncultivated <- sig_table_uncultivated_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = PhylaColors) #+scale_color_brewer(palette="Dark2")
#Diff_Abund_Plot_uncultivated

#Graph together & save a PDF of your file. Final edits can be made in Inkscape.
Diff_Abund_Plot_cultivated + Diff_Abund_Plot_uncultivated
ggsave("Figure3-DESeq-cultivated.pdf", device="pdf", width=10.5, height=7)



####################################################################################
###DESeq Comparison [C1 (left) vs C3 (right), both soils together] -- did not use in manuscript
###Since the soils have an influence on each other, it is better to analyze them as a whole unit
####################################################################################

#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
Crop_microbe <- filter_taxa(Crop_microbe, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the Crop_microbe phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
Crop_microbe <- phyloseq((otu_table(Crop_microbe)+1), sample_data(Crop_microbe), tax_table(Crop_microbe))

#Convert phyloseq to DESeq Data Set object (dds)
Crop_microbe_dds <- phyloseq_to_deseq2(Crop_microbe, ~CropCycle)

#Determine which level should the dataset be set as the REFERENCE sample class
Crop_microbe_dds$CropCycle <- relevel(Crop_microbe_dds$CropCycle, "C1")

#Perform the contrast using standard and recognized parameters and tests
Crop_microbe_dds = DESeq(Crop_microbe_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(Crop_microbe_dds)

#Performing the final calculations and extracting the points
Crop_microbe_dds_results = results(Crop_microbe_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#Crop_microbe_dds_results = lfcShrink(dds=Crop_microbe_dds, contrast = c("SampleType","Soil","Mat"), res=Crop_microbe_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = Crop_microbe_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table = Crop_microbe_dds_results[which(Crop_microbe_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table = cbind(as(sig_table, "data.frame"), as(tax_table(Crop_microbe)[rownames(sig_table), ], "matrix"))

#View(sig_table_C1)
write.csv(sig_table, "16S/ sig_table.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_C1 = Crop_microbe_C1dds_results[which(Crop_microbe_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_C1 = cbind(as(notsig_table_C1, "data.frame"), as(tax_table(Crop_microbe_C1)[rownames(notsig_table_C1), ], "matrix"))
#head(notsig_table_C1)
# write.csv(notsig_table_C1, "notsig_table_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_Phyla <- unique(sig_table$Phylum)
sig_table_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_sub <- subset(sig_table, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("_", " ", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_sub$Genus)

#Plot the logfold changes
sig_table_subp <- ggplot(sig_table_sub, aes(x=log2FoldChange, y=Genus, color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Cycle 1                   Cycle 3") + 
  theme(plot.title = element_text(hjust = 0.5, size=11)) + 
  labs(x="Log2 Fold Change", y ="", tag="A") + 
  theme(plot.tag = element_text(size=20))

#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot <- sig_table_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = PhylaColors) #+scale_color_brewer(palette="Dark2")
Diff_Abund_Plot

#Save a PDF of your file. Final edits can be made in Inkscape.
#ggsave("Figure3-DESeq-cultivated.pdf", device="pdf", width=8, height=9.5)






####################################################################################
### IMPORT AND PREPARE ITS DATASET
####################################################################################
#Import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
Crop_microbe <- import_biom("ITS/table-LKM-with-taxonomy.biom")

#Rename columns
colnames(tax_table(Crop_microbe)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Read in metadata; standard QIIME2 metadata table can be used here
metadata <- read.csv("ITS/metadata-ITS-LKM.csv")

#Create another row for "SampleID"; create a dataframe for the metadata; add metadata to main phyloseq object
row.names(metadata) <- metadata$SampleID
metadata$SampleID <- metadata$SampleID
sample_data(Crop_microbe) <- metadata

#De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_Crop_microbe_0 <- as.data.frame(tax_table(Crop_microbe))

#Make a copy of the table
tax_table_Crop_microbe <- tax_table_Crop_microbe_0

#Remove taxonomic notations
tax_table_Crop_microbe$Kingdom<- gsub("d__", "", tax_table_Crop_microbe$Kingdom)#sometimes "k" is replaced by "d" so make sure that is is properly removed.
tax_table_Crop_microbe$Phylum <- gsub("p__", "", tax_table_Crop_microbe$Phylum)
tax_table_Crop_microbe$Class <- gsub("c__", "", tax_table_Crop_microbe$Class)
tax_table_Crop_microbe$Order <- gsub("o__", "", tax_table_Crop_microbe$Order)
tax_table_Crop_microbe$Family <- gsub("f__", "", tax_table_Crop_microbe$Family)
tax_table_Crop_microbe$Genus <- gsub("g__", "", tax_table_Crop_microbe$Genus)
tax_table_Crop_microbe$Species <- gsub("s__", "", tax_table_Crop_microbe$Species)

row.names(tax_table_Crop_microbe) <- row.names(tax_table_Crop_microbe_0)
tax_table(Crop_microbe) <- as.matrix(tax_table_Crop_microbe)

Crop_microbe_first_last <-subset_samples(Crop_microbe, CycleSoil %in% c("C1_uncultivated", "C3_cultivated"))#keep the first and last sampling points

####################################################################################
###DESeq Comparison [C1, uncultivated (left) vs C3, cultivated (right)] -- USED IN MANUSCRIPT
###This comparison is to show how the very first sample compares to the very last sample, 
###without the noise of everything in between since it does seem that the community shifts in progression
####################################################################################
#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
Crop_microbe_first_last <- filter_taxa(Crop_microbe_first_last, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the Crop_microbe phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
Crop_microbe_first_last <- phyloseq((otu_table(Crop_microbe_first_last)+1), sample_data(Crop_microbe_first_last), tax_table(Crop_microbe_first_last))

#Convert phyloseq to DESeq Data Set object (dds)
Crop_microbe_dds <- phyloseq_to_deseq2(Crop_microbe_first_last, ~CropCycle)

#Determine which level should the dataset be set as the REFERENCE sample class
Crop_microbe_dds$CropCycle <- relevel(Crop_microbe_dds$CropCycle, "C1")

#Perform the contrast using standard and recognized parameters and tests
Crop_microbe_dds = DESeq(Crop_microbe_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(Crop_microbe_dds)

#Performing the final calculations and extracting the points
Crop_microbe_dds_results = results(Crop_microbe_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#Crop_microbe_dds_results = lfcShrink(dds=Crop_microbe_dds, contrast = c("SampleType","Soil","Mat"), res=Crop_microbe_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = Crop_microbe_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table = Crop_microbe_dds_results[which(Crop_microbe_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table = cbind(as(sig_table, "data.frame"), as(tax_table(Crop_microbe)[rownames(sig_table), ], "matrix"))

#View(sig_table_C1)
write.csv(sig_table, "ITS/DESeq_sig_table_C1_C3.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_C1 = Crop_microbe_C1dds_results[which(Crop_microbe_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_C1 = cbind(as(notsig_table_C1, "data.frame"), as(tax_table(Crop_microbe_C1)[rownames(notsig_table_C1), ], "matrix"))
#head(notsig_table_C1)
# write.csv(notsig_table_C1, "notsig_table_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_Phyla <- unique(sig_table$Phylum)
sig_table_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_sub <- subset(sig_table, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("_", " ", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_sub$Genus)

#Plot the logfold changes
sig_table_subp <- ggplot(sig_table_sub, aes(x=log2FoldChange, y=Genus, color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Cycle 1            Cycle 3") + 
  theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  labs(x="Log2 Fold Change", y ="", tag="B") + 
  theme(plot.tag = element_text(size=20))

#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot_ITS <- sig_table_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = FungalColors)

Diff_Abund_Plot_ITS

#Graph together & save a PDF of your file. Final edits can be made in Inkscape.
Diff_Abund_Plot_16S + Diff_Abund_Plot_ITS
ggsave("Figure4-DESeq-first-last-cycles.pdf", device="pdf", width=8, height=6.5)
