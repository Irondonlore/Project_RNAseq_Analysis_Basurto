################################################################################
##### SCRIPT TO PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS WITH DESEQ2 ######
#####                  ON THE AC-45_RNAseq-FFPE DATA                      ######
################################################################################

#| This script is written to perfom DESeq2 over the counts data collected from the
#| RNASeq of 238 patients diagnosed with PCa and BPH. Moreover, the analysis is 
#| done to the raw counts data (with 58k genes) and a counts matrix (with 27k) 
#| previously filtered by counts per million of reads. The filtering was performed
#| considering 50% of patients with PCa having at least approx 12 to 16 reads
#| per gene. The same was done to the BPH patients. Then the total number of genes
#| of both condition was combined and saved into an unique matrix.

#|    Notes from the DESeq2 documentation:

#| DEALING WITH CONTINUES VALUES IN DESEQ2:
#|    1) Convert them to numeric in the design to improve GLM convergence.
#|    2) Centralized them (you can use z-score for example)
#|    3) When having NA/NaN/Inf it is recommended to: 
#|          * Filling the empty cells with:
#|          * The mean value of the variable; or
#|          * the average of the two adjacent cells; or
#|          * a value calculated from a regression with another variable

#| What to do with <NA> values in categorical variables?
#| Recommendation Michale Love: Set it to another value. DESeq2 it is not able to
#| deal with NAN values.

#| To apply contrasts for continuous variables use "name" instead of "contrast"
#| in the function "results"

#| The alpha parameter in the results function it is a significance cutoff used 
#| for optimizing the independent filtering (by default 0.1). If the adjusted p-value 
#| cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

#| INFORMATION: res is a dataframe containing the information about some parameters 
#| required for the DEGs. By plotting the results and giving biological meaning, 
#| we can find the genes that are up or down regulated depending on the 2 conditions. 
#| The columns: 
#|    1) baseMean: the average of the normalized counts taking over all the samples
#|    2) log2FoldChanges: is the log2fold of the value of the gene in the treat 
#|    condition when comparing with the untreated. So, the positive values are the 
#|    upregulated genes in the treated condition and the negative values are the 
#|    downregulated genes in the treated condition
#|    3) lfcSE: estimates the standard error for the log2foldchange
#|    4) stat: are the wold test values for the genes.
#|    5) p-value: it is the p-value of the t-statistics for the gene

#|              Multiple design:
#| The idea behind multiple designs with DESeq2 is to control that variable 
#| contributing to the variance in our data to explore in more detail the gene 
#| expression difference over the conditions one is interested in. For instance, 
#| DESeq2 takes as "covariate" those variable that are added to the model except 
#| for the last one (whose the wald-test is applied)

#| MOREOVER: When doing complex design be careful about false positive that emerge
#| when data is not filtered by raw counts 

################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(pheatmap))
suppressMessages(library(ashr))
suppressMessages(library(calibrate))
suppressMessages(library(vsn))
suppressMessages(library(genefilter))
suppressMessages(library(gplots))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(EnhancedVolcano))
################################################################################


################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PTEN_loss_vs_intact/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE.txt"
counts.file_filtered <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE_cpm_filtered_PCa_BPH.txt"
setwd(dir.proj)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################



################################## DATA ########################################
#| Read table Full counts. NOTE: it is not necessary to do any normalization 
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)
counts_data <- subset(counts_data, select = -c(AC205, AC233))
counts_data_filtered <- read.table(counts.file_filtered, sep = "\t", stringsAsFactors = TRUE)

#| Read table Sample information
sample_info <- read.table(info.file, sep ="\t")

#| Other parameters:
FDR <- 0.05   #| False Discovery Rate Threshold
FC <- 1.5       #| FC change Threshold (=> |log2(FC)| >= 1)
theme_set(theme_classic()) #| For graphic design
################################################################################


## FOR EXPLORING GENES DIFFERENTIALLY EXPRESSED IN PTEN loss VS intact PROTEIN ##
#| Creating the design:
#|    Using an unique variable when changing the design
#|    Defining contrasts (for categorical) and names (for continuous) variables
counts_data <- counts_data[,sample_info$AC.basurto[which(sample_info$Diagnostico == "PCa" & !is.na(sample_info$H_score_cut_0))]]
counts_data_filtered <- counts_data_filtered[,sample_info$AC.basurto[which(sample_info$Diagnostico == "PCa" & !is.na(sample_info$H_score_cut_0))]]
sample_info <- sample_info[which(sample_info$Diagnostico == "PCa" & !is.na(sample_info$H_score_cut_0)), ]

sample_info$H_score_cut_0[which(sample_info$H_score_cut_0 =="PTEN loss")] <- "1"
sample_info$H_score_cut_0[which(sample_info$H_score_cut_0 =="PTEN intact")] <- "0"

sample_info$purity_Zscore <- (sample_info$purity - mean(sample_info$purity))/sd(sample_info$purity)

design <- ~  DV200_Zscore  + Edad_Zscore + H_score_cut_0
tag <- "DV200_Edad_Purity_H_score_cut_0"
contrast <- c("H_score_cut_0", "1", "0")
tag_contrast <- "PTEN_loss_vs_intact"
################################################################################


############################## DESeq2 ANALYSIS #################################

#| Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = design)

#| Run DESeq2
dds <- DESeq(dds)

#| Observing the results
resultsNames(dds)

#| Obtaining the results
res <- results(dds, alpha=FDR, contrast = contrast) 

summary(res)

#| Transforming to a data frame for downstream analysis 
res_dataframe <- data.frame(res)

#| Finding the DEGs
res_dataframe$DEGs <- NA
res_dataframe$DEGs[which((res_dataframe$padj <= FDR) & (!is.na(res_dataframe$padj)) & ((res_dataframe$log2FoldChange>log2(FC)) | (res_dataframe$log2FoldChange<(-log2(FC)))))] <- "Yes"
res_dataframe$DEGs[which(is.na(res_dataframe$DEGs))] <- "No"

#| Assigning gene names
res_dataframe$GeneID <- rownames(res_dataframe)

#| Merging with the genome file
res_dataframe <- merge(res_dataframe, genome_GRCh39.94, by ="GeneID")

#| Saving the table
final_results <- cbind(counts_data,res_dataframe)

#| Saving the table
write.table(final_results, paste("Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_", tag, "_Contrast", tag_contrast, ".txt", sep=""), sep ="\t", row.names =T)
################################################################################

res_dataframe <- read.table("Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_DV200_Edad_H_score_cut_0_ContrastPTEN_loss_vs_intact.txt")

############################## VOLCANO PLOT ####################################


#| To observe the volcano with the name of the top DEGs
res_dataframe$delabel <- NA
res_dataframe$delabel[which((res_dataframe$log2FoldChange < (-log2(2.5)) | (res_dataframe$log2FoldChange > (log2(2)))) & (res_dataframe$padj < FDR))] <- 
  res_dataframe$gene_name[which((res_dataframe$log2FoldChange < (-log2(2.5)) | (res_dataframe$log2FoldChange > log2(2))) & (res_dataframe$padj < FDR))]

#| Assigning direction to the DEGs
res_dataframe$DEGs_direction <- NA
res_dataframe$DEGs_direction[which((res_dataframe$log2FoldChange > log2(1.5)) & (res_dataframe$padj < 0.05))] <- "UP"
res_dataframe$DEGs_direction[which((res_dataframe$log2FoldChange < (-log2(1.5))) & (res_dataframe$padj < 0.05))] <- "DOWN"
res_dataframe$DEGs_direction[which(is.na(res_dataframe$DEGs_direction))] <- "No DEGs"

#| Plotting and saving
ggplot(res_dataframe, aes(x = log2FoldChange, y = -log10(padj), color =DEGs_direction, label=delabel )) + 
  geom_point(size =4,alpha = 0.5) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_color_manual(values = c("midnightblue","darkgrey", "mediumspringgreen") ,guide = "legend") +
  labs(color ="Direction") +
  xlab("Log2 Fold Change") +
  ggtitle("Volcano plot |FC| >= 1.5 and FDR < 0.05 \n Total DEGs PTEN genomic loss vs intact") +
  geom_text_repel(family ="serif") +
  xlim(-10,10)
ggsave("Results/design ~ DV200 + Edad + H-score/VolcanoPlots/Volcano_DGEs_PTEN_loss_vs_intact_FC_1-5_FDR_0-005_nofiltered.pdf", height = 6, width = 6.8)
################################################################################


############################ WATERFALL PLOT ####################################
DEGs <- res_dataframe[which(res_dataframe$DEGs =="Yes"),]
DEGs$Direction <- DEGs$log2FoldChange
DEGs$Direction[which(DEGs$Direction < 0)] <- "Down"
DEGs$Direction[which(DEGs$Direction >= 0 & DEGs$Direction != "Down")] <- "Up"

data3 <- DEGs
data3 <- data3[which(!(duplicated(data3$gene_name))),]
data3$gene_name <- factor(data3$gene_name,levels = data3$gene_name[order(data3$log2FoldChange, decreasing = FALSE)])

ggplot(data3, aes(x = gene_name, y = log2FoldChange, fill =Direction)) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_bar( stat = "identity") +
  xlab("Gene name") +
  ylab("Log2 Fold Change") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values =c("midnightblue", "mediumspringgreen"))
ggsave("Results/design ~ DV200 + Edad + H-score/Log2FC/Log2FC_DEGs_PTEN_loss_vs_intact.pdf", height = 7, width = 8)

#| Top up and down
data3 <- data3[order(data3$log2FoldChange),]
data3$Direction[1:20]
data3$Direction[(dim(data3)[1] -20):dim(data3)[1]]

top_up_down <- rbind(data3[1:20,], data3[(dim(data3)[1] -11):dim(data3)[1],])
ggplot(top_up_down, aes(x = log2FoldChange, y = gene_name, fill =Direction)) +
  theme(text=element_text(size=13,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #coord_cartesian(ylim = c(-5,5)) +
  geom_bar( stat = "identity") +
  xlab("Log2 Fold Change") +
  ylab("Gene name") +
  scale_fill_manual(values =c("midnightblue", "mediumspringgreen"))
ggsave("Results/design ~ DV200 + Edad + H-score/Log2FC/Log2FC_DEGs_PTEN_loss_vs_intact_TOP.pdf", height = 7, width = 8)
################################################################################

length(DEGs$Direction[which(DEGs$Direction == "Down")])
length(DEGs$Direction[which(DEGs$Direction == "Up")])
