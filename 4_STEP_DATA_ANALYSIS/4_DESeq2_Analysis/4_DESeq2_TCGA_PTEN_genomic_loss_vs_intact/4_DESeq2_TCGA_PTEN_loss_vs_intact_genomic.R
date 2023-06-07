################################################################################
####          DESeq2 IN TCGA DATA TO FIND DEGs BY COMPARING               ######
####             PTEN LOSS VS INTACT AT THE GENOMIC LEVEL                 ######
################################################################################


#| 


##########################  LIBRARIES AND DATA #################################
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(variancePartition))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(corrr))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(recount))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(gprofiler2))
suppressMessages(library(survminer))

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Other parameters:
FDR <- 0.05   #| False Discovery Rate Threshold
FC <- 1.5      #| FC change Threshold (>= |log2(FC)| >= )

#| Data directory
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_TCGA_PTEN_genomic_loss_vs_intact/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_merge.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/RawCounts/TCGA_raw_counts.txt"
counts.file_filtered <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/RawCounts/TCGA_raw_counts_filtered.txt"
basurto_tcga.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PTEN_loss_vs_intact/Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_DV200_Edad_H_score_cut_0_ContrastPTEN_loss_vs_intact_filtered.txt"
wgcna_basurto.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors.txt"
setwd(dir.proj)

#| Obtaining data
counts_data <- read.table(counts.file, sep ="\t", header=T)
sample_info <- read.table(info.file, sep ="\t", header=T)
degs_basurto <- read.table(basurto_tcga.file, sep ="\t", header =T)
wgcna_basurto <- read.table(wgcna_basurto.file, sep ="\t", header =T)

#| Changing names in colnames(counts_data) to match with sample_info$PATIENT_ID
colnames(counts_data) <- gsub("\\.", "-", colnames(counts_data))

#| Renaming rownames in the sample_info
rownames(sample_info) <- sample_info$PATIENT_ID

#| Verifying that rownames from sample_info and colnames from counts data match
any(rownames(sample_info) == colnames(counts_data))
################################################################################



################## DESeq2 PTEN LOSS VS INTACT GENOMIC ONLY #####################
#| Filtering by PTEN loss and intact
sample <- sample_info[which(sample_info$PTEN_cna =="0" | sample_info$PTEN_cna =="-2"),]
counts <- counts_data[,sample$PATIENT_ID]

#| Assigning PTEN loss vs intact
sample$PTEN[which(sample$PTEN_cna =="0")] <- "intact"
sample$PTEN[which(sample$PTEN_cna =="-2")] <- "loss"

#| Z-score of AGE
sample$AGE_zcore <- (sample$AGE - mean(sample$AGE))/sd(sample$AGE)

#| Constructing DESeq2 object
dds <- DESeqDataSetFromMatrix(counts, sample, ~ AGE_zcore + PTEN)

#| Run DESeq2
dds <- DESeq(dds)

#| Observing the results 
resultsNames(dds)

#| Obtaining the results
res <- results(dds, alpha=FDR, contrast = c("PTEN", "loss", "intact")) 

#| Transforming res to data.frame
res <- data.frame(res)

#| New variable
res_dataframe <- res

#| Finding the DEGs
res_dataframe$DEGs <- NA
res_dataframe$DEGs[(res_dataframe$padj <= FDR) & (!is.na(res_dataframe$padj)) & ((res_dataframe$log2FoldChange>log2(FC)) | (res_dataframe$log2FoldChange<(-log2(FC))))] <- "Yes"
res_dataframe$DEGs[which(is.na(res_dataframe$DEGs))] <- "No"

#| Assigning gene names
res_dataframe$GeneID <- rownames(res_dataframe)

#| Merging with the genome file
res_dataframe <- merge(res_dataframe, genome_GRCh39.94, by ="GeneID")

#| Saving the table
final_results <- cbind(counts,res_dataframe)

#| Saving table of results
write.table(res_dataframe, "Results/Tables/res_DESeq2_TCGA_PTEN_genomic_loss_vs_intact_design_AGE_PTEN_CNA_no.txt" , sep ="\t")
################################################################################


############################## VOLCANO PLOT ####################################

#| To observe the volcano with the name of the top DEGs
res_dataframe$delabel <- NA
res_dataframe$delabel[which((res_dataframe$log2FoldChange < (-log2(FC)) | (res_dataframe$log2FoldChange > (log2(FC)))) & (res_dataframe$padj < FDR))] <- 
  res_dataframe$gene_name[which((res_dataframe$log2FoldChange < (-log2(FC)) | (res_dataframe$log2FoldChange > log2(FC))) & (res_dataframe$padj < FDR))]

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
  geom_text_repel(family ="serif")
ggsave("Results/VolcanoPlots/Volcano_DGEs_PTEN_loss_vs_intact_FC_1-5_FDR_0-005_nofiltered.pdf", height = 6, width = 6.8)
################################################################################



############################ WATERFALL PLOT ####################################
#| Selecting only DEGs
DEGs <- res_dataframe[which(res_dataframe$DEGs =="Yes"),]
DEGs <- DEGs[which(!(duplicated(DEGs$gene_name))),]
DEGs$gene_name <- factor(DEGs$gene_name ,levels = DEGs$gene_name[order(DEGs$log2FoldChange, decreasing = TRUE)])

#| Plotting and saving
ggplot(DEGs, aes(x = gene_name, y = log2FoldChange, fill =DEGs_direction)) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_bar( stat = "identity",width =0.4) +
  xlab("Gene name") +
  ylab("Log2 Fold Change")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values =c("midnightblue", "mediumspringgreen"))
ggsave("Results/Log2FC/Log2FC_DEGs_PTEN_genomic_loss_vs_intact_WaterFall_nofiltered.pdf", height = 7, width = 8)


#| Top up and down
DEGs <- DEGs[order(DEGs$log2FoldChange),]

top_up_down <- rbind(DEGs[1:25,], DEGs[(dim(DEGs)[1] -20):dim(DEGs)[1],])
ggplot(top_up_down, aes(x = log2FoldChange, y = gene_name, fill =DEGs_direction)) +
  theme(text=element_text(size=13,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  geom_bar( stat = "identity") +
  xlab("Log2 Fold Change") +
  ylab("Gene name") +
  scale_fill_manual(values =c("midnightblue", "mediumspringgreen"))
ggsave("Results/Log2FC/Log2_FC_Gene_name_Waterfall_TOP_no_filtered.pdf", height = 7, width = 8)
################################################################################





################ COMPARING DEGS FROM TCGA WITH DEGS FROM BASURTO ###############
degs_tcga <- res_dataframe
degs_tcga_gene_name <- degs_tcga$GeneID[which(degs_tcga$padj < 0.05 & degs_tcga$log2FoldChange > log2(1.5) | degs_tcga$log2FoldChange < (-log2(1.5)) )]
degs_basurto_genes_name <- degs_basurto$GeneID[which(degs_basurto$padj < 0.05 & degs_basurto$log2FoldChange > log2(1.5)  | degs_basurto$log2FoldChange < (-log2(1.5) ) )]



genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% degs_tcga_gene_name[which(degs_tcga_gene_name %in% degs_basurto_genes_name)])]

#| Make any sense doing the enrichment?
total_gost_DEGs <- gost(list("DEGs PTEN loss vs intact Genomic" = degs_tcga_gene_name[which(degs_tcga_gene_name %in% degs_basurto_genes_name)] ),
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(degs_tcga), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs,capped=FALSE)
total_gost_DEGs <- total_gost_DEGs$result[order(total_gost_DEGs$result$p_value),]
total_gost_DEGs$term_name[1:10]
################################################################################



################################ Log2-Log2 FC ##################################
#| Finding commong genes
tcga <- degs_tcga[which(degs_tcga$GeneID %in% degs_basurto$GeneID ),]
basurto <- degs_basurto[which(degs_basurto$GeneID %in% tcga$GeneID),]

#| Creating a new dataframe
data <- data.frame(Log2FC_TCGA = degs_tcga$log2FoldChange[which(degs_tcga$GeneID %in% basurto$GeneID)],
                   Log2FC_Basurto = basurto$log2FoldChange)


cor <- cor.test(data$Log2FC_TCGA, data$Log2FC_Basurto, method ="pearson")
p <- round(cor$p.value,11)
r <- round(as.numeric(cor$estimate[[1]]),2)

ggplot(data, aes(x= Log2FC_TCGA, y= Log2FC_Basurto)) +
  geom_point(color = "midnightblue") +
  geom_smooth(method = "lm", color="darkgoldenrod3" ) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  xlab("Log2 FC TCGA") +
  ylab("Log2 FC Basurto")+
  annotate(geom="text", x=-2.5, y=-2, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8)
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_TCGA_PTEN_genomic_loss_vs_intact/Results/ScatterPlots/Log2-Log2_FC_Basurto_vs_TCGA_nofiltered.pdf", heigh =4.5, width =6)
################################################################################




######### OVERLAPPING DEGS FROM TCGA WITH WGCNA MODULES FROM BASURTO ###########
module_colors <- unique(wgcna_basurto$moduleColors)
genes_tcga <- res_dataframe$GeneID[which(res_dataframe$DEGs=="Yes")]

#| New dataframe
data_percentage_DGEs <- data.frame(modules = module_colors,
                                   fraction = c(length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="turquoise")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="grey")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="brown")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="red")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="purple")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="yellow")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="blue")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="magenta")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="green")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="black")])])/length(genes_tcga),
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="pink")])])/length(genes_tcga), 
                                                length(genes_tcga[which(genes_tcga %in% rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="greenyellow")])])/length(genes_tcga)),
                                   size = c(length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="grey")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="brown")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="red")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="purple")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="yellow")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="blue")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="magenta")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="green")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="black")]),
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="pink")]), 
                                            length(rownames(wgcna_basurto)[which(wgcna_basurto$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_DGEs, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Fraction of DEGs TCGA (PTEN grenomic loss vs intact) in the WGCNA \n modules of Basurto (PCa samples)")+
  labs(fill = "Modules") +
  labs(size = "Module size") +
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing PTEN genomic loss vs intact TCGA") 
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_TCGA_PTEN_genomic_loss_vs_intact/Results/IntersectionDEGs/BubblePlot_DEGs_TCGA_in_WGCNA_Modules_Basurto_nofiltered.pdf", heigh =6, width =7)
################################################################################
