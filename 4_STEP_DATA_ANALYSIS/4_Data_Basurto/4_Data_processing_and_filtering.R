################################################################################
####    PROCESSING AND FILTERING COUNTS AND SAMPLE INFO BASURTO COHORT     #####
################################################################################

#| This code process the RNAseq data from Basurto cohort containing 240 patients 
#| (198 Prostate Cancer and 40 Benign Hyperplasia). Here I discarded 2 patients 
#| that we diagnosed by "Cancer Incidential". I have applied filtering and normalization
#| to the counts data

#| Notes:

#| NORMALIZATION: DESeq2 uses median of rations, this means that counts are divided
#| by sample-specific size factors determined by median ratio of gene counts relative 
#| to geometric mean per gene. This normalization is usel to perform gene count
#| comparisons between samples and for DE analysis; NOT for within sample comparison

#| VST: We applied this normalization to inpute the value to WGCNA. The description
#| behind the function is the following: The function calculates a variance stabilizing 
#| transformation (VST) from the fitted dispersion-mean relation(s) and then transforms 
#| the count data (normalized by division by the size factors or normalization factors), 
#| yielding a matrix of values which are now approximately homoskedastic (having constant 
#| variance along the range of mean values). The transformation also normalizes with 
#| respect to library size. 

#| Filtering by fragments/counts per million: This filtering takes into account
#| the library size. The threshold must be taking depending on the library size.
#| If we have a library between 50mill and 66mill reads, the appropiate threshold
#| to use is 0.2, because it is taking genes with reads higher that 12 and 15.
#| The advantage of using CPM/FPM filtering rather than counts filtering is that 
#| it accounts for differences in sequencing depth between libraries.

################################ LIBRARIES #####################################
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(edgeR))

#| For plots
theme_set(theme_classic()) 
################################################################################


############################# DATA DIRECTORIES  ################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto"
#info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Basurto_BBDD_AC45_Patients_pheno.xlsx"#"X:/sgarcia/Basurto pheno NDM & AC45/Basurto_BBDD_AC45_Patients_pheno.xlsx"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"#"X:/sgarcia/Basurto pheno NDM & AC45/Basurto_BBDD_AC45_Patients_pheno.xlsx"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE.txt"
dv200.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/DV200_Basurto_238_samples.txt"
genome.file <- "X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt"
setwd(dir.proj)
################################################################################


############################## READING DATA ####################################
#| Read table Full counts 
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Read table Sample information
sample_info <- read.table(info.file, sep ="\t")

#| Sorting sample info by ´AC code´
sample_info <- sample_info[order(sample_info$`AC code`),]

#| Checking they have the same names
any(colnames(counts_data) %in% sample_info$`AC code`)

#| Checking that they are in the same order:
all(colnames(counts_data) == sample_info$`AC code`)

#| Deleting "Cancer Incidental" cases (AC205 y AC233)
counts_data <- subset(counts_data, select = -c(AC205, AC233))
sample_info <- sample_info[which(sample_info$AC.basurto %in% colnames(counts_data)),]

#| Assigning ´AC code´ to rownames
rownames(sample_info) <- sample_info$AC.basurto
names(sample_info)
write.table(sample_info,"Results/Sample_info_table/sample_info_extracted_2.txt", sep ="\t")
################################################################################


################################ PRE-FILTERING #################################
PCa_samples <- rownames(sample_info)[which(sample_info$Diagnostico == "PCa")] # 198
BPH_samples <- rownames(sample_info)[which(sample_info$Diagnostico == "BPH")] # 40

counts_data_PCa <- counts_data[,PCa_samples]
counts_data_BPH <- counts_data[,BPH_samples]

#| To check the average size of the library: Min: 50299454, Max: 66391979, mean: 63088766
library_size <- read.table("X:/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/TableReads.txt", sep ="\t", header =T)
mean(library_size$Mill..reads)

#| FILTERING PER GENE COUNTS: Observing those patients that have genes with counts higher 
#| than 5. Then, counts how many columns these condition was satisfied. If 50% of the sum 
#| of the columns is higher or equal than the sum of those that satisfied the condition 
#| of the counts higher than 5, then select that genes, otherwise, filter them
filter_counts <- rowSums(counts_data>5) >= 0.5*ncol(counts_data)
filter_counts_PCa <- rowSums(counts_data_PCa>5) >= 0.5*ncol(counts_data_PCa)
filter_counts_BPH <- rowSums(counts_data_BPH>5) >= 0.5*ncol(counts_data_BPH)

#| FILTERING BY COUNTS PER MILLION: This filtering follows the same logic as the one above 
#| but instead of filtering by total counts, the CPM is computing, and then those higher than
#| 0.2 (this number will depend on the library size) will be selected, otherwise, filtered.
filter_cpm <- rowSums(cpm(counts_data)>0.2) >= 0.5*ncol(counts_data)
filter_cpm_PCa <- rowSums(cpm(counts_data_PCa)>0.2) >= 0.5*ncol(counts_data_PCa)
filter_cpm_BPH <- rowSums(cpm(counts_data_BPH)>0.2) >= 0.5*ncol(counts_data_BPH)

#| Log2 of counts data
log_counts <- as.data.frame(log2(counts_data + 1))
log_counts_PCa <- as.data.frame(log2(counts_data_PCa + 1))
log_counts_BPH <- as.data.frame(log2(counts_data_BPH + 1))

#| Log2 of filtered counts data
log_counts_filtered <- as.data.frame(log2(counts_data[filter_counts,] + 1))
log_counts_filtered_BPH <- as.data.frame(log2(counts_data_BPH[filter_counts_BPH,] + 1))
log_counts_filtered_PCa <- as.data.frame(log2(counts_data_PCa[filter_counts_PCa,] + 1))

log_counts_filtered_phenotype <- merge(log_counts_filtered_PCa, log_counts_filtered_BPH, by="row.names", all=TRUE)
rownames(log_counts_filtered_phenotype) <- log_counts_filtered_phenotype$Row.names
log_counts_filtered_phenotype <- log_counts_filtered_phenotype[,-c(1)]

#| Log2 of counts data normalized by counts per millions
log_counts_cpm <- cpm(counts_data, log = TRUE)
log_counts_cpm_PCa <- cpm(counts_data_PCa, log = TRUE)
log_counts_cpm_BPH <- cpm(counts_data_BPH, log = TRUE)

#| Log2 of filtered counts data normalized by counts per millions
log_counts_cpm_filtered <- cpm(counts_data, log = TRUE)[filter_cpm,]
log_counts_cpm_filtered_PCa <- cpm(counts_data_PCa, log = TRUE)[filter_cpm_PCa,]
log_counts_cpm_filtered_BPH <- cpm(counts_data_BPH, log = TRUE)[filter_cpm_BPH,]
log_counts_cpm_filtered_phenotype <- merge(log_counts_cpm_filtered_PCa, log_counts_cpm_filtered_BPH, by="row.names", all=TRUE)
rownames(log_counts_cpm_filtered_phenotype) <- log_counts_cpm_filtered_phenotype$Row.names
log_counts_cpm_filtered_phenotype <- log_counts_cpm_filtered_phenotype[,-c(1)]

#| Density plots: I have made this plots to observe the effect of filtering based
#| on filtering by counts and cpm.
pdf("Results/RawDataAnalysis/log2_counts_raw_data.pdf")
plotDensities(log_counts, legend=FALSE, main="Density of log2 counts of raw data")
dev.off()

pdf("Results/RawDataAnalysis/log2_counts_raw_data_CPM.pdf")
plotDensities(log_counts_cpm, legend=FALSE, main="Density of log2 counts of raw data in CPM")
dev.off()

pdf("Results/RawDataAnalysis/log2_counts_raw_data_filtered.pdf")
plotDensities(log_counts_filtered, legend=FALSE, main="Density of log2 counts of raw data in Filtered")
dev.off()

pdf("Results/RawDataAnalysis/log2_counts_raw_data_filtered_Phenotype_PCa_vs_BPH.pdf")
plotDensities(log_counts_filtered_phenotype, legend=FALSE, main="Density of log2 counts of raw data in Filtered")
dev.off()

pdf("Results/RawDataAnalysis/log2_counts_raw_data_CPM_filtered.pdf")
plotDensities(log_counts_cpm_filtered, legend=FALSE, main="Density of log2 counts of raw data in CPM Filtered")
dev.off()

pdf("Results/RawDataAnalysis/log2_counts_raw_data_CPM_filtered_Phenotype_PCa_vs_BPH.pdf")
plotDensities(log_counts_cpm_filtered_phenotype, legend=FALSE, main="Density of log2 counts of raw data in CPM Filtered")
dev.off()

#| Checking dimensionality
dim(log_counts_filtered)                 #| 27503   
dim(log_counts_filtered_phenotype)        #| 29331   
dim(log_counts_cpm_filtered)             #| 24501   
dim(log_counts_cpm_filtered_phenotype)    #| 26467   
dim(log_counts_cpm_filtered_PCa)         #| 24270   
dim(log_counts_cpm_filtered_BPH)         #| 25916  
################################################################################



##########################  DESEQ2 NORMALIZATION  ##############################

#| Selecting the filtering by PCa vs BPH
counts_data_filtered_cpm <- counts_data[rownames(log_counts_cpm_filtered),]
counts_data_filtered_cpm_phenotype <- counts_data[rownames(log_counts_cpm_filtered_phenotype),]
counts_data_filtered_counts_phenotype <- counts_data[rownames(log_counts_filtered_phenotype),]

#| filtering by PTEN loss vs intact
counts_data_filtered_cpm_phenotype_PTEN_loss_intact <- counts_data_filtered_cpm_phenotype[,sample_info$AC.basurto[which(!is.na(sample_info$H_score_cut_0))]]
counts_data_filtered_counts_phenotype_PTEN_loss_intact <- counts_data_filtered_counts_phenotype[,sample_info$AC.basurto[which(!is.na(sample_info$H_score_cut_0))]]

#| Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ 1)
dds_cpm <- DESeqDataSetFromMatrix(countData = counts_data_filtered_cpm, colData = sample_info, design = ~ 1)
dds_cpm_phenotype <- DESeqDataSetFromMatrix(countData = counts_data_filtered_cpm_phenotype, colData = sample_info, design = ~ 1)
dds_counts_phenotype <- DESeqDataSetFromMatrix(countData = counts_data_filtered_counts_phenotype, colData = sample_info, design = ~ 1)

sample_info <- sample_info[which(!is.na(sample_info$H_score_cut_0)),]
dds_cpm_phenotype_PTEN_loss_intact <- DESeqDataSetFromMatrix(countData = counts_data_filtered_cpm_phenotype_PTEN_loss_intact, colData = sample_info, design = ~ 1)
dds_counts_phenotype_PTEN_loss_intact <- DESeqDataSetFromMatrix(countData = counts_data_filtered_counts_phenotype_PTEN_loss_intact, colData = sample_info, design = ~ 1)

#| Estimating the size factors
dds <- estimateSizeFactors(dds)
dds_cpm <- estimateSizeFactors(dds_cpm)
dds_cpm_phenotype <- estimateSizeFactors(dds_cpm_phenotype)
dds_counts_phenotype <- estimateSizeFactors(dds_counts_phenotype)

dds_cpm_phenotype_PTEN_loss_intact <- estimateSizeFactors(dds_cpm_phenotype_PTEN_loss_intact)
dds_counts_phenotype_PTEN_loss_intact <- estimateSizeFactors(dds_counts_phenotype_PTEN_loss_intact)

#| Computing normalization
counts_data_normalized <- counts(dds, normalized=TRUE)
counts_data_filtered_cpm_normalized <- counts(dds_cpm, normalized=TRUE)
counts_data_filtered_cpm_phenotype_normalized <- counts(dds_cpm_phenotype, normalized=TRUE)
counts_data_filtered_counts_phenotype_normalized <- counts(dds_counts_phenotype, normalized=TRUE)

counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact <- counts(dds_cpm_phenotype_PTEN_loss_intact, normalized=TRUE)
counts_data_filtered_counts_phenotype_normalized_PTEN_loss_intact <- counts(dds_counts_phenotype_PTEN_loss_intact, normalized=TRUE)

#| log2 for re-escaling factors
counts_data_normalized <- log(counts_data_normalized + 1, base =2)
counts_data_filtered_cpm_normalized <- log(counts_data_filtered_cpm_normalized + 1, base =2)
counts_data_filtered_cpm_phenotype_normalized <- log(counts_data_filtered_cpm_phenotype_normalized + 1, base =2)
counts_data_filtered_counts_phenotype_normalized <- log(counts_data_filtered_counts_phenotype_normalized + 1, base =2)


counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact <- log(counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact + 1, base =2)
counts_data_filtered_counts_phenotype_normalized_PTEN_loss_intact <- log(counts_data_filtered_counts_phenotype_normalized_PTEN_loss_intact + 1, base =2)

#| Box plot
boxplot(counts_data_normalized)
boxplot(counts_data_filtered_cpm_normalized)
boxplot(counts_data_filtered_cpm_phenotype_normalized)

boxplot(counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact)
boxplot(counts_data_filtered_counts_phenotype_normalized_PTEN_loss_intact)

dim(counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact)
#| Saving counts data without cancer incidental cases
write.table(data.frame(counts_data),"Results/Counts/FullCounts_Basurto.txt", sep ="\t", row.names =T)

#| Saving counts data without cancer incidental cases
write.table(data.frame(counts_data_filtered_cpm_phenotype),"Results/Counts/FullCounts_Basurto_Filtered_CPM_phenotype.txt", sep ="\t", row.names =T)

#| Saving counts data without cancer incidental cases
write.table(data.frame(counts_data_filtered_counts_phenotype),"Results/Counts/FullCounts_Basurto_Filtered_counts_phenotype.txt", sep ="\t", row.names =T)

#| Saving normalized and log2 counts data without cancer incidental cases
write.table(data.frame(counts_data_normalized),"Results/Normalized_Log2_Counts/FullCounts_Basurto_Normalized_Log2.txt", sep ="\t", row.names =T)

#| Saving normalized and log2 counts data filtered by CPM (Phenotype) without cancer incidental cases
write.table(data.frame(counts_data_filtered_cpm_phenotype_normalized),"Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_CPM_phenotype_Normalized_Log2.txt", sep ="\t", row.names =T)

#| Saving normalized and log2 counts data filtered by CPM (Phenotype) without cancer incidental cases
write.table(data.frame(counts_data_filtered_counts_phenotype_normalized),"Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_counts_phenotype_Normalized_Log2.txt", sep ="\t", row.names =T)

#| Saving normalized and log2 counts data filtered by CPM (Phenotype) PTEN loss vs intact without cancer incidental cases
write.table(data.frame(counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact),"Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_CPM_phenotype_Normalized_Log2_PTEN_loss_intact.txt", sep ="\t", row.names =T)

#| Saving normalized and log2 counts data filtered by counts (Phenotype) PTEN loss vs intact without cancer incidental cases
write.table(data.frame(counts_data_filtered_counts_phenotype_normalized_PTEN_loss_intact),"Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_COUNTS_phenotype_Normalized_Log2_PTEN_loss_intact.txt", sep ="\t", row.names =T)

#| Saving counts data with gene name as row names for xCell
gene_names <- genome_GRCh39.94$gene_name[which(rownames(counts_data_filtered_counts_phenotype_normalized) %in% genome_GRCh39.94$GeneID)]
counts_data_filtered_counts_phenotype_normalized <- as.data.frame(counts_data_filtered_counts_phenotype_normalized)
counts_data_filtered_counts_phenotype_normalized_gene_name <- counts_data_filtered_counts_phenotype_normalized
counts_data_filtered_counts_phenotype_normalized_gene_name$GeneID <- rownames(counts_data_filtered_counts_phenotype_normalized_gene_name)
counts_data_filtered_counts_phenotype_normalized_gene_name <- merge(counts_data_filtered_counts_phenotype_normalized_gene_name,genome_GRCh39.94, by ="GeneID" )
counts_data_filtered_counts_phenotype_normalized_gene_name <- counts_data_filtered_counts_phenotype_normalized_gene_name[,c(colnames(counts_data_filtered_counts_phenotype_normalized), "gene_name")]
counts_data_filtered_counts_phenotype_normalized_gene_name <- aggregate(counts_data_filtered_counts_phenotype_normalized_gene_name[ ,colnames(counts_data_filtered_counts_phenotype_normalized)], by = list(counts_data_filtered_counts_phenotype_normalized_gene_name$gene_name), FUN = mean)
rownames(counts_data_filtered_counts_phenotype_normalized_gene_name) <- counts_data_filtered_counts_phenotype_normalized_gene_name$Group.1
counts_data_filtered_counts_phenotype_normalized_gene_name <- counts_data_filtered_counts_phenotype_normalized_gene_name[,-c(1)]

write.table(counts_data_filtered_counts_phenotype_normalized_gene_name,"Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_counts_phenotype_Normalized_Log2_GENE_NAME_ROWS_xCELL_DATA.txt", sep ="\t", row.names =T)

################################################################################

packageVersion("variancePartition")

####################### PROCESSING SAMPLE INFORMATION ##########################
#| Adding PTEN as a clinical feature from the counts data (with z-score transformation):
#|    *) Extracting the PTEN expression from the counts matrix
counts_data_filtered_cpm_phenotype_normalized <- data.frame(counts_data_filtered_cpm_phenotype_normalized)
sample_info$PTEN_Exp <- as.numeric(counts_data_filtered_cpm_phenotype_normalized[c("ENSG00000171862"),])

#| Mean Edad
mean(sample_info$Edad[(sample_info$Diagnostico == "PCa")])  # 62.84
mean(sample_info$Edad[(sample_info$Diagnostico == "BPH")])  # 71.2

#| Z-score of the "Edad"
sample_info$Edad_Zscore <- (sample_info$Edad - mean(sample_info$Edad))/sd(sample_info$Edad)
sample_info$Edad_Zscore <- as.numeric(sample_info$Edad_Zscore)

#| Adding DVD200 (RNA quality)
DV200 <- read.table(dv200.file, sep ="\t")
DV200$`AC code` <- rownames(DV200)
sample_info <- merge(sample_info, DV200, by = "AC code", all.x=TRUE)
sample_info$DV200_Zscore <- (sample_info$DV200 - mean(sample_info$DV200))/sd(sample_info$DV200)
sample_info$DV200 <- as.numeric(sample_info$DV200)

#| Body-Mass-Index
sample_info$BMI_Zscore <- (sample_info$BMI - mean(sample_info$BMI))/sd(sample_info$BMI)

#| Selecting columns of interest
sample_info_extracted <- sample_info[,c("Edad","PTEN_Exp","Diagnostico","DV200", "Edad_Zscore","DV200_Zscore", "DFS.TIME", "DFS.STATUS", "Gleason_score_pieza", "Estadio_pN_nodulospos_pieza", "BMI", "PSA", "MFS.STATUS", "MFS.TIME")]

#| Saving table of results 
write.table(sample_info_extracted, "Results/Sample_info_table/sample_info_extracted_processing_and_filtering.txt")
################################################################################



########################## STABILIZING COUNT VARIANCE ##########################
dim(sample_info)
#| Transforming the data into a DESeq object (Run this line before performing PCA)
dds <- DESeqDataSetFromMatrix(countData = counts_data_filtered_counts_phenotype,
                              colData = sample_info,
                              design = ~ 1)
#| Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

QC =TRUE

if(QC ==TRUE){
  
  #| In order to test for differential expression, we operate on raw counts and 
  #| use discrete distributions. However for other downstream analyses - e.g. for 
  #| visualization or clustering, or WGCNA (see manual) it might be useful to work 
  #| with transformed versions of the count data.
  ##############################################################################
  
  if(nrow(data.frame(counts_data_filtered_cpm_phenotype_normalized_PTEN_loss_intact)) >= 30) { #| For large datasets, apply VST transformation. It 
    #| does not use the design to remove variation in the data.
    #| It uses the design formula to calculate the within-group variability (if blind=FALSE)
    #| or the across-all-samples variability (if blind=TRUE). By default blind=TRUE. FALSE 
    #| should be used to use the data in downstream analysis. If many genes have large 
    #| differences in counts due to the experimental design, it is important to set blind=FALSE 
    #| for downstream analysis.
    
    #| NOTE FOR VST: Previous QC chunk will provide hints on which transformation 
    #| is best to apply to the data (VST with blind=FALSE is recommended by default). 
    #| In order to test for differential expression, we operate on raw counts and 
    #| use discrete distributions. However for other downstream analyses - e.g. 
    #| for visualization or clustering, or WGCNA (see manual) it might be useful 
    #| to work with transformed versions of the count data.
    
    vsd.blind <- vst(dds, blind=TRUE) 
    vsd.noblind <- vst(dds, blind=FALSE)
    
  } else { #| For small datasets (< 30 samples), apply rlog.
    
    vsd.blind <- rlog(dds,blind=TRUE)
    vsd.noblind <- rlog(dds, blind=FALSE)
  }
  
  #Get normalized counts just for comparison with VST results
  dds <- estimateSizeFactors(dds)
  log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
  ntd <- normTransform(dds)
  tag <- "Basurto"
  pdf(paste("Results/VST_Counts/VstPlots/0_CountVariance_Stabilization_logNormalized_", tag,".pdf",sep=''));
  meanSdPlot(log.norm.counts, ranks=FALSE) 
  dev.off()
  pdf(paste("Results/VST_Counts/VstPlots/0_CountVariance_Stabilization_VST_BlindTRUE_", tag,".pdf",sep=''));
  meanSdPlot(assay(vsd.blind), ranks=FALSE)
  dev.off()
  pdf(paste("Results/VST_Counts/VstPlots/0_CountVariance_Stabilization_VST_BlindFALSE_", tag,".pdf",sep=''));
  meanSdPlot(assay(vsd.noblind), ranks=FALSE)
  dev.off()
  
  meanSdPlot(assay(ntd))
  
  plot <- meanSdPlot(assay(ntd))
  plot$gg +
    theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle = 30, hjust=1),axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
    scale_fill_gradient(low = "midnightblue", high = "yellow") 
  ggsave("Results/VST_Counts/VstPlots/Normalized_log2_no_VST.pdf", heigh=5, width = 6)  
  
  plot2 <- meanSdPlot(assay(vsd.blind), ranks=FALSE)
  plot2$gg +
    theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle = 30, hjust=1),axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
    scale_fill_gradient(low = "midnightblue", high = "yellow") 
  ggsave(paste("Results/VST_Counts/VstPlots/0_CountVariance_Stabilization_VST_BlindTRUE_", tag,".pdf",sep=''), heigh=5, width = 6)  
  
  #| NOW PCA with variance stabilization with blind=TRUE and FALSE. Blind =TRUE will
  #|consider the interindividual variance, so that the group information will not be 
  #|taken into account. After evaluation of the variance in the PCA, we will decide if 
  #|we should change blind=FALSE PCA
  #z <- plotPCA(vsd.blind, intgroup = "condition")  #| Specify condition
  #nudge <- position_nudge(y = 1.5)
  #pdf(paste("Results/0_QC_PCA_VST_blindTRUE_",tag,".pdf",sep=''))
  #plot(z + geom_text(aes(label = name), position = nudge))
  #dev.off()
  #z <- plotPCA(vsd.noblind)
  #nudge <- position_nudge(y = 1.5)
  #pdf(paste("Results/0_QC_PCA_VST_blindFALSE_",tag,".pdf", sep=""))
  #plot(z + geom_text(aes(label = name), position = nudge))
  #dev.off()
  
  #Heatmap
  m <- assay(vsd.blind)
  m.cor <- cor(m)
  pdf(paste("Results/VstPlots/0_QC_Heatmap_VST_blindTRUE_",tag,".pdf", sep=""))
  pheatmap(m.cor)
  dev.off()
  m <- assay(vsd.noblind)
  m.cor <- cor(m)
  pdf(paste("Results/VstPlots/0_QC_Heatmap_VST_blindFALSE_",tag,".pdf"))
  pheatmap(m.cor)
  dev.off()
}

#| Saving transformed matrix (Filtered by phenotype):
write.table(assay(vsd.blind),"Results/VST_Counts/FullCounts_Basurto_counts_filtered_phenotype_VST_PTEN_loss_intact.txt", sep ="\t", row.names =T)
################################################################################
