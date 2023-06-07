################################################################################
#                  COX ANALYSIS ON DIFFERENT DATASET                           #
################################################################################

#| This script perform cox regression analysis for accessing the prognosis capability
#| of level of expression of genes on the time of recurrence from different datasets 
#| including Basurto, Taylor, TCGA and Glinsky.

#| The Taylor, TCGA and Glinsky datasets were downloaded from CANCERTOOL database
#| http://genomics.cicbiogune.es/CANCERTOOL/datasetsInfo.php
################################################################################

################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))  
suppressMessages(library(survminer)) 
suppressMessages(library(ggfortify))
suppressMessages(library(viridis))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(biomaRt))
suppressMessages(library(DBI))
suppressMessages(library(RMySQL))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
################################################################################


################################################################################
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/"
setwd(workingDir)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Setting themes for plotting
theme_set(theme_classic())

#| Connecting to MySql database
con <- dbConnect(RMySQL::MySQL(), user='geneticanalyses', password = 'geneticanalyses', host = 'binf-web.cicbiogune.int',port=13306, dbname='geneticanalyses')
dbListTables(con)

#| Threshold of hazard ratio
worse <- 1.5
better <- 0.5

#| padj value threshold (p-value threshold for taylor)
p_threshold <- 0.05
################################################################################


################################################################################
################################# DATASETS #####################################
################################################################################


################################################################################
#######| Basurto 
################################################################################

#| Location of the files
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE_cpm_filtered_PCa_BPH.txt"

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Readtable Full counts. NOTE: it is not necessary to do any normalization 
data_expression_basurto <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Readtable Sample information
data_phenotype_basurto <- read.table(info.file, sep ="\t")

#| Filtering BPH
data_phenotype_basurto <- data_phenotype_basurto[data_phenotype_basurto$Diagnostico == "PCa" & !is.na(data_phenotype_basurto$DFS.STATUS),]
data_expression_basurto <- data_expression_basurto[ ,rownames(data_phenotype_basurto)]

#| DESEQ2 NORMALIZATION: FOR ACCOUNTING WITHIN GENE COMPARISON BETWEEN SAMPLES 

#| Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data_expression_basurto, colData = data_phenotype_basurto, design = ~ 1)

#| Estimating the size factors
dds <- estimateSizeFactors(dds)

#| Computing normalization
data_expression_basurto <- counts(dds, normalized=TRUE)

#| log2 for reescaling factors
data_expression_basurto <- log(data_expression_basurto + 1, base =2)

#| Box plot
boxplot(data_expression_basurto)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_basurto <- c()
coef_gene_basurto <- c()
exp_coef_gene_basurto <- c()

#| For over all the genes
for (i in 1:dim(data_expression_basurto)[1]){
  
  #| Computing the z-score for gene i
  counts_zscore <- as.numeric(data_expression_basurto[i,])
  counts_zscore <- (counts_zscore - mean(counts_zscore))/sd(counts_zscore)
  
  #| Adding the expression value of gene i to data_phenotype_basurto table
  data <- data.frame(data_phenotype_basurto, gene = counts_zscore)
  
  #| Fixing the cox model
  fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ gene + Edad_Zscore + DV200_Zscore, data=data)
  
  #| Summary of the results
  fit_summary <- summary(fit)
  
  #| Computing the p-value
  p_value_basurto <- c(p_value_basurto, fit_summary$coefficients[1, 5])
  
  #| Saving coeff
  coef_gene_basurto <- c(coef_gene_basurto, fit_summary$coefficients[1,1])
  
  #| Extracting and saving the values of the estimated exponential coef
  exp_coef_gene_basurto <- c(exp_coef_gene_basurto, fit_summary$coefficients[1,2])
}

#| Plotting distribution of p-values
data_pvalues_basurto <- data.frame(pvalues = p_value_basurto)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_Basurto.pdf")
ggplot(data_pvalues_basurto, aes(x = p_value_basurto)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from Cox analysis Basurto")
dev.off()

#| False discovery rate
fdr_basurto <- p.adjust(p_value_basurto, method = "fdr")

#| Unique data frame
dataframe_basurto <- data.frame(GeneID = rownames(data_expression_basurto), p_values_basurto = p_value_basurto, padj_basurto = fdr_basurto, coef_basurto = coef_gene_basurto, exp_coef_basurto = exp_coef_gene_basurto)

#| Merging for having gene names
dataframe_basurto <- merge(dataframe_basurto, genome_GRCh39.94, by="GeneID")

#| Sorting by p-values value
dataframe_basurto <- dataframe_basurto[order(dataframe_basurto$padj), ]

#| Exploring the dataframes
head(dataframe_basurto)

#| Saving results from Cox analysis in Basurto
write.table(dataframe_basurto[order(dataframe_basurto$GeneID), ],"Results/Tables/Cox_Analysis_Genes_DFSTIME_Basurto.txt" ,sep ="\t")
################################################################################




################################################################################
#######| Taylor 
################################################################################

#| Extracting the data
data_phenotype_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorphenotype")
data_expression_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorexpression")
data_annotation_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorannotation")

#| Changing colnames of data_annotation_taylor
colnames(data_annotation_taylor) <- c("row_names", "gene_name", "refseq_mrna")

#| Assigning the gene name to rownames
rownames(data_expression_taylor)<- data_expression_taylor$Gene

#| Filtering normal prostate
data_phenotype_taylor <- data_phenotype_taylor[ which(data_phenotype_taylor$tumortype == "Primary Tumor"),]
data_expression_taylor <- data_expression_taylor[, data_phenotype_taylor$Sample] 

#| Checking order and sample ID
length(colnames(data_expression_taylor)) == length(data_phenotype_taylor$Sample)
any(colnames(data_expression_taylor) == data_phenotype_taylor$Sample)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_taylor <- c()
coef_gene_taylor <- c()
exp_coef_gene_taylor <- c()

#| For over all the genes
for (i in 1:dim(data_expression_taylor)[1]){
  
  #| Computing the z-score for gene i
  counts_zscore <- as.numeric(data_expression_taylor[i,])
  counts_zscore <- (counts_zscore - mean(counts_zscore))/sd(counts_zscore)
  
  #| Adding the expression value of gene i to sample_info_S table
  data <- data.frame(data_phenotype_taylor, gene = counts_zscore)
  
  #| Fixing the cox model
  fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ gene, data=data)
  
  #| Summary of the results
  fit_summary <- summary(fit)
  
  #| Computing the p-value
  p_value_taylor <- c(p_value_taylor, fit_summary$coefficients[1, 5])
  
  #| Saving coeff
  coef_gene_taylor <- c(coef_gene_taylor, fit_summary$coefficients[1,1])
  
  #| Extracting and saving the values of the estimated exponential coef
  exp_coef_gene_taylor <- c(exp_coef_gene_taylor, fit_summary$coefficients[1,2])
}

#| Plotting p-values
data_pvalues_taylor <- data.frame(pvalues = p_value_taylor)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_Taylor.pdf")
ggplot(data_pvalues_taylor, aes(x = p_value_taylor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from Cox analysis Taylor")
dev.off()

#| False discovery rate
fdr_taylor <- p.adjust(p_value_taylor, method = "fdr")

#| Unique data frame
dataframe_taylor <- data.frame(refseq_mrna = rownames(data_expression_taylor), p_values_taylor = p_value_taylor, padj_taylor = fdr_taylor, coef_taylor = coef_gene_taylor, exp_coef_taylor = exp_coef_gene_taylor)

#| GeneID
dataframe_taylor <- merge(dataframe_taylor, data_annotation_taylor, by = "refseq_mrna")

#| Sorting by p-values value
dataframe_taylor <- dataframe_taylor[order(dataframe_taylor$p_values), ]

#| Exploring the dataframes
head(dataframe_taylor)
################################################################################



################################################################################
#######| TCGA Fire Horse Legacy data
################################################################################
data_phenotype_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostatephenotype")
data_expression_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostateexpression")

#| Zscore the Age 
data_phenotype_tcga$Age_Zscore <- (as.numeric(data_phenotype_tcga$Age) - mean(as.numeric(data_phenotype_tcga$Age)))/sd(as.numeric(data_phenotype_tcga$Age))

#| Fixing gene name
data_expression_tcga$Gene <- gsub(".{1}$","",data_expression_tcga$Gene)
rownames(data_expression_tcga) <- data_expression_tcga$Gene

#| Deleting Gene column
data_expression_tcga <- data_expression_tcga[,-c(1)]

#| Checking order and sample ID
length(colnames(data_expression_tcga[,1:dim(data_expression_tcga)[2]])) == length(data_phenotype_tcga$Sample)
any(colnames(data_expression_tcga[,1:dim(data_expression_tcga)[2]]) == data_phenotype_tcga$Sample)

#| Finding low counts
low_counts <- rownames(data_expression_tcga)[which(rowSums(data_expression_tcga > 0) > 0.3*dim(data_expression_tcga)[2] )]

#| Filtering low counts counts
data_expression_tcga <- data_expression_tcga[ which(rownames(data_expression_tcga) %in% low_counts),]

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_tcga <- c()
coef_gene_tcga <- c()
exp_coef_gene_tcga <- c()

#| For over all the genes
for (i in 1:dim(data_expression_tcga)[1]){
 
  #| Computing the z-score for gene i
  counts_zscore <- as.numeric(data_expression_tcga[i,])
  counts_zscore <- (counts_zscore - mean(counts_zscore))/sd(counts_zscore) 
  
  #| Adding the expression value of gene i to sample_info_S table
  data <- data.frame(data_phenotype_tcga, gene = counts_zscore)

  #| Fixing the cox model
  fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ gene + Age_Zscore, data=data)
  
  #| Summary of the results
  fit_summary <- summary(fit)
  
  #| Computing the p-value
  p_value_tcga <- c(p_value_tcga, fit_summary$coefficients[1, 5])
  
  #| Saving coeff
  coef_gene_tcga <- c(coef_gene_tcga, fit_summary$coefficients[1,1])
  
  #| Extracting and saving the values of the estimated exponential coef
  exp_coef_gene_tcga <- c(exp_coef_gene_tcga, fit_summary$coefficients[1,2])
}

#| Plotting p-values
data_pvalues_tcga <- data.frame(pvalues = p_value_tcga)
pdf("Results/Images/p_values_histogram_Density_TCGA.pdf")
ggplot(data_pvalues_tcga, aes(x = p_value_tcga)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from Cox analysis TCGA")
dev.off()

#| False discovery rate
fdr_tcga <- p.adjust(p_value_tcga, method = "fdr")

#| Unique data frame
dataframe_tcga <- data.frame(gene_name = rownames(data_expression_tcga), p_values_tcga = p_value_tcga, padj_tcga = fdr_tcga, coef_tcga = coef_gene_tcga, exp_coef_tcga = exp_coef_gene_tcga)

#| Sorting by p-values value
dataframe_tcga <- dataframe_tcga[order(dataframe_tcga$padj), ]

#| Exploring the dataframes
head(dataframe_tcga)
################################################################################


################################################################################
#######| Glinsky
################################################################################
data_phenotype_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyphenotype")
data_expression_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyexpression")

#| Zscore the Age 
data_phenotype_glinsky$Age_Zscore <- (as.numeric(data_phenotype_glinsky$Age) - mean(as.numeric(data_phenotype_glinsky$Age)))/sd(as.numeric(data_phenotype_glinsky$Age))

#| Giving gene name to row names
rownames(data_expression_glinsky) <- data_expression_glinsky$Gene

#| Deleting column "Gene" form gene expression data
data_expression_glinsky <- data_expression_glinsky[,-c(1)]

#| Checking order and sample ID
length(colnames(data_expression_glinsky[,1:dim(data_expression_glinsky)[2]])) == length(data_phenotype_glinsky$Sample)
any(colnames(data_expression_glinsky[,1:dim(data_phenotype_glinsky)[2]]) == data_phenotype_glinsky$Sample)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_glinsky <- c()
coef_gene_glinsky <- c()
exp_coef_gene_glinsky <- c()

#| For over all the genes
for (i in 1:dim(data_expression_glinsky)[1]){
  
  #| Computing the z-score for gene i
  counts_zscore <- as.numeric(data_expression_glinsky[i,])
  counts_zscore <- (counts_zscore - mean(counts_zscore))/sd(counts_zscore) 
  
  #| Adding the expression value of gene i to sample_info_S table
  data <- data.frame(data_phenotype_glinsky, gene = counts_zscore)
  
  #| Fixing the cox model
  fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ gene + Age_Zscore , data=data)
  
  #| Summary of the results
  fit_summary <- summary(fit)
  
  #| Computing the p-value
  p_value_glinsky <- c(p_value_glinsky, fit_summary$coefficients[1, 5])
  
  #| Saving coeff
  coef_gene_glinsky <- c(coef_gene_glinsky, fit_summary$coefficients[1,1])
  
  #| Extracting and saving the values of the estimated exponential coef
  exp_coef_gene_glinsky <- c(exp_coef_gene_glinsky, fit_summary$coefficients[1,2])
}

#| Plotting p-values
data_pvalues_glinsky <- data.frame(pvalues = p_value_glinsky)
pdf("Results/Images/p_values_histogram_Density_Glinsky.pdf")
ggplot(data_pvalues_glinsky, aes(x = p_value_glinsky)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from Cox analysis Glinsky")
dev.off()

#| False discovery rate
fdr_glinsky <- p.adjust(p_value_glinsky, method = "fdr")

#| Unique data frame
dataframe_glinsky <- data.frame(gene_name = rownames(data_expression_glinsky), p_values_glinsky = p_value_glinsky, padj_glinsky = fdr_glinsky, coef_glinsky = coef_gene_glinsky, exp_coef_glinsky = exp_coef_gene_glinsky)

#| Sorting by p-values value
dataframe_glinsky <- dataframe_glinsky[order(dataframe_glinsky$p_values), ]

#| Exploring the dataframes
head(dataframe_glinsky)
################################################################################



################################################################################
#| Saving results results in table
################################################################################

#| Saving all dataframes in an unique list
list_dataframes <- list(dataframe_basurto[,-c(1,6,7,8,9,10)], dataframe_taylor[,-c(1,6)], dataframe_tcga, dataframe_glinsky)

#| Merging them up
final_table <- Reduce(function(x, y) merge(x, y, all=TRUE), list_dataframes) 

nonNAN_datasets <- c()
sig_datasets <- c()
worse_sig_datasets <- c()
better_sig_datasets <- c()
total_coherence_halfdatasets <- c()

for (i in 1:dim(final_table)[1]){

  p_values <- as.numeric(final_table[i,c("p_values_basurto", "p_values_taylor", "p_values_tcga", "p_values_glinsky")])
  exp_coef <- as.numeric(final_table[i,c("exp_coef_basurto", "exp_coef_taylor", "exp_coef_tcga", "exp_coef_glinsky")])
  
  #| non NAN
  nonNAN_datasets <- c(nonNAN_datasets, length(p_values[which(!is.na(p_values))]))
  
  #| In how many datasets the gene is significant
  sig_datasets <- c(sig_datasets, length(p_values[which(p_values < 0.05)]))
  
  #| Worse significant
  worse_sig_datasets <- c(worse_sig_datasets, length(exp_coef[which(p_values < 0.05 & exp_coef > 1)]))
  
  #| Better significant
  better_sig_datasets <- c(better_sig_datasets, length(exp_coef[which(p_values < 0.05 &  exp_coef< 1)]))
  
  #| Total coherence
  if(length(exp_coef[which(p_values < 0.05 & exp_coef > 1)]) >= ceiling(length(p_values[which(!is.na(p_values))])/2)){
      total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "YES_WORSE")
  }else{
    if(length(exp_coef[which(p_values < 0.05 & exp_coef < 1)]) >= ceiling(length(p_values[which(!is.na(p_values))])/2)){
      total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "YES_BETTER")
    }else{
      total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "NO")
    }
  }
}

#| cbinding columns
final_table <- cbind(final_table, nonNAN_datasets, sig_datasets, worse_sig_datasets, better_sig_datasets, total_coherence_halfdatasets)

#| Saving final table
write.table(final_table, "Results/Tables/Final_table_Cox_Analysis_Genes_DFSTIME_Basurto_Taylor_TCGA_Glinsky.txt", sep = "\t")
################################################################################



################################################################################
#| Hazard ratio as a function of the p-values for different datasets
################################################################################

#| Hazard ratio as a function of the p-value BASURTO
hazard_pvalue_basurto  <- data.frame( p_values = dataframe_basurto$padj_basurto, hazard = dataframe_basurto$exp_coef_basurto)
pdf("Results/Images/p_value_Hazard_ratio/p_values_Hazard_ratio_basurto.pdf")
ggplot(hazard_pvalue_basurto, aes(x=p_values, y = hazard)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab("Hazard ratio") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Hazard ratio vs p-adj-values Basurto")
dev.off()

##| Hazard ratio as a function of the p-value TAYLOR
hazard_pvalue_taylor <- data.frame( p_values = dataframe_taylor$p_values_taylor, hazard = dataframe_taylor$exp_coef_taylor)
pdf("Results/Images/p_value_Hazard_ratio/p_values_Hazard_ratio_taylor.pdf")
ggplot(hazard_pvalue_taylor, aes(x=p_values, y = hazard)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-values") +
  ylab("Hazard ratio") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
  ggtitle("Hazard ratio vs p-values Taylor")
dev.off()

#| Hazard ratio as a function of the p-value TCGA
hazard_pvalue_tcga <- data.frame( p_values = dataframe_tcga$padj_tcga, hazard = dataframe_tcga$exp_coef_tcga)
pdf("Results/Images/p_value_Hazard_ratio/p_values_Hazard_ratio_TCGA.pdf")
ggplot(hazard_pvalue_tcga, aes(x=p_values, y = hazard)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab("Hazard ratio") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
  ggtitle("Hazard ratio vs p-adj-values TCGA")
dev.off()

#| Hazard ratio as a function of the p-value GLINSKY
hazard_pvalue_glinsky <- data.frame( p_values = final_table$padj_glinsky, hazard = final_table$exp_coef_glinsky)
pdf("Results/Images/p_value_Hazard_ratio/p_values_Hazard_ratio_Glinsky.pdf")
ggplot(hazard_pvalue_glinsky, aes(x=p_values, y = hazard)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab("Hazard ratio") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Hazard ratio vs p-adj-values Glinsky")
dev.off()
################################################################################



################################################################################
#| Hazard ratio of every dataset
################################################################################

#| Basurto and Taylor
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_Basurto_Taylor.pdf")
dataframe_taylor$gene_name <- dataframe_taylor$hgnc_symbol
Hazard_Basurto_Taylor <- merge(dataframe_basurto, dataframe_taylor, by = "gene_name")
ggplot(Hazard_Basurto_Taylor, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio Basurto") +
  ylab("Hazard ratio Taylor") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

#| Basurto and TCGA
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_Basurto_TCGA.pdf")
dataframe_tcga$gene_name <- dataframe_tcga$GeneID
Hazard_Basurto_TCGA <- merge(dataframe_basurto, dataframe_tcga, by = "gene_name")
ggplot(Hazard_Basurto_TCGA, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio Basurto") +
  ylab("Hazard ratio TCGA") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

#| Basurto and Glinsky
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_Basurto_Glinsky.pdf")
dataframe_glinsky$gene_name <- dataframe_glinsky$GeneID
Hazard_Basurto_Glinsky <- merge(dataframe_basurto, dataframe_glinsky, by = "gene_name")
ggplot(Hazard_Basurto_Glinsky, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio Basurto") +
  ylab("Hazard ratio Glinsky") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

#| Taylor and TCGA
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_Taylor_TCGA.pdf")
Hazard_Taylor_TCGA <- merge(dataframe_taylor, dataframe_tcga, by = "gene_name")
ggplot(Hazard_Taylor_TCGA, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio Taylor") +
  ylab("Hazard ratio TCGA") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

#| Taylor and Glinsky
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_Taylor_Glinsky.pdf")
Hazard_Taylor_Glinsky <- merge(dataframe_taylor, dataframe_glinsky, by = "gene_name")
ggplot(Hazard_Taylor_Glinsky, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio Taylor") +
  ylab("Hazard ratio Glinsky") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

#| TCGA and Glinsky
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Images/hazard_hazard_plots/hazard_hazard_TCGA_Glinsky.pdf")
Hazard_TCGA_Glinsky <- merge(dataframe_tcga, dataframe_glinsky, by = "gene_name")
ggplot(Hazard_TCGA_Glinsky, aes(x=exp_coef.x, y = exp_coef.y)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("Hazard ratio TCGA") +
  ylab("Hazard ratio Glinsky") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()

################################################################################



################################################################################
#######| Intersecting
################################################################################


################################################################################
#| Worse and Better prognosis from every dataset
################################################################################


#| Worse
worse_basurto <- final_table[final_table$exp_coef_basurto> worse & final_table$padj_basurto<p_threshold,]

worse_taylor <- dataframe_taylor[dataframe_taylor$exp_coef >worse & dataframe_taylor$p_values < p_threshold,]
worse_tcga <- dataframe_tcga[dataframe_tcga$exp_coef>worse & dataframe_tcga$padj<p_threshold,]
worse_glinsky <- dataframe_glinsky[dataframe_glinsky$exp_coef>worse & dataframe_glinsky$padj<p_threshold,]

#| Better
better_basurto <- dataframe_basurto[dataframe_basurto$exp_coef<better & dataframe_basurto$padj<p_threshold,]
better_taylor <- dataframe_taylor[dataframe_taylor$exp_coef < better & dataframe_taylor$p_values < p_threshold,]
better_tcga <- dataframe_tcga[dataframe_tcga$exp_coef<better & dataframe_tcga$padj<p_threshold,]
better_glinsky <- dataframe_glinsky[dataframe_glinsky$exp_coef<better & dataframe_glinsky$padj<p_threshold,]



venn.diagram(x = list(worse_basurto$gene_name, worse_taylor$hgnc_symbol[!is.na(worse_taylor$hgnc_symbol)], worse_tcga$GeneID[!is.na(worse_tcga$GeneID)], worse_glinsky$GeneID),
             category.names = c("worse_basurto", "worse_taylor", "worse_tcga", "worse_glinsky"),
             filename = "VennProof_Worse_zscores.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c("green", "black", "yellow", "magenta"),
             fill = c(alpha("green",0.3), alpha("white",), alpha("yellow",0.3),alpha("magenta",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")


venn.diagram(x = list(better_basurto$gene_name, better_taylor$hgnc_symbol[!is.na(better_taylor$hgnc_symbol)],better_tcga$GeneID[!is.na(better_tcga$GeneID)], better_glinsky$GeneID),
             category.names = c("better_basurto", "better_taylor", "better_tcga", "better_glinsky"),
             filename = "VennProof_better_zscores.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c("green", "black", "yellow", "magenta"),
             fill = c(alpha("green",0.3), alpha("white",), alpha("yellow",0.3),alpha("magenta",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")


genes_worse <- intersect(intersect(worse_basurto$gene_name,worse_taylor$hgnc_symbol), worse_tcga$GeneID)
basurto <- c()
taylor <- c()
tcga <- c()

for (i in 1:length(genes_worse)){
  basurto <- c(basurto,worse_basurto$exp_coef[worse_basurto$gene_name == genes_worse[i]])
  taylor <- c(taylor, worse_taylor$exp_coef[worse_taylor$hgnc_symbol == genes_worse[i] & !is.na(worse_taylor$hgnc_symbol)])
  print(taylor)
  print(genes_worse[i])
  tcga <- c(tcga, worse_tcga$exp_coef[worse_tcga$GeneID == genes_worse[i] & !is.na(worse_tcga$GeneID)])
}


join_data <- data.frame( gene = genes_worse, Basurto = basurto, Taylor = taylor[-c(8,12)], TCGA = tcga)


genes_better <- intersect(intersect(better_basurto$gene_name, better_taylor$hgnc_symbol[!is.na(better_taylor$hgnc_symbol)]), better_tcga$GeneID)
basurto <- c()
taylor <- c()
tcga <- c()

for (i in 1:length(genes_better)){
  basurto <- c(basurto,better_basurto$exp_coef[better_basurto$gene_name == genes_better[i]])
  taylor <- c(taylor, better_taylor$exp_coef[better_taylor$hgnc_symbol == genes_better[i] & !is.na(better_taylor$hgnc_symbol)])
  print(taylor)
  print(genes_better[i])
  tcga <- c(tcga, better_tcga$exp_coef[better_tcga$GeneID == genes_better[i]])
}

length(as.numeric(data_expression_taylor[data_expression_taylor==dataframe_taylor$refseq_mrna[which(dataframe_taylor$hgnc_symbol == "SKA3")] & data_phenotype_taylor$,-c(1)]))
length(data_phenotype_taylor$DFS.TIME)


################################################################################
###| Plotting probability of survival by quantiles
################################################################################

i <- 15

#| Worse BASURTO
data_phenotype_basurto$quantile <- cut(as.numeric(data_expression_basurto[which(rownames(data_expression_basurto) == genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name==genes_worse[i])]),]),quantile(as.numeric(data_expression_basurto[which(rownames(data_expression_basurto) == genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name==genes_worse[i])]),])),include.lowest=TRUE,labels=FALSE)
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS)~quantile, data=data_phenotype_basurto)
pdf(paste("Results/Images/survival_plots/Basurto_Survival_plot_quantiles_gene_",genes_worse[i], ".pdf", sep =""))
ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           title=paste("Kaplan-Meier Curve gene ", genes_worse[i], " Basurto", sep =""),# add title to plot
           ggtheme = theme_survminer( base_family="serif"), 
           font.family = "serif"
)
dev.off()

#| Worse TAYLOR
data_phenotype_taylor$quantile <- cut(as.numeric(data_expression_taylor[which(rownames(data_expression_taylor) == dataframe_taylor$refseq_mrna[which(dataframe_taylor$hgnc_symbol == genes_worse[i])]),]),quantile(as.numeric(data_expression_taylor[which(rownames(data_expression_taylor) == dataframe_taylor$refseq_mrna[which(dataframe_taylor$hgnc_symbol == genes_worse[i])]),])),include.lowest=TRUE,labels=FALSE)
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS)~quantile, data=data_phenotype_taylor)
pdf(paste("Results/Images/survival_plots/Taylor_Survival_plot_quantiles_gene_",genes[i], ".pdf", sep =""))
ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           title=paste("Kaplan-Meier Curve gene ", genes_worse[i], " Taylor", sep =""),# add title to plot
           ggtheme = theme_survminer( base_family="serif"), 
           font.family = "serif"
)
dev.off()

#| Worse TCGA
data_phenotype_tcga$quantile <- cut(as.numeric(as.numeric(data_expression_tcga[which(rownames(data_expression_tcga) == genes_worse[i]),])),quantile(as.numeric(data_expression_tcga[which(rownames(data_expression_tcga) == genes_worse[i]),])), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS)~quantile, data=data_phenotype_tcga)
pdf(paste("Results/Images/survival_plots/TCGA_Survival_plot_quantiles_gene_",genes_worse[i], ".pdf", sep =""))
ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           title=paste("Kaplan-Meier Curve gene ", genes_worse[i], " TCGA", sep =""),# add title to plot
           ggtheme = theme_survminer( base_family="serif"), 
           font.family = "serif"
)
dev.off()

