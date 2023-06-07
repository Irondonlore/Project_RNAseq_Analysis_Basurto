################################################################################
#### PCA, HISTOGRAMS, VARIANCE  PARTITION AND CORRELATION ANALYSIS OF      #####
####                         THE AC-45_RNAseq_FFPE                         #####
################################################################################

#| This script is written to perform some statistical analysis such as PCA, histograms
#| correlation and variance partition analysis as a first approach to the study 
#| of the RNAseq data of 240 patients of Basurto's hospital which information is
#| contained in the Big_data folder at X:/DATA_shared/AC-45_RNAseq-FFPE 

#| Moreover, at the beginning was performed an analysis to test the variability of 
#| the variance as a function of the mean. To avoid this dependence, a vst transformation
#| using DESeq2 was applying to to counts matrix.

#| The source for the variance partition analysis can be found in the following link:
#| https://www.bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf

#| Quick definition of variancePartition analysis in R:
#| variancePartition fits a linear (mixed) model that jointly considers the 
#| contribution of all specified variables on the expression of each gene.
#| the results of variancePartition give insight into the expression data at multiple
#| levels to identify which variable is the responsible for the variation in the 
#| expression. Moreover:
#| 1) variancePartition provides a natural interpretation of multiple variables
#| 2) variancePartition quantifies the contribution of each variable
#| 3) variancePartition interprets contribution of each variable to each gene
#| 4) variancePartition can assess contribution of one variable (i.e. Individual) 
#|    separately in subset of the data defined by another variable (i.e. Tissue)

#| When variancePartition is applied to RNA-seq expression data we can combine it
#| with DESeq2.One of the advantages of using first DESeq2 is that it deals with a
#| first processing on the data. 
#| the data

#| What to do with missing values? (NA/NaN/Inf)
#| variancePartition fits a regression model for each gene and drops samples that 
#| have NA/NaN/Inf values in each model fit. 

#| Which variables should be included?
#| It is useful to include all variables in the first analysis and then drop variables 
#| that have minimal effect. But be careful with those variables who are highly  
#| correlated among them and also be careful about the batch effects (technical issues).

#| Assessing the correlation between pairs of variables:
#| Evaluating the correlation between variables in a important part in interpreting
#| variancePartition results. When comparing two continuous variables, Pearson
#| correlation is widely used. But variancePartition includes categorical variables 
#| in the model as well. In order to accommodate the correlation between a continuous 
#| and a categorical variable, or two categorical variables we used canonical
#| correlation analysis.

#| Canonical Correlation Analysis (CCA) is similar to correlation between two vectors, 
#| except that CCA can accommodate matricies as well. For a pair of variables, 
#| canCorPairs assesses the degree to which they co-vary and contain the same 
#| information. Variables in the formula can be a continuous variable or a discrete 
#| variable expanded to a matrix

#| TMA  
################################################################################


################################ LIBRARIES #####################################
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

#| For plots
theme_set(theme_classic()) 
################################################################################


############################## DATA DIRECTORIES ################################
dir.proj <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/"
info.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/variancePartitioninput.txt"
counts.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_counts_phenotype_Normalized_Log2.txt"
foxo.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TTGTTT_V_FOXO4_01.txt"
foxo_pathway.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/FOXO_Pathway.txt"
pi3k_akt_mtor.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt"
pten_loss.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Signature_PTEN_loss_Paper_Imada_et_al_BMC_cancer.xlsx"
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("W:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################



############################### READING DATA ###################################
#| Counts data (normalized)
counts_data <- read.table(counts.file, sep ="\t", header=T)

#| Sample information
sample_info <- read.table(info.file, sep ="\t", header=T)
sample_info_2 <- read_xlsx("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Basurto_BBDD_AC45_Patients_pheno.xlsx")
################################################################################




############################# PI3K-PTEN ACTIVITY ###############################
#| AKT levels and Diagnosis
pdf("Results/Box_plots/AKT1_expression_vs_PCa_and_BPH.pdf")
ggplot(sample_info_extracted, aes(x=Diagnostico, y=AKT1_Exp_log2, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  
  ggtitle("AKT1 Expression vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + 
  xlab("Diagnosis") + 
  labs(fill ="Diagnosis") + 
  ylab("AKT1 Expression (Normalized log2 counts)")
dev.off()

#| mTOR levels and Diagnosis
pdf("Results/Box_plots/mTOR_expression_vs_PCa_and_BPH.pdf")
ggplot(sample_info_extracted, aes(x=Diagnostico, y=mTOR_Exp_log2, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  
  ggtitle("mTOR Expression vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + 
  xlab("Diagnosis") + 
  labs(fill ="Diagnosis") + 
  ylab("mTOR Expression (Normalized log2 counts)")
dev.off()

#| FOXO1 levels and Diagnosis
pdf("Results/Box_plots/FOXO1_expression_vs_PCa_and_BPH.pdf")
ggplot(sample_info_extracted, aes(x=Diagnostico, y=FOXO1_Exp_log2, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  
  ggtitle("FOXO1 Expression vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + 
  xlab("Diagnosis") + 
  labs(fill ="Diagnosis") + 
  ylab("FOXO1 Expression (Normalized log2 counts)")
dev.off()

#| AKT1 and PTEN expression 
cor <- cor.test(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,18)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_vs_AKT1_Expression_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= AKT1_Exp_log2)) + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN vs AKT1 Expression\n 198 PCa samples") +  
  ylab("AKT1 Expression (Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=10.5, y=8.9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| AKT1 and PTEN expression (Tumor purity)
cor <- cor.test(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,18)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_vs_AKT1_Expression_198_PCa_samples_purity.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= AKT1_Exp_log2, color = purity)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN vs AKT1 Expression\n 198 PCa samples") +  
  ylab("AKT1 Expression (Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=10.5, y=9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D") +
  labs(color="Purity")
dev.off()

#| AKT1 and PTEN expression (Mean expression FOXO signature)
cor <- cor.test(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,18)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_vs_AKT1_Expression_198_PCa_samples_Mean_expression_FOXO_signature.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= AKT1_Exp_log2, color = `Mean expression FOXO`)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN vs AKT1 Expression\n 198 PCa samples") +  
  ylab("AKT1 Expression (Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=10.5, y=9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D") +
  labs(color="FOXO Signature")
dev.off()

#| AKT1 and FOXO1 expression 
cor <- cor.test(sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$FOXO1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/AKT1_vs_FOXO1_Expression_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= AKT1_Exp_log2, y= FOXO1_Exp_log2)) + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 vs FOXO1 Expression\n 198 PCa samples") +  
  ylab("FOXO1 Expression (Normalized log2 counts)")+ 
  xlab("AKT1 Expression (Normalized log2 counts)")+
  annotate(geom="text", x=9.1, y=8.7, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| AKT1 and FOXO1 expressionn (Tumor purity)
cor <- cor.test(sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$FOXO1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/AKT1_vs_FOXO1_Expression_198_PCa_samples_purity.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= AKT1_Exp_log2, y= FOXO1_Exp_log2, color = purity)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 vs FOXO1 Expression\n 198 PCa samples") +  
  ylab("FOXO1 Expression (Normalized log2 counts)")+ 
  xlab("AKT1 Expression (Normalized log2 counts)")+
  annotate(geom="text", x=9.1, y=8.7, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  scale_color_viridis(discrete=F, option="D") +
  labs(color="Purity")
dev.off()


#| AKT1 and PTEN expression (Mean expression FOXO signature)
cor <- cor.test(sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$FOXO1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/AKT1_vs_FOXO1_Expression_198_PCa_samples_Mean_expression_FOXO_signature.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= AKT1_Exp_log2, y= FOXO1_Exp_log2, color = `Mean expression FOXO`)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 vs FOXO1 Expression\n 198 PCa samples") +  
  ylab("FOXO1 Expression (Normalized log2 counts)")+ 
  xlab("AKT1 Expression (Normalized log2 counts)")+
  annotate(geom="text", x=9.1, y=8.7, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  scale_color_viridis(discrete=F, option="D") +
  labs(color="FOXO Signature")
dev.off()


#| PTEN and FOXO1 expression 
cor <- cor.test(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$FOXO1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,13)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_vs_FOXO1_Expression_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= FOXO1_Exp_log2)) + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN vs FOXO1 Expression\n 198 PCa samples") +  
  ylab("FOXO1 Expression (Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=12, y=8.6, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()



#| AR and FOXO1 expression 
cor <- cor.test(sample_info_extracted$AR_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$FOXO1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,13)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/AR_vs_FOXO1_Expression_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= AR_Exp_log2, y= FOXO1_Exp_log2)) + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AR vs FOXO1 Expression\n 198 PCa samples") +  
  ylab("FOXO1 Expression (Normalized log2 counts)")+ 
  xlab("AR Expression (Normalized log2 counts)")+
  annotate(geom="text", x=11, y=8.6, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| AR and AKT
cor <- cor.test(sample_info_extracted$AR_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$AKT1_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/AR_vs_AKT1_Expression_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= AR_Exp_log2, y= AKT1_Exp_log2)) + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AR vs AKT1 Expression\n 198 PCa samples") +  
  ylab("AKT1 Expression (Normalized log2 counts)")+ 
  xlab("AR Expression (Normalized log2 counts)")+
  annotate(geom="text", x=11, y=8.6, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()


summary(sample_info_extracted$purity)

plot(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico =="PCa")], sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico =="PCa")])


ggplot(sample_info_extracted, aes(x= AR_Exp_log2, y= AKT1_Exp_log2, color =Diagnostico)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="green") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AR vs AKT1 Expression\n 198 PCa samples") +  
  ylab("AKT1 Expression (Normalized log2 counts)")+ 
  xlab("AR Expression (Normalized log2 counts)")+
  annotate(geom="text", x=11, y=8.6, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
################################################################################





################################################################################
#| Sample correlation
################################################################################
data <- assay(vsd.blind)

dmat <- as.matrix(cor(data,method="spearman"))

#| PCa and BPH
pdf("Results/Sample_Correlation/SampleCorrelation_PCa_and_BPH_counts_data_filtered_phenotype_VST.pdf")
pheatmap(dmat,border_color=NA,annotation_legend=T, fontsize = 6, fontsize_col = 1.5, fontsize_row = 1.5)
dev.off()

#| PCa
PCa_samples_dat <- subset(data, select = PCa_samples)
dmat_PCa <- as.matrix(cor(PCa_samples_dat,method="spearman"))

pdf("Results/Sample_Correlation/SampleCorrelation_PCa_counts_filtered_phenotype_VST.pdf")
pheatmap(dmat_PCa,border_color=NA,annotation_legend=T,fontsize = 6, fontsize_col = 1.5, fontsize_row = 1.5)
dev.off()

#| BPH
BPH_samples_dat <- subset(data, select = BPH_samples)
dmat_BPH <- as.matrix(cor(BPH_samples_dat,method="spearman"))

pdf("Results/Sample_Correlation/SampleCorrelation_BPH_counts_data_filtered_phenotype_VST.pdf")
pheatmap(dmat_BPH,border_color=NA,annotation_legend=T,fontsize = 6, fontsize_col = 4, fontsize_row = 4)
dev.off()
################################################################################




################################################################################
#| Gene dispersion
################################################################################
d <- DESeq2::estimateDispersions(dds)
pdf("Results/Dispersions/DispersionPlot_data_counts_CPM_filtered_phenotype_PCa_BPH.pdf")
plotDispEsts(d)
dev.off()
################################################################################




################################################################################
#| Clustering
################################################################################



################################################################################
#| PCA: By sample groups
################################################################################

#| Tranforming the data (using the same variable in case the lines above we not run)
vsd.blind <- vst(dds, blind=TRUE)

#| To fix theme
theme_set(theme_classic()) 

#| 1) intgroup = "PCa"
pcaData <- plotPCA(vsd.blind, intgroup= "Diagnostico", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#png("Results/PCA/PCA_PCa.png", width = 780, height = 580,pointsize = 15)
pdf("Results/PCA/PCA_Diagnosis.pdf")
ggplot(pcaData, aes(PC1, PC2, color=Diagnostico))+ geom_point(size=2.5) + 
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(color='Diagnosis')  +
  coord_fixed()  +  ggtitle("PCA: PCa vs BPH") + scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) 
dev.off()

#| intgroup = c("PTEN")
pcaData <- plotPCA(vsd.blind, intgroup= c("PTEN_Exp_Zscore"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_PTEN.pdf")
ggplot(pcaData, aes(PC1, PC2, color=PTEN_Exp_Zscore)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: PTEN Expression")  + scale_color_viridis(discrete=F, option="A")
dev.off()

#| intgroup = c("Diagnostico","PTEN_Exp")
pcaData <- plotPCA(vsd.blind, intgroup= c("Diagnostico","PTEN_Exp_Zscore"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_Diagnostico_PTEN.pdf")
ggplot(pcaData, aes(PC1, PC2, color=PTEN_Exp_Zscore, shape= Diagnostico)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: PTEN Expression and Diagnostic")  + scale_color_viridis(discrete=F, option="A") +
  labs(color="PTEN Expression", fill=  )
dev.off()

#| intgroup = c("Diagnostico","Edad")
pcaData <- plotPCA(vsd.blind, intgroup= c("Diagnostico","Edad"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_Diagnostico_Edad.pdf")
ggplot(pcaData, aes(PC1, PC2, color=Edad, shape= Diagnostico)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: Edad and Diagnostic")  + scale_color_viridis(discrete=F, option="A")
dev.off()

#| intgroup = c("DV200")
pcaData <- plotPCA(vsd.blind, intgroup= c("DV200.value"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_DV200.pdf")
ggplot(pcaData, aes(PC1, PC2, color=DV200.value)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: DV200.value")  + scale_color_viridis(discrete=F, option="A")
dev.off()

#| intgroup = c("BMI")
pcaData <- plotPCA(vsd.blind, intgroup= c("BMI"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_BMI.pdf")
ggplot(pcaData, aes(PC1, PC2, color=BMI)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: BMI")  + scale_color_viridis(discrete=F, option="A")
dev.off()

#| intgroup = c("PSA")
pcaData <- plotPCA(vsd.blind, intgroup= c("PSA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("Results/PCA/PCA_PSA.pdf")
ggplot(pcaData, aes(PC1, PC2, color=PSA)) +
  geom_point(size=2.5) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()  +  ggtitle("PCA: Prostate Specific Antigen (PSA)")  + scale_color_viridis(discrete=F, option="A")
dev.off()
################################################################################




################################################################################
#| Histogram
################################################################################
pdf("Results/Histograms/Density_Edad_Diagnostico.pdf")
ggplot(sample_info_extracted, aes(x=Edad, fill =Diagnostico)) + geom_density(alpha=0.8)+
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Density plot: Age and Diagnostic") + xlab("Age") + ylab("Density") +  labs(fill='Diagnosis')
#+geom_text(x=50, y=0.06, label="62", family = "serif",size=3)
dev.off()

pdf("Results/Histograms/Density_PTEN_Diagnostico.pdf")
ggplot(sample_info_extracted, aes(x=PTEN_Exp_Zscore, fill=Diagnostico)) + geom_density(alpha=0.8)+
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Density plot: PTEN expression and Diagnosis") + xlab("PTEN Expression") +  ylab("Density") +labs(fill='Diagnosis')
dev.off()

pdf("Results/Histograms/Density_DV200.pdf")
ggplot(sample_info_extracted, aes(x=DV200.value)) + geom_density(fill="midnightblue", alpha=0.9)+
  scale_fill_manual(values = "aquamarine2")+theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Density plot: DV200") + xlab("DV200") +  ylab("Density") 
dev.off()





################################################################################
#| Variance Partition Analysis
################################################################################
#| Defining the formula
sample_info_extracted$DFS.TIME <- as.numeric(sample_info_extracted$DFS.TIME)
formula <- ~ DV200.value + Edad+ DFS.TIME
varPart <- fitExtractVarPartModel( as.data.frame(counts_data_filtered_cpm_phenotype), formula, sample_info_extracted)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )
pdf("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Variance_Partition/plotVarViolin_DV200_Edad_DFSTIME_DFSTIME.pdf")
plotVarPart( vp , main = "Variance Partition Analysis")
dev.off()
sample_info_extracted$Diagnostico
################################################################################

#######################  VARIANCE PARTITION ANALYSIS  ##########################
sample_info_2 <- sample_info_2[which(sample_info_2$`AC code` %in% colnames(counts_data)),]
sample_info$ICC <- sample_info_2$ICC

formula <- ~ purity + stromal + immune + (1|Diagnostico)
varPart <- fitExtractVarPartModel( counts_data, formula, sample_info)
vp <- sortCols( varPart )
colnames(vp) <- c("Purity", "Immune", "Stromal", "Diagnosis", "Residuals")
plotPercentBars( vp[1:10,] )
data <- stack(vp)

data
ggplot(data, aes(x =ind, y =values, fill =ind)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), ,legend.position = "none", axis.text.x = element_text(angle = 30, hjust=1),axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Variance Partition Basurto") +  
  scale_fill_manual(values = viridis(6)) + 
  ylab("Variance explained (%)") + 
  xlab("Variables")
ggsave("Results/Variance_Partition/VariancePartition_Basurto_Diagnosis_Purity_Immune_Stroma.pdf", heigh= 5, width=6)
################################################################################


################################################################################
#| Correlation Analysis
################################################################################

#| PTEN and Edad
pdf("Results/Correlation_plots/Correlation_PTEN_Edad.pdf")
ggplot(sample_info_extracted, aes(x= PTEN_Exp_Zscore, y= Edad)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Age vs PTEN Expression (r = 0.33)") + xlab("PTEN Expression") + ylab("Age")
dev.off()
cor(sample_info_extracted$PTEN_Exp, sample_info_extracted$Edad, method = 'pearson')

#| PTEN and Edad Classification PCa
pdf("Results/Correlation_plots/Correlation_PTEN_Edad_Diagnostico.pdf")
ggplot(sample_info_extracted, aes(x= PTEN_Exp_Zscore, y= Edad, color = Diagnostico)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"), legend.key.size = unit(1, 'cm')) + 
  ggtitle("Age vs PTEN Expression") + ylab("Age")+ xlab("PTEN Expression") + scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + labs(color ="Diagnosis")
dev.off()
cor(sample_info_extracted$PTEN_Exp[sample_info_extracted$Diagnostico == "PCa"], sample_info_extracted$Edad[sample_info_extracted$Diagnostico == "PCa"], method = 'pearson')  #| 0.02709694
cor(sample_info_extracted$PTEN_Exp[sample_info_extracted$Diagnostico == "BPH"], sample_info_extracted$Edad[sample_info_extracted$Diagnostico == "BPH"], method = 'pearson')  #| 0.262104

#| PTEN and DV200
pdf("Results/Correlation_plots/Correlation_PTEN_DV200.pdf")
ggplot(sample_info_extracted, aes(x= PTEN_Exp_Zscore, y= DV200.value)) + geom_point(size=2.5, color="midnightblue") + 
  geom_smooth(method=lm, color="green") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("DV200 vs PTEN Expression") + ylab("DV200") + xlab("PTEN Expression")
dev.off()

#| Edad and DFS.TIME
pdf("Results/Correlation_plots/Correlation_Edad_DFS.TIME.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$DFS.TIME !=0 & !is.na(sample_info_extracted$DFS.TIME)),], aes(x= DFS.TIME, y= Edad, color = Diagnostico)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Age vs DFS.TIME") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + ylab("Age")
dev.off()


#| Edad and BMI
pdf("Results/Correlation_plots/Correlation_Edad_BMI.pdf")
ggplot(sample_info_extracted, aes(x= BMI, y= Edad, color = Diagnostico)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Age vs BMI") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + ylab("Age") +labs(color ="Diagnosis")
dev.off()
cor(sample_info_extracted$BMI[sample_info_extracted$Diagnostico == "PCa"], sample_info_extracted$Edad[sample_info_extracted$Diagnostico == "PCa"], method = 'pearson')  #| 0.09568011
cor(sample_info_extracted$BMI[sample_info_extracted$Diagnostico == "BPH"], sample_info_extracted$Edad[sample_info_extracted$Diagnostico == "BPH"], method = 'pearson')  #| -0.07629329

#| PTEN and BMI
pdf("Results/Correlation_Plots/Correlation_PTEN_BMI.pdf")
ggplot(sample_info_extracted, aes(x= PTEN_Exp, y= BMI, color = Diagnostico)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs BMI") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + xlab("PTEN Expression") + labs(color ="Diagnosis")
dev.off()
cor(sample_info_extracted$PTEN_Exp[sample_info_extracted$Diagnostico == "PCa"], sample_info_extracted$BMI[sample_info_extracted$Diagnostico == "PCa"], method = 'pearson')  #| -0.065
cor(sample_info_extracted$PTEN_Exp[sample_info_extracted$Diagnostico == "BPH"], sample_info_extracted$BMI[sample_info_extracted$Diagnostico == "BPH"], method = 'pearson')  #| 0.22

#| Converting H-scores to numeric
sample_info_extracted$H_score <- as.numeric(sample_info_extracted$H_score)
sample_info_extracted$DFS.TIME <- as.numeric(sample_info_extracted$DFS.TIME)

#| Edad and H_score
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_Edad_H_score.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= Edad)) + geom_point(size=2.5) + 
  geom_smooth(method=lm, color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  xlab("H-score")+
  ggtitle("Age vs H_score PTEN Assessment") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + ylab("Age") + labs(color ="DFS STATUS")
dev.off()

#| PTEN and H_score 
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= PTEN_Exp_log2)) + geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  xlab("H-score") +
  annotate(geom="text", x=150, y=10.2, label=paste("r = ",format(round(cor(sample_info_extracted$PTEN_Exp_Zscore[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))]),2), nsmall =2),sep=""), color="black", family ="serif", size =6) 
dev.off()

#| PTEN and H_score wihout zeros
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_No_Zeros.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score) & sample_info_extracted$H_score != 0),], aes(x= H_score, y= PTEN_Exp_log2)) + geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=150, y=10.2, label=paste("r = ", format(round(cor(sample_info_extracted$PTEN_Exp_Zscore[which(!is.na(sample_info_extracted$H_score)& sample_info_extracted$H_score != 0)], sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score)& sample_info_extracted$H_score != 0)], method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6) +
  labs(color ="DFS STATUS")
dev.off()


#| DFS == 1
sample <- sample_info_extracted[which(sample_info_extracted$DFS.STATUS == 1 & !is.na(sample_info_extracted$DFS.STATUS) & !is.na(sample_info_extracted$H_score) ),]
data <- counts_data_filtered_cpm_phenotype_normalized[,sample$`AC basurto`]

#| PTEN and H_score DFS == 1
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_Recurrence.pdf")
ggplot(sample, aes(x= H_score, y= PTEN_Exp_log2)) + geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment \n Only Recurrence") +  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=150, y=10.2, label=paste("r = ", format(round(cor(sample_info_extracted$PTEN_Exp_Zscore[which(!is.na(sample_info_extracted$H_score)& sample_info_extracted$H_score != 0)], sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score)& sample_info_extracted$H_score != 0)], method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6) +
  labs(color ="DFS STATUS")

dev.off()

#| Edad and H_score
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_Edad_H_score_Recurrence.pdf")
ggplot(sample, aes(x= H_score, y= Edad)) + geom_point(size=2.5) + 
  geom_smooth(method=lm, color ="midnightblue") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Age vs H_score PTEN\n Recurrence") + ylab("Age") + xlab("H-score")+
  annotate(geom="text", x=150, y=50, label=paste("r = ", format(round(cor(sample$Edad,sample$H_score, method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6) 
dev.off()


#| DFS.TIME and H_score all zero filtered
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_DFS_TIME_H_score_Recurrence.pdf")
ggplot(sample, aes(x= DFS.TIME, y= H_score)) + geom_point(size=2.5) + 
  geom_smooth(method=lm) + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("DFS TIME vs H_score PTEN Assessment \n Recurrence") +  
  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  xlab("Time of recurrence") + labs(color ="DFS STATUS") +
  ylab("H-score")+
  annotate(geom="text", x=60, y=170, label=paste("r = ", format(round(cor(sample$Edad,sample$H_score, method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6) 
dev.off()
################################################################################


################################################################################
#| BoxPlot
################################################################################

#| PTEN and Diagnostico
pdf("Results/Box_plots/BoxPlot_PTEN_PCa_vs_BPH.pdf")
ggplot(sample_info, aes(x=Diagnostico, y=PTEN_Exp_Zscore, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  ggtitle("PTEN Expression vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + xlab("Diagnosis") + labs(fill ="Diagnosis") + ylab("PTEN Expression")
dev.off()

#| Edad and Diagnostico
pdf("Results/Box_plots/BoxPlot_Edad_PCa_vs_BPH.pdf")
ggplot(sample_info_extracted, aes(x=Diagnostico, y=Edad, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  ggtitle("Age vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + labs(fill ="Diagnosis") + xlab("Diagnosis") + ylab("Age")
dev.off()

#| DV200 and Diagnostico
pdf("Results/Box_plots/BoxPlot_DV200_PCa_vs_BPH.pdf")
ggplot(sample_info_extracted, aes(x=Diagnostico, y=DV200.value, fill = Diagnostico)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  ggtitle("DV200 vs PCa and BPH")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) + labs(fill ="Diagnosis") + xlab("Diagnosis") +
  ylab("DV200")
dev.off()

#| H-score and Gleason score
sample_info_extracted$Gleason_score_pieza<- as.character(sample_info_extracted$Gleason_score_pieza)

my_comparisons <- list( c("6", "7"), c("7", "8"), c("8", "9"), c("6", "9") )
pdf("Results/Box_plots/PTEN_Assessment_Boxplot_Gleason_score_H_score")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$Gleason_score_pieza)&!is.na(sample_info_extracted$H_score) & !is.na(sample_info_extracted$DFS.STATUS)),], aes(x= Gleason_score_pieza, y= H_score, fill =Gleason_score_pieza)) + 
  geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Gleason score vs H-score PTEN Assessment") +  
  scale_fill_manual(values = viridis(4)) + 
  xlab("Gleason score") + ylab("H-score")+
  labs(fill ="Gleason score") + stat_compare_means(comparisons = my_comparisons, family = "Times New Roman")
dev.off()
################################################################################

