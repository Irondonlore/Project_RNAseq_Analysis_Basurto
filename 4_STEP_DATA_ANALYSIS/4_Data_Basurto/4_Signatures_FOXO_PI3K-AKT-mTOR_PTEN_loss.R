#############  ANALYZING DIFFERENT SIGNATURES IN BASURTO's COHORT  #############

#| I have extracted and analyse different signatures that can be useful to explore
#| PTEN loss.

#| FOXO signature: https://www2.stat.duke.edu/~sayan/genesets/Jan2006/cards/C3/TTGTTT_V_FOXO4_01.html

#| PI3K-AKT-mTOR signature: https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html

#| FOXO pathway signature: https://www.gsea-msigdb.org/gsea/msigdb/cards/PID_FOXO_PATHWAY

#| PTEN loss signature: https://bmccancer.biomedcentral.com/articles/10.1186/s12885-021-08593-y

################################################################################


###############################  LIBRARIES  ####################################
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(variancePartition))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(corrr))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
################################################################################


############################## DATA DIRECTORIES ################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted_processing_and_filtering.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_CPM_phenotype_Normalized_Log2.txt"
foxo.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TTGTTT_V_FOXO4_01.txt"
foxo_pathway.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/FOXO_Pathway.txt"
pi3k_akt_mtor.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt"
pten_loss.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Signature_PTEN_loss_Paper_Imada_et_al_BMC_cancer.xlsx"
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################


############################### READING DATA ###################################
#| Counts data (normalized)
counts_data <- read.table(counts.file, sep ="\t", header=T)

#| Sample information
sample_info <- read.table(info.file, sep ="\t", header=T)

#| FOXO signature
foxo_signature <- read.table(foxo.file, sep ="\t", quote="",header=T)  #dim= 1262 4

#| PI3K/AKT/MTOR signature
pi3k_akt_mtor_signature <- read.table(pi3k_akt_mtor.file, sep =",", quote="",header=T)

#| PTEN loss
signature_PTEN_Loss <-  read_xlsx(pten_loss.file, sheet= "S1_BHM_signature")

#| FOXO pathway GSEA
foxo_pathway_signature <- read.table(foxo_pathway.file, sep =",", quote="",header=T)
################################################################################


####################### COMPUTING Z-SCORES OF NORMALIZED DATA ################## 
counts_zscore <- scale(t(as.data.frame(counts_data)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)

#| Boxplot zscore counts
boxplot(counts_zscore, col = c("orange"), font.family ="serif")
################################################################################


############################  FOXO SIGNATURE  ##################################
#| Changing name columns
colnames(foxo_signature) <- c("gene_name", "STATUS", "LOOKED_UP_SYMBOL", "TITLE")

#| Merging gene names to the foxo signature
foxo_signature <- merge(foxo_signature,genome_GRCh39.94, by ="gene_name" )  #dim = 842 10

#| Finding gene of the FOXO signature inside the gene expression data
counts_foxo <- counts_zscore[which(rownames(counts_zscore) %in% foxo_signature$GeneID),]

#| Computing the mean over all the columns (all samples)
samples_media_expression <- colMeans(counts_foxo, na.rm=TRUE)

#| Saving to the sample_info table
sample_info$`Mean expression FOXO` <- samples_media_expression
  
#| Mean expression FOXO vs Diagnosis
pdf("Results/Box_plots/Mean_expression_FOXO_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info, aes(x =Diagnostico, y =`Mean expression FOXO`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Mean expression FOXO signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("Mean expression FOXO signature (Z score of normalized counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression FOXO vs PTEN_expression
cor <- cor.test(sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
p <- round(cor$p.value,16)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_Expression_Mean_FOXO_Signature_198_PCa_samples_Mean_expression_FOXO_signature.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= `Mean expression FOXO`)) + 
  geom_point(size=2.5, color ="darkgoldenrod3") + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN expression vs FOXO signature \n 198 PCa samples") +  
  ylab("Mean expression FOXO (Z score of Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=10.5, y=0.5, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| Mean expression FOXO and purity
cor_PCa <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
cor_BPH <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "BPH")], sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "BPH")], method ="pearson", exact = FALSE)

p_PCa <- round(cor_PCa$p.value,8)
r_PCa <- round(as.numeric(cor_PCa$estimate[[1]]),2)

p_BPH <- round(cor_BPH$p.value,1)
r_BPH <- round(as.numeric(cor_BPH$estimate[[1]]),2)

pdf("Results/Correlation_plots/Purity_Mean_FOXO_Signature_198_PCa_samples.pdf")
ggplot(sample_info_extracted, aes(x= purity, y= `Mean expression FOXO`, color = Diagnostico, fill =Diagnostico)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("") +  
  ylab("Mean expression FOXO (Z scores Normalized log2 counts)")+ 
  xlab("Purity") +
  scale_color_manual(values=c("midnightblue", "green")) +
  scale_fill_manual(values=c("midnightblue", "green"))+
  annotate(geom="text", x=0.75, y=0.45, label=paste("r = ",r_BPH,", p < ", p_BPH, sep=""), color="black", family ="serif", size =4.5)+
  annotate(geom="text", x=0.75, y=-0.09, label=paste("r = ",r_PCa,", p < ", p_PCa, sep=""), color="black", family ="serif", size =4.5)
dev.off()
################################################################################



#########################  PI3K/AKT/MTOR SIGNATURE  ############################

#| Converting it to data.frame
pi3k_akt_mtor_signature <- data.frame(gene_name = colnames(pi3k_akt_mtor_signature))

#| Merging gene names to the foxo signature
pi3k_akt_mtor_signature <- merge(pi3k_akt_mtor_signature,genome_GRCh39.94, by ="gene_name" )  

#| Finding gene of the pi3k_akt_mtor signature inside the gene expression data
counts_pi3k_akt_mtor <- counts_zscore[which(rownames(counts_zscore) %in% pi3k_akt_mtor_signature$GeneID),]

#| Computing the mean over all the columns (all samples)
samples_media_expression <- colMeans(counts_pi3k_akt_mtor, na.rm=TRUE)

#| Saving to the sample_info table
sample_info_extracted$`Mean expression PI3K-AKT-mTOR` <- samples_media_expression


#| Mean expression PI3K-AKT-mTOR vs Diagnosis
pdf("Results/Box_plots/Mean_expression_PI3K-AKT-mTOR_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info_extracted, aes(x =Diagnostico, y =`Mean expression PI3K-AKT-mTOR`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Mean expression PI3K-AKT-mTOR signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("Mean expression PI3K-AKT-mTOR signature (Z scores Normalized log2 counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression PI3K-AKT-mTOR vs Diagnosis
pdf("Results/Box_plots/PC1_PI3K-AKT-mTOR_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info_extracted, aes(x =Diagnostico, y =`PC1 PI3K-AKT-mTOR`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PC1 PI3K-AKT-mTOR signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("PC1 PI3K-AKT-mTOR signature (Z scores Normalized log2 counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression PI3K-AKT-mTOR vs PTEN_expression
cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,7)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_Expression_Mean_PI3K-AKT-mTOR_Signature_198_PCa_samples_Mean_expression_FOXO_signature.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= `Mean expression PI3K-AKT-mTOR`)) + 
  geom_point(size=2.5, color ="darkgoldenrod3") + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN expression vs PI3K-AKT-mTOR signature \n 198 PCa samples") +  
  ylab("Mean expression PI3K-AKT-mTOR (Z score Normalized log2 counts)")+ 
  xlab("PTEN Expression (Z score Normalized log2 counts)")+
  annotate(geom="text", x=10.5, y=-0.25, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()


#| Mean expression PI3K-AKT-mTOR vs Mean expression FOXO 
cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,3)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_PI3K_mTOR_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PI3K-AKT-mTOR`, y =`Mean expression FOXO`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=0.5, y=0.46, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression FOXO (Z score Normalized log2 counts)")+ 
  xlab("Mean expression PI3K-AKT-mTOR (Z score Normalized log2 counts)")+
  ggtitle("Correlation between FOXO and PI3K-AKT-mTOR signature\n 198 PCa samples")
dev.off()


#| Mean expression PI3K-AKT-mTOR vs Mean expression FOXO pathway
cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,9)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_pathway_PI3K_mTOR_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PI3K-AKT-mTOR`, y =`Mean expression FOXO pathway`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.2, y=0.9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression FOXO pathway (Z scores Normalized log2 counts)")+ 
  xlab("Mean expression PI3K-AKT-mTOR (Z scores Normalized log2 counts)")+
  ggtitle("Correlation between FOXO pathway and PI3K-AKT-mTOR signature\n 198 PCa samples")
dev.off()

#| Mean expression Pi3k-AKT-mTOR
cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,12)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_PTEN_Loss_PI3K_mTOR_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PI3K-AKT-mTOR`, y =`Mean expression PTEN_loss`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.2, y=0.7, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression PTEN loss (Z scores Normalized log2 counts)")+ 
  xlab("Mean expression PI3K-AKT-mTOR (Z scores Normalized log2 counts)")+
  ggtitle("Correlation between FOXO pathway and PI3K-AKT-mTOR signature\n 198 PCa samples")
dev.off()

#| Mean expression PI3K-AKT-mTOR and purity
cor_PCa <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
cor_BPH <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "BPH")], sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "BPH")], method ="pearson", exact = FALSE)

p_PCa <- round(cor_PCa$p.value,3)
r_PCa <- round(as.numeric(cor_PCa$estimate[[1]]),2)

p_BPH <- round(cor_BPH$p.value,3)
r_BPH <- round(as.numeric(cor_BPH$estimate[[1]]),2)

pdf("Results/Correlation_plots/Purity_Mean_PI3K-AKT-mTOR_Signature_198_PCa_samples.pdf")
ggplot(sample_info_extracted, aes(x= purity, y= `Mean expression PI3K-AKT-mTOR`, color = Diagnostico, fill =Diagnostico)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("") +  
  ylab("Mean expression PI3K-AKT-mTOR pathway (Z scores Normalized log2 counts)")+ 
  xlab("Purity") +
  scale_color_manual(values=c("midnightblue", "green")) +
  scale_fill_manual(values=c("midnightblue", "green"))+
  annotate(geom="text", x=0.95, y=-0.5, label=paste("r = ",r_BPH,", p < ", p_BPH, sep=""), color="black", family ="serif", size =4.5)+
  annotate(geom="text", x=0.75, y=0.35, label=paste("r = ",r_PCa,", p < ", p_PCa, sep=""), color="black", family ="serif", size =4.5)
dev.off()
################################################################################



########################### SIGNATURE PTEN LOSS ################################
#| Data
signature_PTEN_Loss <- signature_PTEN_Loss$X1[which(signature_PTEN_Loss$concordant >= 0.95)]
signature_PTEN_Loss <- data.frame(gene_name =signature_PTEN_Loss)

#| Merging gene names to the foxo signature
signature_PTEN_Loss <- merge(signature_PTEN_Loss,genome_GRCh39.94, by ="gene_name" )  

#| Finding gene of the FOXO signature inside the gene expression data
counts_signatue_PTEN_loss <- counts_zscore[which(rownames(counts_zscore) %in% signature_PTEN_Loss$GeneID),]

#| Computing the mean over all the columns (all samples)
samples_media_expression <- colMeans(counts_signatue_PTEN_loss, na.rm=TRUE)

#| Saving to the sample_info table
sample_info_extracted$`Mean expression PTEN_loss` <- samples_media_expression

#| Mean expression PTEN_loss vs Diagnosis
pdf("Results/Box_plots/Mean_expression_PTEN_loss_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info_extracted, aes(x =Diagnostico, y =`Mean expression PTEN_loss`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Mean expression PTEN_loss signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("Mean expression PTEN_loss signature (Z score of normalized counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression PTEN_loss vs PTEN_expression
cor <- cor.test(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
p <- round(cor$p.value,3)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_Expression_Mean_PTEN_loss_Signature_198_PCa_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= `Mean expression PTEN_loss`)) + 
  geom_point(size=2.5, color ="darkgoldenrod3") + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN expression vs PTEN loss signature \n 198 PCa samples") +  
  ylab("Mean expression PTEN_loss (Z scores Normalized log2 counts)")+ 
  xlab("PTEN Expression (Normalized log2 counts)")+
  annotate(geom="text", x=12, y=0.9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| Mean expression PTEN_loss and purity
cor_PCa <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
cor_BPH <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "BPH")], sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "BPH")], method ="pearson", exact = FALSE)

p_PCa <- round(cor_PCa$p.value,8)
r_PCa <- round(as.numeric(cor_PCa$estimate[[1]]),2)

p_BPH <- round(cor_BPH$p.value,1)
r_BPH <- round(as.numeric(cor_BPH$estimate[[1]]),2)

pdf("Results/Correlation_plots/Purity_Mean_PTEN_loss_Signature_198_PCa_samples.pdf")
ggplot(sample_info_extracted, aes(x= purity, y= `Mean expression PTEN_loss`, color = Diagnostico, fill =Diagnostico)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("") +  
  ylab("Mean expression PTEN loss (Z scores Normalized log2 counts)")+ 
  xlab("Purity") +
  scale_color_manual(values=c("midnightblue", "green")) +
  scale_fill_manual(values=c("midnightblue", "green"))+
  annotate(geom="text", x=0.95, y=-0.5, label=paste("r = ",r_BPH,", p < ", p_BPH, sep=""), color="black", family ="serif", size =4.5)+
  annotate(geom="text", x=0.75, y=0.1, label=paste("r = ",r_PCa,", p < ", p_PCa, sep=""), color="black", family ="serif", size =4.5)
dev.off()


#| Mean expression PTEN_loss vs Mean expression FOXO
cor <- cor.test(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,1)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_PTEN_loss_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PTEN_loss`, y =`Mean expression FOXO`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.2, y=0.5, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression FOXO (Z scores Normalized log2 counts)")+ 
  xlab("Mean expression PTEN loss (Z scores Normalized log2 counts)")+
  ggtitle("Correlation between FOXO and PTEN loss signature\n 198 PCa samples")
dev.off()

#| Mean expression PTEN_loss vs Mean expression FOXO pathway
cor <- cor.test(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,21)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_pathway_PTEN_loss_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PTEN_loss`, y =`Mean expression FOXO pathway`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.2, y=0.5, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression FOXO pathway (Z scores Normalized log2 counts)")+ 
  xlab("Mean expression PTEN loss (Z scores Normalized log2 counts)")+
  ggtitle("Correlation between FOXO pathway and PTEN loss signature\n 198 PCa samples")
dev.off()

#| Mean expression PTEN_loss vs Mean expression PI3K-AKT-mTOR
cor <- cor.test(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,10)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_PI3K-AKT-mTOR_PTEN_loss_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression PTEN_loss`, y =`Mean expression PI3K-AKT-mTOR`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.2, y=0.5, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression PI3k-AKT-mTOR (Z scores Normalized log2 counts)")+ 
  xlab("Mean expression PTEN loss (Z scores Normalized log2 counts)")+
  ggtitle("Correlation between PI3k-AKT-mTOR and PTEN loss signature\n 198 PCa samples")
dev.off()
################################################################################



############################  FOXO PATHWAY  ####################################
#| Converting it to data.frame
foxo_pathway_signature <- data.frame(gene_name = colnames(foxo_pathway_signature))

#| Merging gene names to the foxo signature
foxo_pathway_signature <- merge(foxo_pathway_signature,genome_GRCh39.94, by ="gene_name" )  #dim = 842 10

#| Finding gene of the pi3k_akt_mtor signature inside the gene expression data
counts_foxo_pathway <- counts_zscore[which(rownames(counts_zscore) %in% foxo_pathway_signature$GeneID),]

#| Computing the mean over all the columns (all samples)
samples_media_expression <- colMeans(counts_foxo_pathway, na.rm=TRUE)

#| Saving to the sample_info table
sample_info_extracted$`Mean expression FOXO pathway` <- samples_media_expression

#| Mean expression FOXO vs Diagnosis
pdf("Results/Box_plots/Mean_expression_FOXO_pathway_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info_extracted, aes(x =Diagnostico, y =`Mean expression FOXO pathway`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Mean expression FOXO pathway signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("Mean expression FOXO pathway signature (Z scores Normalized log2 counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression FOXO pathway vs Diagnosis
pdf("Results/Box_plots/PC1_FOXO_pathway_Diagnosis_Basurto_Cohort_238_samples.pdf")
ggplot(sample_info_extracted, aes(x =Diagnostico, y =`PC1 FOXO pathway`, fill =Diagnostico)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PC1 FOXO pathway signature vs Diagnosis") +  
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue")) + 
  ylab("PC1 FOXO pathway signature (Z scores Normalized log2 counts)") + 
  labs(fill ="Diagnosis") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

#| Mean expression FOXO vs PTEN_expression
cor <- cor.test(sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/PTEN_Expression_Mean_FOXO_pathway_Signature_198_PCa_samples_Mean_expression_FOXO_signature.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x= PTEN_Exp_log2, y= `Mean expression FOXO pathway`)) + 
  geom_point(size=2.5, color ="darkgoldenrod3") + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN expression vs FOXO pathway signature \n 198 PCa samples") +  
  ylab("Mean expression FOXO pathway (Z scores Normalized log2 counts)")+ 
  xlab("PTEN Expression (Z scores Normalized log2 counts)")+
  annotate(geom="text", x=12, y=0.9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6)
dev.off()

#| Mean expression FOXO vs Mean expression FOXO pathway
cor <- cor.test(sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson")
p <- round(cor$p.value,5)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_pathway_and_FOXO_signatures_PCa_samples_Basurto.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x =`Mean expression FOXO`, y =`Mean expression FOXO pathway`)) +
  geom_point(size =2.3, color = "deeppink4")+ 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))  +
  geom_smooth(method = "lm", color ="lightblue4") +
  annotate(geom="text", x=-0.3, y=0.9, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =6) +
  ylab("Mean expression FOXO pathway (Z scores normalized log2 counts)")+ 
  xlab("Mean expression FOXO (Z scores normalized log2 counts)")+
  ggtitle("Correlation between FOXO pathway and PI3K-AKT-mTOR signature\n 198 PCa samples")
dev.off()

#| Mean expression FOXO pathway and purity
cor_PCa <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact = FALSE)
cor_BPH <- cor.test(sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "BPH")], sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$Diagnostico == "BPH")], method ="pearson", exact = FALSE)

p_PCa <- round(cor_PCa$p.value,1)
r_PCa <- round(as.numeric(cor_PCa$estimate[[1]]),2)

p_BPH <- round(cor_BPH$p.value,1)
r_BPH <- round(as.numeric(cor_BPH$estimate[[1]]),2)

pdf("Results/Correlation_plots/Purity_Mean_FOXO_pathway_Signature_198_PCa_samples.pdf")
ggplot(sample_info_extracted, aes(x= purity, y= `Mean expression FOXO pathway`, color = Diagnostico, fill =Diagnostico)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("") +  
  ylab("Mean expression FOXO pathway (Z scores Normalized log2 counts)")+ 
  xlab("Purity") +
  scale_color_manual(values=c("midnightblue", "green")) +
  scale_fill_manual(values=c("midnightblue", "green"))+
  annotate(geom="text", x=0.95, y=-0.5, label=paste("r = ",r_BPH,", p < ", p_BPH, sep=""), color="black", family ="serif", size =4.5)+
  annotate(geom="text", x=0.75, y=0.2, label=paste("r = ",r_PCa,", p < ", p_PCa, sep=""), color="black", family ="serif", size =4.5)
dev.off()
################################################################################
