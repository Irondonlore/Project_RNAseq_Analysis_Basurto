################################################################################
########## ESTIMATING TUMOR PURITY, STROMAL AND IMMUNE SCORES ##################
################################################################################

#| ESTIMATE provides researchers with scores for tumor purity, the level of stromal 
#| cells present, and the infiltration level of immune cells in tumor tissues based 
#| on expression data. Normal cells in tumor tissue not only influence the tumor signal
#| in molecular studies but also play an important role in cancer biology. The estimate 
#| package predicts the presence of stromal and immune cells in tumor tissue using gene 
#| expression data
################################################################################

################################ LIBRARIES #####################################
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(corrr))


#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("W:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Set working directory
setwd("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/")
################################################################################


############################## READING DATA ####################################
counts_data <- read.table("Results/Counts/FullCounts_Basurto_Filtered_counts_phenotype.txt", sep ="\t")
sample_info <- read.table("Results/Sample_info_table/sample_info_extracted.txt", sep ="\t")
################################################################################



################################ ESTIMATE ######################################

#| Assining a new variable
counts_data_gene_name <- counts_data
counts_data_gene_name$GeneID <- rownames(counts_data_gene_name)

#| Merging the genome_GRCh39.94 with the counts data by GeneID
counts_data_gene_name <- merge(counts_data_gene_name, genome_GRCh39.94, by ="GeneID")

#| Ignoring duplicated gene_names
counts_data_gene_name <- counts_data_gene_name[-which(duplicated(counts_data_gene_name$gene_name)),]

#| Assigning gene_names to rownames of counts data
rownames(counts_data_gene_name) <- counts_data_gene_name$gene_name
counts_data_gene_name <- counts_data_gene_name[,-c(240:245)]
counts_data_gene_name <- counts_data_gene_name[,-c(1)]

#| Applying estrimate
est <- counts_data_gene_name |> 
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)

#| Adding Stromal, Immune, estimate and purity values to sample_info_extracted table
names(sample_info)
sample_info <- sample_info[,-c(23:26)]
sample_info<- cbind(sample_info, est[,-c(1)])

#| How is related the purity of the tumor with the PTEN expression at the transcriptional level?
names(sample_info)
#| 1) Correlation  plot
cor <- cor.test(sample_info$PTEN_Exp_log2, sample_info$purity, method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 19)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Tumor_purity_as_function_PTEN_normalized_log2_counts_and_Diagnosis_PCa_vs_BPH.pdf")
ggplot(sample_info, aes(x= PTEN_Exp_log2, y = purity, color =Diagnostico)) + geom_point(size =3) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Tumor purity as a function of PTEN transcriptional levels")+
  xlab("PTEN expression (Normalized Log2 counts)") + 
  ylab("Tumor purity") +geom_smooth(method=lm, color ="green") +
  annotate(geom="text", x=10.6, y=0.72, label=paste("r = ", r, ", p = ", p, sep =""), color="black", family ="serif", size =5) +
  scale_color_manual(values = c("goldenrod3", "midnightblue"))
dev.off()

#| 2) BoxPlot: Tumor purity vs Diagnosis
ggplot(sample_info, aes(x= Diagnostico, y = purity, fill =Diagnostico)) + geom_boxplot() +
  theme(text=element_text(size=16,  family="serif"),legend.position = "none", legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #ggtitle("Tumor purity and Diagnosis (BPH and PCa)")+
  xlab("Diagnosis") + 
  ylab("Tumor purity")  +
  scale_fill_manual(values = c("lightseagreen", "midnightblue")) +
  labs(fill='Diagnosis') +
  stat_compare_means(label ="p.format", family="serif", size =5)
ggsave("Results/Box_plots/Tumor_purity_and_Diagnosis.pdf", height = 5, width = 5.3)

#| 3) BoxPlot: Stromal Score vs Diagnosis
ggplot(sample_info, aes(x= Diagnostico, y = stromal, fill =Diagnostico)) + geom_boxplot() +
  theme(text=element_text(size=16,  family="serif"),legend.position = "none", legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #ggtitle("Tumor purity and Diagnosis (BPH and PCa)")+
  xlab("Diagnosis") + 
  ylab("Stromal score")  +
  scale_fill_manual(values = c("lightseagreen", "midnightblue")) +
  labs(fill='Diagnosis') +
  stat_compare_means(label ="p.format", family="serif", size =5)
ggsave("Results/Box_plots/StromalScore_and_Diagnosis.pdf", height = 5, width = 5.3)

#| 4) BoxPlot: Immune Score vs Diagnosis
ggplot(sample_info, aes(x= Diagnostico, y = immune, fill =Diagnostico)) + geom_boxplot() +
  theme(text=element_text(size=16,  family="serif"),legend.position = "none", legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #ggtitle("Tumor purity and Diagnosis (BPH and PCa)")+
  xlab("Diagnosis") + 
  ylab("Immune score")  +
  scale_fill_manual(values = c("lightseagreen", "midnightblue")) +
  labs(fill='Diagnosis') +
  stat_compare_means(label ="p.format", family="serif", size =5)
ggsave("Results/Box_plots/ImmuneScore_and_Diagnosis.pdf", height = 5, width = 5.3)
################################################################################
