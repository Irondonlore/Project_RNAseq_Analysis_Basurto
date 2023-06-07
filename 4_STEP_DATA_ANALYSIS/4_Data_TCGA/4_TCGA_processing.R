#########################  TCGA PROCESSING DATA  ###############################

#| Here, I have processed the TCGA data, following the info from this paper: 
#| https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10050165/ (a pan-cancer genomic 
#| analysis demostratic that PTEN loss tumors are associated with immune evasion 
#| and poor ourcome from TCGA. 

#| Raw read counts for TCGA-tumors were downloaded using FANTOM-CAT/recount2 
#| (https://jhubiostatistics.shinyapps.io/recount/) in the "Download Data with R".
#| Particularly, this webside have processed the RNAseq from TCGA for obtaining
#| matrix of row counts.

################################################################################


###############################  LIBRARIES  ####################################
suppressMessages(library(DESeq2))
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
suppressMessages(library(recount))
suppressMessages(library(schoolmath))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(gprofiler2))
suppressMessages(library(xCell))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
################################################################################


############################## DATA DIRECTORIES ################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/data_clinical_patient.txt"
pten_cna.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/TCGA PTEN info/cna.txt"
pten_protein.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/TCGA PTEN info/Protein expression (RPPA).txt"
pten_mrna.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/TCGA PTEN info/mRNA expression (RNA Seq V2 RSEM).txt"
pi3k_akt_mtor.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt"
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################


################################################################################
#| Load the data
load(file.path('TCGA', 'rse_gene.Rdata'))

#| Scale counts. rse is a RangedSummarizedExperiment object. Here it is contained the information for many cancers
rse <- scale_counts(rse_gene)

#| Filtering Prostate cancer patients (Classified by Prostate Adenocarcinoma)
prostate_rse <- rse[,which(rse$gdc_cases.project.name == "Prostate Adenocarcinoma")]

#| Filtering primary tumors
prostate_rse <- prostate_rse[,which(prostate_rse$gdc_cases.samples.sample_type == "Primary Tumor")]

#| Ignoring duplicated samples
prostate_rse <- prostate_rse[,which(!duplicated(prostate_rse$gdc_cases.submitter_id))]

#| Extracting sample ID
sample_id <- prostate_rse$gdc_cases.submitter_id

#| Counts data
counts_data <- assay(prostate_rse)
counts_data <- data.frame(counts_data)

#| Changing Colnames from the counts data
colnames(counts_data) <- sample_id

#| Assigning correct ensembl ID
rownames(counts_data) <- genes$group_name

#| Verifying that are counts and not RSEM
decimals <- c()
for (i in 1:dim(counts_data)[2]){
  decimals <- c(decimals, sum(floor(as.numeric(counts_data[,c(i)])) - as.numeric(counts_data[,c(i)])) )
}
decimals  #| Vector with 0s

#| Finding low counts
low_counts <- rownames(counts_data)[which(rowSums(counts_data > 5) > 0.5*dim(counts_data)[2] )]

#| Filtering low counts counts
counts_data <- counts_data[ which(rownames(counts_data) %in% low_counts),]

#| Correct ensembl ID
genes <- rowData(prostate_rse)
genes <- data.frame(genes)
genes <- genes[which(genes$gene_id %in% rownames(counts_data)),]
symbol <- as.character(genes$symbol)

#| Assigning a column with gene symbols
counts_data$gene_name <- symbol
counts_data <- counts_data[which(!is.na(counts_data$gene_name)),]

#| Merging to obtain symbol names with the genome file
counts_data <- merge(counts_data, genome_GRCh39.94, by ="gene_name")

#| Aggregating dataframe
counts_data <- aggregate(counts_data, by=list(counts_data$GeneID), mean)

#| Naming rownames as ENSEMBL ID
rownames(counts_data) <- counts_data$Group.1

#| Deleting the columns created
counts_data <- counts_data[,-c(1,2,500, 501, 502, 503, 504, 505)]

#| Round values
counts_data <- round(counts_data)
################################################################################


################################################################################
#| PTEN CNA 
pten_cna <- read.table(pten_cna.file, sep = "\t", header =T)
colnames(pten_cna) <- c("STUDY_ID","PATIENT_ID", "PTEN_cna")
pten_cna$PATIENT_ID <- sub("-01","",pten_cna$PATIENT_ID )
pten_cna <- pten_cna[which(pten_cna$PATIENT_ID %in% colnames(counts_data)),]

#| PTEN protein
pten_protein <- read.table(pten_protein.file, sep ="\t", header =T)
colnames(pten_protein) <- c("STUDY_ID","PATIENT_ID", "PTEN_protein")
pten_protein$PATIENT_ID <- sub("-01","", pten_protein$PATIENT_ID )
pten_protein <- pten_protein[which(pten_protein$PATIENT_ID %in% colnames(counts_data)),]

#| Read table Sample information
sample_info <- read.table(info.file, sep = "\t", header = T)
sample_info <- sample_info[which(sample_info$PATIENT_ID %in% colnames(counts_data)),]

#| Adding CNA, mRNA, and protein info to sample_info
df_list <- list(sample_info, pten_cna,pten_protein)      

#merge all data frames together
sample_info <- Reduce(function(x, y) merge(x, y, by ="PATIENT_ID"), df_list) 

#| Order columns and rows
sample_info <- sample_info[order(sample_info$PATIENT_ID),]
counts_data <- counts_data[,order(colnames(counts_data))]

#| VERIFYING: Are all the PATIENT_ID of sample_info equal to the columns of counts_data?
any(sample_info$PATIENT_ID == colnames(counts_data))
################################################################################


############################ DATA NORMALIZATION ################################
#| Normalize with DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ 1)

#| Estimating the size factors
dds <- estimateSizeFactors(dds)

#| Computing normalization
counts_data_normalized <- counts(dds, normalized=TRUE)

#| log2 for re-escaling factors
counts_data_normalized <- log(counts_data_normalized + 1, base =2)
counts_data_normalized <- data.frame(counts_data_normalized)

#| Box plot
#boxplot(counts_data_normalized)
#boxplot(counts_data)

#| Saving row counts and normalized
write.table(counts_data, "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/RawCounts/TCGA_raw_counts_filtered.txt", sep ="\t")
write.table(counts_data_normalized, "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/NormalizedCounts/TCGA_raw_counts_filtered_normalized_log2.txt", sep ="\t")

#| Adding PTEN mRNA info from normalized data
sample_info$PTEN_mrna <- as.numeric(counts_data_normalized[which(rownames(counts_data_normalized) == "ENSG00000171862"),])


#| PI3K/AKT/MTOR signature
pi3k_akt_mtor_signature <- read.table(pi3k_akt_mtor.file, sep =",", quote="",header=T)
pi3k_akt_mtor_signature <- data.frame(gene_name = colnames(pi3k_akt_mtor_signature))

#| Computing z-scores
counts_zscore <- scale(t(as.data.frame(counts_data_normalized)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)

counts_pi3k_akt_mtor <- counts_zscore[which(rownames(counts_zscore) %in% genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% pi3k_akt_mtor_signature$gene_name)]),]
samples_media_expression <- colMeans(counts_pi3k_akt_mtor, na.rm=TRUE)
samples_media_expression <- as.data.frame(samples_media_expression)
sample_info$`PI3K-AKT-mTOR` <- samples_media_expression$samples_media_expression
################################################################################


#######################  VARIANCE PARTITION ANALYSIS  ##########################
formula <- ~  (1|PTEN_cna) + `PI3K-AKT-mTOR` + PTEN_mrna + PTEN_protein + Purity + AGE

varPart <- fitExtractVarPartModel( counts_data, formula, sample_info)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )

data <- stack(vp)

ggplot(data, aes(x =ind, y =values, fill =ind)) +
  geom_boxplot()+ theme(text=element_text(size=16,  family="serif"), ,legend.position = "none", axis.text.x = element_text(angle = 30, hjust=1),axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Variance Partition TCGA") +  
  scale_fill_manual(values = viridis(7)) + 
  ylab("Variance explained (%)") + 
  xlab("Variables")
ggsave("VariancePartition/VariancePartition_TCGA_PTEN_CNA_mRNA_Protein_PI3k-AKT-mTOR_AGE.pdf", heigh= 6, width=6)
################################################################################




################################## xCell #######################################

#| To apply xCell, we need to fix our counts_data. Rownames must be gene symbol
expression_data <- counts_data

#| Creating a new column with gene symbol IDs
expression_data$gene_name <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(counts_data))]

#| Given there are repetitive genes (with the same counts), we apply aggregate
expression_data <- aggregate(expression_data, by=list(expression_data$gene_name), mean)

#| Asigning new rownames
rownames(expression_data) <- expression_data$Group.1

#| Deleting undesired columns
expression_data <- expression_data[,-c(1,499)]

#| Applying xCell 
data <- xCellAnalysis(expression_data)
data <- t(data)
data <- data.frame(data)
data <- data[sample_info$PATIENT_ID[which(sample_info$PTEN_cna == "-2" | sample_info$PTEN_cna == "0")],]
data$PTEN_Status <- sample_info$PTEN_cna[which(sample_info$PTEN_cna == "-2" | sample_info$PTEN_cna == "0")]

data$PTEN_Status[which(data$PTEN_Status=="-2")] <- "loss"
data$PTEN_Status[which(data$PTEN_Status=="0")] <- "intact"

#| Saving results
write.table(data, "xCell/Table_results_TCGA.txt", sep ="\t")

#| Which one has the most significant difference?
p_value <-  c()
for (i in 1:(dim(data)[2] - 1) ){
  
  wilcox_test <- wilcox.test(data[which(data$PTEN_Status == "loss"), i], data[which(data$PTEN_Status == "intact"), i])
  p_value <- c(p_value, wilcox_test$p.value )
}

data_plot <- data.frame( datasets = colnames(data[,-c(68)]),
                         p_value =p_value )

data_plot <- data_plot[order(data_plot$p_value, decreasing = T),]
data_plot$log <- -log(data_plot$p_value,2)

ggplot(data_plot[48:67,], aes(x = log, y = reorder(datasets,log))) +
  geom_bar(stat = "identity", fill ="lightseagreen", alpha =0.7) + 
  scale_fill_manual(values = c("midnightblue", "darkgoldenrod3")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Datasets") +
  xlab("-log2(p-value)")
ggsave("xCell/xCell_Results_PTEN_genomic_loss_vs_intact_wilxon_test.pdf", heigh= 6, width =6)


ggplot(data, aes(x= PTEN_Status, y = StromaScore, fill =PTEN_Status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Stroma Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "PTEN status")
ggsave("xCell/BoxPlots/StromaScore_PTEN_genomic_loss_vs_intact.pdf", heigh= 5, width =6)

ggplot(data, aes(x= PTEN_Status, y = ImmuneScore, fill =PTEN_Status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Immune Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "PTEN status")+
  ylim(0,0.1)
ggsave("xCell/BoxPlots/ImmuneScore_PTEN_genomic_loss_vs_intact.pdf", heigh= 5, width =6)

ggplot(data, aes(x= PTEN_Status, y = Fibroblasts, fill =PTEN_Status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Fibroblasts") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "PTEN status")
ggsave("xCell/BoxPlots/Fibroblasts_PTEN_genomic_loss_vs_intact.pdf", heigh= 5, width =6)

ggplot(data, aes(x= PTEN_Status, y = Epithelial.cells, fill =PTEN_Status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Epithelial cells") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "PTEN status") +
  ylim(0,0.0025)
ggsave("xCell/BoxPlots/EpithelialCells_PTEN_genomic_loss_vs_intact.pdf", heigh= 5, width =6)

ggplot(data, aes(x= PTEN_Status, y = MicroenvironmentScore, fill =PTEN_Status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Microenvironment Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "PTEN status")+
  ylim(0,0.1)
ggsave("xCell/BoxPlots/MicroenvironmentScore_PTEN_genomic_loss_vs_intact.pdf", heigh= 5, width =6)
################################################################################


######################### ESTIMATING TUMOR PURITY #############################
est <- filter_common_genes(expression_data, id = "hgnc_symbol", tidy = FALSE, tell_missing = TRUE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)

#| Changing names
colnames(est) <- c("PATIENT_ID", "Stromal", "Immune", "Estimate", "Purity")

#| Merging dataframes
sample_info <- merge(sample_info, est, by ="PATIENT_ID")

#| Changin rownames of sample_info
rownames(sample_info) <- sample_info$PATIENT_ID
################################################################################


################################ PLOTS #########################################
#| PTEN protein distribution
ggplot(sample_info[ which(!is.na(sample_info$PTEN_protein)),], aes(x= PTEN_protein)) +
  geom_density(color = "blue", fill ="darkblue") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein (RPPA)") + 
  ylab("Density")
ggsave("DensityPlots/PTEN_protein_Density.pdf", heigh=6, width = 6)

#| Plotting PTEN protein vs PTEN CNA
my_comparisons <- list( c("-1", "0"), c("-2", "0"), c("0", "1") )
ggplot(sample_info[which(sample_info$PTEN_cna !="NP"),], aes(x =PTEN_cna, y= PTEN_protein, fill =PTEN_cna)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN CNA") + 
  ylab("PTEN protein (RPPA)")  +
  scale_fill_manual(values = viridis(5)) +
  labs(fill='PTEN CNA') +
  stat_compare_means(comparisons = my_comparisons,label ="p.format",size =5, label.family ="serif")
ggsave("BoxPlots/PTEN_protein_PTEN_CNA.pdf", heigh=6, width = 6)

#| Plotting PTEN protein vs PTEN CNA only complete loss and intact at genomic level
ggplot(sample_info[which(sample_info$PTEN_cna !="NP" & sample_info$PTEN_cna !="1"& sample_info$PTEN_cna !="2" & sample_info$PTEN_cna !="-1"),], aes(x =PTEN_cna, y= PTEN_protein, fill =PTEN_cna)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN CNA") + 
  ylab("PTEN protein (RPPA)")  +
  scale_fill_manual(values = viridis(5)) +
  labs(fill='PTEN CNA') +
  stat_compare_means(label ="p.format",family="serif", size =5)
ggsave("BoxPlots/PTEN_protein_PTEN_CNA_0_-2.pdf", heigh=6, width = 6)

#| PTEN mRNA PTEN CNA
my_comparisons <- list( c("-1", "0"), c("-2", "0"), c("0", "1") )
ggplot(sample_info[which(sample_info$PTEN_cna !="NP"),], aes(x =PTEN_cna, y= PTEN_mrna, fill =PTEN_cna)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN CNA") + 
  ylab("PTEN mRNA normalized log2 counts")  +
  scale_fill_manual(values = viridis(5)) +
  labs(fill='PTEN CNA') +
  stat_compare_means(comparisons= my_comparisons,label ="p.format",family="serif", size =5)
ggsave("BoxPlots/PTEN_protein_mRNA_CNA.pdf", heigh=6, width = 6)

#| PTEN mRNA PTEN CNA only -2 and 0
ggplot(sample_info[which(sample_info$PTEN_cna !="NP" & sample_info$PTEN_cna !="1"& sample_info$PTEN_cna !="2" & sample_info$PTEN_cna !="-1"),], aes(x =PTEN_cna, y= PTEN_mrna, fill =PTEN_cna)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN CNA") + 
  ylab("PTEN mRNA normalized log2 counts")  +
  scale_fill_manual(values = viridis(5)) +
  labs(fill='PTEN CNA') +
  stat_compare_means(label ="p.format",family="serif", size =5)
ggsave("BoxPlots/PTEN_protein_mRNA_CNA_0_-2.pdf", heigh=6, width = 6)

#| PTEN mRNA PTEN protein and CNA
data <- sample_info[which(sample_info$PTEN_cna =="-2" | sample_info$PTEN_cna =="0" ),]
data$`PTEN CNA` <- data$PTEN_cna
data <- data[which(!is.na(data$PTEN_mrna)),]
cor1 <- cor.test(data$PTEN_mrna[which(data$PTEN_cna =="0")], data$PTEN_protein[which(data$PTEN_cna =="0")], method ="pearson")
p1 <- round(cor1$p.value,2)
r1 <- round(as.numeric(cor1$estimate[[1]]),2)
cor2 <- cor.test(data$PTEN_mrna[which(data$PTEN_cna =="-2")], data$PTEN_protein[which(data$PTEN_cna =="-2")], method ="pearson")
p2 <- round(cor2$p.value,2)
r2<- round(as.numeric(cor2$estimate[[1]]),2)
ggplot(data, aes(x = PTEN_mrna, y=PTEN_protein, color=`PTEN CNA`, fill=`PTEN CNA`)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN mRNA normalized log2 counts") + 
  ylab("PTEN protein (RPPA)") +
  scale_color_manual(values =viridis(2)) +
  annotate(geom="text", x=10, y=-1.3, label=paste("r = ",r1,", p < ", p1, sep=""), color="black", family ="serif", size =3.8)+
  annotate(geom="text", x=11, y=1, label=paste("r = ",r2,", p < ", p2, sep=""), color="black", family ="serif", size =3.8)
ggsave("ScatterPlots/PTEN_protein_mRNA_Protein_CNA.pdf", heigh = 6, width = 6.5)

#| PTEN mRNA and protein with all CNA values
data <- sample_info[which(sample_info$PTEN_cna !="NP"),]
data$`PTEN CNA` <- data$PTEN_cna
data <- data[which(!is.na(data$PTEN_mrna)),]
ggplot(data, aes(x = PTEN_mrna, y=PTEN_protein, color=`PTEN CNA` , fill=`PTEN CNA` )) + 
  geom_point(size=2.5) + 
  #geom_smooth(method = "lm") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN mRNA normalized log2 counts") + 
  ylab("PTEN protein (RPPA)") +
  scale_color_manual(values =viridis(5)) #+
ggsave("ScatterPlots/PTEN_protein_mRNA_Protein_CNA_all.pdf", heigh = 6, width = 6.5)


#| PTEN mRNA PTEN protein only with no loss at genomic level
data <- sample_info[which(sample_info$PTEN_cna =="0" & !is.na(sample_info$PTEN_mrna)& !is.na(sample_info$`PI3K-AKT-mTOR`) ),]
cor <- cor.test(data$PTEN_mrna, data$PTEN_protein, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(x =PTEN_mrna, y= PTEN_protein, color =`PI3K-AKT-mTOR`)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="limegreen") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("CNA = 0")+
  xlab("PTEN mRNA normalized log2 counts") + 
  ylab("PTEN protein (RPPA)")  +
  scale_color_viridis(option="A") +
  annotate(geom="text", x=12.3, y=-0.3, label=paste("r = ",r2,", p < ", p2, sep=""), color="black", family ="serif", size =3.8)
ggsave("ScatterPlots/PTEN_protein_mRNA_Protein_PI3K-AKT-mTOR.pdf", heigh = 6, width = 7)

#| PTEN protein and PI3K-AKT-mTOR signature
data <- sample_info[which(!is.na(sample_info$`PI3K-AKT-mTOR`) & sample_info$PTEN_cna != "NP"),]
data$`PTEN CNA` <- data$PTEN_cna
cor <- cor.test(data$PTEN_protein, data$`PI3K-AKT-mTOR`, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(x =`PI3K-AKT-mTOR`, y= PTEN_protein, color =`PTEN CNA` )) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="limegreen") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  ylab("PTEN protein (RPPA)") +
  annotate(geom="text", x=0.3, y=-1.3, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8) +
  scale_color_manual(values =viridis(5))
ggsave("ScatterPlots/PTEN_protein_Protein_PI3K-AKT-mTOR.pdf", heigh = 6, width = 6)


#| PTEN mRNA and PI3K-AKT-mTOR signature
ggplot(sample_info, aes(x =`PI3K-AKT-mTOR`, y= PTEN_mrna)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="limegreen") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  ylab("PTEN mRNA normalized log2 counts") 
ggsave("ScatterPlots/PTEN_mRNA_PI3K-AKT-mTOR.pdf", heigh = 6, width = 6)


#| PTEN mRNA and PI3K-AKT-mTOR signature with CNA
data <- sample_info[which(sample_info$PTEN_cna != "NP"),]
cor <- cor.test(data$PTEN_mrna, data$`PI3K-AKT-mTOR`, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(x =PTEN_mrna , y= `PI3K-AKT-mTOR`, color =PTEN_cna)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="limegreen") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ylab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  xlab("PTEN mRNA normalized log2 counts") +
  labs(color ="PTEN CNA") + scale_color_manual(values =viridis(6)) +
  annotate(geom="text", x=8.5, y=-0.3, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8)
ggsave("ScatterPlots/PTEN_mRNA_PI3K-AKT-mTOR_CNA.pdf", heigh = 6, width = 6)


#| PTEN mRNA and PI3K-AKT-mTOR signature with CNA = 0
data <- sample_info[which(sample_info$PTEN_cna =="0"),]
cor <- cor.test(data$PTEN_mrna, data$`PI3K-AKT-mTOR`, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(x =PTEN_mrna, y= `PI3K-AKT-mTOR`,  color =PTEN_cna)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  ylab("PTEN mRNA normalized log2 counts") +
  labs(color ="PTEN CNA") +
  scale_color_manual(values= "darkgoldenrod3") +
  annotate(geom="text", x=12, y=-0.3, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8)
#  scale_color_viridis(option="A")
ggsave("ScatterPlots/PTEN_mRNA_PI3K-AKT-mTOR_CNA_0.pdf", heigh = 6, width = 6)

#| PTEN mRNA and PI3K-AKT-mTOR signature with CNA = -1
data <- sample_info[which(sample_info$PTEN_cna =="-1"),]
cor <- cor.test(data$PTEN_mrna, data$`PI3K-AKT-mTOR`, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(y =`PI3K-AKT-mTOR`, x= PTEN_mrna, color =PTEN_cna)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  ylab("PTEN mRNA normalized log2 counts") +
  labs(color ="PTEN CNA") +
  scale_color_manual(values= "darkgoldenrod3") +
  annotate(geom="text", x=9.5, y=0.4, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8)
ggsave("ScatterPlots/PTEN_mRNA_PI3K-AKT-mTOR_CNA_-1.pdf", heigh = 6, width = 6)

#| PTEN mRNA and PI3K-AKT-mTOR signature with CNA = -2
data <- sample_info[which(sample_info$PTEN_cna =="-2"),]
cor <- cor.test(data$PTEN_mrna, data$`PI3K-AKT-mTOR`, method ="pearson")
p <- round(cor$p.value,2)
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(data, aes(y =`PI3K-AKT-mTOR`, x= PTEN_mrna, color =PTEN_cna)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="midnightblue") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR z scores of normalized log2 counts") + 
  ylab("PTEN mRNA normalized log2 counts") +
  labs(color ="PTEN CNA") +
  scale_color_manual(values= "darkgoldenrod3")+
  annotate(geom="text", x=8, y=0.4, label=paste("r = ",r,", p < ", p, sep=""), color="black", family ="serif", size =3.8)
ggsave("ScatterPlots/PTEN_mRNA_PI3K-AKT-mTOR_CNA_-2.pdf", heigh = 6, width = 6)


#| PTEN mRNA  and PI3K-AKT-mTOR signature color PTEN protein
ggplot(sample_info[which(!is.na(sample_info$`PI3K-AKT-mTOR`)),], aes(x =`PI3K-AKT-mTOR`, y= PTEN_mrna, color =PTEN_protein)) + 
  geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="limegreen") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PI3K-AKT-mTOR signature") + 
  ylab("PTEN mRNA normalized log2 counts") 

################################################################################


################################################################################
#| Saving sample information
write.table(sample_info, "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_merge.txt", sep ="\t")
################################################################################









########################## STABILIZING COUNT VARIANCE ##########################

#| Transforming the data into a DESeq object (Run this line before performing PCA)
dds <- DESeqDataSetFromMatrix(countData = counts_data,
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
  
  if(nrow(counts_data) >= 30) { #| For large datasets, apply VST transformation. It 
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
  tag <- "TCGA"
  pdf(paste("VST/0_CountVariance_Stabilization_logNormalized_", tag,".pdf",sep=''));
  meanSdPlot(log.norm.counts, ranks=FALSE) 
  dev.off()
  pdf(paste("VST/0_CountVariance_Stabilization_VST_BlindTRUE_", tag,".pdf",sep=''));
  meanSdPlot(assay(vsd.blind), ranks=FALSE)
  dev.off()
  pdf(paste("VST/0_CountVariance_Stabilization_VST_BlindFALSE_", tag,".pdf",sep=''));
  meanSdPlot(assay(vsd.noblind), ranks=FALSE)
  dev.off()

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
write.table(assay(vsd.blind),"VST/FullCounts_TCGA_filtered_VST.txt", sep ="\t", row.names =T)
################################################################################
