################################################################################
#                                       #
################################################################################

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
suppressMessages(library(DBI))
suppressMessages(library(RMySQL))
suppressMessages(library(viridis))
suppressMessages(library(VennDiagram))
################################################################################


################################################################################
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Correlation_Gene_DFSTIME/"
setwd(workingDir)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Setting themes for plotting
theme_set(theme_classic())

#| Connecting to MySql database
con <- dbConnect(RMySQL::MySQL(), user='geneticanalyses', password = 'geneticanalyses', host = 'binf-web.cicbiogune.int',port=13306, dbname='geneticanalyses')
dbListTables(con)

#| padj value threshold 
p_threshold <- 0.05
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

#| Filtering BPH and non recurrence cases
data_phenotype_basurto <- data_phenotype_basurto[which(data_phenotype_basurto$Diagnostico == "PCa" & !is.na(data_phenotype_basurto$DFS.STATUS) & data_phenotype_basurto$DFS.STATUS == 1),]
data_expression_basurto <- data_expression_basurto[ ,rownames(data_phenotype_basurto)]


#| To detect outliers ENSG00000251562
data_expression_basurto[which(rowSums(data_expression_basurto>100000) >= 1),]


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

#| Empty vectors to save the p values
p_value_basurto <- c()
cor_basurto <- c()

#| For over all the genes
for (i in 1:dim(data_expression_basurto)[1]){
  
  cor <- cor.test(data_phenotype_basurto$DFS.TIME, as.numeric(data_expression_basurto[i,]), method ="spearman")
  
  #| Correlation value
  cor_basurto <- c(cor_basurto, as.numeric(cor[[4]]))
  
  #| P-value
  p_value_basurto <- c(p_value_basurto, cor[[3]])
  
  }

#| Plotting distribution of Correlation
data_cor_basurto <- data.frame(cor = cor_basurto)
pdf("Results/Images/cor_histogram/correlation_histogram_Density_Basurto_Correlation.pdf")
ggplot(data_cor_basurto, aes(x = cor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(expression(rho)) + ylab("Density")+ ggtitle("correlation values Basurto")
dev.off()

#| Plotting distribution of p-values
data_pvalues_basurto <- data.frame(pvalues = p_value_basurto)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_Basurto_Correlation.pdf")
ggplot(data_pvalues_basurto, aes(x = pvalues)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from correlation analysis Basurto")
dev.off()


#| False discovery rate
fdr_basurto <- p.adjust(p_value_basurto, method = "fdr")

#| Unique data frame
dataframe_basurto <- data.frame(GeneID = rownames(data_expression_basurto), p_values_basurto = p_value_basurto, padj_basurto = fdr_basurto, cor_basurto = cor_basurto)

#| Merging for having gene names
dataframe_basurto <- merge(dataframe_basurto, genome_GRCh39.94, by="GeneID")

#| Sorting by p-values value
dataframe_basurto <- dataframe_basurto[order(dataframe_basurto$padj), ]

#| Exploring the dataframes
head(dataframe_basurto)

#| Saving results from Cox analysis in Basurto
write.table(dataframe_basurto[order(dataframe_basurto$GeneID), ],"Results/Tables/Correlation_Analysis_Genes_DFSTIME_Basurto.txt" ,sep ="\t")

##| cor as a function of the p-value BASURTO
cor_pvalue_basurto <- data.frame( p_values = dataframe_basurto$padj, cor = dataframe_basurto$cor)
pdf("Results/Images/cor_p_value_plot/correlation_p_values_basurto.pdf")
ggplot(cor_pvalue_basurto, aes(x=p_values, y = cor)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab(expression(rho)) +
  ggtitle("Correlation vs p-adjusted-values Basurto")+
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()
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
data_phenotype_taylor <- data_phenotype_taylor[ which(data_phenotype_taylor$tumortype == "Primary Tumor" & data_phenotype_taylor$DFS.STATUS==1),]
data_expression_taylor <- data_expression_taylor[, data_phenotype_taylor$Sample] 

#| Checking order and sample ID
length(colnames(data_expression_taylor)) == length(data_phenotype_taylor$Sample)
any(colnames(data_expression_taylor) == data_phenotype_taylor$Sample)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_taylor <- c()
cor_taylor <- c()

#| For over all the genes
for (i in 1:dim(data_expression_taylor)[1]){
  
  cor <- cor.test(data_phenotype_taylor$DFS.TIME, as.numeric(data_expression_taylor[i,]), method ="spearman")
  
  #| Correlation value
  cor_taylor <- c(cor_taylor, as.numeric(cor[[4]]))
  
  #| P-value
  p_value_taylor <- c(p_value_taylor, cor[[3]])
  
}

#| Plotting distribution of Correlation
data_cor_taylor <- data.frame(cor = cor_taylor)
pdf("Results/Images/cor_histogram/correlation_histogram_Density_Taylor_Correlation.pdf")
ggplot(data_cor_taylor, aes(x = cor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(expression(rho)) + ylab("Density")+ ggtitle("correlation values Taylor")
dev.off()

#| Plotting p-values
data_pvalues_taylor <- data.frame(pvalues = p_value_taylor)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_Taylor_Correlation.pdf")
ggplot(data_pvalues_taylor, aes(x = pvalues)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ 
  ggtitle("p-values from correlation analysis Taylor")
dev.off()

#| False discovery rate
fdr_taylor <- p.adjust(p_value_taylor, method = "fdr")

#| Unique data frame
dataframe_taylor <- data.frame(refseq_mrna = rownames(data_expression_taylor), p_values_taylor = p_value_taylor, padj_taylor = fdr_taylor, cor_taylor = cor_taylor)

#| GeneID
dataframe_taylor <- merge(dataframe_taylor, data_annotation_taylor, by = "refseq_mrna")

#| Sorting by p-values value
dataframe_taylor <- dataframe_taylor[order(dataframe_taylor$p_values), ]

#| Exploring the dataframes
head(dataframe_taylor)

##| cor as a function of the p-value TAYLOR
cor_pvalue_taylor <- data.frame( p_values = dataframe_taylor$padj, cor = dataframe_taylor$cor)
pdf("Results/Images/cor_p_value_plot/correlation_p_values_taylor.pdf")
ggplot(cor_pvalue_taylor, aes(x=p_values, y = cor)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-values") +
  ylab(expression(rho)) +
  ggtitle("Correlation vs p-adjusted-values Taylor")+
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()
################################################################################




################################################################################
#######| TCGA Fire Horse Legacy data
################################################################################
data_phenotype_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostatephenotype")
data_expression_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostateexpression")

#| Fixing gene name
data_expression_tcga$Gene <- gsub(".{1}$","",data_expression_tcga$Gene)
rownames(data_expression_tcga) <- data_expression_tcga$Gene

#| Deleting Gene column
data_expression_tcga <- data_expression_tcga[,-c(1)]

#| Filtering non recurrence samples
data_phenotype_tcga <- data_phenotype_tcga[which(data_phenotype_tcga$DFS.STATUS == 1),]
data_expression_tcga <- data_expression_tcga[,data_phenotype_tcga$Sample]

#| Checking order and sample ID
length(colnames(data_expression_tcga[,1:dim(data_expression_tcga)[2]])) == length(data_phenotype_tcga$Sample)
any(colnames(data_expression_tcga[,1:dim(data_expression_tcga)[2]]) == data_phenotype_tcga$Sample)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_tcga <- c()
cor_tcga <- c()

#| For over all the genes
for (i in 1:dim(data_expression_tcga)[1]){
  
  cor <- cor.test(data_phenotype_tcga$DFS.TIME, as.numeric(data_expression_tcga[i,]), method ="spearman")
  
  #| Correlation value
  cor_tcga <- c(cor_tcga, as.numeric(cor[[4]]))
  
  #| P-value
  p_value_tcga <- c(p_value_tcga, cor[[3]])
}


#| Plotting distribution of Correlation
data_cor_tcga <- data.frame(cor = cor_tcga)
pdf("Results/Images/cor_histogram/correlation_histogram_Density_TCGA_Correlation.pdf")
ggplot(data_cor_tcga, aes(x = cor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(expression(rho)) + ylab("Density")+ ggtitle("correlation values TCGA")
dev.off()

#| Plotting p-values
data_pvalues_tcga <- data.frame(pvalues = p_value_tcga)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_TCGA_Correlation.pdf")
ggplot(data_pvalues_tcga, aes(x = pvalues)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from correlation analysis TCGA")
dev.off()

#| False discovery rate
fdr_tcga <- p.adjust(p_value_tcga, method = "fdr")

#| Unique data frame
dataframe_tcga <- data.frame(gene_name = rownames(data_expression_tcga), p_values_tcga = p_value_tcga, padj_tcga = fdr_tcga, cor_tcga = cor_tcga)

#| Sorting by p-values value
dataframe_tcga <- dataframe_tcga[order(dataframe_tcga$padj), ]

#| Exploring the dataframes
head(dataframe_tcga)

##| cor as a function of the p-value TCGA
cor_pvalue_tcga <- data.frame( p_values = dataframe_tcga$p_values, cor = dataframe_tcga$cor)
pdf("Results/Images/cor_p_value_plot/correlation_p_values_TCGA.pdf")
ggplot(cor_pvalue_tcga, aes(x=p_values, y = cor)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab(expression(rho)) +
  ggtitle("Correlation vs p-adjusted-values TCGA")+
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()
################################################################################



################################################################################
#######| Glinsky
################################################################################
data_phenotype_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyphenotype")
data_expression_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyexpression")

#| Giving gene name to row names
rownames(data_expression_glinsky) <- data_expression_glinsky$Gene

#| Deleting column "Gene" form gene expression data
data_expression_glinsky <- data_expression_glinsky[,-c(1)]

#| Filtering out non recurrence cases
data_phenotype_glinsky <-data_phenotype_glinsky[which(data_phenotype_glinsky$DFS.STATUS == 1),]
data_expression_glinsky <- data_expression_glinsky[,data_phenotype_glinsky$Sample]

#| Checking order and sample ID
length(colnames(data_expression_glinsky[,1:dim(data_expression_glinsky)[2]])) == length(data_phenotype_glinsky$Sample)
any(colnames(data_expression_glinsky[,1:dim(data_phenotype_glinsky)[2]]) == data_phenotype_glinsky$Sample)

#| Empty vectors to save the p values and information of the hazard ratio from the cox analysis
p_value_glinsky <- c()
cor_glinsky <- c()

#| For over all the genes
for (i in 1:dim(data_expression_glinsky)[1]){
  
  cor <- cor.test(data_phenotype_glinsky$DFS.TIME, as.numeric(data_expression_glinsky[i,]), method ="spearman")
  
  #| Correlation value
  cor_glinsky <- c(cor_glinsky, as.numeric(cor[[4]]))
  
  #| P-value
  p_value_glinsky <- c(p_value_glinsky, cor[[3]])
}

#| Plotting distribution of Correlation
data_cor_glinsky <- data.frame(cor = cor_glinsky)
pdf("Results/Images/cor_histogram/correlation_histogram_Density_Glinsky_Correlation.pdf")
ggplot(data_cor_glinsky, aes(x = cor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(expression(rho)) + ylab("Density")+ ggtitle("correlation values Glinsky")
dev.off()

#| Plotting p-values
data_pvalues_glinsky <- data.frame(pvalues = p_value_glinsky)
pdf("Results/Images/p_values_histogram/p_values_histogram_Density_Glinsky_Correlation.pdf")
ggplot(data_pvalues_glinsky, aes(x = pvalues)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)  +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + ylab("Density")+ ggtitle("p-values from Cox analysis Glinsky")
dev.off()

#| False discovery rate
fdr_glinsky <- p.adjust(p_value_glinsky, method = "fdr")

#| Unique data frame
dataframe_glinsky <- data.frame(gene_name = rownames(data_expression_glinsky), p_values_glinsky = p_value_glinsky, padj_glinsky = fdr_glinsky, cor_glinsky= cor_glinsky)

#| Sorting by p-values value
dataframe_glinsky <- dataframe_glinsky[order(dataframe_glinsky$p_values), ]

#| Exploring the dataframes
head(dataframe_glinsky)

##| cor as a function of the p-value GLINSKY
cor_pvalue_glinsky <- data.frame( p_values = dataframe_glinsky$p_values, cor = dataframe_glinsky$cor)
pdf("Results/Images/cor_p_value_plot/correlation_p_values_Glinsky.pdf")
ggplot(cor_pvalue_glinsky, aes(x=p_values, y = cor)) + geom_point( shape =21,colour = "black", size =2) +
  xlab("P-adjusted-values") +
  ylab(expression(rho)) +
  ggtitle("Correlation vs p-adjusted-values Glinsky")+
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
dev.off()
################################################################################



################################################################################
#| Saving results in table
################################################################################

#| Saving all dataframes in an unique list
list_dataframes <- list(dataframe_basurto[,-c(1,5,6,7,8,9)], dataframe_taylor[,-c(1,5)], dataframe_tcga, dataframe_glinsky)

#| Merging them up
final_table <- Reduce(function(x, y) merge(x, y, all=TRUE), list_dataframes) 

nonNAN_datasets <- c()
sig_datasets <- c()
up_sig_datasets <- c()
down_sig_datasets <- c()
total_coherence_halfdatasets <- c()

for (i in 1:dim(final_table)[1]){
  
  p_values <- as.numeric(final_table[i,c("p_values_basurto", "p_values_taylor", "p_values_tcga", "p_values_glinsky")])
  cor <- as.numeric(final_table[i,c("cor_basurto", "cor_taylor", "cor_tcga", "cor_glinsky")])
  
  #| non NAN
  nonNAN_datasets <- c(nonNAN_datasets, length(p_values[which(!is.na(p_values))]))
  
  #| In how many datasets the gene is significant
  sig_datasets <- c(sig_datasets, length(p_values[which(p_values < p_threshold)]))
  
  #| Worse significant
  up_sig_datasets <- c(up_sig_datasets, length(cor[which(p_values < p_threshold & cor > 0)]))
  
  #| Better significant
  down_sig_datasets <- c(down_sig_datasets, length(cor[which(p_values < p_threshold &  cor< 0)]))
  
  #| Total coherence
  if(length(cor[which(p_values < p_threshold & cor > 0)]) >= ceiling(length(p_values[which(!is.na(p_values))])/2)){
    total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "YES_UP")
  }else{
    if(length(cor[which(p_values < p_threshold & cor < 0)]) >= ceiling(length(p_values[which(!is.na(p_values))])/2)){
      total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "YES_DOWN")
    }else{
      total_coherence_halfdatasets <- c(total_coherence_halfdatasets, "NO")
    }
  }
}

#| cbinding columns
final_table <- cbind(final_table, nonNAN_datasets, sig_datasets, up_sig_datasets, down_sig_datasets, total_coherence_halfdatasets)

#| Saving final table
write.table(final_table, "Results/Tables/Final_table_Correlation_Analysis_Genes_DFSTIME_Basurto_Taylor_TCGA_Glinsky.txt", sep = "\t")
################################################################################





################################################################################
#| Setting a threshold and finding common genes from the correlation analysis
################################################################################


#venn.diagram(x = list(genes_basurto_negative, genes_taylor_negative, genes_tcga_negative, genes_glinsky_negative),
#             category.names = c("genes_basurto_negative", "genes_taylor_negative", "genes_tcga_negative", "genes_glinsky_negative"),
#             filename = "VennProof_Correlation_Negative.png",
#             output=FALSE,
#             imagetype="png" ,
#             height = 980 , 
#             width = 980 , 
#             resolution = 400,
#             compression = "lzw",
#             # Circles
#             lwd = 2,
#             col = c("green", "black", "yellow", "magenta"),
#             fill = c(alpha("green",0.3), alpha("white",), alpha("yellow",0.3),alpha("magenta",0.3)),
#             
#             # Numbers
#             cex = .6,
#             fontface = "bold",
#             fontfamily = "sans",
#             
#             # Set names
#             cat.cex = 0.5,
#             cat.fontface = "bold",
#             cat.default.pos = "outer",
#             cat.fontfamily = "serif")


#venn.diagram(x = list(genes_basurto_positive, genes_taylor_positive, genes_tcga_positive, genes_glinsky_positive),
#             category.names = c("genes_basurto_positive", "genes_taylor_positive", "genes_tcga_positive", "genes_glinsky_positive"),
#             filename = "VennProof_Correlation_Positive.png",
#             output=FALSE,
#             imagetype="png" ,
#             height = 980 , 
#             width = 980 , 
#             resolution = 400,
#             compression = "lzw",
#             # Circles
#             lwd = 2,
#             col = c("green", "black", "yellow", "magenta"),
#             fill = c(alpha("green",0.3), alpha("white",), alpha("yellow",0.3),alpha("magenta",0.3)),
#             
#             # Numbers
#             cex = .6,
#             fontface = "bold",
#             fontfamily = "sans",
#             
#             # Set names
#             cat.cex = 0.5,
#             cat.fontface = "bold",
#             cat.default.pos = "outer",
#             cat.fontfamily = "serif")
#
