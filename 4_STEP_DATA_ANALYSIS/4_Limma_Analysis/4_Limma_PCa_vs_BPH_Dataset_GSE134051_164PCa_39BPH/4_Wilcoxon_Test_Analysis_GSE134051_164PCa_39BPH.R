################################################################################
####           WILCOXON TEST TO COMPARE WITH DESEQ2 RESULTS                   ##
################################################################################

#| DESeq2 is a tool based on fitting the data to a Negative Binomial Distribution
#| to then test for differences in gene expression by applying a wald test. However,
#| there is recent publication where they suggest that DESeq2 exaggerates the number
#| of false positive when the sample size is large and that an appropiate test to use 
#| is wilcoxon test instead of wald test. 
#| https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4

#| This script is aim to perform Wilcoxon test to the GSE134051 data using
#| the code employed in the publication.
################################################################################



################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(DESeq2))
suppressWarnings(library(edgeR, quietly = T))
suppressMessages(library(viridis))
suppressMessages(library(ggplot2))
suppressMessages(library(GEOquery))
################################################################################

dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/4_Limma_PCa_vs_BPH_Dataset_GSE134051_164PCa_39BPH"
setwd(dir.proj)


#| Loading series and platform data from the GEO  
gse_GSE134051 <- getGEO("GSE134051", GSEMatrix =T, AnnotGPL = F)

#| Extracting sample information 
sample_info_GSE134051 <- pData(gse_GSE134051$GSE134051_series_matrix.txt.gz)

#| Gene annotation
geneAnnotation_GSE134051 <- fData(gse_GSE134051$GSE134051_series_matrix.txt.gz)

#| Gene expression data
geneExprs_GSE134051 <- exprs(gse_GSE134051$GSE134051_series_matrix.txt.gz) 
geneExprs_GSE134051 <- as.data.frame(geneExprs_GSE134051)
geneExprs_GSE134051$ENS_ID <- geneAnnotation_GSE134051$SPOT_ID

#| Filtering the NoEntry 
geneExprs_GSE134051 <- geneExprs_GSE134051[geneExprs_GSE134051$ENS_ID != "NoEntry",]

#| Aggregating repeated values: Computing the mean of the expression for repetitive genes
geneExprs_GSE134051 <- aggregate(geneExprs_GSE134051, by =list(geneExprs_GSE134051$ENS_ID), FUN =mean)

#| Assigning the row names 
rownames(geneExprs_GSE134051) <- geneExprs_GSE134051$Group.1

#| Deleting non-interesting columns
geneExprs_GSE134051 <- subset(geneExprs_GSE134051, select = -c(Group.1, ENS_ID))
geneExprs_GSE134051



################################################################################
#############                      WILCOXON TEST                    ############ 

conditions <- factor(t(sample_info_GSE134051["diagnosis:ch1"]))

#| Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(geneExprs_GSE134051),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(geneExprs_GSE134051[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues,method = "fdr")

#| Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- geneExprs_GSE134051[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- geneExprs_GSE134051[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))

#| Saving Results
output <- data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(output) <- rownames(geneExprs_GSE134051)
head(output)
write.table(output_cpm, "Results/Wilcoxon_results/wilcoxon_output_GSE134051_data.txt", sep ="\t")



################################################################################
###                    COMPARING WITH AC-45_RNAseq-FFPE DATA                 ###     
output_AC_45_RNAseq_FFPE <- read.table("../Results/Wilcoxon_results/wilcoxon_output_AC-45_RNAseq_FFPE_data.txt", sep ="\t")

data_frame_combined <- merge(output, output_AC_45_RNAseq_FFPE, by = "row.names")
head(data_frame_combined)

theme_set(theme_classic())

# Log2FC-Log2FC plot
ggplot(data_frame_combined, aes(x=log2foldChange.x, y =log2foldChange.y)) +
  geom_point() +theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Log2 FC GSE134051") +
  ylab("Log2 FC AC-45_RNAseq-FFPE data") +
  ggtitle("Log2FC-Log2FC plot by hand. Correlation = 0.54")+geom_smooth(method=lm, colour= "green")
cor(data_frame_combined$log2foldChange.x, data_frame_combined$log2foldChange.y)


# Histograms p-values
pvalues_GSE <-data.frame(pValues = data_frame_combined$pValues.x) 
pvalues_AC_45 <-data.frame(pValues=data_frame_combined$pValues.y)

pvalues_GSE$t <- "GSE134051" 
pvalues_AC_45$t <- "AC-45_RNAseq-FFPE"

hist_data <- rbind(pvalues_GSE, pvalues_AC_45)
head(hist_data)
ggplot(hist_data, aes(x=pValues, fill =t)) + geom_density(alpha = 0.88) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("p-values") + scale_fill_manual(values = c("aquamarine2", "midnightblue")) +
  ylab("Density")


#| DEGs

DEGs_up <- data_AC45_RNAseq_FFPE[(data_AC45_RNAseq_FFPE$padj <= 0.05) & !(is.na(data_AC45_RNAseq_FFPE$padj)) & (data_AC45_RNAseq_FFPE$log2FoldChange >=2.5) ,]
DEGs_down <- data_AC45_RNAseq_FFPE[(data_AC45_RNAseq_FFPE$padj <= 0.05) & !(is.na(data_AC45_RNAseq_FFPE$padj)) & (data_AC45_RNAseq_FFPE$log2FoldChange <= -2.5) ,]


uno <- data_frame_combined[data_frame_combined$FDR.x < 0.05 & data_frame_combined$log2foldChange.x>= 0.2,]
uno

dos <- data_frame_combined[data_frame_combined$FDR.y < 0.05 & data_frame_combined$log2foldChange.y >= 2,]
dos

tres <- data.frame(GeneID = unique(c(uno$Row.names, dos$Row.names)))
tres <- merge(tres,genome_GRCh39.94 , by ="GeneID")
length(tres$gene_name)


uno <- data_frame_combined[data_frame_combined$FDR.x < 0.05 & data_frame_combined$log2foldChange.x < (-0.1),]
dim(uno)

dos <- data_frame_combined[data_frame_combined$FDR.y < 0.05 & data_frame_combined$log2foldChange.y <= (-2),]
dos

l <- uno$Row.names[match(uno$Row.names, dos$Row.names)]

l[!is.na(l)]
tres <- data.frame(GeneID = l[!is.na(l)])
tres <- merge(tres,genome_GRCh39.94 , by ="GeneID")
length(tres$gene_name)

min(data_frame_combined$log2foldChange.y )
