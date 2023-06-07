################################################################################
####      FROM GEO TO DEGS AND ANALYSIS OF THE GSE134051 DATASET            ####
################################################################################
#| Description:
#| The GSE134051 contains the total RNA of RPE-derived freshly frozen tissue of 
#| 164 PCa and 39 BPH patients was assessed for gene expression profiling by custom 
#| microarrays, including the 48 samples analyzed previously by RNA-seq.

#| The paper related to the data GSE consisted on the testing of new drugs for PCa 
#| treatment. They found there is a potential lncRNA that can influence the expression 
#| of the tumor suppressor protein p53 which is involved in cell cycle processes. 
#| They explored the knockout of this gene on cell linages and created a potential 
#| drug which improve the treatment of xenograph mice with PCa what does not have 
#| too many effects in the neuronal system.  

#| The principal objective in this script is to analyse the GSE134051 database
#| and to compare it with the basurto's data (AC-45_RNAseq-FFPE):
#|    1) Comparing DEGs with Design: PCa vs BPH
################################################################################



################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(GEOquery))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(variancePartition))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(ashr))
suppressMessages(library(calibrate))
suppressMessages(library(vsn))
suppressMessages(library(genefilter))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(stringr))
suppressMessages(library(VennDiagram))
suppressMessages(library(gprofiler2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
################################################################################



################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/4_Limma_PCa_vs_BPH_Dataset_GSE134051_164PCa_39BPH"
setwd(dir.proj)
#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################



################################################################################
#| FROM GEO TO DIFFERENTIALLY EXPRESSED GENES
################################################################################

#| Loading series and platform data from the GEO  
gse_GSE134051 <- getGEO("GSE134051", GSEMatrix =T, AnnotGPL = F)

#| Extracting sample information 
sample_info_GSE134051 <- pData(gse_GSE134051$GSE134051_series_matrix.txt.gz)
head(sample_info_GSE134051)
#| Changing name of follow-up time in months:ch1 by DFS.Time and convert it to numeric
#| for variance partition analysis
sample_info_GSE134051$DFS.TIME <- sample_info_GSE134051$`follow-up time in months:ch1`
sample_info_GSE134051$DFS.TIME[is.na(sample_info_GSE134051$DFS.TIME)]<- 0
sample_info_GSE134051$DFS.TIME <- as.numeric(sample_info_GSE134051$DFS.TIME)

#| Gene annotation
geneAnnotation_GSE134051 <- fData(gse_GSE134051$GSE134051_series_matrix.txt.gz)

#| Gene expression data
geneExprs_GSE134051 <- exprs(gse_GSE134051$GSE134051_series_matrix.txt.gz) 
geneExprs_GSE134051 <- as.data.frame(geneExprs_GSE134051)
geneExprs_GSE134051$ENS_ID <- geneAnnotation_GSE134051$SPOT_ID

#| A checkpoint of the expression data (NOT NEEDED FOR THIS DATASET)
#quantile_geneExprs_GSE134051 <- as.numeric(quantile(geneExprs_GSE134051, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#LogC <- (quantile_geneExprs_GSE134051[5] > 100) || (quantile_geneExprs_GSE134051[6] - quantile_geneExprs_GSE134051[1] > 50 && quantile_geneExprs_GSE134051[2] > 0)

#| Filtering the NoEntry 
geneExprs_GSE134051 <- geneExprs_GSE134051[geneExprs_GSE134051$ENS_ID != "NoEntry",]

#| Aggregating repeated values: Computing the mean of the expression for repetitive genes
geneExprs_GSE134051 <- aggregate(geneExprs_GSE134051, by =list(geneExprs_GSE134051$ENS_ID), FUN =mean)

#| Assigning the row names 
rownames(geneExprs_GSE134051) <- geneExprs_GSE134051$Group.1

#| Deleting non-interesting columns
geneExprs_GSE134051 <- subset(geneExprs_GSE134051, select = -c(Group.1, ENS_ID))

#| Extracting PTEN expression
sample_info_GSE134051$PTEN <- unlist(geneExprs_GSE134051[rownames(geneExprs_GSE134051) == "ENSG00000171862",])

#| Saving sample data and counts data (This step may present some problems)
write.table(geneExprs_GSE134051, "Data/counts_data_GSE134051.txt", sep = "\t", row.names = T)
write.table(sample_info_GSE134051, "Data/sample_info_GSE134051.txt", sep = "\t", row.names =T)
################################################################################




################################################################################
#| DEGs with limma
################################################################################

#| Creating the  classification of PCa vs BPH: Group
sample_info_GSE134051$Group <- "PCa"
sample_info_GSE134051$Group[sample_info_GSE134051$`diagnosis:ch1` == "benign prostate hyperplasia"] <- "BPH"

#| Design 1: PCa vs BPH  #######################################################
design <- model.matrix(~Group + 0, sample_info_GSE134051)
colnames(design) <- c("BPH", "PCa")
head(design)

#| Linear fitting 
fit <- lmFit(geneExprs_GSE134051, design)  

#| Set up contrasts of interest (BPH vs PCa) and recalculate model coefficients
cts <- paste("PCa", "BPH", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

#| Compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
fit2

#| topTable (tT)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(fit2)[1])
head(tT)

#| Extracting only genes that matched once
str_genes <- str_extract(rownames(tT), "ENSG[0-9]{11};")
tT_filtered <-tT[-which(!is.na(str_genes)),]
tT_filtered


#| Transforming to gene symbol  ################################################
new_data_frame <- data.frame(GeneID = rownames(tT_filtered))

new_data_frame_merge <- merge(new_data_frame, genome_GRCh39.94, by ="GeneID")

rownames(genome_GRCh39.94) <- genome_GRCh39.94$GeneID
  
rownames(new_data_frame_merge) <- new_data_frame_merge$GeneID

tT_filtered_f <- tT_filtered[rownames(new_data_frame_merge),]

tT_filtered_f$gene_name <- new_data_frame_merge$gene_name

new_data_frame_merge <- new_data_frame_merge[order(new_data_frame_merge$GeneID),]


#| SAVING THE RESULTS  #########################################################
write.table(tT_filtered_f, "Data/tT_filtered_GSE134051.txt", sep = "\t")


#| VISUALIZATION  ##############################################################
dim(tT_filtered_f)
dim(tT_filtered)

#| Volcano plot
pdf("Results/VolcanoPlots/VolcanoPlot_GSE134051.pdf")
volcano <- EnhancedVolcano(tT_filtered,
                           lab = rownames(tT_filtered),
                           x = 'logFC',
                           y = 'adj.P.Val', selectLab = 0,
                           subtitle = "Differential expression PCa vs BPH",
                           col=c("#000004FF", "lightgoldenrod1","#721F81FF", "#F1605DFF" ))
volcano + theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=18, hjust = 0.5, face ="bold"))
dev.off()

#| Comparing PCa vs BPH with AC-45_RNAseq-FFPE data  ###########################

#| Loading the data
data_AC45_RNAseq_FFPE <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/Results/Tables/Design ~ 1 + DV200 + Diagnostico/resDESeq2_dds_DV200_Diagnostico_ContrastPCa_vs_BPH.txt", sep ="\t")
data_AC45_RNAseq_FFPE_filtered <- data_AC45_RNAseq_FFPE[rownames(tT_filtered_f),]
any(rownames(data_AC45_RNAseq_FFPE[rownames(tT_filtered_f),]) == rownames(tT_filtered_f))

#| LFC - LFC plot
theme_set(theme_classic())

data <- data.frame(data_AC45_RNAseq_FFPE_filtered$log2FoldChange)
colnames(data) <- "log2FoldChange"
data$logFC <- tT_filtered_f$logFC

#| Correlation
cor <- cor.test(data$log2FoldChange,data$logFC )

pdf("Results/log2FC/Log2FCGSE134051_AC-45_RNAseq-FFPE.pdf")
ggplot(data, aes(x=log2FoldChange, y = logFC)) + geom_point(size=2.5, colour = "midnightblue") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Log2FC AC-45_RNAseq-FFPE") +
  ylab("Log2FC GSE134051") +geom_smooth(method=lm, colour= "green") +
  ggtitle(paste("Correlation: ", round(cor$estimate, digits = 3), sep =""))
dev.off()


#| Resigning new threshold for log2foldchange
DEGs_up <- data_AC45_RNAseq_FFPE[(data_AC45_RNAseq_FFPE$padj <= 0.05) & !(is.na(data_AC45_RNAseq_FFPE$padj)) & (data_AC45_RNAseq_FFPE$log2FoldChange >=2.5) ,]
DEGs_down <- data_AC45_RNAseq_FFPE[(data_AC45_RNAseq_FFPE$padj <= 0.05) & !(is.na(data_AC45_RNAseq_FFPE$padj)) & (data_AC45_RNAseq_FFPE$log2FoldChange <= -2.5) ,]

DEGs_up$gene_name
DEGs_down$gene_name

DEGs_up$gene_name[order(DEGs_up$padj)]
DEGs_up$gene_name[order(DEGs_up$log2FoldChange, decreasing = TRUE)]
?order
#| DEGs up and down from tT_filtered
DEGs_up_tT <- tT_filtered_f[(tT_filtered_f$adj.P.Val <= 0.05) & !(is.na(tT_filtered_f$adj.P.Val)) & (tT_filtered_f$logFC >=1) ,]
DEGs_down_tT <- tT_filtered_f[(tT_filtered_f$adj.P.Val <= 0.05) & !(is.na(tT_filtered_f$adj.P.Val)) & (tT_filtered_f$logFC <= -1) ,]

DEGs_up_tT$gene_name
DEGs_down_tT$gene_name

dim(DEGs_up_tT)
dim(DEGs_down_tT)

#| Matching common genes: UP
common_genes_PCa_vs_BPH_up <- c(DEGs_up_tT$gene_name, DEGs_up$gene_name) 
common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))]
write.table(common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))], "Results/up-regulated_intersected_genes_GSE134051_AC-45_RNAseq-FFPE.txt", sep = "\t", col.names ="Genes")

#| Matching common genes: DOWN
common_genes_PCa_vs_BPH_down <- c(DEGs_down_tT$gene_name, DEGs_down$gene_name) 
common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))]
write.table(common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))], "Results/down-regulated_intersected_genes_GSE134051_AC-45_RNAseq-FFPE.txt", sep = "\t", col.names ="Genes")

#| VENN Diagram: UP
venn.diagram(x = list(DEGs_up_tT$gene_name, DEGs_up$gene_name),
             category.names = "VennDiagram_up",
             filename = "Results/VennDiagrams/Venn_diagramm_DEGs_Module_design_up.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c("midnightblue", "black"),
             fill = c(alpha("midnightblue",0.3), alpha("white",)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")

#| VENN Diagram: DOWN
venn.diagram(x = list(DEGs_down_tT$gene_name, DEGs_down$gene_name),
             category.names = "VennDiagram_down",
             filename = "Results/VennDiagrams/Venn_diagramm_DEGs_Module_design_down.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c("midnightblue", "black"),
             fill = c(alpha("midnightblue",0.3), alpha("white",)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")

DEGs_up$gene_name
#| ENRICHMENT OF COMMON GENES: gProfiler
#|    1) Up-regulated:
total_DEGs_gost <- gost(list("UP-REGULATED GENES Comparison PCa vs BPH GSE134051 and AC-45_RNAseq-FFPE data" = common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))]), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "custom", custom_bg = DEGs_up$gene_name, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)


pdf("Results/Enrichment_GSE134051_AC-45_RNAseq-FFPE/gProfiler_EnrichmentPlot_Up-regulated_genes_intersection_GSE134051.pdf")
gostplot(total_DEGs_gost, interactive = FALSE)
dev.off()
total_DEGs_gost$result$term_name[1:10]
l <- total_DEGs_gost$result$term_name[order(total_DEGs_gost$result$p_value)]
l[1:20]
0################################################################################


################################################################################
#| Checking normalization  (The GEO give the data normalized by Z-Score)
################################################################################
theme_set(theme_classic())
data_boxplot <- stack(geneExprs_GSE134051[,-c(1)])

extracted_geo <- sample_info_GSE134051[, c("geo_accession","Group")]
colnames(extracted_geo) <- c("ind", "Group")

merge_data_boxplot <- merge(data_boxplot, extracted_geo, by = "ind")
merge_data_boxplot <- merge_data_boxplot[order(merge_data_boxplot$Group),]

pdf("Results/BoxPlots/Normalization_BoxPlot_GSE134051.pdf")
ggplot(merge_data, aes(x =ind, y = values, color =Group, fill = Group)) + geom_boxplot(outlier.size=0.09)+ 
  theme(text=element_text(size=16,  family="serif"),axis.text.x=element_blank(), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Samples") + scale_colour_manual(values = c(alpha("darkgoldenrod3",0.3), alpha("midnightblue",0.3))) +
  scale_fill_manual(values = c("darkgoldenrod3", "midnightblue"))
dev.off()
################################################################################



################################################################################
#| Histogram of counts
################################################################################
pdf("Results/Histogram/Histograms_counts_GSE134051_data.pdf")
ggplot(data_boxplot, aes(x = values)) + geom_histogram(binwidth = 0.25, colour ="black", fill=alpha("midnightblue",0.9) ) +
  scale_fill_manual(values = c(alpha("darkgoldenrod3",0.3))) +
  theme(text=element_text(size=16,  family="serif"),axis.text.x=element_blank(), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Values") +
  ylab("Count") +
  ggtitle("Histogram of counts for GSE134051 dataset")
dev.off()
################################################################################


################################################################################
#| Sample correlation
################################################################################
dmat <- as.matrix(cor(geneExprs_GSE134051, method="spearman"))

#| PCa and BPH
#pdf("Results/SampleCorrelation/SampleCorrelation_PCa_and_BPH_counts_data.pdf")
pheatmap(dmat,border_color=NA,annotation_legend=T, fontsize = 6, fontsize_col = 1.5, fontsize_row = 1.5)
#dev.off()


################################################################################
#| PCA: By sample groups
################################################################################

#| To fix theme
theme_set(theme_classic()) 

pca <- prcomp(t(geneExprs_GSE134051))
percentVar <- 100*(pca$sdev^2/sum(pca$sdev^2))

#| PCA of PCa (Yes, No)
pdf("Results/PCA/PCA_PCa_GSE134051.pdf")
cbind(sample_info_GSE134051, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, color=Group)) + geom_point(size=2.5) + 
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",round(percentVar[1], digits=2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2], digits=2),"% variance")) +
  coord_fixed() + ggtitle("PCA Plot: PCa vs BPH") + 
  scale_colour_manual(values = c("darkgoldenrod3", "midnightblue"))
dev.off()

#| PTEN
pdf("Results/PCA/PCA_PTEN_GSE134051.pdf")
cbind(sample_info_GSE134051, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, color=PTEN)) + geom_point(size=2.5) + 
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste0("PC1: ",round(percentVar[1], digits=2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2], digits=2),"% variance")) +
  coord_fixed() + ggtitle("PCA: PTEN expression") + scale_color_viridis(discrete=F, option="A")
dev.off()


################################################################################
#| Variance Partition Analysis
################################################################################
form <- ~  (1|`diagnosis:ch1`) + DFS.TIME 
varPart <- fitExtractVarPartModel( geneExprs_GSE134051, form, sample_info_GSE134051)
vp <- sortCols( varPart )
plotPercentBars( vp[1:10,] )
pdf("Results/ViolinPlots/plotVarViolin_Diagnosis_GSE132051.pdf")
plotVarPart( vp, main = "Variance Partition Analysis")
dev.off()


################################################################################
#| BoxPlot
################################################################################

#| PTEN and Diagnostico
pdf("Results/BoxPlots/BoxPlot_PTEN_PCa_vs_BPH_GSE134051.pdf")
ggplot(sample_info_GSE134051, aes(x=Group, y=PTEN, fill = Group)) + geom_boxplot() + 
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("aquamarine2", "midnightblue")) +  ggtitle("PTEN Expression vs PCa and BPH Patients")+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)
dev.off()

