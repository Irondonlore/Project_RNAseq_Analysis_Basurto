################################################################################
####           WILCOXON TEST TO COMPARE WITH DESEQ2 RESULTS                   ##
################################################################################

#| DESeq2 is a tool based on fitting the data to a Negative Binomial Distribution
#| to then test for differences in gene expression by applying a wald test. However,
#| there is recent publication where they suggest that DESeq2 exaggerates the number
#| of false positive when the sample size is large and that an appropiate test to use 
#| is wilcoxon test instead of wald test. 
#| https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4

#| This script is aim to perform Wilcoxon test to the AC-RNAseq-FFPE data using
#| the code employed in the publication.
################################################################################



################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(DESeq2))
suppressWarnings(library(edgeR, quietly = T))
suppressMessages(library(viridis))
suppressMessages(library(ggplot2))
################################################################################



################################################################################
################################## DATA ########################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/Results/Tables/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE.txt"
counts.file_filtered <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE_cpm_filtered_PCa_BPH.txt"
setwd(dir.proj)

#| Readtable Full counts. NOTE: it is not necessary to do any normalization 
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)
counts_data_filtered <- read.table(counts.file_filtered, sep = "\t", stringsAsFactors = TRUE)
hist(unlist(counts_data_filtered), breaks =10090,xlim = c(0,5000))
#| Readtable Sample information
sample_info <- read.table(info.file, sep ="\t")
conditions <- factor(t(sample_info["Diagnostico"]))

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################


mean(sample_info$Edad[sample_info$Diagnostico=="BPH"])
mean(sample_info$Edad[sample_info$Diagnostico=="PCa"])

################################################################################
#############   WILCOXON TEST FOR NORMALIZED COUNTS DATA WITH CPM   ############ 

#| Counts matrix processed with edgeR package
y_cpm <- DGEList(counts=counts_data, group=conditions)

#| Remove rows consistently have zero or very low counts (min counts = 10)
keep <- filterByExpr(y_cpm)
y_cpm <- y_cpm[keep,keep.lib.sizes=FALSE]

#| Perform TMM normalization and transfer to CPM (Counts Per Million). This step
#| it is recommended by researches because wilcoxon it is not a regression based 
#| model and thus it cannot adjust for possible cofunding factors
y_cpm <- calcNormFactors(y_cpm,method="TMM")
count_norm_cpm <-cpm(y_cpm)
count_norm_cpm<-as.data.frame(count_norm_cpm)

#| Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm_cpm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm_cpm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues,method = "fdr")

#| Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm_cpm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm_cpm[,c(which(conditions==conditionsLevel[2]))]
foldChanges_cpm <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))

#| Saving Results
output_cpm <- data.frame(log2foldChange=foldChanges_cpm, pValues=pvalues, FDR=fdr)
rownames(output_cpm) <- rownames(count_norm_cpm)
write.table(output_cpm, "Results/Wilcoxon_results/wilcoxon_output_AC-45_RNAseq_FFPE_data.txt", sep ="\t")




##########################   Visualization   ###################################

par(mfrow = c(2,1))
hist(as.numeric(output_cpm$pValue))
hist(as.numeric(output_vst$pValue))


par(mfrow = c(1,1))
hist(as.numeric(output_cpm$FDR))
hist(as.numeric(output_vst$FDR))

m <- merge( output_cpm, output_vst, by = "row.names")
plot(m$FDR.x, m$FDR.y)


############################ COMPARING WITH DESEQ2 #############################

compare <- data.frame(DESeq2_WF = 10232, 
                      DESeq_F = 28576, 
                      Wilcoxon_NF = 8629, 
                      Wilcoxon_F_CPM = 7337, 
                      WIlcoxon_F_vst = 41)

compare <- data.frame(DEGs = c(6329, 4095, 3338 , 3301),
                      Method = c("DESeq2 without F", "DESeq F", "Wilcoxon CPM F", "Wilcoxon CPM F Phe"))

theme_set(theme_classic()) 

pdf("Results/Wilcoxon_vs_DESeq2.pdf")
ggplot(data=compare, aes(x = Method, y = DEGs, fill = Method)) + geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=60,hjust=1), text=element_text( size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.3, face ="bold")) +
  scale_fill_viridis(discrete=T, option="A") + ggtitle("Differentially Expressed genes with DESeq2 and Wilcoxon")
dev.off()



#| Venn Diagram 

venn.diagram(x = list(modules$gene_name, DEGs_down$gene_name, DEGs_down_tT$gene_name),
             category.names = c(sprintf("Module %s", color), "DOWN DEGs Basurto", "DOWN DEGs GSE134051"),
             filename = "Results/VennDiagram/Diagnostico/Venn_Turquoise_DOWN_DEGs_Basurto_and_GSE134051.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c(color, "black", "yellow"),
             fill = c(alpha(color,0.3), alpha("white",), alpha("yellow",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")