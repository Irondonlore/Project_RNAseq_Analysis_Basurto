################################################################################
################################ G PROFILER ####################################
################################################################################

#| This code is written to use gProfiler on the DEGs found on different designs
#| with DESeq2 to explore the differences between PTEN genomic loss vs intact in
#| TCGA for prostate cancer patients with primary tumor

#| In concrete: gProfileR is a tool for the interpretation of large gene lists which
#|can be run using a web interface or through R. The core tool takes a gene list 
#|as input and performs statistical enrichment analysis using hypergeometric testing
#|similar to clusterProfiler. Source: https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/gProfileR_REVIGO.html#:~:text=gProfileR%20is%20a%20tool%20for,hypergeometric%20testing%20similar%20to%20clusterProfiler.

#| You can use gProfiler for a wide selection of organisms, and the tool accepts 
#|your gene list as input. If your gene list is ordered (e.g. by padj. values), 
#|then gProfiler will take the order of the genes into account when outputting 
#|enriched terms or pathways.

#| The colors in the outputs of gProfiler represents the quality of the evidence, 
#|or how confidence was the functional annotation (weaker is depicted in Blue)

#|g:Profiler's best known functionality is the over-representation analysis to 
#|identify significantly enriched biological functions and pathways obtained from
#|well established data sources which include, among others, Gene Ontology (GO)
#| KEGG and Reactome.

#| The main type of data is Emsemble.

#|There are four main API wrapper functions in gprofiler2:
#|  1) gost for functional enrichment analysis
#|  2) gconvert for mapping gene identifiers between different namespaces
#|  3) gorth for mapping orthologous genes across species
#|  3) gsnpense for mapping SNP rs-IDs to chromosome positions, genes and variant effects

#| Definition of the parameters used in the gost function:
#|    *) organanism:
#|    *) ordered_query:
#|    *) multi_query:
#|    *) significant:
#|    *) exclude_iea:
#|    *) measure_underrepresentation:
#|    *) evcodes: 
#|    *) user_threshold:
#|    *) correction_method:
#|    *) domain_scope: 
#|    *) custom_bg:
#|    *) numeric_ns:
#|    *) sources:
#|    *) as_short_link:

#| Paper: https://f1000research.com/articles/9-709
#| Another source: https://biit.cs.ut.ee/gprofiler/page/r
################################################################################


############################ LIBRARIES AND DATA ################################
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(dplyr))
suppressMessages(library(viridis))

#| For plots
theme_set(theme_classic()) 

#| Data directories
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Enrichment_Analysis/gProfiler PTEN genomic loss vs intact TCGA/"
degs_tcga.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_TCGA_PTEN_genomic_loss_vs_intact/Results/Tables/res_DESeq2_TCGA_PTEN_genomic_loss_vs_intact_design_AGE_PTEN_CNA_no.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/NormalizedCounts/TCGA_raw_counts_filtered_normalized_log2.txt"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_merge.txt"
setwd(dir.proj)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Obtaining data
counts_data <- read.table(counts.file, sep ="\t", header=T)
sample_info <- read.table(info.file, sep ="\t", header=T)
degs_tcga <- read.table(degs_tcga.file, sep ="\t", header =T)

#| Changing names in colnames(counts_data) to match with sample_info$PATIENT_ID
colnames(counts_data) <- gsub("\\.", "-", colnames(counts_data))

#| Renaming rownames in the sample_info
rownames(sample_info) <- sample_info$PATIENT_ID

#| Verifying that rownames from sample_info and colnames from counts data match
any(rownames(sample_info) == colnames(counts_data))
################################################################################

dim(degs_tcga[which(degs_tcga$DEGs =="Yes"),])

########################### ENRICHMENT ANALYSIS ################################
total_gost_DEGs_tcga <- gost(list("DEGs PTEN genomic loss vs intact " = degs_tcga$GeneID[which(degs_tcga$DEGs=="Yes")]),
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(degs_tcga), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs_tcga,capped=FALSE)
ggsave("Results/EnrichmentPlots/gost_DEGs_TCGA_PTEN_genomic_loss_vs_intact_nofiltered.pdf", heigh =6, width =7)

total_gost_DEGs_tcga <- total_gost_DEGs_tcga$result[order(total_gost_DEGs_tcga$result$p_value),]
total_gost_DEGs_tcga$term_name[1:70]
total_gost_DEGs_tcga$`Term name` <- paste(total_gost_DEGs_tcga$term_name, "\n (N = ",total_gost_DEGs_tcga$term_size, ")",sep ="")

#| Top most significant
ggplot(total_gost_DEGs_tcga[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=12,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(size= "Intersection size")+
  ggtitle("Enrichment all DGEs (1527) PTEN genomic loss vs intact \n |FC| > 1.5 & FDR < 0.05") +
  xlab("Data")
ggsave("Results/EnrichmentPlots/Enrichment_All_DEGs_PTEN_genomic_loss_vs_intact_nofiltered.pdf", height = 6, width = 7)


#| Up regulated (NO RESULTS FOUND)
total_gost_DEGs <- gost(list("DEGs PTEN loss vs intact Genomic" = rownames(degs_tcga)[which(degs_tcga$DEGs=="Yes" & degs_tcga$log2FoldChange > 0 )]),
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(degs_tcga), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs,capped=FALSE)
total_gost_DEGs <- total_gost_DEGs$result[order(total_gost_DEGs$result$p_value),]
total_gost_DEGs$term_name[1:30]


#| Down regulated(NO RESULTS FOUND)
total_gost_DEGs <- gost(list("DEGs PTEN loss vs intact Genomic" = rownames(degs_tcga)[which(degs_tcga$DEGs=="Yes" & degs_tcga$log2FoldChange < 0 )]),
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(degs_tcga), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs,capped=FALSE)
total_gost_DEGs <- total_gost_DEGs$result[order(total_gost_DEGs$result$p_value),]
total_gost_DEGs$term_name[1:30]
################################################################################


