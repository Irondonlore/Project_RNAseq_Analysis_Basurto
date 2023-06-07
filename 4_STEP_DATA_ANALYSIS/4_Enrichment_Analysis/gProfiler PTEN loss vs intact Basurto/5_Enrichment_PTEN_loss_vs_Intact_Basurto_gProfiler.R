################################################################################
################################ G PROFILER ####################################
################################################################################

#| This code is written to use gProfiler on the DEGs found on different designs
#| with DESeq2 to explore the differences between PCa and BPH

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

#| Setting working dir
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Enrichment_Analysis/gProfiler PTEN loss vs intact Basurto/"
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Other parameters:
FDR <- 0.05   #| False Discovery Rate Threshold
FC <- 1.5      #| FC change Threshold (>= |log2(FC)| >= )
################################################################################



#################################### DATA ######################################
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PTEN_loss_vs_intact/Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_DV200_Edad_H_score_cut_0_ContrastPTEN_loss_vs_intact.txt"

#| DESeq2 DATA:
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")

#| Filtering those genes who satisfied the condition for DEGs classified by "Yes"
DEGs <- DEGs_genes[which(!is.na(DEGs_genes$padj) & (DEGs_genes$padj < FDR) & ((DEGs_genes$log2FoldChange > log2(FC)) | (DEGs_genes$log2FoldChange < (-log2(FC))))),]
DEGs_genes_up <- DEGs[which( DEGs$log2FoldChange >=0),]
DEGs_genes_down <- DEGs[which(DEGs$log2FoldChange <0),]
################################################################################


################################# ENRICHMENT ###################################
total_gost_DEGs <- gost(list("DEGs" = rownames(DEGs)), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(DEGs_genes), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs, interactive = FALSE, capped =FALSE) 
ggsave("Results/gost_PTEN_loss_vs_intact_DEGs.pdf", height = 5, width = 8.3)

results_DEGs <- total_gost_DEGs$result[order(total_gost_DEGs$result$p_value),]

results_DEGs$`Term name` <-  paste(results_DEGs$term_name, "\n (N = ",results_DEGs$term_size, ")",sep ="")

#| Top 10 most significant
ggplot(results_DEGs[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment DGEs |FC| >= 1.5 and FDR < 0.05") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Total_DEGs_Enrichment_PTEN_loss_vs_intact.pdf", height = 6, width = 6)

#| Top with GO:BP
dat1_filtered <- results_DEGs[which(results_DEGs$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=10, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment |FC| >= 1.5 and FDR < 0.05 GO:BP") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/GO;BP_Total_DEGs_PTEN_loss_vs_intact_Enrichment.pdf", height = 4.5, width = 4.3)

#| Top with GO:CC
dat1_filtered <- results_DEGs[which(results_DEGs$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment |FC| >= 1.5 and FDR < 0.05 GO:CC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/GO;CC_Total_DEGs_PTEN_loss_vs_intact_Enrichment.pdf", height = 4.5, width = 4)

#| Top with TF
dat1_filtered <- results_DEGs[which(results_DEGs$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment |FC| >= 1.5 and FDR < 0.05 TF") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/TF_Total_DEGs_PTEN_loss_vs_intact_Enrichment.pdf", height = 4.5, width = 5)

################################################################################




############################ ENRICHMENT FILTERED ###############################
total_gost_DEGs_filtered <- gost(list("DEGs" = rownames(DEGs_filtered)), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(DEGs_genes_filtered), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs_filtered, interactive = FALSE, capped =FALSE) 
ggsave("Results/gost_PTEN_loss_vs_intact_DEGs_filtered.pdf", height = 6, width = 7)

results_DEGs_filtered <- total_gost_DEGs_filtered$result[order(total_gost_DEGs_filtered$result$p_value),]

results_DEGs_filtered$`Term name` <-  paste(results_DEGs_filtered$term_name, "\n (N = ",results_DEGs_filtered$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_DEGs_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment DGEs |FC| >= 1.5 and FDR < 0.05") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Total_DEGs_Enrichment_PTEN_loss_vs_intact_filtered.pdf", height = 6, width = 6)

#| Top with GO:CC
dat1_filtered <- results_DEGs_filtered[which(results_DEGs_filtered$source == "GO:CC"),]
ggplot(dat1_filtered[1:5,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment |FC| >= 1.5 and FDR < 0.05 GO:CC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/GO;CC_Total_DEGs_PTEN_loss_vs_intact_Enrichment_filtered.pdf", height = 6, width = 5)


#| Volcano plot
DEGs_genes_filtered$delabel <- NA
DEGs_genes_filtered$delabel[which((DEGs_genes_filtered$log2FoldChange < (-log2(2.5)) | DEGs_genes_filtered$log2FoldChange > (log2(2.5))) & DEGs_genes_filtered$padj < FDR)] <- 
  DEGs_genes_filtered$gene_name[which((DEGs_genes_filtered$log2FoldChange < (-log2(2.5)) | DEGs_genes_filtered$log2FoldChange > (log2(2.5))) & DEGs_genes_filtered$padj < FDR)]

DEGs_genes_filtered$DEGs <- NA
DEGs_genes_filtered$DEGs[which((DEGs_genes_filtered$log2FoldChange > log2(FC)) & DEGs_genes_filtered$padj < FDR)] <- "UP"
DEGs_genes_filtered$DEGs[which((DEGs_genes_filtered$log2FoldChange < (-log2(FC))) & DEGs_genes_filtered$padj < FDR)] <- "DOWN"
DEGs_genes_filtered$DEGs[which(is.na(DEGs_genes_filtered$DEGs))] <- "No DEGs"

ggplot(DEGs_genes_filtered, aes(x = log2FoldChange, y = -log10(padj), color =DEGs, label=delabel )) + 
  geom_point(size =4,alpha = 0.5) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_color_manual(values = c("midnightblue","darkgrey", "mediumspringgreen") ,guide = "legend") +
  labs(color ="Module color") +
  xlab("Log2 Fold Change") +
  ggtitle("Volcano plot |FC| >= 1.5 and FDR < 0.05 \n Total DEGs PTEN loss vs intact filtered") +
  geom_text_repel(family ="serif") +xlim(-10,10)
ggsave("Results/Volcano_DGEs_PTEN_loss_vs_intact_filtered.pdf", height = 6, width = 6)
################################################################################

