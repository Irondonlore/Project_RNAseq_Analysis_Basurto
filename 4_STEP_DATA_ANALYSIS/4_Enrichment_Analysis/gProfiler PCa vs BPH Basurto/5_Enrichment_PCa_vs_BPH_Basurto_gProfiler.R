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


################################################################################
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(dplyr))
suppressMessages(library(viridis))

dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Enrichment_Analysis/gProfiler PCa vs BPH Basurto/"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PCa_vs_BPH/Results/design ~ DV200 + Edad + Diagnostico/Tables/resDESeq2_dds_DV200_Edad_Diagnostico_ContrastPCa_vs_BPH.txt"

setwd(dir.proj)

theme_set(theme_classic()) 
################################################################################



#################################### DATA ######################################
#| DEGs
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")

#| DEGs up
DEGs_genes_up <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange >=0),]

#| DEGs down
DEGs_genes_down <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange <0),]

#| DEGs all
DEGs_all <- DEGs_genes[which(DEGs_genes$DEGs == "Yes"),]
################################################################################



##############################   gProfiler  DEGs ###############################
#| Enrichment of all DGEs 
total_gost_DEGs_all <- gost(list("DEGs Design <- ~  1 + DV200 + Edad + Diagnostico. Contrast: PCa vs BPH" = DEGs_all$GeneID),
                            organism = "hsapiens", ordered_query = FALSE, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                            measure_underrepresentation = FALSE, evcodes = TRUE,
                            user_threshold = 0.05, correction_method = "fdr",
                            domain_scope = "custom", custom_bg = DEGs_genes$GeneID, 
                            numeric_ns = "", sources = NULL, as_short_link = FALSE)

gostplot(total_gost_DEGs_all, interactive = FALSE, capped=FALSE)
ggsave("Results/Plots/gost_plot_gprofiler_PCa_vs_BPH_DGEs_All.pdf", height = 5, width = 8.1)

results_DGEs_all <- total_gost_DEGs_all$result[order(total_gost_DEGs_all$result$p_value),]
results_DGEs_all$`Term name` <-  paste(results_DGEs_all$term_name, "\n (N = ",results_DGEs_all$term_size, ")",sep ="")

write.table(results_DGEs_all[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")],"Results/Tables/Table_results_gprofiler_PCa_vs_BPH_regulated_genes.txt", sep = "\t", row.names = T)

#results_DGEs_all <- read.table("Results/Tables/")
#| Top 15 most significant
ggplot(results_DGEs_all[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  ggtitle("Enrichment DGEs |FC| >= 2 and FDR < 0.05") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/Total_DEGs_Enrichment_PCa_vs_BPH.pdf", height = 6, width = 6)

#| Top with REAC
dat1_filtered <- results_DGEs_all[which(results_DGEs_all$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/REAC_Total_DEGs_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:BP
dat1_filtered <- results_DGEs_all[which(results_DGEs_all$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;BP_Total_DEGs_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:CC
dat1_filtered <- results_DGEs_all[which(results_DGEs_all$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;CC_Total_DEGs_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:MF
dat1_filtered <- results_DGEs_all[which(results_DGEs_all$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;MF_Total_DEGs_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with TF
dat1_filtered <- results_DGEs_all[which(results_DGEs_all$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/TF_Total_DEGs_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)
################################################################################




##############################  UP-REGULATED   #################################
total_DEGs_gost_up <- gost(list("UP-REGULATED. Design <- ~  1 + DV200 + Edad + Diagnostico. PCa vs BPH" = DEGs_genes_up$GeneID),
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = DEGs_genes$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Results/Plots/gost_plot_gprofiler_PCa_vs_BPH_Up_regulated_genes.pdf")
gostplot(total_DEGs_gost_up, interactive = FALSE, capped=FALSE)
dev.off()



results_DGEs_up <- total_DEGs_gost_up$result[order(total_DEGs_gost_up$result$p_value),]
results_DGEs_up$`Term name` <-  paste(results_DGEs_up$term_name, "\n (N = ",results_DGEs_up$term_size, ")",sep ="")

write.table(results_DGEs_up$result[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")],"Results/Tables/Table_results_gprofiler_PCa_vs_BPH_Up_regulated_genes.txt")

#| Top 15 most significant
ggplot(results_DGEs_up[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment (Top 15 most significant) all Up DGEs ") +
  xlab("Data")
ggsave("Results/Plots/Total_DEGs_UP_Enrichment_PCa_vs_BPH.pdf", height = 6, width = 8)

#| Top with REAC
dat1_filtered <- results_DGEs_up[which(results_DGEs_up$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/REAC_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:BP
dat1_filtered <- results_DGEs_up[which(results_DGEs_up$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;BP_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:CC
dat1_filtered <- results_DGEs_up[which(results_DGEs_up$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;CC_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)

#| Top with GO:MF
dat1_filtered <- results_DGEs_up[which(results_DGEs_up$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/GO;MF_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 5.5)


#| Top with MIR
dat1_filtered <- results_DGEs_up[which(results_DGEs_up$source == "MIRNA"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/MIRNA_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 3.5)
################################################################################




#############################  DOWN-REGULATED   ################################
total_DEGs_gost_down <- gost(list("DOWN-REGULATED. Design <- ~  1 + DV200 + Edad + Diagnostico. PCa vs BPH" = DEGs_genes_down$GeneID),
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = DEGs_genes$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Results/Plots/gost_plot_gprofiler_PCa_vs_BPH_Down_regulated_genes.pdf")
gostplot(total_DEGs_gost_down, interactive = FALSE, capped=FALSE)
dev.off()

results_DGEs_down <- total_DEGs_gost_down$result[order(total_DEGs_gost_down$result$p_value),]
results_DGEs_down$`Term name` <-  paste(results_DGEs_down$term_name, "\n (N = ",results_DGEs_down$term_size, ")",sep ="")

write.table(results_DGEs_down$result[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")],"Results/Tables/Table_results_gprofiler_PCa_vs_BPH_Down_regulated_genes.txt")

#| Top 15 most significant
ggplot(results_DGEs_down[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment (Top 15 most significant) all Up DGEs ") +
  xlab("Data")
ggsave("Results/Plots/Total_DEGs_DOWN_Enrichment_PCa_vs_BPH.pdf", height = 6, width = 8)


#| Top with REAC
dat1_filtered <- results_DGEs_down[which(results_DGEs_down$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  #ggtitle("Enrichment |FC| >= 2 and FDR < 0.05 REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Results/Plots/REAC_Total_DEGs_UP_PCa_vs_BPH_Enrichment.pdf", height = 4.5, width = 4.7)
################################################################################