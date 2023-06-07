################################################################################
#####             ENRICHMENT ANALYSIS OF MODULES DE FROM WGCNA             #####                        
################################################################################

#| WGCNA detected 12 modules in the data of the 147 PCa patients with PTEN loss 
#| vs intact detected by IHC. In this code, I have analyse the most correlated
#| genes with our traits of interest. I have perform Enrichment analysis with
#| gprofiler and I have analyzed the enrichment on each categories. 

#| I have also intersect the DEGs from the comparison PTEN loss vs intact to each
#| of our modules to observe the fraction of DEGs in every module. The DEGs are 
#| those from Basurto and TCGA (in the case of TCGA the comparison was at the genomic
#| level)

#| Moreover, I have selected a list of non-coding genes that are known in literature
#| to regulate PTEN to find if those non-codings are present in the modules that
#| were enriched in MIRNAs

################################################################################


############################  LIBRARIES AND DATA  ##############################
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(viridis))
suppressMessages(library(readxl))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(readxl))

#| For plots
theme_set(theme_classic())

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("W:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Pathways directories
workingDir <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/"
data.file_WGCNA <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_new.txt"
data.file_DEGs <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PTEN_loss_vs_intact/Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_DV200_Edad_H_score_cut_0_ContrastPTEN_loss_vs_intact.txt"
ptenloss_signature.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Signature_PTEN_loss_Paper_Imada_et_al_BMC_cancer.xlsx"
non_coding.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/NonCodingGenes.xlsx"
setwd(workingDir)

#| WGCNA data
wgcna_Modules <- read.table(data.file_WGCNA, sep = "\t")
module_colors <- unique(wgcna_Modules$moduleColors)

#| PTEN loss signature from Imada et al.
PTEN_loss_signatue <- read_xlsx(ptenloss_signature.file, sheet= "S5_TCGA_signature_ERGfusion")

#| Data contaning non coding genes in the human genome
nonCoding_genes <- read_xlsx(non_coding.file) 

#| DEGs by comparing PTEN loss vs intact Basurto
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")
DEGs <- DEGs_genes[which((DEGs_genes$log2FoldChange > log2(1.5) | DEGs_genes$log2FoldChange < (-log2(1.5))) & DEGs_genes$padj < 0.05 & !is.na(DEGs_genes$padj < 0.05)),]
################################################################################


#################################### Green #####################################
genes_green <- rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "green"]
total_gost_green <- gost(list("Module Green WGCNA" = genes_green), 
                         organism = "hsapiens", ordered_query = FALSE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = TRUE,
                         user_threshold = 0.05, correction_method = "fdr",
                         domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots new/ggplot_results_WGCNA_Module_green.pdf")
gostplot(total_gost_green, interactive = FALSE, capped =FALSE)
dev.off()

results_green <- total_gost_green$result[order(total_gost_green$result$p_value),]

write.table(results_green[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment tables new/ggprofiler_results_Enrichment_WGCNA_module_green.txt", sep = "\t", row.names = T)

#results_green <- read.table("Tables/Enrichment tables/ggprofiler_results_Enrichment_WGCNA_module_green_new.txt", sep = "\t", header = T)

results_green$`Term name` <- paste(results_green$term_name, "\n (N = ",results_green$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_green[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/Enrichment_module_Green.pdf", height = 6, width = 6)


#| Top with MIRNA 
dat1_filtered <- results_green[which(results_green$source == "MIRNA"),]
ggplot(dat1_filtered[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green MIRNA") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/MIRNA_Enrichment_module_Green.pdf", height = 5, width = 3.5)

#| Top with GO:BP 
dat1_filtered <- results_green[which(results_green$source == "GO:BP"),]
ggplot(dat1_filtered[1:14,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green GO:BP") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;BP_Enrichment_module_Green.pdf", height = 6, width = 4.6)


#| Top with GO:CC 
dat1_filtered <- results_green[which(results_green$source == "GO:CC"),]
ggplot(dat1_filtered[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green GO:CC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;CC_Enrichment_module_Green.pdf", height = 6, width = 4.6)


#| Top with GO:MF 
dat1_filtered <- results_green[which(results_green$source == "GO:MF"),]
ggplot(dat1_filtered[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green GO:MF") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;MF_Enrichment_module_Green.pdf", height = 6, width = 4.6)

#| Top with TF 
dat1_filtered <- results_green[which(results_green$source == "TF"),]
ggplot(dat1_filtered[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green TF") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/TF_Enrichment_module_Green.pdf", height = 6, width = 5)

#| Top with REAC
dat1_filtered <- results_green[which(results_green$source == "REAC"),]
ggplot(dat1_filtered[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Green REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/REAC_Enrichment_module_Green.pdf", height = 6, width = 5)
################################################################################


################################## Magenta #####################################
genes_magenta <- rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "magenta"]
total_gost_magenta <- gost(list("Module magenta WGCNA" = genes_magenta), 
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_magenta, interactive = FALSE, capped =FALSE)
ggsave("Images/Enrichment plots new/ggplot_results_WGCNA_Module_magenta.pdf", height = 5, width = 8)


results_magenta <- total_gost_magenta$result[order(total_gost_magenta$result$p_value),]

write.table(results_magenta[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment tables new/ggprofiler_results_Enrichment_WGCNA_module_magenta.txt", sep = "\t", row.names = T)

#results_magenta <- read.table("Tables/Enrichment tables/ggprofiler_results_Enrichment_WGCNA_module_turquoise.txt", sep ="\t",  header = T)

results_magenta$`Term name` <- paste(results_magenta$term_name, "\n (N = ",results_magenta$term_size, ")",sep ="")


#| Top 15 most significant
ggplot(results_magenta[1:11,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module magenta") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/Enrichment_module_magenta.pdf", height = 5, width = 5)

#| Top with TF
dat1_filtered <- results_magenta[which(results_magenta$source == "TF"),]
ggplot(dat1_filtered[1:11,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta TF") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/TF_Enrichment_module_magenta.pdf", height = 4.5, width = 4.8)

#| Top with MIRNA
dat1_filtered <- results_magenta[which(results_magenta$source == "MIRNA"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta MIRNA") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/MIRNA_Enrichment_module_magenta.pdf", height = 4.5, width = 3.3)

#| Top with HPA
dat1_filtered <- results_magenta[which(results_magenta$source == "HPA"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta HPA") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/HPA_Enrichment_module_magenta.pdf", height = 4.5, width = 4.5)


#| Top with REAC
dat1_filtered <- results_magenta[which(results_magenta$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/REAC_Enrichment_module_magenta.pdf", height = 4.5, width = 5)

#| Top with GO:CC
dat1_filtered <- results_magenta[which(results_magenta$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta GO:CC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;CC_Enrichment_module_magenta.pdf", height = 4.5, width = 4.4)
################################################################################



################################### Red #####################################
genes_red <- rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "red"]
total_gost_red <- gost(list("Module Red WGCNA" = genes_red), 
                          organism = "hsapiens", ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = TRUE,
                          user_threshold = 0.05, correction_method = "fdr",
                          domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)

gostplot(total_gost_red, interactive = FALSE, capped =FALSE)
ggsave("Images/Enrichment plots new/ggplot_results_WGCNA_Module_red.pdf", height = 5, width = 8)


results_red <- total_gost_red$result[order(total_gost_red$result$p_value),]

write.table(results_red[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment tables/ggprofiler_results_Enrichment_WGCNA_module_red_new.txt", sep = "\t", row.names = T)

results_red$`Term name` <- paste(results_red$term_name, "\n (N = ",results_red$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_red[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Red") +
  xlab("Data")
ggsave("Images/Enrichment plots new/Enrichment_module_red.pdf", height = 4.5, width = 6)

#| Top with REAC
dat1_filtered <- results_red[which(results_red$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/REAC_Enrichment_module_red.pdf", height = 4.3, width = 4.5)

#| Top with TF
dat1_filtered <- results_red[which(results_red$source == "TF"),]
ggplot(dat1_filtered[1:12,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/TF_Enrichment_module_red.pdf", height = 4.5, width = 5)

#| Top with GO:CC
dat1_filtered <- results_red[which(results_red$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;CC_Enrichment_module_red.pdf", height = 4.5, width = 4.6)

#| Top with GO:BP
dat1_filtered <- results_red[which(results_red$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/GO;BP_Enrichment_module_red.pdf", height = 4.5, width = 4.4)

#| Top with MIRNA
dat1_filtered <- results_red[which(results_red$source == "MIRNA"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta MIRNA") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/MIRNA_Enrichment_module_red.pdf", height = 4.5, width = 3.3)

################################################################################


##################################### Darkgrey ######################################
genes_darkgrey <- rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "darkred"]
total_gost_darkgrey <- gost(list("Module darkgrey WGCNA" = genes_darkgrey), 
                       organism = "hsapiens", ordered_query = FALSE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                       measure_underrepresentation = FALSE, evcodes = TRUE,
                       user_threshold = 0.05, correction_method = "fdr",
                       domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots new/ggplot_results_WGCNA_Module_darkgrey.pdf")
gostplot(total_gost_darkgrey, interactive = FALSE, capped =FALSE)
dev.off()

results_darkgrey <- total_gost_darkgrey$result[order(total_gost_darkgrey$result$p_value),]

write.table(results_greenyellow[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment tables/ggprofiler_results_Enrichment_WGCNA_module_greenyellow_new.txt", sep = "\t", row.names = T)


#results_greenyellow <- read.table("Tables/Enrichment tables/ggprofiler_results_Enrichment_WGCNA_module_greenyellow.txt", sep = "\t", header = T)

results_darkgrey$`Term name` <-  paste(results_darkgrey$term_name, "\n (N = ",results_darkgrey$term_size, ")",sep ="")


#| Top most significant
ggplot(results_darkgrey[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module darkgrey") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots new/Enrichment_module_darkgrey.pdf", height = 6, width = 6)
################################################################################



############################ PLOT: SIZE OF MODULES #############################
#| Size
size <- c()
for (i in 1:length(module_colors)){
  size <- c(size,length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors ==module_colors[i])]) )
}

length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == "grey")])

genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == "green")])]
#| New dataframe
data_size <- data.frame(modules = module_colors,
                                   size = size)
colors <- sort(unique(moduleColors))
ggplot(data_size, aes(x= size, y = reorder(modules,size), fill =modules)) +
  geom_bar(stat="identity", color ="black") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.position = "none",legend.key.size = unit(0.2, 'cm'),plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors) +
  xlab("Module size") +
  ylab("Module by color name")
ggsave("Images/5_Module_size_new.pdf", height = 5, width = 6.6)
################################################################################



####################### VOLCANO WITH MODULE COLORS #############################
#| Enrichment over all the modules
modules_color <- sort(module_colors)
dat_Modules_DEGs <- merge(wgcna_Modules, DEGs_genes, by ="row.names") 

dat_Modules_DEGs$delabel <- NA
dat_Modules_DEGs$delabel[which(dat_Modules_DEGs$log2FoldChange < (-log2(2.5)) | dat_Modules_DEGs$log2FoldChange > (log2(1.8)) & dat_Modules_DEGs$padj <= 0.05)] <- 
  dat_Modules_DEGs$gene_name[which(dat_Modules_DEGs$log2FoldChange < (-log2(2.5)) | dat_Modules_DEGs$log2FoldChange > (log2(1.8)) & dat_Modules_DEGs$padj <= 0.05)]

ggplot(dat_Modules_DEGs, aes(x = log2FoldChange, y = -log10(padj), color=moduleColors, label=delabel )) + 
  geom_point(size =4,alpha = 0.5) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_color_manual(values = modules_color ,guide = "legend") +
  labs(color ="Module color") +
  xlab("Log2 Fold Change") +
  geom_text_repel()
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_PTEN_loss_vs_intact/Results/design ~ DV200 + Edad + H-score/VolcanoPlots/Volcano_DEGs_Module_colors_new_2.pdf", height = 5.5, width = 7)
################################################################################



################# % OF DEGS CONTAINED IN THE DIFFERENT MODULES #################

#| FC == 1.5
dges_modules <- wgcna_Modules[which(rownames(wgcna_Modules) %in% DEGs$GeneID),]
colors <- unique(wgcna_Modules$moduleColors[which(rownames(wgcna_Modules) %in% DEGs$GeneID)])

data_percentage_DGEs <- data.frame(modules = colors,
                                   fraction = c(length(rownames(dges_modules[which(dges_modules$moduleColors == colors[1]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[2]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[3]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[4]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[5]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[6]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[7]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[8]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[9]),]))/length(DEGs$GeneID),
                                                length(rownames(dges_modules[which(dges_modules$moduleColors == colors[10]),]))/length(DEGs$GeneID)),
                                   size = c(length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[1])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[2])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[3])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[4])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[5])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[6])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[7])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[8])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[9])]),
                                            length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors[10])])))
#| Plotting the results
ggplot(data_percentage_DGEs, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors[order(colors)]) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Fraction of DEGs (PTEN loss vs intact)\n in the WGCNA modules  (PCa samples)")+
  labs(fill = "Modules") +
  labs(size = "Module size") +
  # ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing PTEN loss vs intact") +
  ylim(0,1)
ggsave("Images/5_Fraction_DEGs_in_WGCNA_modules_PTEN_loss_vs_intact_bubble_plot_new.pdf", height=4, width =5)

length(rownames(dges_modules[which(dges_modules$moduleColors == "yellow"),]))
################################################################################




######################### NON CODINGS REGULATING PTEN ##########################

#| Finding the non-coding genes
NonCodings <- nonCoding_genes$`Approved symbol`[which(nonCoding_genes$`Locus type` =="RNA, long non-coding" | nonCoding_genes$`Locus type` =="RNA, micro")]

#| For every module
black <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "black"])]
green <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "green"])]
red <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "red"])]
turquoise <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "turquoise"])]
grey <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "grey"])]
brown <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "brown"])]
blue <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "blue"])]
pink <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "pink"])]
greenyellow <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "greenyellow"])]
magenta <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "magenta"])]
purple <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "purple"])]
yellow <- genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% rownames(wgcna_Modules)[wgcna_Modules$moduleColors == "yellow"])]

#| DEGs
deg <- DEGs$gene_name[which((DEGs$log2FoldChange > log2(1.5) | DEGs$log2FoldChange < (-log2(1.5))) & DEGs$padj < 0.05)]

black[which(black=="MIR383")]
green[which(green=="MALAT1")]
red[which(red=="NBAT1")]
black[grep("WNT", black)]
red[grep("WNT", red)]

#| non-coding that regulate PTEN
list_nonCoding_RegulatingPTEN <- c("MIR21", "MIR130", "MIR451", "MIR221", "MIR222", "MIR301A", "MIR214", "MIR494", "MIR29A", "MIR29A",
                                   "MIR1555P", "MIR130B", "MIR185", "MIR29", "MIR101",  "MIR29B", "MIR17", "MIR19", "MIR181A",
                                   "MIR136", "MIR106A", "MIR103A", "MIR217", "MIR93", "MIR9", "MIR106B", "MIR510", "MIR1443P",
                                   "MIR19A", "MIR19A3P", "MIR19B", "MIR18A", "MIR26A1", "MIR26B","MIR302A", "MIR506", "MIR4310",
                                   "MIR216A", "MIR193A", "MIR33A", "MIR144", "MIR23A", "MIR20A","MIR22",  "MIR202", "MIR499A", 
                                   "MIR1908", "MIR544", "MIR26A", "MIR206", "MIR486", "MIR718", "MIR589" , "MIR200A" , "MIR552" ,
                                   "MIR186", "MIR371","MIR548","MIR130A","MIR298", "MIR561", "MIR4262", "MIR543", "MIR155",
                                   "MIR205", "MIR32", "MIR126", "MIR10A", "MIR152","MIR4310","MIR22HG",
                                   "GAS5", "Linc-USP16", "BGL3","lincRNA-p21", "PTENP1", "MIR148A" , "MALAT1",
                                   "CASC2", "MEG3","XIST", "NBAT1", "FER1L4","NEAT1", "FER1L4",  "RP11-79H23.3", "LINC00470",
                                   "SLC25A5-AS1", "Gas5", "TP73-AS1", "CNOT6L","ZEB2", "DKK", "VCAN", "VAPA", "LINC00702", "LINC00689")

non_coding_ENS <- c("ENSG00000284536","ENSG00000284485","ENSG00000284416","ENSG00000284286","ENSG00000284219","ENSG00000284204","ENSG00000284190","ENSG00000284157","ENSG00000284038","ENSG00000284032","ENSG00000283904","ENSG00000283871","ENSG00000283844","ENSG00000283824","ENSG00000283819","ENSG00000283815","ENSG00000283762","ENSG00000265172","ENSG00000264850","ENSG00000264850","ENSG00000260455","ENSG00000251562","ENSG00000245532","ENSG00000237984","ENSG00000234741","ENSG00000234741","ENSG00000233117","ENSG00000229807","ENSG00000227372","ENSG00000224281","ENSG00000216031","ENSG00000214548","ENSG00000212040","ENSG00000208036","ENSG00000208023","ENSG00000208009","ENSG00000207996","ENSG00000207980","ENSG00000207973","ENSG00000207951","ENSG00000207947","ENSG00000207942","ENSG00000207941","ENSG00000207932","ENSG00000207927","ENSG00000207870","ENSG00000207798","ENSG00000207757","ENSG00000207731","ENSG00000207725","ENSG00000207721","ENSG00000207698","ENSG00000207641","ENSG00000207635","ENSG00000207614","ENSG00000207607","ENSG00000207604","ENSG00000207548","ENSG00000199161","ENSG00000199121","ENSG00000199085","ENSG00000199075","ENSG00000194717","ENSG00000186594","ENSG00000177640","ENSG00000169554","ENSG00000138767","ENSG00000132204","ENSG00000101558","ENSG00000088340","ENSG00000088340","ENSG00000038427")
list_nonCoding_RegulatingPTEN <- non_coding_ENS
#| New dataframe
data_percentage_DGEs <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% list_nonCoding_RegulatingPTEN)]),
                                                length(grey[which(grey %in% list_nonCoding_RegulatingPTEN)]),
                                                length(brown[which(brown %in% list_nonCoding_RegulatingPTEN)]),
                                                length(red[which(red %in% list_nonCoding_RegulatingPTEN)]),
                                                length(purple[which(purple %in% list_nonCoding_RegulatingPTEN)]),
                                                length(yellow[which(yellow %in% list_nonCoding_RegulatingPTEN)]),
                                                length(blue[which(blue %in% list_nonCoding_RegulatingPTEN)]),
                                                length(magenta[which(magenta %in% list_nonCoding_RegulatingPTEN)]),
                                                length(green[which(green %in% list_nonCoding_RegulatingPTEN)]),
                                                length(black[which(black %in% list_nonCoding_RegulatingPTEN)]),
                                                length(pink[which( pink %in% list_nonCoding_RegulatingPTEN)]),
                                                length(greenyellow[which(greenyellow %in% list_nonCoding_RegulatingPTEN)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_DGEs, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of Non-coding genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_NONCODING_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many non-codings in the DEGs?
deg[which(deg %in% list_nonCoding_RegulatingPTEN)]

ACo_genes <- deg[grep("AC0",deg)]
DEGs_genes$GeneID[which(DEGs_genes$gene_name %in% ACo_genes)]

#| How many EMT genes are in each module?
EMT_genes <- read.table("../EMT-genes.txt", sep =",")
EMT_genes <- as.character(EMT_genes)

#| How many EMT in the DEGs?
deg[which(deg %in% EMT_genes)]

#| New dataframe
data_percentage_EMT <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% EMT_genes)]),
                                                length(grey[which(grey %in% EMT_genes)]),
                                                length(brown[which(brown %in% EMT_genes)]),
                                                length(red[which(red %in% EMT_genes)]),
                                                length(purple[which(purple %in% EMT_genes)]),
                                                length(yellow[which(yellow %in% EMT_genes)]),
                                                length(blue[which(blue %in% EMT_genes)]),
                                                length(magenta[which(magenta %in% EMT_genes)]),
                                                length(green[which(green %in% EMT_genes)]),
                                                length(black[which(black %in% EMT_genes)]),
                                                length(pink[which( pink %in% EMT_genes)]),
                                                length(greenyellow[which(greenyellow %in% EMT_genes)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_EMT, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of EMT genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_EMT_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)



#| How many NOTCH genes are in each module?
NOTCH_genes <- read.table("../NOTCH-genes.txt", sep =",")
NOTCH_genes <- as.character(NOTCH_genes)

#| How many NOTCH in the DEGs?
deg[which(deg %in% NOTCH_genes)]

#| New dataframe
data_percentage_NOTCH <- data.frame(modules = module_colors,
                                  fraction = c(length(turquoise[which(turquoise %in% NOTCH_genes)]),
                                               length(grey[which(grey %in% NOTCH_genes)]),
                                               length(brown[which(brown %in% NOTCH_genes)]),
                                               length(red[which(red %in% NOTCH_genes)]),
                                               length(purple[which(purple %in% NOTCH_genes)]),
                                               length(yellow[which(yellow %in% NOTCH_genes)]),
                                               length(blue[which(blue %in% NOTCH_genes)]),
                                               length(magenta[which(magenta %in% NOTCH_genes)]),
                                               length(green[which(green %in% NOTCH_genes)]),
                                               length(black[which(black %in% NOTCH_genes)]),
                                               length(pink[which( pink %in% NOTCH_genes)]),
                                               length(greenyellow[which(greenyellow %in% NOTCH_genes)])),
                                  size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_NOTCH, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of NOTCH genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_NOTCH_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)



#| How many TNFA genes are in each module?
TNFA_genes <- read.table("../TNFA-genes.txt", sep =",")
TNFA_genes <- as.character(TNFA_genes)

#| How many TNFA in the DEGs?
deg[which(deg %in% TNFA_genes)]

#| New dataframe
data_percentage_TNFA <- data.frame(modules = module_colors,
                                    fraction = c(length(turquoise[which(turquoise %in% TNFA_genes)]),
                                                 length(grey[which(grey %in% TNFA_genes)]),
                                                 length(brown[which(brown %in% TNFA_genes)]),
                                                 length(red[which(red %in% TNFA_genes)]),
                                                 length(purple[which(purple %in% TNFA_genes)]),
                                                 length(yellow[which(yellow %in% TNFA_genes)]),
                                                 length(blue[which(blue %in% TNFA_genes)]),
                                                 length(magenta[which(magenta %in% TNFA_genes)]),
                                                 length(green[which(green %in% TNFA_genes)]),
                                                 length(black[which(black %in% TNFA_genes)]),
                                                 length(pink[which( pink %in% TNFA_genes)]),
                                                 length(greenyellow[which(greenyellow %in% TNFA_genes)])),
                                    size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_TNFA, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of TNFA genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_TNFA_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)



#| How many PI3K genes are in each module?
PI3K_genes <- read.table("../PI3K-genes.txt", sep =",")
PI3K_genes <- as.character(PI3K_genes)

#| How many TNFA in the DEGs?
deg[which(deg %in% PI3K_genes)]

#| New dataframe
data_percentage_PI3K <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% PI3K_genes)]),
                                                length(grey[which(grey %in% PI3K_genes)]),
                                                length(brown[which(brown %in% PI3K_genes)]),
                                                length(red[which(red %in% PI3K_genes)]),
                                                length(purple[which(purple %in% PI3K_genes)]),
                                                length(yellow[which(yellow %in% PI3K_genes)]),
                                                length(blue[which(blue %in% PI3K_genes)]),
                                                length(magenta[which(magenta %in% PI3K_genes)]),
                                                length(green[which(green %in% PI3K_genes)]),
                                                length(black[which(black %in% PI3K_genes)]),
                                                length(pink[which( pink %in% PI3K_genes)]),
                                                length(greenyellow[which(greenyellow %in% PI3K_genes)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_PI3K, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of PI3K genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_PI3K_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many TGFB genes are in each module?
TGFB_genes <- read.table("../TGFB-genes.txt", sep =",")
TGFB_genes <- as.character(TGFB_genes)

#| How many TGFB in the DEGs?
deg[which(deg %in% TGFB_genes)]

#| New dataframe
data_percentage_TGFB <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% TGFB_genes)]),
                                                length(grey[which(grey %in% TGFB_genes)]),
                                                length(brown[which(brown %in% TGFB_genes)]),
                                                length(red[which(red %in% TGFB_genes)]),
                                                length(purple[which(purple %in% TGFB_genes)]),
                                                length(yellow[which(yellow %in% TGFB_genes)]),
                                                length(blue[which(blue %in% TGFB_genes)]),
                                                length(magenta[which(magenta %in% TGFB_genes)]),
                                                length(green[which(green %in% TGFB_genes)]),
                                                length(black[which(black %in% TGFB_genes)]),
                                                length(pink[which( pink %in% TGFB_genes)]),
                                                length(greenyellow[which(greenyellow %in% TGFB_genes)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_TGFB, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of NOTCH genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_TGFB_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many TGFB genes are in each module?
WNT_genes <- read.table("../WNT-genes.txt", sep =",")
WNT_genes <- as.character(WNT_genes)

#| New dataframe
data_percentage_WNT <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% WNT_genes)]),
                                                length(grey[which(grey %in% WNT_genes)]),
                                                length(brown[which(brown %in% WNT_genes)]),
                                                length(red[which(red %in% WNT_genes)]),
                                                length(purple[which(purple %in% WNT_genes)]),
                                                length(yellow[which(yellow %in% WNT_genes)]),
                                                length(blue[which(blue %in% WNT_genes)]),
                                                length(magenta[which(magenta %in% WNT_genes)]),
                                                length(green[which(green %in% WNT_genes)]),
                                                length(black[which(black %in% WNT_genes)]),
                                                length(pink[which( pink %in% WNT_genes)]),
                                                length(greenyellow[which(greenyellow %in% WNT_genes)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_WNT, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of NOTCH genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_WNT_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)



#| How many TGFB genes are in each module?
HEDGEHOG_genes <- read.table("../HEDGEHOG-genes.txt", sep =",")
HEDGEHOG_genes <- as.character(HEDGEHOG_genes)

#| How many TGFB in the DEGs?
deg[which(deg %in% HEDGEHOG_genes)]

#| New dataframe
data_percentage_HEDGEHOG <- data.frame(modules = module_colors,
                                  fraction = c(length(turquoise[which(turquoise %in% HEDGEHOG_genes)]),
                                               length(grey[which(grey %in% HEDGEHOG_genes)]),
                                               length(brown[which(brown %in% HEDGEHOG_genes)]),
                                               length(red[which(red %in% HEDGEHOG_genes)]),
                                               length(purple[which(purple %in% HEDGEHOG_genes)]),
                                               length(yellow[which(yellow %in% HEDGEHOG_genes)]),
                                               length(blue[which(blue %in% HEDGEHOG_genes)]),
                                               length(magenta[which(magenta %in% HEDGEHOG_genes)]),
                                               length(green[which(green %in% HEDGEHOG_genes)]),
                                               length(black[which(black %in% HEDGEHOG_genes)]),
                                               length(pink[which( pink %in% HEDGEHOG_genes)]),
                                               length(greenyellow[which(greenyellow %in% HEDGEHOG_genes)])),
                                  size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                           length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_HEDGEHOG, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of NOTCH genes that regulate PTEN \n in each module")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_HEDGEHOG_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many 200B_200C genes are in each module?
MIR200B_200C_genes <- read.table("../MIR200B_200C_genes.txt", sep =",")
MIR200B_200C_genes <- as.character(MIR200B_200C_genes)

#| How many 200B_200C in the DEGs?
deg[which(deg %in% MIR200B_200C_genes)]

#| New dataframe
data_percentage_MIR200B_200C <- data.frame(modules = module_colors,
                                       fraction = c(length(turquoise[which(turquoise %in% MIR200B_200C_genes)]),
                                                    length(grey[which(grey %in% MIR200B_200C_genes)]),
                                                    length(brown[which(brown %in% MIR200B_200C_genes)]),
                                                    length(red[which(red %in% MIR200B_200C_genes)]),
                                                    length(purple[which(purple %in% MIR200B_200C_genes)]),
                                                    length(yellow[which(yellow %in% MIR200B_200C_genes)]),
                                                    length(blue[which(blue %in% MIR200B_200C_genes)]),
                                                    length(magenta[which(magenta %in% MIR200B_200C_genes)]),
                                                    length(green[which(green %in% MIR200B_200C_genes)]),
                                                    length(black[which(black %in% MIR200B_200C_genes)]),
                                                    length(pink[which( pink %in% MIR200B_200C_genes)]),
                                                    length(greenyellow[which(greenyellow %in% MIR200B_200C_genes)])),
                                       size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                                length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_MIR200B_200C, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of genes in the WGCNA modules beloging ro the MIR200C/200B pathway")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_MIR200B_200C_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many 200B_200C genes are in each module?
MIR9_genes <- read.table("../MIR9_genes.txt", sep =",")
MIR9_genes <- as.character(MIR9_genes)

#| How many 200B_200C in the DEGs?
deg[which(deg %in% MIR9_genes)]

#| New dataframe
data_percentage_MIR9 <- data.frame(modules = module_colors,
                                           fraction = c(length(turquoise[which(turquoise %in% MIR9_genes)]),
                                                        length(grey[which(grey %in% MIR9_genes)]),
                                                        length(brown[which(brown %in% MIR9_genes)]),
                                                        length(red[which(red %in% MIR9_genes)]),
                                                        length(purple[which(purple %in% MIR9_genes)]),
                                                        length(yellow[which(yellow %in% MIR9_genes)]),
                                                        length(blue[which(blue %in% MIR9_genes)]),
                                                        length(magenta[which(magenta %in% MIR9_genes)]),
                                                        length(green[which(green %in% MIR9_genes)]),
                                                        length(black[which(black %in% MIR9_genes)]),
                                                        length(pink[which( pink %in% MIR9_genes)]),
                                                        length(greenyellow[which(greenyellow %in% MIR9_genes)])),
                                           size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                                    length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_MIR9, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of genes in the WGCNA modules beloging ro the MIR9 pathway")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_MIR9_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)




#| How many 200B_200C genes are in each module?
MIR21_genes <- read.table("../MIR21_genes.txt", sep =",")
MIR21_genes <- as.character(MIR21_genes)

#| How many 200B_200C in the DEGs?
deg[which(deg %in% MIR21_genes)]

#| New dataframe
data_percentage_MIR21 <- data.frame(modules = module_colors,
                                   fraction = c(length(turquoise[which(turquoise %in% MIR21_genes)]),
                                                length(grey[which(grey %in% MIR21_genes)]),
                                                length(brown[which(brown %in% MIR21_genes)]),
                                                length(red[which(red %in% MIR21_genes)]),
                                                length(purple[which(purple %in% MIR21_genes)]),
                                                length(yellow[which(yellow %in% MIR21_genes)]),
                                                length(blue[which(blue %in% MIR21_genes)]),
                                                length(magenta[which(magenta %in% MIR21_genes)]),
                                                length(green[which(green %in% MIR21_genes)]),
                                                length(black[which(black %in% MIR21_genes)]),
                                                length(pink[which( pink %in% MIR21_genes)]),
                                                length(greenyellow[which(greenyellow %in% MIR21_genes)])),
                                   size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                            length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_MIR21, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of genes in the WGCNA modules beloging ro the MIR21 pathway")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_MIR21_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)



#| How many 200B_200C genes are in each module?
MIR202_genes <- read.table("../MIR202_genes.txt", sep =",")
MIR202_genes <- as.character(MIR202_genes)

#| How many 200B_200C in the DEGs?
deg[which(deg %in% MIR202_genes)]

#| New dataframe
data_percentage_MIR202 <- data.frame(modules = module_colors,
                                    fraction = c(length(turquoise[which(turquoise %in% MIR202_genes)]),
                                                 length(grey[which(grey %in% MIR202_genes)]),
                                                 length(brown[which(brown %in% MIR202_genes)]),
                                                 length(red[which(red %in% MIR202_genes)]),
                                                 length(purple[which(purple %in% MIR202_genes)]),
                                                 length(yellow[which(yellow %in% MIR202_genes)]),
                                                 length(blue[which(blue %in% MIR202_genes)]),
                                                 length(magenta[which(magenta %in% MIR202_genes)]),
                                                 length(green[which(green %in% MIR202_genes)]),
                                                 length(black[which(black %in% MIR202_genes)]),
                                                 length(pink[which( pink %in% MIR202_genes)]),
                                                 length(greenyellow[which(greenyellow %in% MIR202_genes)])),
                                    size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                             length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_MIR202, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of genes in the WGCNA modules beloging ro the MIR202 pathway")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_MIR202_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)


#| How many 200B_200C genes are in each module?
MIR498_genes <- read.table("../MIR498_genes.txt", sep =",")
MIR498_genes <- as.character(MIR498_genes)

#| How many 200B_200C in the DEGs?
deg[which(deg %in% MIR498_genes)]

#| New dataframe
data_percentage_MIR498 <- data.frame(modules = module_colors,
                                     fraction = c(length(turquoise[which(turquoise %in% MIR498_genes)]),
                                                  length(grey[which(grey %in% MIR498_genes)]),
                                                  length(brown[which(brown %in% MIR498_genes)]),
                                                  length(red[which(red %in% MIR498_genes)]),
                                                  length(purple[which(purple %in% MIR498_genes)]),
                                                  length(yellow[which(yellow %in% MIR498_genes)]),
                                                  length(blue[which(blue %in% MIR498_genes)]),
                                                  length(magenta[which(magenta %in% MIR498_genes)]),
                                                  length(green[which(green %in% MIR498_genes)]),
                                                  length(black[which(black %in% MIR498_genes)]),
                                                  length(pink[which( pink %in% MIR498_genes)]),
                                                  length(greenyellow[which(greenyellow %in% MIR498_genes)])),
                                     size = c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="turquoise")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="grey")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="red")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColor =="purple")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="brown")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="yellow")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="magenta")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="green")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="black")]),
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="pink")]), 
                                              length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors =="greenyellow")])))

#| Plotting and saving
ggplot(data_percentage_MIR498, aes(x =modules, y = fraction, fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=12,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("black" ,"blue","brown", "green","greenyellow","grey","magenta","pink"  ,"purple", "red","turquoise","yellow")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Number of genes in the WGCNA modules beloging ro the MIR498 pathway")+
  labs(fill = "Modules") +
  labs(size = "Module size") #+
ggsave("Images/5_MIR498_genes_in_MODULES_bubble_plot.pdf", height=6, width =7)

################################################################################




