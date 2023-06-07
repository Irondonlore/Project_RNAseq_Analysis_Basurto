################################################################################
#####             ENRICHMENT ANALYSIS OF MODULES DE FROM WGCNA             #####                        
################################################################################

#| 

################################################################################



################################################################################
###############################  LIBRARIES  ####################################
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(viridis))

options(stringsAsFactors = FALSE)
################################################################################

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)


workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/"
setwd(workingDir)

data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + Diagnostico/Tables/resDESeq2_dds_DV200_Edad_Diagnostico_ContrastPCa_vs_BPH.txt"
data.file_DEGs2 <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + DFS-STATUS/Tables/resDESeq2_dds_DV200_Edad_DFS-LA;LI_ContrastLA_vs_LI.txt"

data.file_DEGs_filtered <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + Diagnostico/Tables/resDESeq2_dds_DV200_Edad_Diagnostico_ContrastPCa_vs_BPH_filtered.txt"
data.file_DEGs2_filtered <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + DFS-STATUS/Tables/resDESeq2_dds_DV200_Edad_DFS-LA;LI_ContrastLA_vs_LI_filtered.txt"

#| DESeq2 DATA:
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")
DEGs_genes2 <- read.table(data.file_DEGs2, sep = "\t")

DEGs_genes_filtered <- read.table(data.file_DEGs_filtered, sep = "\t")
DEGs_genes2_filtered <- read.table(data.file_DEGs2_filtered, sep = "\t")

#| Filtering those genes who satisfied the condition for DEGs classified by "Yes"
DEGs_genes_up <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange >=0),]
DEGs_genes_down <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange <0),]

DEGs_genes_up2 <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes" & DEGs_genes2$log2FoldChange >=0),]
DEGs_genes_down2 <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes" & DEGs_genes2$log2FoldChange <0),]

DEGs_genes_up_filtered <- DEGs_genes_filtered[which(DEGs_genes_filtered$DEGs == "Yes" & DEGs_genes_filtered$log2FoldChange >=0),]
DEGs_genes_down_filtered <- DEGs_genes_filtered[which(DEGs_genes_filtered$DEGs == "Yes" & DEGs_genes_filtered$log2FoldChange <0),]

DEGs_genes_up2_filtered <- DEGs_genes2[which(DEGs_genes2_filtered$DEGs == "Yes" & DEGs_genes2_filtered$log2FoldChange >=0),]
DEGs_genes_down2_filtered <- DEGs_genes2[which(DEGs_genes2_filtered$DEGs == "Yes" & DEGs_genes2_filtered$log2FoldChange <0),]

DEGs_all <- DEGs_genes[which(DEGs_genes$DEGs == "Yes"),]
DEGs2_all <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes"),]

DEGs_all_filtered <- DEGs_genes_filtered[which(DEGs_genes_filtered$DEGs == "Yes"),]
DEGs2_all_filtered <- DEGs_genes2_filtered[which(DEGs_genes2_filtered$DEGs == "Yes"),]

################################################################################


################################################################################
geneInfo <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors.txt",  sep = "\t")
modules <- unique(geneInfo$moduleColors)


#| Green
genes_green <- rownames(geneInfo)[geneInfo$moduleColors == "green"]
total_gost_green <- gost(list("Module Green WGCNA" = genes_green), 
                              organism = "hsapiens", ordered_query = FALSE, 
                              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                              measure_underrepresentation = FALSE, evcodes = TRUE,
                              user_threshold = 0.05, correction_method = "fdr",
                              domain_scope = "custom", custom_bg = rownames(geneInfo), 
                              numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_green.pdf")
gostplot(total_gost_green, interactive = FALSE)
dev.off()

results_green <- total_gost_green$result[order(total_gost_green$result$p_value),]
results_green$term_name[1:50]

gp_mod = results_green[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]

gp_mod$GeneRatio = paste0(gp_mod$intersection_size, "/", gp_mod$query_size)

gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)

names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size","geneID", "GeneRatio", "BgRatio")

gp_mod$geneID = gsub(",", "/", gp_mod$geneID)

row.names(gp_mod) = gp_mod$ID

# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
# define as enrichResult object
gp_mod_enrich = new("enrichResult", result = gp_mod)

enrichplot::dotplot(gp_mod_cluster)

barplot(gp_mod_enrich, showCategory = 40, font.size = 16) +
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")

write.table(results_green[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_green.txt", sep = "\t")


#| Turquoise
genes_turquoise <- rownames(geneInfo)[geneInfo$moduleColors == "turquoise"]
total_gost_turquoise <- gost(list("Module Turquoise WGCNA" = genes_turquoise), 
                                  organism = "hsapiens", ordered_query = FALSE,
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                                  measure_underrepresentation = FALSE, evcodes = TRUE,
                                  user_threshold = 0.05, correction_method = "fdr",
                                  domain_scope = "custom", custom_bg = rownames(geneInfo), 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_turquoise.pdf")
gostplot(total_gost_turquoise, interactive = FALSE, capped = FALSE)
dev.off()

results_turquoise <- total_gost_turquoise$result[order(total_gost_turquoise$result$p_value),]
results_turquoise$term_name[1:20]

write.table(results_turquoise[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_turquoise.txt", sep = "\t", row.names = T)
d <- results_turquoise[,c("term_id", "term_name")]
unique(results_turquoise$source)

#| Blue
genes_blue <- rownames(geneInfo)[geneInfo$moduleColors == "blue"]
total_gost_blue <- gost(list("Module Blue WGCNA" = genes_blue),
                             organism = "hsapiens", ordered_query = FALSE,
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = TRUE,
                             user_threshold = 0.05, correction_method = "fdr",
                             domain_scope = "custom", custom_bg = rownames(geneInfo),
                             numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_blue.pdf")
gostplot(total_gost_blue, interactive = FALSE)
dev.off()
results_blue <- total_gost_blue$result[order(total_gost_blue$result$p_value),]
results_blue$term_name[1:50]
write.table(results_blue[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_blue.txt", sep = "\t", row.names = T)

results_blue <- read.table("Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_blue.txt", sep = "\t", header = T)

results_blue$`Term name` <-  paste(results_blue$term_name, "\n (N = ",results_blue$term_size, ")",sep ="")

#| Top most significant
ggplot(results_blue[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Blue") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/Enrichment_module_blue.pdf", height = 6, width = 6)


#| Top with GO:BP 
dat1_filtered <- results_blue[which(results_blue$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Blue GO:BP") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/GO;BP_Enrichment_module_blue.pdf", height = 6, width = 4)


#| Top with MIRNA 
dat1_filtered <- results_blue[which(results_blue$source == "MIRNA"),]
ggplot(dat1_filtered[1:8,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Blue MIRNA") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/MIRNA_Enrichment_module_blue.pdf", height = 6, width = 4)


#| Top with GO:CC
dat1_filtered <- results_blue[which(results_blue$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Blue GO:CC") +
  labs(size= "Intersection size") +
  xlab("Data")
ggsave("Images/Enrichment plots/GO;CC_Enrichment_module_Blue.pdf", height = 6, width = 4)



#| Magenta
genes_magenta <- rownames(geneInfo)[geneInfo$moduleColors == "magenta"]
total_gost_magenta <- gost(list("Module Magenta WGCNA" = genes_magenta), 
                             organism = "hsapiens", ordered_query = FALSE, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = TRUE,
                             user_threshold = 0.05, correction_method = "fdr",
                             domain_scope = "custom", custom_bg = rownames(geneInfo), 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_Magenta.pdf")
gostplot(total_gost_magenta, interactive = FALSE)
dev.off()
results_magenta <- total_gost_magenta$result[order(total_gost_magenta$result$p_value),]
results_magenta$term_name[1:50]
write.table(results_magenta[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_magenta.txt", sep = "\t", row.names = T)
results_magenta <- read.table("Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_magenta.txt", sep = "\t", header = T)

results_magenta$`Term name` <-  paste(results_magenta$term_name, "\n (N = ",results_magenta$term_size, ")",sep ="")

#| Top most significant
ggplot(results_magenta[1:15,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Magenta") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/Enrichment_module_magenta.pdf", height = 6, width = 6)

#| Top with GO:BP 
dat1_filtered <- results_magenta[which(results_magenta$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Magenta GO:BP") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/GO;BP_Enrichment_module_magenta.pdf", height = 6, width = 4)


#| Top with GO:MF 
dat1_filtered <- results_magenta[which(results_magenta$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Magenta GO:MF") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/GO;MF_Enrichment_module_magenta.pdf", height = 6, width = 4)

#| Top with GO:CC
dat1_filtered <- results_magenta[which(results_magenta$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  ggtitle("Enrichment Module Magenta GO:CC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("Images/Enrichment plots/GO;CC_Enrichment_module_magenta.pdf", height = 6, width = 4)





#| greenyellow
genes_greenyellow <- rownames(geneInfo)[geneInfo$moduleColors == "greenyellow"]
total_gost_greenyellow <- gost(list("Module greenyellow WGCNA" = genes_greenyellow), 
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = rownames(geneInfo), 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_greenyellow.pdf")
gostplot(total_gost_greenyellow, interactive = FALSE, capped=FALSE)
dev.off()
results_greenyellow <- total_gost_greenyellow$result[order(total_gost_greenyellow$result$p_value),]
results_greenyellow$term_name[1:50]
write.table(results_greenyellow[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_greenyellow.txt", sep = "\t", row.names = T)

#| Green
genes_green<- rownames(geneInfo)[geneInfo$moduleColors == "green"]
total_gost_green <- gost(list("Module Green WGCNA" = genes_green), 
                               organism = "hsapiens", ordered_query = FALSE, 
                               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                               measure_underrepresentation = FALSE, evcodes = TRUE,
                               user_threshold = 0.05, correction_method = "fdr",
                               domain_scope = "custom", custom_bg = rownames(geneInfo), 
                               numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_green.pdf")
gostplot(total_gost_green, interactive = FALSE)
dev.off()
results_green <- total_gost_green$result[order(total_gost_green$result$p_value),]
results_green$term_name[1:50]
write.table(results_green[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_green.txt", sep = "\t", row.names = T)

length(genes_blue)


#| Pink
genes_pink <- rownames(geneInfo)[geneInfo$moduleColors == "pink"]
total_gost_pink<- gost(list("Module Pink WGCNA" = genes_pink), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(geneInfo), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_pink.pdf")
gostplot(total_gost_pink, interactive = FALSE)
dev.off()
results_pink <- total_gost_blue$result[order(total_gost_pink$result$p_value),]
results_pink$term_name[1:50]
write.table(results_blue[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_pink.txt", sep = "\t", row.names = T)



#| Salmon
genes_salmon <- rownames(geneInfo)[geneInfo$moduleColors == "salmon"]
total_gost_salmon <- gost(list("Module Salmon WGCNA" = genes_salmon), 
                       organism = "hsapiens", ordered_query = FALSE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                       measure_underrepresentation = FALSE, evcodes = TRUE,
                       user_threshold = 0.05, correction_method = "fdr",
                       domain_scope = "custom", custom_bg = rownames(geneInfo), 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_salmon.pdf")
gostplot(total_gost_salmon, interactive = FALSE)
dev.off()
results_salmon <- total_gost_salmon$result[order(total_gost_salmon$result$p_value),]
results_salmon$term_name[1:100]
write.table(results_blue[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_pink.txt", sep = "\t", row.names = T)



#| cyan
genes_cyan <- rownames(geneInfo)[geneInfo$moduleColors == "cyan"]
total_gost_cyan <- gost(list("Module Cyan WGCNA" = genes_cyan), 
                          organism = "hsapiens", ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = TRUE,
                          user_threshold = 0.05, correction_method = "fdr",
                          domain_scope = "custom", custom_bg = rownames(geneInfo), 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_cyan.pdf")
gostplot(total_gost_cyan, interactive = FALSE, capped=FALSE)
dev.off()
results_cyan <- total_gost_cyan$result[order(total_gost_cyan$result$p_value),]
results_cyan$term_name[1:100]
write.table(results_blue[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_cyan.txt", sep = "\t", row.names = T)


#| tan
genes_tan <- rownames(geneInfo)[geneInfo$moduleColors == "tan"]
total_gost_tan <- gost(list("Module Tan WGCNA" = genes_tan), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(geneInfo), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
pdf("Images/Enrichment plots/ggplot_results_WGCNA_Module_cyan.pdf")
gostplot(total_gost_cyan, interactive = FALSE, capped=FALSE)
dev.off()
results_cyan <- total_gost_cyan$result[order(total_gost_cyan$result$p_value),]
results_cyan$term_name[1:100]
write.table(results_blue[,c("p_value", "term_name")], "Tables/Enrichment table/Enrichment_ggprofiler_results_WGCNA_Module_cyan.txt", sep = "\t", row.names = T)





################### PCa vs BPH inteception WGCNA and DEGS ######################

#| Turquoise
genes_down <- DEGs_genes_down$GeneID[which(DEGs_genes_down$GeneID %in% rownames(geneInfo)[geneInfo$moduleColors == "turquoise"])]
total_gost <- gost(list("Module  WGCNA" = genes_down), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(geneInfo), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped=FALSE)
total_gost$result[order(total_gost$result$p_value), ]$term_name[1:50]

#| Green
genes <- DEGs_genes$GeneID[which(DEGs_genes$GeneID %in% rownames(geneInfo)[geneInfo$moduleColors == "green"])]
total_gost <- gost(list("Module  WGCNA" = genes), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(geneInfo), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped=FALSE)
total_gost$result[order(total_gost$result$p_value), ]$term_name[1:50]
total_gost$result$p_value    
################################################################################




################### PCa vs BPH inteception WGCNA and DEGS ######################

#| Turquoise
genes_down <- DEGs_genes_down2$GeneID[which(DEGs_genes_down2$GeneID %in% rownames(geneInfo)[geneInfo$moduleColors == "turquoise"])]
total_gost <- gost(list("Module  WGCNA" = genes_down), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(geneInfo), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped=FALSE)
total_gost$result[order(total_gost$result$p_value), ]$term_name[1:50]

#| Magenta
genes <- DEGs2$GeneID[which(DEGs2$GeneID %in% rownames(geneInfo)[geneInfo$moduleColors == "magenta"])]
total_gost <- gost(list("Module  WGCNA" = genes), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(geneInfo), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped=FALSE)
total_gost$result[order(total_gost$result$p_value), ]$term_name[1:50]
total_gost$result$p_value    


#| Green
genes <- DEGs2$GeneID[which(DEGs2$GeneID %in% rownames(geneInfo)[geneInfo$moduleColors == "green"])]
total_gost <- gost(list("Module  WGCNA" = genes), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(geneInfo), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped=FALSE)
total_gost$result[order(total_gost$result$p_value), ]$term_name[1:50]
total_gost$result$p_value    
################################################################################






################################################################################
#| Which module are located the DGEs from DESeq2 comparing PTEN loss vs intact
module_eigengenes
DGEs <- c("SI","CRISP3","MFNG","CHGA","NOL4","LGI1","UNC13A","SLC44A5","VWA5B1","HFM1","GPR85","MS4A8","CST1","SCG2",
          "FAR2P1","POTEF","POTEI","POTEKP","RNU1-3","LINC01087","MTND1P23","CEACAM22P","AC002511.1","CYP4F62P","FAR2P3",
          "MIR2052HG","AC023790.2","AC060766.5","CEACAM20","RNA5-8SN2","RNA5-8SN2","AL157396.1")

genes <- genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% DGEs)]

modules_DGEs <- geneInfo$moduleColors[which(rownames(geneInfo) %in% genes)]

modules <- unique(modules_DGEs)

length(geneInfo$moduleColors[which(geneInfo$moduleColors =="green")])
data_DGEs_modules <- data.frame(module = modules,
                                fraction = c(length(genes[which(modules_DGEs == "grey")])/length(genes),
                                             length(genes[which(modules_DGEs == "green")])/length(genes),
                                             length(genes[which(modules_DGEs == "greenyellow")])/length(genes),
                                             length(genes[which(modules_DGEs == "blue")])/length(genes),
                                             length(genes[which(modules_DGEs == "magenta")])/length(genes)))
ggplot(data_DGEs_modules, aes(x =module, y = fraction, fill = module)) + 
  geom_bar(stat="identity", position ="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = c("blue", "green", "greenyellow", "grey","magenta")) +
  xlab("Modules WGCNA (PCa samples)") +
  ylab("Fraction of DEGs (PTEN loss vs intact) in the WGCNA modules")+
  labs(fill = "Modules") +
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing PTEN loss vs intact") +
  ylim(0,1)
  

DGEs[which(modules_DGEs == "green")]

DEGs_H_score <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_PTEN_loss_intact_H-score/Results/design ~ DV200 + Edad + H-score/Tables/resDESeq2_dds_DV200_Edad_H_score_cut_0_ContrastPTEN_loss_vs_intact_filtered.txt", sep = "\t")
DEGs_H_score$gene_name[which(DEGs_H_score$DEGs =="Yes" & DEGs_H_score$log2FoldChange <0)]
