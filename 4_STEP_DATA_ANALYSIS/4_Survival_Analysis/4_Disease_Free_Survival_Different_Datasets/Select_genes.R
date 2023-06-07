



################################################################################
library(gprofiler2)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(viridis)
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets"
setwd(workingDir)
################################################################################


data_DFS_Cox <- read.table("Results/Tables/Final_table_Cox_Analysis_Genes_DFSTIME_Basurto_Taylor_TCGA_Glinsky.txt", sep = "\t", header= T)

better_all <- data_DFS_Cox$gene_name[which(data_DFS_Cox$total_coherence_halfdatasets == "YES_BETTER" &  data_DFS_Cox$sig_datasets>=3)]
worse_all <- data_DFS_Cox$gene_name[which(data_DFS_Cox$total_coherence_halfdatasets == "YES_WORSE" &  data_DFS_Cox$sig_datasets>=3)]
worse_all
better_all
any(is.na(data_DFS_Cox$gene_name))


better_basurto$gene_name
#| Enrichment

gost_better <- gost(list("Better prognosis" = better_basurto$GeneID),
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                    measure_underrepresentation = FALSE, evcodes = TRUE,
                    user_threshold = 0.05, correction_method = "fdr",
                    domain_scope = "custom", custom_bg = dataframe_basurto$GeneID, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_better, interactive = FALSE)
gost_better$result$term_name[1:30]  

gost_worse <- gost(list("Worse prognosis" = worse_basurto$GeneID),
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = dataframe_basurto$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_worse, interactive = FALSE)
gost_worse$result$term_name[1:30]  
