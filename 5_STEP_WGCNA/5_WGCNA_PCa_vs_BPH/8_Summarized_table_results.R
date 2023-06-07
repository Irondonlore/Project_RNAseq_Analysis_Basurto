################################################################################
####    SUMMARIZED TABLE OF RESULTS
################################################################################

#| Results from 
#|    1) Cox with DFS.TIME and gene expression
#|    2) Correlation with DFS.TIME and gene expression
#|    3) Correlation with FOXO mean vector and gene expression
#|    4) WGCNA results

################################################################################

dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/"
setwd(dir.proj)

################################################################################
###     DATA 
################################################################################

#| 1)
cox_DFSTIME <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Disease_Free_Survival_Different_Datasets/Results/Tables/Cox_Analysis_Genes_DFSTIME_Basurto.txt", sep = "\t")
cox_DFSTIME <- cox_DFSTIME[,c("GeneID","p_values_basurto", "padj_basurto", "coef_basurto", "exp_coef_basurto")]
colnames(cox_DFSTIME) <- c("GeneID","p_values_basurto_cox_DFS", "padj_basurto_cox_DFS", "coef_basurto_cox_DFS", "exp_coef_basurto_cox_DFS")

#| 2)
correlation_DFSTIME <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Correlation_Gene_DFSTIME/Results/Tables/Correlation_Analysis_Genes_DFSTIME_Basurto.txt", sep = "\t")
correlation_DFSTIME <- correlation_DFSTIME[,c("GeneID", "p_values_basurto", "padj_basurto", "cor_basurto")]
colnames(correlation_DFSTIME) <- c("GeneID", "p_values_basurto_cor_DFS", "padj_basurto_cor_DFS", "cor_basurto_DFS")

#| 3)
correlation_FOXO <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Tables/Correlation_Analysis_Genes_Expression_FOXO_Basurto.txt", sep = "\t")
correlation_FOXO <- correlation_FOXO[, c("GeneID", "p_values_basurto", "padj_basurto", "cor_basurto")]
colnames(correlation_FOXO) <- c("GeneID", "p_values_basurto_cor_FOXO", "padj_basurto_cor_FOXO", "cor_basurto_FOXO")

#| 4)
wgcna <- read.table("Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors.txt", sep ="\t")
wgcna$GeneID <- rownames(wgcna)

#| Checking dimensionality
dim(cox_DFSTIME)
dim(correlation_DFSTIME)
dim(correlation_FOXO)
dim(wgcna)

#| Joining datasets in an unique data frame
joined_table <- merge()

joined_table <- list(wgcna, cox_DFSTIME, correlation_DFSTIME, correlation_FOXO)  
Reduce(function(x, y) merge(x, y, by="GeneID"), joined_table)  

