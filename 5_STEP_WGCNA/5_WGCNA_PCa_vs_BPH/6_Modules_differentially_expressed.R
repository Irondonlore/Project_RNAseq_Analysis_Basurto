################################################################################
#####   FINDING THE MODULES HAVING THE MOST SIGNIFICANT DIFFERENT AMONG    #####
#####                            CONDITIONS                                #####                        
################################################################################

#| This code uses the Limma package to fix linear models on the Eigengene of every
#| module of WGCNA to analyze differences across one or multiple variables indicated
#| in the design matrix (which are associated to clinical, histological or biochemical
#| parameters). 
#| https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html

################################################################################


################################################################################
###############################  LIBRARIES  ####################################
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(limma))
options(stringsAsFactors = FALSE)
################################################################################



################################################################################
#################################  DATA  #######################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/"
tag <- "AC-45_RNAseq-FFPE"

setwd(dir.proj)

data1 <- load(file = paste(results.file, "Data/", "2_blockwiseModules_AC_45_RNAseq_FFPE.RData", sep =""))
data2 <- load(file = paste(results.file, "Data/", tag, "_dataInput.RData", sep =""))
wgcna_Modules <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/Tables/3_WGCNA_Summarized_Table.txt",  sep = "\t")

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################



################################################################################
#| Module eigengenes
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
module_eigengenes <- MEs

#| Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~datTraits$Edad  + datTraits$DV200.value + datTraits$H_score_cut_0)
module_eigengenes <- module_eigengenes[which(!is.na(datTraits$H_score_cut_0)),]
dim(des_mat)
#Run linear model on each module. Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

#Apply multiple testing correction and obtain stats in a data frame.

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

## Removing intercept from test coefficients
#Lets take a look at the results. They are sorted with the most significant results at the top.

head(stats_df)
stats_df


greenyellow, magenta, blue, green, pink