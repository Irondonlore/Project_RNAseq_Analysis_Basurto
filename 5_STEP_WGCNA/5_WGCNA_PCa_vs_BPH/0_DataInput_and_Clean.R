################################################################################
#####             DATA INPUT AND CLEAN FOR WGCNA ANALYSIS                    ###
################################################################################

#| * Include a description * 

################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(skimr))
suppressMessages(library(DataExplorer))

options(stringsAsFactors = FALSE)
################################################################################


################################################################################
################################    DATA    ####################################
####################################### #########################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/VST_Counts/FullCounts_Basurto_cpm_filtered_phenotype_VST.txt"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/"
tag <- "AC-45_RNAseq-FFPE"

counts_data <- read.table(counts.file, header=T) 
datTraits <- read.table(info.file, sep ="\t", header =T)
setwd(dir.proj)

#################    DEALING WITH CATEGORICAL VARIABLES    #####################

#| Assigning 0 to BPH and 1 to PCa
datTraits$Diagnostico[datTraits$Diagnostico == "BPH"] <- 0
datTraits$Diagnostico[datTraits$Diagnostico == "PCa"] <- 1

#| Assigning 0 to PTEN loss and 1 to PTEN intact CUT 0
datTraits$H_score_cut_0[datTraits$H_score_cut_0 == "PTEN loss"] <- 1
datTraits$H_score_cut_0[datTraits$H_score_cut_0 == "PTEN intact"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN intact CUT 10
datTraits$H_score_cut_10[datTraits$H_score_cut_10 == "PTEN loss"] <- 1
datTraits$H_score_cut_10[datTraits$H_score_cut_10 == "PTEN intact"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN intact CUT 20
datTraits$H_score_cut_20[datTraits$H_score_cut_20 == "PTEN loss"] <- 1
datTraits$H_score_cut_20[datTraits$H_score_cut_20 == "PTEN intact"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN intact CUT 30
datTraits$H_score_cut_30[datTraits$H_score_cut_30 == "PTEN loss"] <- 1
datTraits$H_score_cut_30[datTraits$H_score_cut_30 == "PTEN intact"] <- 0

#| Assigning 0 to No info of pAKT and 1 to positive for pAKT 
datTraits$pAKT_positive[which(datTraits$pAKT_positive == "No pAKT")] <- 0
datTraits$pAKT_positive[which(datTraits$pAKT_positive == "pos")] <- 1

#| In the case of the variable group (or those which have more than two levels) 
#| it is better to separate then into groups of categorical variable using the 
#| WGCNA function called binarizeCategoricalVariable()
#binarize_out = binarizeCategoricalVariable(datTraits$Group, includePairwise = TRUE)
#datTraits <- cbind(datTraits, binarize_out)

#| Converting to numeric values
datTraits$Diagnostico <- as.numeric(datTraits$Diagnostico)
datTraits$H_score_cut_0 <-  as.numeric(datTraits$H_score_cut_0)
datTraits$H_score_cut_10 <- as.numeric(datTraits$H_score_cut_10)
datTraits$H_score_cut_20 <- as.numeric(datTraits$H_score_cut_20)
datTraits$H_score_cut_30 <- as.numeric(datTraits$H_score_cut_30)
datTraits$pAKT_positive <- as.numeric(datTraits$pAKT_positive)
names(datTraits)
#| Selecting the columns of interest
datTraits <- subset( datTraits, select = c( Diagnostico, 
                                            Edad, 
                                            PTEN_Exp_log2, 
                                            DV200.value, 
                                            DFS.TIME,
                                            DFS.STATUS,
                                            purity,
                                            Mean.expression.FOXO,
                                            Mean.expression.FOXO.pathway,
                                            Mean.expression.PTEN_loss,
                                            Mean.expression.PI3K.AKT.mTOR,
                                            H_score_cut_0,
                                            H_score_cut_10,
                                            H_score_cut_20,
                                            H_score_cut_30,
                                            pAKT_positive))

################################################################################


################################################################################
##############################    DATA CLEAN   #################################
################################################################################
if (dim(counts_data)[2] == dim(datTraits)[1]){
  
  #| We need to transpose the expression data for further analysis.
  datExpr <- as.data.frame(t(counts_data))
  
  #| Checking data for excessive missing values and identification of oulier microarray 
  #| samples
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  gsg$allOK
  
  #| IF the last statement is TRUE, all genes have passed the cuts. If not, we remove 
  #| the offending genes and samples from the data 
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  #| Next, we cluster the samples (which is not the same as the clustering genes 
  #| that will come later) to see if the are any obvious outliers
  sampleTree <- hclust(dist(datExpr), method ="average")
  
  #| Plot the sample tree: Open a graphic output window of size 12 by 9 inches the 
  #| user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  pdf(file = paste(results.file,"Images/0_sampleClustering.pdf", sep =""))
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
 
}

#| IS THERE ANY OUTLIERS?
outlier = FALSE

if (outlier){
  #| IF TRUE: We can remove it by hand or by choosing a height cut that will remove 
  #| the offending sample. We can draw a line and determine the cluster under the 
  
  h <- 90
  color <- "red"
  abline(h=h, col=color)
  clust <- cutreeStatic(sampleTree, cutHeight = 90, minSize = 60)
  table(clust)
  
  #| IF clust 1 contains the samples we want to keep we do:
  keepSamples <- (clust == 1) 
  datExpr <- datExpr[keepSamples,]
}

#| Before continuing with network construction and module detection, we visualize 
#| how the clinical traits relate to the sample dendrogram. But first we need to 
#| convert traits to a color representation: white means low, read high, grey 
#| means missing entry
if (dim(counts_data)[2] == dim(datTraits)[1]){

  collectGarbage()
  
  #| Associating numbers to colors
  traitColors = numbers2colors(datTraits, signed = FALSE)
  
  #| Plot the sample dendrogram and the colors underneath
  pdf(paste(results.file, "Images/0_sample_dendrogram_and_trait_heatmap.pdf", sep =""))
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap", cex.colorLabels = 0.5, cex.dendroLabels = 0.2, family ="serif")
  dev.off()
  
  save(datExpr, datTraits, file = paste(results.file,"Data/", tag, "_dataInput.RData", sep =""))
  
}




