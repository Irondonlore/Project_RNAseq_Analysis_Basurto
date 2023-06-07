################################################################################
##### RELATING MODULES TO EXTERNAL INFORMATION WITH DATA AC-45_RNAseq-FFPE #####
################################################################################

#| This script related the module found with blockwise method with the clinical
#| traits of the AC-45_RNAseq-FFPE data. The idea is to find those modules having
#| high correlation with the traits of interest

################################################################################



###############################  LIBRARIES  ####################################
suppressMessages(library(WGCNA))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(RColorBrewer))
options(stringsAsFactors = FALSE)
################################################################################


################################################################################
#################################  DATA  #######################################
################################################################################
dir.proj <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/"
results.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/"
tag <- "AC-45_RNAseq-FFPE"

setwd(dir.proj)

load(file = paste(results.file, "Data/", "2_blockwiseModules_AC_45_RNAseq_FFPE_new_2.RData", sep =""))
load(file = paste(results.file, "Data/", tag, "_dataInput.RData", sep =""))

data_loaded = TRUE
names(datTraits)

################################################################################


if(data_loaded){
  
  #| Relabel blockwise modules
  bwLabels = matchLabels(net$colors, moduleLabels)
  
  #| Number of modules identified
  nModules <- length(table(bwLabels))
  
  #| Convert labels to colors for plotting
  bwModuleColors = labels2colors(bwLabels)
  
  #| Plotting dendogram and module colors
  sizeGrWindow(12,9)
  pdf(paste(results.file,"Images/3_SingleBlock_geneDendogram_and_ModuleColors_new.pdf", sep =""))
  plotDendroAndColors(geneTree,
                      cbind(bwModuleColors),
                      c("Modules"),
                      main = "Single block gene dendrogram and module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  
  ################    Quantifying module-trait associations    #################   
  
  #| Define numbers of genes and samples
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  #| Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  #| Plotting Module-Trait relationships
  sizeGrWindow(8,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 15, 2, 0.5), family ="serif");
  
  # Display the correlation values within a heatmap plot
  
  #| 1) Edad, DV200, Purity, DFS.TIME, H_score_cut_0, Signatures (FOXO, PTEN loss, PI3K-AKT-mTOR)
  datTraits_1 <- subset(datTraits, select=c(Edad, DV200.value, purity, DFS.STATUS, DFS.TIME))
  moduleTraitCor_1 <- cor(MEs, datTraits_1, use = "p")
  moduleTraitPvalue_1 <- corPvalueStudent(moduleTraitCor_1, nSamples)
  textMatrix_1 = paste(signif(moduleTraitCor_1, 2), "\n(",
                     signif(moduleTraitPvalue_1, 1), ")", sep = "")
  pdf(paste(results.file, "Images/3_Module-trait relationships_Edad_DV200_purity_DFS-TIME_DFS-STATUS_new.pdf", sep =""))
  labeledHeatmap(Matrix = moduleTraitCor_1,
                 xLabels = names(datTraits_1),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix_1,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 cex.lab.y= 0.5,
                 cex.lab.x= 0.7,
                 cex.lab= 0.6,
                 main = paste("Module-trait relationships"))
  dev.off()
  
  #| 2) Edad, DV200, Purity, DFS.TIME, H_score_cut_0, Signatures (FOXO, PTEN loss, PI3K-AKT-mTOR), pAKT_positive
  datTraits_2 <- subset(datTraits, select=c(Edad, DV200.value, purity, stromal,immune, H_score_cut_0, Mean.expression.PI3K.AKT.mTOR))
  colnames(datTraits_2) <- c("Age", "DV200", "Purity", "Stromal", "Immune", "PTEN_status", "PI3K-AKT-mTOR")
  moduleTraitCor_2 <- cor(MEs, datTraits_2, use = "p")
  moduleTraitPvalue_2 <- corPvalueStudent(moduleTraitCor_2, nSamples)
  textMatrix_2 = paste(signif(moduleTraitCor_2, 2), "\n(",
                       signif(moduleTraitPvalue_2, 1), ")", sep = "")
  sizeGrWindow(12,11)
  pdf(paste(results.file, "Images/3_Module-trait relationships_Edad_DV200_purity_H_score_Signatures_PTENLOSS_FOXO_PI3K-AKT-mTOR_stromal_immune_2.pdf", sep =""))
 
  labeledHeatmap(Matrix = moduleTraitCor_2,
                 xLabels = names(datTraits_2),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 #colors = greenWhiteRed(50),
                 textMatrix = textMatrix_2,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 cex.lab.y= 0.45,
                 cex.lab.x= 0.4,
                 cex.lab= 0.55,
                 colors  = colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                 main = paste("Module-trait relationships"))
  dev.off()#ggsave("Results/Images/3_Module-trait relationships_Edad_DV200_purity_H_score_Signatures_PTENLOSS_FOXO_PI3K-AKT-mTOR_stromal_immune.pdf", height= 5, width = 6)
  
  #| Eigengene heatmap
  pdf("Results/Images/3_Eigengene_heatmap_2.pdf")
  plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90,
                        heatmapColors =colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                        letterSubPlots=T ) 
  dev.off()
  
  ###############  GENE SIGNIFICANCE AND MODULE MEMBERSHIP  ####################
  
  #|       Module Membership
  
  #| Name of the color modules
  modNames <- substring(names(MEs), 3)
  
  #| Gene Module Membership and MM p-value
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use ="p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  #|       Gene significance
  
  #| Name of the traits
  traitNames <- names(datTraits)
  
  #| Gene Trait significance
  geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
  names(GSPvalue) = paste("p.GS.", traitNames, sep="")
  
  #| Saving the information
  geneInfo <- data.frame(t(datExpr),geneModuleMembership, 
                         geneTraitSignificance, moduleColors, 
                         row.names=colnames(datExpr))
  
  #| Summarized table 
  write.table(geneInfo, paste(results.file, "Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_new.txt", sep =""), row.names = T, sep = "\t")

}

#| Ploting GS vs MM: Green module
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
pdf("Results/Images/3_GS_MM_PTEN_loss_Green.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 15]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#| Ploting GS vs MM: Magenta module
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
pdf("Results/Images/3_GS_MM_PTEN_loss_Magenta.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 15]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
