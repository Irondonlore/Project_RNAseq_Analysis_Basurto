################################################################################
##### RELATING MODULES TO EXTERNAL INFORMATION WITH DATA AC-45_RNAseq-FFPE #####
################################################################################

#| This script related the module found with blockwise method with the clinical
#| traits of the AC-45_RNAseq-FFPE data. The idea is to find those modules having
#| high correlation with the traits of interest

################################################################################


################################################################################
###############################  LIBRARIES  ####################################
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))

options(stringsAsFactors = FALSE)
################################################################################


################################################################################
#################################  DATA  #######################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/"
tag <- "AC-45_RNAseq-FFPE"

setwd(dir.proj)

load(file = paste(results.file, "Data/", "2_blockwiseModules_AC_45_RNAseq_FFPE.RData", sep =""))
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
  pdf(paste(results.file,"Images/3_SingleBlock_geneDendogram_and_ModuleColors.pdf", sep =""))
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
  
  #| 1) Diagnosis, Edad, DV200, Purity
  datTraits_1 <- subset(datTraits, select=c(Diagnostico, Edad, DV200.value, purity, DFS.STATUS))
  moduleTraitCor_1 <- cor(MEs, datTraits_1, use = "p")
  moduleTraitPvalue_1 <- corPvalueStudent(moduleTraitCor_1, nSamples)
  textMatrix_1 = paste(signif(moduleTraitCor_1, 2), "\n(",
                     signif(moduleTraitPvalue_1, 1), ")", sep = "")
  pdf(paste(results.file, "Images/3_Module-trait relationships_Diagnosis_Edad_DV200_purity_DFS-STATUS.pdf", sep =""))
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
  
  #| 2) Diagnostico, Mean.expression.PI3K.AKT.mTOR, H_score_cut_0, H_score_cut_10, H_score_cut_20, H_score_cut_30, pAKT_positive
  datTraits_2 <- subset(datTraits, select=c(Diagnostico,PTEN_Exp_log2, H_score_cut_0, H_score_cut_10, H_score_cut_20, H_score_cut_30, pAKT_positive))
  moduleTraitCor_2 <- cor(MEs, datTraits_2, use = "p")
  moduleTraitPvalue_2 <- corPvalueStudent(moduleTraitCor_2, nSamples)
  textMatrix_2 = paste(signif(moduleTraitCor_2, 2), "\n(",
                       signif(moduleTraitPvalue_2, 1), ")", sep = "")
  pdf(paste(results.file, "Images/3_Module-trait relationships_Diagnosis_PTEN_Exp_H_score_cut_0_10_20_30_pAKT_positive.pdf", sep =""))
  labeledHeatmap(Matrix = moduleTraitCor_2,
                 xLabels = names(datTraits_2),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix_2,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 cex.lab.y= 0.5,
                 cex.lab.x= 0.7,
                 cex.lab= 0.6,
                 main = paste("Module-trait relationships"))
  dev.off()
  
  #| 3) Diagnostico, signatures
  datTraits_3 <- subset(datTraits, select=c(Mean.expression.PI3K.AKT.mTOR, Mean.expression.FOXO, Mean.expression.FOXO.pathway, Mean.expression.PTEN_loss))
  moduleTraitCor_3 <- cor(MEs, datTraits_3, use = "p")
  moduleTraitPvalue_3 <- corPvalueStudent(moduleTraitCor_3, nSamples)
  textMatrix_3 = paste(signif(moduleTraitCor_3, 2), "\n(",
                       signif(moduleTraitPvalue_3, 1), ")", sep = "")
  pdf(paste(results.file, "Images/3_Module-trait relationships_Diagnosis_Signatures.pdf", sep =""))
  labeledHeatmap(Matrix = moduleTraitCor_3,
                 xLabels = names(datTraits_3),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix_3,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 cex.lab.y= 0.5,
                 cex.lab.x= 0.7,
                 cex.lab= 0.6,
                 main = paste("Module-trait relationships"))
  dev.off()

  
  #| 4) Extracting the most important variables
  datTraits_4 <- subset(datTraits, select=c(Diagnostico, purity, Mean.expression.PI3K.AKT.mTOR, Mean.expression.FOXO, Mean.expression.FOXO.pathway, Mean.expression.PTEN_loss, H_score_cut_0, H_score_cut_20, pAKT_positive))
  moduleTraitCor_4 <- cor(MEs, datTraits_4, use = "p")
  moduleTraitPvalue_4 <- corPvalueStudent(moduleTraitCor_4, nSamples)
  textMatrix_4 = paste(signif(moduleTraitCor_4, 2), "\n(",
                       signif(moduleTraitPvalue_4, 1), ")", sep = "")
  pdf(paste(results.file, "Images/3_Module-trait relationships_Diagnosis_Signatures_H-score_pAKT_activity.pdf", sep =""))
  labeledHeatmap(Matrix = moduleTraitCor_4,
                 xLabels = names(datTraits_4),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix_4,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 cex.lab.y= 0.5,
                 cex.lab.x= 0.4,
                 cex.lab= 0.6,
                 main = paste("Module-trait relationships"))
  dev.off()
  
  names(datTraits)
  
  
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
  write.table(geneInfo, paste(results.file, "Tables/3_geneInfo_counts_geneMM_geneS_moduleColors.txt", sep =""), row.names = T, sep = "\t")

}


#| Intramodular analysis: identifying genes with high GS and MM
#| Interesting modules: Green, turquoise, blue, magenta
#| Using the GS and MM measures, we can identify genes that have a high significance 
#| for a given trait as well as high module membership in interesting modules

module = "magenta"
trait ="Mean_expression_FOXO"
n_trait = 7
column = match(module, modNames)
moduleGenes = moduleColors==module


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
pdf(paste("Results/Images/3_MM_vs_GS_Module",module,"_",trait, ".pdf",sep =""))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 7]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, 
                   family="serif")
dev.off()


###############################################################################
#################         MODULES EXPLORATION        ###########################
################################################################################

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| WGCNA Modules per gene
wgcna_Modules <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_WGCNA/Results/Tables/3_WGCNA_Summarized_Table.txt", sep = "\t")

#| Module: green
color <- "green"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")
modules$gene_name[1:500]
total_DEGs_gost <- gost(list("Enrichment Module Green" = modules$GeneID), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(modules), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_DEGs_gost, interactive = FALSE)
results <- total_DEGs_gost$result[order(total_DEGs_gost$result$p_value),]
results$term_name[1:20]

#| Module: blue
color <- "blue"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")
modules$gene_name[1:20]
total_DEGs_gost <- gost(list("Enrichment Module Blue" = modules$GeneID), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(modules), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_DEGs_gost, interactive = FALSE)
results <- total_DEGs_gost$result[order(total_DEGs_gost$result$p_value),]
results$term_name[1:20]

#| Module: Magenta
color <- "magenta"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")
modules$gene_name[1:20]
total_DEGs_gost <- gost(list("Enrichment Module Magenta" = modules$GeneID), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(modules), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_DEGs_gost, interactive = FALSE)
results <- total_DEGs_gost$result[order(total_DEGs_gost$result$p_value),]
results$term_name[1:30]


#| Module: Turquoise
color <- "turquoise"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")
modules$gene_name[1:20]
total_DEGs_gost <- gost(list("Enrichment Module Magenta" = modules$GeneID), 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = rownames(modules), 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_DEGs_gost, interactive = FALSE)
results <- total_DEGs_gost$result[order(total_DEGs_gost$result$p_value),]
results$term_name[1:20]



################################################################################
#|  Relating genes modules with genes identified in the intersection of the two dataset

#| Module: magenta
color <- "magenta"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")

#| Venn Diagram

#| UP 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))])
common[which(duplicated(common))]
# [1] "HPN"     "NPY"     "GDF15"   "ACSM1"   "GJB1"    "TRGC2"   "RF00017"
venn.diagram(x = list(modules$gene_name, DEGs_up$gene_name, DEGs_up_tT$gene_name),
             category.names = c(sprintf("Module %s", color), "UP DEGs Basurto", "UP DEGs GSE134051"),
             filename = "Results/VennDiagram/Diagnostico/Venn_Magenta_UP_DEGs_Basurto_and_GSE134051.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c(color, "black", "yellow"),
             fill = c(alpha(color,0.3), alpha("white",), alpha("yellow",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")

#| DOWN 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))])
common[which(duplicated(common))]
#| None


#| Module: green
color <- "green"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")

#| Venn Diagram

#| UP 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))])
common[which(duplicated(common))]
#| [1] "FOLH1"      "AGR2"       "GOLM1"      "CLDN8"      "ERG"        "OR51E2"     "SLC38A11"  
#| [8] "GCNT1"      "PLA2G2A"    "OR51C1P"    "AC144450.1" "TRGJP2"     "TRGC1"      "TRGJP"     
#| [15] "TRGV9"      "PCA3"       "AC141930.1" "AMACR"      "RPL7P16"  

venn.diagram(x = list(modules$gene_name, DEGs_up$gene_name, DEGs_up_tT$gene_name),
             category.names = c(sprintf("Module %s", color), "UP DEGs Basurto", "UP DEGs GSE134051"),
             filename = "Results/VennDiagram/Diagnostico/Venn_Green_UP_DEGs_Basurto_and_GSE134051.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c(color, "black", "yellow"),
             fill = c(alpha(color,0.3), alpha("white",), alpha("yellow",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")

#| DOWN 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))])
common[which(duplicated(common))]
#| None


#| Module: blue
color <- "blue"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")

#| Venn Diagram

#| UP 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))])
common[which(duplicated(common))]
#| [1] RF00017 

#| DOWN 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))])
common[which(duplicated(common))]
#| "BMP5"  "NELL2"

venn.diagram(x = list(modules$gene_name, DEGs_down$gene_name, DEGs_down_tT$gene_name),
             category.names = c(sprintf("Module %s", color), "DOWN DEGs Basurto", "DOWN DEGs GSE134051"),
             filename = "Results/VennDiagram/Diagnostico/Venn_Blue_DOWN_DEGs_Basurto_and_GSE134051.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c(color, "black", "yellow"),
             fill = c(alpha(color,0.3), alpha("white",), alpha("yellow",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")


#| Module: Turquoise
color <- "turquoise"
modules <- rownames(geneInfo)[geneInfo$moduleColors == color]
modules <- modules[!is.na(modules)]
modules <- data.frame(GeneID = modules)
modules <- merge(modules ,genome_GRCh39.94, by = "GeneID")

#| Venn Diagram

#| UP 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_up[which(duplicated(common_genes_PCa_vs_BPH_up))])
common[which(duplicated(common))]
#| [1] RF00017 

#| DOWN 
common <- c(modules$gene_name[which(!duplicated(modules$gene_name))], common_genes_PCa_vs_BPH_down[which(duplicated(common_genes_PCa_vs_BPH_down))])
common[which(duplicated(common))]
#| "CFD"

venn.diagram(x = list(modules$gene_name, DEGs_down$gene_name, DEGs_down_tT$gene_name),
             category.names = c(sprintf("Module %s", color), "DOWN DEGs Basurto", "DOWN DEGs GSE134051"),
             filename = "Results/VennDiagram/Diagnostico/Venn_Turquoise_DOWN_DEGs_Basurto_and_GSE134051.png",
             output=FALSE,
             imagetype="png" ,
             height = 980 , 
             width = 980 , 
             resolution = 400,
             compression = "lzw",
             # Circles
             lwd = 2,
             col = c(color, "black", "yellow"),
             fill = c(alpha(color,0.3), alpha("white",), alpha("yellow",0.3)),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "serif")
