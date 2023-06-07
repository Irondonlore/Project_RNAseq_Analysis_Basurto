################################################################################
###############  INTERSECTION OF THE DEGS AND WGCNA MODULES  ###################
###############        FOR THE DATA OF NDM VS LOC            ###################
################################################################################

# Previously, we perform an analysis over the data (joint) with WGCNA following
#the respective workflow. Additionally, we did the differential gene expression
#analysis with DESeq2 over the data NDM vs LILA considering the two conditions
#Now, we would like to intersect some of the modules found in WGCNA (the ones
#from my trait of interest) with the differential expressed genes detected
#by the analysis with DESeq. To do so, here it is implemented a VennDiagram
#to detect the genes (intersection) that appears in both analysis. The next step
#is to carry out an enrichment analysis to decipher what is the biological signi-
#cance of the genes detected.


# NOTE: I need to fix the problem of the ensemble ID because.

################################################################################
################################## LIBRARIES ###################################
################################################################################
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))

theme_set(theme_classic()) #| For graphic desig
################################################################################

################################################################################
#################################### DATA ######################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/"
data.file_WGCNA <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_vs_BPH/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors.txt"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + Diagnostico/Tables/resDESeq2_dds_DV200_Edad_Diagnostico_ContrastPCa_vs_BPH_filtered.txt"
data.file_DEGs2 <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_DESeq2_Analysis/4_DESeq2_Basurto_Diagnostico/Results/design ~ DV200 + Edad + DFS-STATUS/Tables/resDESeq2_dds_DV200_Edad_DFS-LA;LI_ContrastLA_vs_LI_filtered.txt"

setwd(dir.proj)

#| WGCNA DATA:
wgcna_Modules <- read.table(data.file_WGCNA, sep = "\t")

#| DESeq2 DATA:
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")
DEGs_genes2 <- read.table(data.file_DEGs2, sep = "\t")

#| Filtering those genes who satisfied the condition for DEGs classified by "Yes"
DEGs_genes <- DEGs_genes[which(DEGs_genes$DEGs == "Yes"),]
DEGs_genes_up <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange >=0),]
DEGs_genes_down <- DEGs_genes[which(DEGs_genes$DEGs == "Yes" & DEGs_genes$log2FoldChange <0),]

DEGs_genes2 <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes"),]
DEGs_genes_up2 <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes" & DEGs_genes2$log2FoldChange >=0),]
DEGs_genes_down2 <- DEGs_genes2[which(DEGs_genes2$DEGs == "Yes" & DEGs_genes2$log2FoldChange <0),]

DEGs_genes_GeneID <- rownames(DEGs_genes)
DEGs_genes_up_GeneID <- rownames(DEGs_genes_up)
DEGs_genes_down_GeneID <- rownames(DEGs_genes_down)

DEGs_genes_GeneID2 <- rownames(DEGs_genes2)
DEGs_genes_up_GeneID2 <- rownames(DEGs_genes_up2)
DEGs_genes_down_GeneID2 <- rownames(DEGs_genes_down2)


################################################################################


##################### FRACTION OF DEGS IN WGCNA MODULES ########################
############################# PCa VS BPH #######################################
modules <- unique(wgcna_Modules$moduleColors)
fraction <- c()
mod <- c()
name <- c()
size <- c()

for (i in 1:length(modules)){
  all <- length(DEGs_genes_GeneID[which(DEGs_genes_GeneID %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_GeneID)
  up <- length(DEGs_genes_up_GeneID[which(DEGs_genes_up_GeneID %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_up_GeneID)
  down <- length(DEGs_genes_down_GeneID[which(DEGs_genes_down_GeneID %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_down_GeneID)
  fraction <- c(all, up, down, fraction)
  mod <- c(modules[i], modules[i], modules[i], mod)
  name <- c(paste("All (",length(DEGs_genes_GeneID), ")",sep=""), paste("Up (",length(DEGs_genes_up_GeneID), ")",sep=""), paste("Down (", length(DEGs_genes_down_GeneID),")",sep=""), name)
  size <- c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), size)
}

mod2<- c()
for (i in 1:length(mod)){
  mod2 <- c(paste(mod[length(mod)-i+1],"\n (",size[length(mod)-i+1],")", sep=''), mod2)
}
dat <- data.frame(
  group = name,
  Modules = mod,
  y = fraction,
  size = size
)

#| Barplot
pdf("Results/Images/5_Fraction_DEGs_in_WGCNA_modules_PCa_vs_BPH.pdf")
ggplot(dat, aes(x=Modules, y=y, fill=group)) + 
  geom_bar(stat="identity", position ="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1, size =6),text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values =viridis(2)) +
  xlab("Modules WGCNA") +
  ylab("Fraction of DEGs in the WGCNA modules")+
  labs(fill = "DEGs") +
  ylim(0,1) +
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing PCa vs BPH") 
dev.off()

#| Bubble plot
pdf("Results/Images/5_Fraction_DEGs_in_WGCNA_modules_PCa_vs_BPH_bubble_plot.pdf")
ggplot(dat, aes(x = group, y = Modules)) + 
  geom_point(aes(size = size, fill = y), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "Fraction of DEGs \n in Modules") +
  labs(color ="Module size") +
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing PCa vs BPH") + 
  xlab("Up and Down regulated DEGs (PCa vs BPH)")
dev.off()
################################################################################


##################### FRACTION OF DEGS IN WGCNA MODULES ########################
############################# LA VS LI #######################################
modules <- unique(wgcna_Modules$moduleColors)
fraction <- c()
mod <- c()
name <- c()
size <- c()
for (i in 1:length(modules)){
  
  all <- length(DEGs_genes_GeneID2[which(DEGs_genes_GeneID2 %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_GeneID2)
  up <- length(DEGs_genes_up_GeneID2[which(DEGs_genes_up_GeneID2 %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_up_GeneID2)
  down <- length(DEGs_genes_down_GeneID2[which(DEGs_genes_down_GeneID2 %in% rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])])])/length(DEGs_genes_down_GeneID2)
  fraction <- c(all, up, down, fraction)
  mod <- c(modules[i], modules[i], modules[i], mod)
  name <- c(paste("All (",length(DEGs_genes_GeneID2), ")",sep=""), paste("Up (",length(DEGs_genes_up_GeneID2), ")",sep=""), paste("Down (", length(DEGs_genes_down_GeneID2),")",sep=""), name)
  size <- c(length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == modules [i])]), size)
}

mod2<- c()
for (i in 1:length(mod)){
  mod2 <- c(paste(mod[length(mod)-i+1],"\n (",size[length(mod)-i+1],")", sep=''), mod2)
}

dat <- data.frame(
  group = name,
  Modules = mod,
  y = fraction,
  size = size
)

pdf("Results/Images/5_Fraction_DEGs_in_WGCNA_modules_LA_vs_LI.pdf")
ggplot(dat, aes(x=Modules, y=y, fill=group)) + 
  geom_bar(stat="identity", position ="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1, size =6),text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values =viridis(2)) +
  xlab("Modules WGCNA") +
  ylab("Fraction of DEGs in the WGCNA modules")+
  labs(fill = "DEGs") +
  ylim(0,1)+
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing LA vs LI")
dev.off()

#| Bubble plot
pdf("Results/Images/5_Fraction_DEGs_in_WGCNA_modules_LA_vs_LI_bubble_plot.pdf")
ggplot(dat, aes(x = group, y = Modules)) + 
  geom_point(aes(size = size, fill = y), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "Fraction of DEGs \n in Modules") +
  labs(color ="Module size") +
  ggtitle("Evaluating fraction of DEGs in each module of WGCNA \nDEGs obtained by comparing LA vs LI") +
  xlab("Up and Down regulated DEGs (LA vs LI)")
dev.off()

################################################################################





###############################  GRID PLOT  ####################################






############################# INTERSECTED GENES ################################
################################################################################
myCol <- brewer.pal(0, "Pastel2")

for (color in unique(wgcna_Modules$moduleColors)){
  
  #| Extracting the genes from the modules of the wgcna data
  modules <- rownames(wgcna_Modules)[wgcna_Modules$moduleColors == color]
  modules <- modules[!is.na(modules)]
  
  #| Matching the DEGs and the wgcna genes
  genes <- match(modules, DEGs_genes)
  genes <- genes[!is.na(genes)]
  genes <- DEGs_genes[genes]
  genes <- as.data.frame(genes)
  
  #Creating a new directory
  dir.create(paste(results.file,tag, sep =""))
  
  #| Saving the intersected genes
  write.table(genes, paste(results.file, tag, "Intersected_genes_DEGs_WGCNA_design_", tag, "_Module",color,".txt", sep =""), sep ="\t")
  
  #| Visualization: VENN DIAGRAM
  venn.diagram(x = list(modules, DEGs_genes),
               category.names = c(sprintf("Module %s", color), "DEGs"),
               filename = sprintf("%s%s/Venn_diagramm_DEGs_Module%s_design_%s.png",results.file, tag, color, tag),
               output=FALSE,
               imagetype="png" ,
               height = 980 , 
               width = 980 , 
               resolution = 400,
               compression = "lzw",
               # Circles
               lwd = 2,
               col = c(color, "black"),
               fill = c(alpha(color,0.3), alpha("white",)),
               
               # Numbers
               cex = .6,
               fontface = "bold",
               fontfamily = "sans",
               
               # Set names
               cat.cex = 0.5,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.fontfamily = "serif")
}

