################################################################################
########## VISUALIZATION OF THE NETWORK TOPOLOGY FROM THE TOM  #################
################################################################################

#| Visualization of the TOM from the WGCNA results of the RNAseq analysis for the 
#| AC-45_RNA-FFPE data.

#|  In theory, The topological overlap of two nodes reflects their similarity in 
#| terms of the commonality of the nodes they connect to. This matrix is used
#| to visualize the network topology. It provides a better visualization than the
#| one obtained with the adjacency matrix (confirm this quote)
################################################################################



################################################################################
###############################  LIBRARIES  ####################################
################################################################################
library(WGCNA)
library(viridis)
options(stringsAsFactors = FALSE)
################################################################################


################################################################################
#################################  DATA  #######################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_WGCNA/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_WGCNA/Results/"
setwd(dir.proj)

# Loading the data 
load(file = paste(results.file, "Data/", "2_blockwiseModules_AC_45_RNAseq_FFPE.RData", sep =""))
load(file = paste(results.file, "Data/", "2_blockwiseModules_TOM_AC_45_RNAseq_FFPE-block.1.RData", sep =""))
load(file = paste(results.file, "Data/", "AC-45_RNAseq-FFPE_dataInput.RData", sep =""))

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
################################################################################



################################################################################
#################    VISUALIZING THE NETWORK WITH THE TOM    ###################

#| The TOM contained in this space comes in a dist format. It is necessary to trans-
#| form it into a matrix object. In R, the dist() function is used to compute a 
#| distance matrix. But the result you get back isn't really a matrix, it's a "dist"
#| object. Under the hood, the "dist" object is stored as a simple vector. 
TOM = as.matrix(TOM)

#| Computing the dissimilarity TOM (who is useful for separating genes and construct
#| the hieralchical clustering)
dissTOM = 1-TOM

#| Transform dissTOM with a power to make moderately strong connections more 
#| visible in the heatmap (power of..?)
plotTOM = dissTOM^10

#| Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

#| Open a new windows
sizeGrWindow(9,9)

#| This plot will contain the hieralchical clustering for ALL GENES. This takes
#| too much time
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

TOMplot
#| NOTE: With many genes this could mean a problem when running in the personal
#| computer. An alternative approach for visualizing the network quickly is restrict 
#| the space of possible genes in our network:

#| Select a particular amount of genes
nSelect =5040

#| For reproducibility, we set the random seed. This is because the selection of
#| the subset of genes is random.
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

#| Comment from the manual: There's no simple way of restricting a clustering tree  
#| to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

#| Open a graphical window
sizeGrWindow(9,9)

#| Taking the dissimilarity to a power, let's say 10, makes the plot more informative 
#| by effectively changing the color palette; setting the diagonal to NA also 
#| improves the clarity of the plot
plotDiss = selectTOM^10

#| To have all the squares of the diagonal in white.
diag(plotDiss) = NA

#pdf(file = "Images/4_NetworkHeatmapSelectedGenes.pdf")
# Plotting the restricted-in-space TOM :)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
#dev.off()
################################################################################

#| Plot the dendrogram
par(cex = 1.0, family= "serif")
pdf("Results/Images/4_Eigengene_dendrogram.pdf")
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

par(cex = 1.0, family= "serif")
pdf("Results/Images/4_Eigengene_adjaccency_heatmap.pdf")
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, heatmapColors = rocket(30))
dev.off()
