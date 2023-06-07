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
library(igraph)

options(stringsAsFactors = FALSE)
################################################################################


#################################  DATA  #######################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/"
tag <- "AC-45_RNAseq-FFPE"

#| Setting working dir
setwd(dir.proj)

#| Loading the data
load(file = "2_blockwiseModules_TOM_AC_45_RNAseq_FFPE_new_2-block.1.RData")

load(file = paste(results.file, "Data/", "2_blockwiseModules_AC_45_RNAseq_FFPE_new_2.RData", sep =""))
load(file = paste(results.file, "Data/", tag, "_dataInput.RData", sep =""))

#| Number of genes and number of samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Read in the annotation file
annot <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Module information
module_info <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_new.txt")

#| TOM as matrix
TOM <- as.matrix(TOM)
rownames(TOM) <- names(datExpr)
colnames(TOM) <- names(datExpr)
################################################################################


####################   EXPORTING NETWORK TO CYTOSCAPE    #######################
modules <- "green"

#| Filtering genes by GS and MM
genes <- rownames(module_info)[which(module_info$GS.H_score_cut_0>0.2 & module_info$MMmagenta>0.6 & module_info$moduleColors == "green")]

#| Select module probes
probes <- names(datExpr[,genes])

inModule = is.finite(match(names(datExpr),genes))
modProbes = genes
modGenes = annot$gene_name[which(annot$GeneID %in%genes)];
head(TOM)
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.09,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
################################################################################


################################################################################
####################   EXPORTING NETWORK TO CYTOSCAPE    #######################

# Select modules to explore
modules <- c("black", "darkred", "green")

#| Selecting those genes with the most biological significance
genes <- rownames(module_info)[which(module_info$GS.H_score_cut_0>0.09)]
module_info$moduleColors[which(rownames(module_info) %in% genes)]
#| Select module probes
probes <- names(datExpr[])
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[which(inModule)]
modGenes <- annot$GeneID[which(modProbes%in%annot$GeneID)]

#| Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

colnames(modTOM) <- annot$gene_name[match(colnames(modTOM), annot$GeneID)]
rownames(modTOM) <- annot$gene_name[match(rownames(modTOM), annot$GeneID)]

# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])




g <- graph.adjacency(
  modTOM,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
?graph_from_adjacency_matrix


modTOM[modTOM >0.1]

genes

# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "blue"#module_info$moduleColors[which(genes %in% as.character(V(g)$name))]




V(g)$name <- annot$gene_name[match(V(g)$name, annot$GeneID)]
 
# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(datExpr[,genes], 1, mean)) + 1.0) * 3

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 1.0

#Convert the graph adjacency object into a minimum spanning tree based on Prims algorithm
mst <- mst(g, weights=edgeweights)

plot(
  mst,
  layout=layout.fruchterman.reingold(g),
  edge.curved=TRUE,
  vertex.label.dist=0.4,
  vertex.label.color="black",
  vertex.size=vSizes,
  asp=FALSE,
  vertex.label.cex=0.5,
  edge.arrow.mode=0.4,
  main="My first graph",
  edge.width=edgeweights,
  vertex.label.degree	=0.1
  )
?plot

modTOM
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
