################################################################################
###############  COMPARING BASURTO WITH OTHER DATASETS #########################
################################################################################

#| Comparing PTEN expression over different datasets containing Normal prostate
#| and Primary tumor. 

#| Datasets:
#|  *) Basurto (N= 40, PT= 198)
#|  *) Taylor (N=29, PT=131)
#|  *) Lapointe (N=9, PT=13)
#|  *) Tomlins (N= 23, PT= 32)
#|  *) Varambally (N= 6, PT=7)

#| Datasets were extracted from CANCERTOOL http://genomics.cicbiogune.es/CANCERTOOL/
################################################################################


################################ LIBRARIES #####################################
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(biomaRt))
suppressMessages(library(DBI))
suppressMessages(library(RMySQL))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
################################################################################



################################################################################
workingDir <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/"
setwd(workingDir)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("W:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Setting themes for plotting
theme_set(theme_classic())

#| Connecting to MySql database
con <- dbConnect(RMySQL::MySQL(), user='geneticanalyses', password = 'geneticanalyses', host = 'binf-web.cicbiogune.int',port=13306, dbname='geneticanalyses')
dbListTables(con)

#| For plots
theme_set(theme_classic()) 
################################################################################



################################### BASURTO #################################### 
#| Location of the files
info.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Normalized_Log2.txt"

#| Readtable Full counts. NOTE: it is not necessary to do any normalization 
data_expression_basurto <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Readtable Sample information
data_phenotype_basurto <- read.table(info.file, sep ="\t")

#| PTEN Basurto
#PTEN_basurto <- as.numeric(data_expression_basurto[which(rownames(data_expression_basurto) == "ENSG00000171862"),data_phenotype_basurto$AC.basurto[which(data_phenotype_basurto$Diagnostico=="PCa")]])
PTEN_basurto <- as.numeric(data_expression_basurto[which(rownames(data_expression_basurto) == "ENSG00000171862"),])

#| Adding PTEN infor to sample data
data_phenotype_basurto$PTEN_mrna <- PTEN_basurto
################################################################################


#################################### TAYLOR ####################################
#| Extracting the data
data_phenotype_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorphenotype")
data_expression_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorexpression")
data_annotation_taylor <- dbGetQuery(conn = con, statement = "SELECT * from taylorannotation")

#| Changing colnames of data_annotation_taylor
colnames(data_annotation_taylor) <- c("row_names", "gene_name", "refseq_mrna")

#| Assigning the gene name to rownames
rownames(data_expression_taylor)<- data_expression_taylor$Gene

#| PTEN
data_annotation_taylor$refseq_mrna[which(data_annotation_taylor$gene_name=="PTEN")] #| "NM_000314"

#| PTEN expression
#PTEN_taylor <- as.numeric(data_expression_taylor[rownames(data_expression_taylor) == "NM_000314",data_phenotype_taylor$Sample[ which(data_phenotype_taylor$tumortype == "Primary Tumor")]])
samples <- data_phenotype_taylor$Sample[ which(data_phenotype_taylor$tumortype == "Primary Tumor"  | data_phenotype_taylor$tumortype == "Normal")]
PTEN_taylor <- as.numeric(data_expression_taylor[rownames(data_expression_taylor) == "NM_000314", samples])

data_phenotype_taylor <- data_phenotype_taylor[which(data_phenotype_taylor$Sample %in% samples),]
data_phenotype_taylor$PTEN_mrna <- PTEN_taylor
################################################################################


#################################### GRASSO ####################################
#| Extracting the data
data_phenotype_grasso <- dbGetQuery(conn = con, statement = "SELECT * from grassophenotype")
data_expression_grasso <- dbGetQuery(conn = con, statement = "SELECT * from grassoexpression")

rownames(data_expression_grasso) <- data_expression_grasso$Gene

data_phenotype_grasso <- data_phenotype_grasso[which(data_phenotype_grasso$tumortype == "Normal" | data_phenotype_grasso$tumortype == "Primary Tumor"),]

PTEN_grasso <- as.numeric(data_expression_grasso[which(rownames(data_expression_grasso) == "PTEN"),data_phenotype_grasso$Sample])

data_phenotype_grasso$PTEN_mrna <- PTEN_grasso
################################################################################


#################################### LAPOINTE ####################################
#| Extracting the data
data_phenotype_lapointe <- dbGetQuery(conn = con, statement = "SELECT * from lapointephenotype")
data_expression_lapointe <- dbGetQuery(conn = con, statement = "SELECT * from lapointeexpression")

rownames(data_expression_lapointe) <- data_expression_lapointe$Gene

data_phenotype_lapointe <- data_phenotype_lapointe[which(data_phenotype_lapointe$tumortype == "Normal" | data_phenotype_lapointe$tumortype == "Primary Tumor"),]

PTEN_lapointe <- as.numeric(data_expression_lapointe[which(rownames(data_expression_lapointe) == "PTEN"),data_phenotype_lapointe$Sample])

colnames(data_expression_lapointe[,-c(1)]) == data_phenotype_lapointe$Sample

data_phenotype_lapointe$PTEN_mrna <- PTEN_lapointe
################################################################################


#################################### TOMLINS ####################################
#| Extracting the data
data_phenotype_tomlins <- dbGetQuery(conn = con, statement = "SELECT * from tomlinsphenotype")
data_expression_tomlins <- dbGetQuery(conn = con, statement = "SELECT * from tomlinsexpression")

rownames(data_expression_tomlins) <- data_expression_tomlins$Gene

colnames(data_expression_tomlins[,-c(1)]) == data_phenotype_tomlins$Sample

data_phenotype_tomlins <- data_phenotype_tomlins[which(data_phenotype_tomlins$tumortype == "Normal" | data_phenotype_tomlins$tumortype == "Primary Tumor"),]

PTEN_tomlins <- as.numeric(data_expression_tomlins[which(rownames(data_expression_tomlins) == "PTEN"),data_phenotype_tomlins$Sample])

data_phenotype_tomlins$PTEN_mrna <- PTEN_tomlins
################################################################################


#################################### VARAMBALLY ################################
#| Extracting the data
data_phenotype_varambally <- dbGetQuery(conn = con, statement = "SELECT * from varamballyphenotype")
data_expression_varambally <- dbGetQuery(conn = con, statement = "SELECT * from varamballyexpression")

rownames(data_expression_varambally) <- data_expression_varambally$Gene

colnames(data_expression_varambally[,-c(1)]) == data_phenotype_varambally$Sample

data_phenotype_varambally <- data_phenotype_varambally[which(data_phenotype_varambally$tumortype == "Normal" | data_phenotype_varambally$tumortype == "Primary Tumor"),]

PTEN_varambally <- as.numeric(data_expression_varambally[which(rownames(data_expression_varambally) == "PTEN"),data_phenotype_varambally$Sample])

data_phenotype_varambally$PTEN_mrna <- PTEN_varambally
################################################################################


######################## TCGA Fire Horse Legacy data ###########################
data_phenotype_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostatephenotype")
data_expression_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostateexpression")

#| Fixing gene name
data_expression_tcga$Gene <- gsub(".{1}$","",data_expression_tcga$Gene)
rownames(data_expression_tcga) <- data_expression_tcga$Gene

#| Deleting Gene column
data_expression_tcga <- data_expression_tcga[,-c(1)]

#| PTEN TCGA
PTEN_tcga <- as.numeric(data_expression_tcga[which(rownames(data_expression_tcga)=="PTEN"),])
################################################################################


################################# Glinsky ######################################
data_phenotype_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyphenotype")
data_expression_glinsky <- dbGetQuery(conn = con, statement = "SELECT * from glinskyexpression")

#| Giving gene name to row names
rownames(data_expression_glinsky) <- data_expression_glinsky$Gene

#| Deleting column "Gene" form gene expression data
data_expression_glinsky <- data_expression_glinsky[,-c(1)]

#| PTEN glinksy
PTEN_glinsky <- as.numeric(data_expression_glinsky[which(rownames(data_expression_glinsky) == "PTEN"),])
################################################################################


############ TESTING PTEN EXPRESSION ON NORMAL VS PRIMARY TUMOR ################
data_phenotype_basurto$status_diagnosis <- NA
data_phenotype_basurto$status_diagnosis[which(data_phenotype_basurto$Diagnostico == "PCa")] <- "Primary Tumor"
data_phenotype_basurto$status_diagnosis[which(data_phenotype_basurto$Diagnostico == "BPH")] <- "Normal"

data <- list(Basurto = data.frame( PTEN_mrna = data_phenotype_basurto$PTEN_mrna, tumortype =data_phenotype_basurto$status_diagnosis),
             Grasso = data.frame(PTEN_mrna= data_phenotype_grasso$PTEN_mrna, tumortype=data_phenotype_grasso$tumortype),
             Lapointe = data.frame(PTEN_mrna= data_phenotype_lapointe$PTEN_mrna, tumortype=data_phenotype_lapointe$tumortype),
             Taylor = data.frame(PTEN_mrna= data_phenotype_taylor$PTEN_mrna, tumortype=data_phenotype_taylor$tumortype),
             Tomlins = data.frame(PTEN_mrna=data_phenotype_tomlins$PTEN_mrna, tumortype=data_phenotype_tomlins$tumortype),
             Varambally = data.frame(PTEN_mrna=data_phenotype_varambally$PTEN_mrna, tumortype=data_phenotype_varambally$tumortype))

data <- data.frame(PTEN_mrna = c(data_phenotype_basurto$PTEN_mrna, data_phenotype_grasso$PTEN_mrna, data_phenotype_lapointe$PTEN_mrna, data_phenotype_taylor$PTEN_mrna,data_phenotype_tomlins$PTEN_mrna,data_phenotype_varambally$PTEN_mrna),
                   tumortype = c(data_phenotype_basurto$status_diagnosis, data_phenotype_grasso$tumortype, data_phenotype_lapointe$tumortype, data_phenotype_taylor$tumortype,data_phenotype_tomlins$tumortype,data_phenotype_varambally$tumortype),
                   dataset = c(rep("Basurto", times=length(data_phenotype_basurto$PTEN_mrna)), rep("Grasso", times=length(data_phenotype_grasso$PTEN_mrna)),rep("Lapointe", times=length(data_phenotype_lapointe$PTEN_mrna)),rep("Taylor", times=length(data_phenotype_taylor$PTEN_mrna)),rep("Tomlins", times=length(data_phenotype_tomlins$PTEN_mrna)),rep("Varambally", times=length(data_phenotype_varambally$PTEN_mrna)))
)

data <- data[which(!is.na(data$PTEN_mrna)),]
data$PTEN_mrna <- (data$PTEN_mrna -mean(data$PTEN_mrna))/sd(data$PTEN_mrna)

ggplot(data, aes(x = dataset, y = PTEN_mrna, fill =tumortype)) +
  geom_boxplot(outlier.shape=NA) +
  theme(text=element_text(size=11,  family="serif"), axis.text.x = element_text(angle = 20, v=0.8,size=11 ),axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
  scale_fill_manual(values =c("aquamarine3", "midnightblue")) +
  ylab("PTEN expression scaled by z scores") +
  xlab("Datasets") +
  stat_compare_means(label ="p.format",family="serif", size =2.9, label.y= c(1.1,-1.2,-1.2,0.8,-1.2, 0.8)) +
  labs(fill="Tumor Type")
ggsave("Results/Box_plots/PTEN_expression_Basurto_Grasso_Lapinte_Taylor_Tomlins_Varambally_Nomal_vs_Primary_tumor.pdf", height = 3.7, width = 5)
################################################################################




######################### EXPRESSION BY DATA SET ###############################

#| BASURTO
ggplot(data_phenotype_basurto, aes(x= status_diagnosis, y = PTEN_mrna, fill =status_diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white",dotsize =0.5, alpha =0.8) +
  stat_compare_means(label ="p.format",vjust=-0.5,hjust =-0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Prostate tissue")+
  ylab("PTEN mRNA (Normalized log2 counts)") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("Results/Box_plots/TumorType_Basurto.pdf", heigh=4.5, width = 5)

#| TAYLOR
ggplot(data_phenotype_taylor, aes(x= tumortype, y = PTEN_mrna, fill =tumortype)) +
  geom_boxplot(outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white",dotsize =0.5, alpha =0.8) +
  stat_compare_means(label ="p.format", vjust=-0.5,hjust =-0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Prostate tissue")+
  ylab("PTEN mRNA (Normalized log2 counts)") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("Results/Box_plots/TumorType_taylor.pdf", heigh=4.5, width = 5)

#| GRASSO
ggplot(data_phenotype_grasso, aes(x= tumortype, y = PTEN_mrna, fill =tumortype)) +
  geom_boxplot(outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white", dotsize =0.5, alpha =0.8) +
  stat_compare_means(label ="p.format", vjust=-0.5,hjust =-0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Prostate tissue")+
  ylab("PTEN mRNA (Normalized log2 counts)") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("Results/Box_plots/TumorType_grasso.pdf", heigh=4.5, width = 5)

#| LAPOINTE
ggplot(data_phenotype_lapointe, aes(x= tumortype, y = PTEN_mrna, fill =tumortype)) +
  geom_boxplot(outlier.colour = NA, ) +
  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white",dotsize =0.5, alpha =0.8) +
  stat_compare_means(label ="p.format", vjust=-0.5,hjust =-0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Prostate tissue")+
  ylab("PTEN mRNA (Normalized log2 counts)") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("Results/Box_plots/TumorType_lapointe.pdf", heigh=4.5, width = 5)

#| VARAMBALLY
ggplot(data_phenotype_varambally, aes(x= tumortype, y = PTEN_mrna, fill =tumortype)) +
  geom_boxplot(outlier.colour = NA, ) +
  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white",dotsize =0.5, alpha =0.8) +
  stat_compare_means(label ="p.format", vjust=-0.5,hjust =-0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Prostate tissue")+
  ylab("PTEN mRNA (Normalized log2 counts)") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("Results/Box_plots/TumorType_varambally.pdf", heigh=4.5, width = 5)
################################################################################

