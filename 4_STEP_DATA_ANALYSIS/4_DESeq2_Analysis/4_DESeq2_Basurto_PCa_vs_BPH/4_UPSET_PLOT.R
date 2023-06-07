################################################################################
#####                             UPSET PLOT                              ######
################################################################################

#| Upset plot for the results of the DEGs with different designs.

#| Packages to install:
#| install_github("jokergoo/ComplexHeatmap")
#| install.packages("UpSetR")

suppressMessages(library(ggplot2))
suppressMessages(library(UpSetR))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(circlize))
suppressMessages(library(devtools))
suppressMessages(library(viridis))

################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/"
setwd(dir.proj)
################################################################################


################################################################################
################################## DATA ########################################
################################################################################
set1 <- read.table("Results/Tables/resDESeq2_dds_DV200_Diagnostico_ContrastPCa_vs_BPH.txt", sep = "\t")
set2 <- read.table("Results/Tables/resDESeq2_dds_DV200_Edad_NameEdad.txt", sep = "\t")
set3 <- read.table("Results/Tables/resDESeq2_dds_DV200_Diagnostico_DiagnosticoEdad_ContrastPCa_vs_BPH.txt", sep = "\t")
set4 <- read.table("Results/Tables/resDESeq2_dds_DV200_Diagnostico_DiagnosticoEdad_NameDiagnosticoPCa.Edad.txt", sep = "\t")
set5 <- read.table("Results/Tables/resDESeq2_dds_DV200_Diagnostico_DiagnosticoEdad_NameDiagnosticoBPH.Edad_Zscore.txt", sep = "\t")
set6 <- read.table("Results/Tables/resDESeq2_dds_DV200_Group_ContrastLA_vs_LI.txt", sep = "\t")
set7 <- read.table("Results/Tables/resDESeq2_dds_DV200_DiagnosticoEdad_NameDiagnosticoBPH.Edad_Zscore.txt", sep = "\t")
set8 <- read.table("Results/Tables/resDESeq2_dds_DV200_DiagnosticoEdad_NameDiagnosticoPCa.Edad.txt", sep = "\t")
set9 <- read.table("Results/Tables/resDESeq2_dds_DV200_GroupEdad_ContrastLA_vs_LI.txt", sep = "\t")
set10 <- read.table("Results/Tables/resDESeq2_dds_DV200_GroupEdad_NameGroupLA.Edad.txt", sep = "\t")
set11 <- read.table("Results/Tables/resDESeq2_dds_DV200_GroupEdad_NameGroupLI.Edad.txt", sep = "\t")

set1  <- rownames(set1)[set1$DEGs == "Yes"]
set2  <- rownames(set2)[set2$DEGs == "Yes"]
set3  <- rownames(set3)[set3$DEGs == "Yes"]
set4  <- rownames(set4)[set4$DEGs == "Yes"]
set5  <- rownames(set5)[set5$DEGs == "Yes"]
set6  <- rownames(set6)[set6$DEGs == "Yes"]
set7  <- rownames(set7)[set7$DEGs == "Yes"]
set8  <- rownames(set8)[set8$DEGs == "Yes"]
set9  <- rownames(set9)[set9$DEGs == "Yes"]
set10 <- rownames(set10)[set10$DEGs == "Yes"]
set11 <- rownames(set11)[set11$DEGs == "Yes"]

#| Saving the DEGs in a list
lt <- list(set2, set4, set5, set6, set7, set8, set9, set10, set11)

#| Converting the list to a binary matrix
mt <- list_to_matrix(lt)
colnames(mt) <- c("~ DV200 + Edad\n
                  Name: Edad", 
                  "~ DV200 + Diagnostico + Diagnostico:Edad\n
                  Name: DiagnosticoPCa.Edad",
                  "~ DV200 + Diagnostico + Diagnostico:Edad\n
                  Name: DiagnosticoBPH.Edad",
                  "~ DV200 + Group\n
                  Contrast: LA vs LI",
                  "~ DV200 + Diagnpstoco:Edad\n
                  Name: DiagnosticoBPH.Edad",
                  "~ DV200 + Diagnpstoco:Edad\n
                  Name: DiagnosticoPCa.Edad",
                  "~ DV200 + Group:Edad\n
                  Contrast: LA_vs_LI",
                  "~ DV200 + Group:Edad\n
                  Name: GroupLA.Edad",
                  "~ DV200 + Group:Edad\n
                  Name: GroupLI.Edad")

#| Combination matrix to compute the size of the sets and the combination sets
mc <- make_comb_mat(mt)

#| Creating the plot
pdf("Upset_plot.pdf")
UpSet(mc, comb_col = "#440154FF",pt_size = unit(1, "mm"), lwd = 0.5,
      left_annotation = upset_left_annotation(mc, width = unit(4,"cm"), annotation_name_gp = gpar(fontsize=10, fontfamily ="serif", fontface ="bold")), column_names_gp = gpar(fontsize = 1, fontface ="bold"), 
      row_names_gp = gpar(fontsize = 5,  fontfamily ="serif",fontface ="bold"),
      column_dend_gp = gpar(fontsize = 6),top_annotation = upset_top_annotation(mc, annotation_name_gp = gpar(fontsize=10, fontfamily ="serif",fontface ="bold") ))
dev.off()



# Selected designs

################################################################################
################################## DATA ########################################
###############################################################################
set1 <- read.table("Results/Tables/resDESeq2_dds_DV200_GroupEdad_NameGroupLA.Edad.txt", sep = "\t")
set2 <- read.table("Results/resDESeq2_dds_DV200_GroupDFS.TIME_Zscore_NameGroupLA.DFS.TIME_Zscore.txt", sep = "\t")
set3 <- read.table("Results/Tables/resDESeq2_dds_DV200_Edad_GroupDFS.TIME_Zscore_NameGroupLA.DFS.TIME_Zscor.txt", sep = "\t")

set1  <- rownames(set1)[set1$DEGs == "Yes"]
set2  <- rownames(set2)[set2$DEGs == "Yes"]
set3  <- rownames(set3)[set3$DEGs == "Yes"]

lt <- list(set1, set2, set3)

mt <- list_to_matrix(lt)

colnames(mt) <- c("~ DV200 + Group:Edad\n
                  Name: GroupLA.Edad", 
                  "~ DV200 + Group:DFS.TIME\n
                  Name: GroupLA.DFS.TIME",
                  "~ DV200 + Edad + Group:DFS.TIME\n
                  Name: GroupLA.DFS.TIME")

mc <- make_comb_mat(mt)

pdf("Upset_plot_selected_designs.pdf")
UpSet(mc, comb_col = "#440154FF",pt_size = unit(4, "mm"), lwd = 0.8,
      left_annotation = upset_left_annotation(mc, width = unit(4,"cm"), annotation_name_gp = gpar(fontsize=14, fontfamily ="serif", fontface ="bold")), column_names_gp = gpar(fontsize = 1, fontface ="bold"), 
      row_names_gp = gpar(fontsize = 10,  fontfamily ="serif",fontface ="bold"),
      column_dend_gp = gpar(fontsize = 6),top_annotation = upset_top_annotation(mc, annotation_name_gp = gpar(fontsize=14, fontfamily ="serif",fontface ="bold") ))
dev.off()
