################################### xCELL ######################################

#| Code for inference the score of cell type using xCell.

################################################################################

########################### LIBRARIES AND DATA #################################
library(xCell)
library(ggplot2)

theme_set(theme_classic())

exprMatrix <- read.table("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_counts_phenotype_Normalized_Log2_GENE_NAME_ROWS_xCELL_DATA.txt",header=TRUE,row.names=1, as.is=TRUE)
data_trait <- read.table("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt", sep ="\t", header=T)
###############################################################################

#| Running xCell
data <-xCellAnalysis(exprMatrix)


############################### PCa VS BPH #####################################
h_data <- t(data)
h_data <- as.data.frame(h_data)
h_data$Diagnosis <- data_trait$Diagnostico

data_frame <- data.frame(Diagnosis = data_trait$Diagnostico,
                         StromaScore = as.numeric(data[which(rownames(data) == "StromaScore"),]),
                         ImmuneScore = as.numeric(data[which(rownames(data) == "ImmuneScore"),]),
                         MicroenvironmentScore = as.numeric(data[which(rownames(data) == "MicroenvironmentScore"),]),
                         Epithelial_cells = as.numeric(data[which(rownames(data) == "Epithelial cells"),]),
                         Smooth_muscle = as.numeric(data[which(rownames(data) == "Smooth muscle"),]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),]),
                         Epithelial_cells = as.numeric(data[which(rownames(data) == "Epithelial cells"),]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),]))

ggplot(data_frame, aes(x= Diagnosis, y = StromaScore, fill =Diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Diagnosis")+
  ylab("Stroma Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosiss")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Stroma_Score_Results_Diagnosis_PCa_vs_BPH.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= Diagnosis, y = ImmuneScore, fill =Diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Diagnosis")+
  ylab("Immune Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosis")+ 
  ylim(0,0.025)
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Immune_Score_Results_Diagnosis_PCa_vs_BPH.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= Diagnosis, y = MicroenvironmentScore, fill =Diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Diagnosis")+
  ylab("Microenvironment Score") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosis")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Microenvironment_score_Results_Diagnosis_PCa_vs_BPH.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= Diagnosis, y = Fibroblasts, fill =Diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),  legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Diagnosis")+
  ylab("Fibroblasts") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosis")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Fibroblasts_Results_Diagnosis_PCa_vs_BPH.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= Diagnosis, y = Epithelial_cells, fill =Diagnosis)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Diagnosis")+
  ylab("Epithelial") +
  scale_fill_manual(values = c("lightseagreen", "midnightblue"))+
  labs(fill = "Diagnosis")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Epithelial_Results_Diagnosis_PCa_vs_BPH.pdf", heigh= 4.5, width =4.8)

#| Which one has the most significant difference?

p_value <-  c()
difference <- c()

for (i in 1:(dim(h_data)[2] - 1) ){
  
  wilcox_test <- wilcox.test(h_data[which(h_data$Diagnosis == "PCa"), i], h_data[which(h_data$Diagnosis == "BPH"), i])
  p_value <- c(p_value, wilcox_test$p.value )
  difference <- c(difference, mean(h_data[which(h_data$Diagnosis == "PCa"), i])- mean(h_data[which(h_data$Diagnosis == "BPH"), i]))
}

h_k <- data.frame( datasets = colnames(h_data[,-c(68)]),
                   p_value = p_value,
                   difference = difference)


h_k <- h_k[order(h_k$p_value, decreasing = T),]
h_k$log <- -log(h_k$p_value,2)
h_k$difference[which(h_k$difference >= 0)] <- "Positive"
h_k$difference[which(h_k$difference < 0)] <- "Negative"

ggplot(h_k[56:67,], aes(x = log, y = reorder(datasets,log), fill =difference)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("black", "azure3")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Datasets") +
  xlab("-log2(p-value)")+
  labs(fill ="PCa - BPH")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/xCell_Results_Diagnosis_PCa_Vs_BPH_wilxon_test.pdf", heigh= 5.5, width =7.5)

library(viridis)
viridis(5)
################################################################################



########################### PTEN LOSS VS INTACT ################################
h_data <- t(data)
h_data <- as.data.frame(h_data)
h_data <- h_data[data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")],]
h_data$PTEN_Status <- data_trait$H_score_cut_0[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]


data_frame <- data.frame(PTEN_status = data_trait$H_score_cut_0[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")],
                         StromaScore = as.numeric(data[which(rownames(data) == "StromaScore"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         ImmuneScore = as.numeric(data[which(rownames(data) == "ImmuneScore"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         MicroenvironmentScore = as.numeric(data[which(rownames(data) == "MicroenvironmentScore"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Epithelial_cells = as.numeric(data[which(rownames(data) == "Epithelial cells"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Smooth_muscle = as.numeric(data[which(rownames(data) == "Smooth muscle"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]),
                         Fibroblasts = as.numeric(data[which(rownames(data) == "Fibroblasts"),data_trait$AC.basurto[which(data_trait$H_score_cut_0 == "PTEN loss" | data_trait$H_score_cut_0 == "PTEN intact")]]))



ggplot(data_frame, aes(x= PTEN_status, y = StromaScore, fill =PTEN_status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5, method="wilcox")+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'),legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Stroma Score") +
  scale_fill_manual(values =c("midnightblue", "orange"))+
  labs(fill = "PTEN status")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Stroma_Score_Results_PTEN_loss_vs_intact.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= PTEN_status, y = ImmuneScore, fill =PTEN_status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Immune Score") +
  scale_fill_manual(values =c("midnightblue", "orange"))+
  labs(fill = "PTEN status")+ 
  ylim(0,0.025)
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Immune_Score_Results_PTEN_loss_vs_intact.pdf", heigh=4.5, width =4.8)

ggplot(data_frame, aes(x= PTEN_status, y = MicroenvironmentScore, fill =PTEN_status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Microenvironment Score") +
  scale_fill_manual(values =c("midnightblue", "orange"))+
  labs(fill = "PTEN status")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Microenvironment_score_Results_PTEN_loss_vs_intact.pdf", heigh= 4.5, width =4.8)

ggplot(data_frame, aes(x= PTEN_status, y = Fibroblasts, fill =PTEN_status)) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(label ="p.format", vjust=0.5,family="serif", size =5)+
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN Status")+
  ylab("Fibroblasts") +
  scale_fill_manual(values = c("midnightblue", "orange"))+
  labs(fill = "PTEN status")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Fibroblasts_Results_PTEN_loss_vs_intact.pdf", heigh= 4.5, width =4.8)

#| Which one has the most significant difference?

p_value <-  c()
difference <- c()

for (i in 1:(dim(h_data)[2] - 1) ){

  wilcox_test <- wilcox.test(h_data[which(h_data$PTEN_Status == "PTEN loss"), i], h_data[which(h_data$PTEN_Status == "PTEN intact"), i])
  p_value <- c(p_value, wilcox_test$p.value )
  difference <- c(difference, mean(h_data[which(h_data$PTEN_Status == "PTEN loss"), i]) - mean(h_data[which(h_data$PTEN_Status == "PTEN intact"), i]))
  
}

h_k <- data.frame( datasets = colnames(h_data[,-c(68)]),
                   p_value =p_value,
                   difference =difference)


h_k <- h_k[order(h_k$p_value, decreasing = T),]
h_k$log <- -log(h_k$p_value,2)
h_k$difference[which(h_k$difference >= 0)] <- "Positive"
h_k$difference[which(h_k$difference < 0)] <- "Negative"

ggplot(h_k[56:67,], aes(x = log, y = reorder(datasets,log), fill =difference)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("midnightblue", "darkgoldenrod3")) +
  scale_fill_manual(values = c("black", "azure3")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Datasets") +
  xlab("-log2(p-value)") +
  labs(fill = "PTEN loss - intact")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/xCell_Results_PTEN_loss_vs_intact_wilxon_test.pdf", heigh= 5.5, width =7.5)
###############################################################################