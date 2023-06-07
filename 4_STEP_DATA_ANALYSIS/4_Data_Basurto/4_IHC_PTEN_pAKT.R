####################  PTEN AND pAKT ASSESSMENT BY IHC  #########################

#| The way PTEN was assessed in the Basurto's cohort was by IHC. The patholist measured
#| different intensities for the appearance of PTEN in the tumor tissue of 198 
#| patients. Here I have processed this data and I have added it to the sample_info 
#| table

#| Moreover, given that the assessment of PTEN represents a continues value (H-score)
#| we explore a threshold of this value, to binary classify PTEN as loss vs intact.
#| How the H-score was computed? see: http://e-immunohistochemistry.info/web/H-score.htm

#| The assessment of pAKT was not very reliable

#| NOTE: For the assessment of pAKT, values like %Tumor, control in the stroma 
#| (neg or pos) are the same to the ones for PTEN. The explanation of why this is 
#| the same is because the assessment was done in adjacent laminal tissues (almost 
#| the same region)
#| In the code above I have saved the values of the assessment of pAKT as vector of
#| length 560, to then bind into the dataframe containing the PTEN assessment along 
#| with the other parameters in common. 

################################################################################


###############################  LIBRARIES  ####################################
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(corrr))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
################################################################################


############################## DATA DIRECTORIES ################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted_processing_and_filtering.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_counts_phenotype_Normalized_Log2.txt"
foxo.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TTGTTT_V_FOXO4_01.txt"
foxo_pathway.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/FOXO_Pathway.txt"
pi3k_akt_mtor.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt"
pten_loss.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/Signature_PTEN_loss_Paper_Imada_et_al_BMC_cancer.xlsx"
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################



############################### READING DATA ###################################
#| Counts data (normalized)
counts_data <- read.table(counts.file, sep ="\t", header=T)

#| Sample information
sample_info <- read.table(info.file, sep ="\t", header=T)
################################################################################



################################ pAKT ASSESSMENT ###############################
Assessment.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TMA assessment PTEN and pAKT.xlsx"
pAKT_TMA.sheet <- "TMA pAKT"

#| Reading Matrix TMA
pAKT_TMA <- read_excel(Assessment.file, sheet = pAKT_TMA.sheet)

#| Finding TMA IDs of the pAKT
TMA_ID_pAKT <- pAKT_TMA$`1`[grepl("TMA", pAKT_TMA$`1`)]
pAKT_vector <- c()

#| Over all the data
for (i in 1:length(TMA_ID_pAKT)){
  
  if (i < length(TMA_ID_pAKT)){
    
    df_pAKT <- t(pAKT_TMA[(which(pAKT_TMA$`1` == TMA_ID_pAKT[i])+1):(which(pAKT_TMA$`1` == TMA_ID_pAKT[i+1])-1),])
    
  } else {
    
    df_pAKT <- t(pAKT_TMA[(which(pAKT_TMA$`1` == TMA_ID_pAKT[i])+1):(dim(pAKT_TMA)[1]+3),])
    
  }
  
  for (j in 0:(dim(df_pAKT)[2]-1)){
    value_pAKT <- df_pAKT[1:dim(df_pAKT)[1],j+1]
    pAKT_vector <- c(pAKT_vector, as.character(value_pAKT))
    
    
  }
  
}
################################################################################




################################ PTEN ASSESSMENT ###############################
Assessment.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TMA assessment PTEN and pAKT.xlsx"
TMA_matrix_ID <-"Matriz TMA"
PTEN_Assessment.sheet <-"PTEN"

#| Reading Matrix TMA
data_TMA_matrix_ID <- read_excel(Assessment.file, sheet = TMA_matrix_ID)
data_TMA_matrix_ID <- data_TMA_matrix_ID[, 1:8]

#| Reading Assessment data
data_PTEN_Assessment <- read_excel(Assessment.file, sheet =PTEN_Assessment.sheet)

#| Renaming first column
colnames(data_TMA_matrix_ID)[1] <- "TMA ID"

#| Finding TMA IDs
TMA_ID <- data_TMA_matrix_ID$`TMA ID`[grepl("TMA", data_TMA_matrix_ID$`TMA ID`)]
AC_ID <- c()

#| Over all the data
for (i in 1:length(TMA_ID)){
  
  if (i < length(TMA_ID)){
    df <- data_TMA_matrix_ID[(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i])+1):(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i+1])-3),]
  } else {
    df <- data_TMA_matrix_ID[(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i])+1):dim(data_TMA_matrix_ID)[1],]
  }
  
  for (j in 0:(dim(df)[2]-1)){
    value <- df[1:dim(df)[1],dim(df)[2]-j]
    AC_ID <- c(AC_ID, value[[1]])
  }
  
}

#| Adding a column with AC ID information to the assessment data
data_PTEN_Assessment$`AC TMA` <- AC_ID

#| Changing AC ID to AC Basurto's ID
data_PTEN_Assessment$`AC basurto` <- gsub("AC-","AC",data_PTEN_Assessment$`AC TMA`)

#| H-score = [(0 x % negative cells) + (1 x %weakly positive cells) + (2 x %moderately positive cells) + (3 x %strongly positive cells)] 
data_PTEN_Assessment$H_score <- (data_PTEN_Assessment$`% negativo`)*0 + (data_PTEN_Assessment$`% débil`)*1 + (data_PTEN_Assessment$`% moderado`)*2 + (data_PTEN_Assessment$`% intenso`)*3
data_PTEN_Assessment$`% tumor numeric` <- as.numeric(data_PTEN_Assessment$`% tumor`)

#| Assigning NA to those value where there was tumor but negative for PTEN in the stroma
data_PTEN_Assessment$H_score[which(data_PTEN_Assessment$estroma == "neg")] <- NA

#| Adding pAKT assessement to the data frame of PTEN
data_PTEN_Assessment$pAKT <- pAKT_vector

#| Aggregating dataframes
H_score_dataframe <- aggregate( H_score ~ `AC basurto`, data_PTEN_Assessment, mean )
Percentage_tumor <- aggregate( `% tumor numeric` ~ `AC basurto`, data_PTEN_Assessment, mean )

#| Merging dataframes
sample_info_extracted <- merge(sample_info_extracted, H_score_dataframe, by.y="AC basurto", all.x =TRUE)
sample_info_extracted <- merge(sample_info_extracted, Percentage_tumor, by.y="AC basurto", all.x =TRUE)

rownames(sample_info_extracted) <- sample_info_extracted$`AC basurto`

#| Dealing with results of pAKT
pAKT_pos <- data_PTEN_Assessment[which((data_PTEN_Assessment$`AC basurto` %in% sample_info_extracted$`AC basurto`[which(sample_info_extracted$Diagnostico == "PCa")] & data_PTEN_Assessment$estroma == "pos" & data_PTEN_Assessment$pAKT=="pos" )),]
pAKT_pos_ID <- unique(pAKT_pos$`AC basurto`)

#| Adding to the dataframe those values where pAKT was positive
sample_info_extracted$pAKT_positive <- "No pAKT"
sample_info_extracted$pAKT_positive[which(sample_info_extracted$`AC basurto` %in% pAKT_pos_ID)] <- "pos"
################################################################################



####################  RESULTS PTEN AND pAKT WHOLE TMA  #########################
#| Percentage of neg, pos, and NA in the estroma WHOLE TMA
stroma <- data_PTEN_Assessment$estroma[which(data_PTEN_Assessment$`AC basurto` %in% 
                                               sample_info_extracted$`AC basurto`[which(sample_info_extracted$Diagnostico == "PCa")])]

data_stroma <- data.frame(
  group=c("% Positive ","% Negative ","% NA ") , 
  value=c(length(stroma[which(stroma == "pos")])/length(stroma),
          length(stroma[which(stroma == "neg")])/length(stroma),
          length(stroma[which(is.na(stroma))])/length(stroma)),
  value2=c("Stroma", "Stroma", "Stroma")
)

pdf("Results/Bar_plots/Quality_control_stroma_whole_TMA_pAKT_PTEN_Assessment_PCA_samples.pdf")
ggplot(data_stroma, aes(x =value2, y = value, fill = group)) + geom_bar(stat = "identity") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Quality control in the stroma from the IHC of \n PTEN and pAKT (Considering the whole TMA and PCa samples)") +
  ylab("Percentage (%) of PCa tissues of the whole TMA") +
  scale_fill_manual(values = viridis(3)) +
  labs(fill="Percentage") +
  ylim(0,1)
dev.off()

pos_pie_chart <- cumsum(data_stroma$value) - 0.5*data_stroma$value
data_stroma$value <- round(data_stroma$value,2)

pdf("Results/Pie_chart/Quality_control_stroma_whole_TMA_pAKT_PTEN_Assessment_PCA_samples_circle.pdf")
ggplot(data_stroma, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + geom_text(aes(y = pos_pie_chart, label = value), color = "black", size =4, family ="serif") +
  labs(fill="Percentage") + theme_void() +scale_fill_manual(values = viridis(3)) +
  ggtitle("Quality control in the stroma from the IHC of \n PTEN and pAKT (Considering the whole TMA and PCa samples)") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 1, face ="bold")) 
dev.off()  

#| Histogram percentage tumor, debil, moderado e intenso PTEN WHOLE TMA and stroma positive
PTEN_tumor <- data.frame(tumor = data_PTEN_Assessment$`% tumor numeric`[which(data_PTEN_Assessment$`AC basurto` %in% sample_info_extracted$`AC basurto`[which(sample_info_extracted$Diagnostico == "PCa" & data_PTEN_Assessment$estroma =="pos")])])

pdf("Results/Histograms/Density_Percentage_tumor_tissue_PTEN_and_pAKT_whole_TMA_PCa_samples.pdf")
ggplot(PTEN_tumor, aes(x=tumor)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7)+theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=12.5, hjust = 0.5, face ="bold")) +
  xlab("% Tumor in the whole TMA") + ylab("Density")+ ggtitle(paste("Density histogram of the % of the tumor tissue for PTEN and pAKT assessment\n ( ",dim(PTEN_tumor)[1], " possitions in the TMA (pos stroma))",sep ="" )) + scale_x_continuous(n.breaks=10)
dev.off()

#| Percentage of neg, pos, tumor, debil, moderado e intenso PTEN WHOLE TMA MEAN
PTEN_data <- data_PTEN_Assessment[which(data_PTEN_Assessment$`AC basurto` %in% sample_info_extracted$`AC basurto`[which(sample_info_extracted$Diagnostico == "PCa" & data_PTEN_Assessment$estroma =="pos")]),]

data_PTEN_neg_debil_mod_i <- data.frame(
  group=c("% Negative ","% Weak ","% Moderate ", "% Intense") , 
  value=c(mean(PTEN_data$`% negativo`[which(!is.na(PTEN_data$`% negativo`))]),
          mean(PTEN_data$`% débil`[which(!is.na(PTEN_data$`% débil`))]),
          mean(PTEN_data$`% moderado`[which(!is.na(PTEN_data$`% moderado`))]),
          mean(PTEN_data$`% intenso`[which(!is.na(PTEN_data$`% intenso`))]))
)

pdf("Results/Bar_plots/Quality_control_PTEN_ASSESSMENT_whole_TMA_PCA_samples.pdf")
ggplot(data_PTEN_neg_debil_mod_i, aes(x =group, y = value, fill = group)) + geom_bar(stat = "identity") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle(paste("PTEN assessment by IHC for the whole TMA considering PCa samples\n",
                dim(PTEN_tumor)[1], " possitions in the TMA (pos stroma)", sep ="")) +
  xlab("PTEN assessment parameters") + 
  ylab("Mean %") +
  scale_fill_manual(values = viridis(4)) +
  labs(fill="Percentage") 
dev.off()
################################################################################




######################  RESULTS PTEN AND pAKT BASUTO  ##########################
#| H-score histogram
ggplot(sample_info_extracted, aes(H_score)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill ="blue", linewidth=0.7) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("H-score PTEN") + ylab("Density") +
  ggtitle("H???score density histogram PTEN")
ggsave("Results/Histograms/Density_H-score_PTEN_NA_stroma.pdf", height = 7, width = 7) 


#| % of samples with H-score different from NA and % of samples with pAKT and % of overlaping between non-NA H-score and pos pAKT
data_H_score_PTEN_pos_pAKT <- data.frame(group = c("% H-score \n >=0", "% pos pAKT"),
                                         value = c(length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa" & !is.na(sample_info_extracted$H_score))])/length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa")]),
                                                   length(sample_info_extracted$pAKT_positive[which(sample_info_extracted$Diagnostico == "PCa" & sample_info_extracted$pAKT_positive == "pos")])/length(sample_info_extracted$pAKT_positive[which(sample_info_extracted$Diagnostico == "PCa")])))
pdf("Results/Bar_plots/Percentage_confirmed_PTEN_pAKT_assessment_by_IHC.pdf")
ggplot(data_H_score_PTEN_pos_pAKT, aes(x =group, y = value, fill = group)) + geom_bar(stat = "identity") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  ggtitle("Evaluating PTEN and pAKT assessment over \nthe Basurto cohort (198 PCa samples)")+
  xlab("Assesment parameters") + 
  ylab("Percentage (%) for the PCa sample") +
  scale_fill_manual(values = c("darkslategrey", "darkgoldenrod2")) +
  labs(fill="Percentage") + ylim(0,1)
dev.off()

#| Percentage of PCa samples with H-score 
data_PCa_H_score <- data.frame(group =c("PCa", "PCa"),
                               value =c(length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa" & is.na(sample_info_extracted$H_score))])/length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa")]),
                                        length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa" & !is.na(sample_info_extracted$H_score))])/length(sample_info_extracted$H_score[which(sample_info_extracted$Diagnostico == "PCa")])),
                               value2=c("H-score = NA", "H-score >= 0"))
pdf("Results/Bar_plots/Evaluating_percentage_H-scores_from_PTEN_assessment_by_IHC_Basurto_PCa_samples.pdf")
ggplot(data_PCa_H_score, aes(x =group, y = value, fill = value2)) + geom_bar(stat = "identity") +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(2, 'cm'), plot.title=element_text(size=14.5, hjust = 0.5, face ="bold")) +
  ggtitle("Evaluating percentage of H-scores from PTEN assessment \n by IHC in the Basurto's cohort (198 PCa)")+
  xlab("") + 
  ylab("Percentage (%) of PCa samples") +
  scale_fill_manual(values = c("cyan4", "darkgoldenrod2")) +
  ylim(0,1) +
  labs(fill = "H-score values")
dev.off()

#| PTEN Expression vs H-score (and pAKT activity)
cor <- cor.test(sample_info_extracted$PTEN_Exp_Zscore[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))], method ="spearman", exact=FALSE)
cor$p.value
round(as.numeric(cor$estimate[[1]]),2)

pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_pAKT_Activity.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= PTEN_Exp_log2, color = pAKT_positive)) + geom_point(size=2.5) + 
  geom_smooth(color ="firebrick4") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  xlab("H-score") +
  labs(color ="pAKT activity") +
  scale_color_manual(values = c("midnightblue", "darkgoldenrod3")) +
  annotate(geom="text", x=150, y=10.2, label=expression(paste(rho, " = ", "-0.03", ", "," p < 0.75", sep ="" )), color="black", family ="serif", size =6)
dev.off()

#| PTEN Expression vs H-score (and purity)
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_purity.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= PTEN_Exp_log2, color = purity)) + geom_point(size=2.5) + 
  geom_smooth(color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  xlab("H-score") +
  labs(color ="Purity") +
  annotate(geom="text", x=150, y=10.2, label=expression(paste(rho, " = ", "-0.03", ", "," p < 0.75", sep ="" )), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D")
dev.off()

#| Mean expression PI3K-AKT-mTOR vs H-score
cor <- cor.test(sample_info_extracted$PTEN_Exp_Zscore[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(!is.na(sample_info_extracted$H_score))], method ="spearman", exact=FALSE)
cor$p.value
round(as.numeric(cor$estimate[[1]]),2)

ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= `Mean expression PI3K-AKT-mTOR`)) + geom_point(size=2.5) + 
  geom_smooth(color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("Mean expression PI3K-AKT-mTOR (Z score of normalized log2 counts)")+ 
  xlab("H-score") +
  annotate(geom="text", x=50, y=0.5, label=expression(paste(rho, " = ", "-0.1", ", "," p < 0.24", sep ="" )), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D")



#| PTEN assesment vs %tumor 
pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_percentage_tumor.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= PTEN_Exp_log2, color = `% tumor numeric`)) + geom_point(size=2.5) + 
  geom_smooth(color ="black") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  xlab("H-score") +
  labs(color ="mean % Tumor") +
  annotate(geom="text", x=150, y=10.2, label=expression(paste(rho, " = ", "-0.03", ", "," p < 0.75", sep ="" )), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D")
dev.off()

#|  % tumor vs H-score 
pdf("Results/Correlation_plots/Correlation_Percentage_tumor_vs_H_score_pAKT.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= `% tumor numeric`, color = pAKT_positive)) + geom_point(size=2.5) + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("% Tumor vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("% Tumor PTEN Assessment")+ 
  xlab("H-score") +
  labs(color ="mean % Tumor") +
  annotate(geom="text", x=150, y=10.2, label=expression(paste(rho, " = ", "-0.03", ", "," p < 0.75", sep ="" )), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=T, option="D")
dev.off()

#| Mean expression FOXO vs H-score (and pAKT activity)
cor <- cor.test(sample_info_extracted$`Mean expression FOXO`[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))], method ="spearman", exact=FALSE)
cor$p.value
round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Mean_expression_FOXO_vs_H_score_pAKT_activity.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= `Mean expression FOXO`, color =pAKT_positive)) + geom_point(size=2.5) + 
  geom_smooth(color ="firebrick4") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("Mean Expression FOXO signature vs H-score PTEN \nAssessment 147 PCa samples") +  
  ylab("Mean expression FOXO signature")+ 
  xlab("H-score") +
  labs(color ="pAKT activity") +
  scale_color_manual(values = c("midnightblue", "darkgoldenrod3")) +
  annotate(geom="text", x=150, y=7.6, label=expression(paste(rho, " = ", "-0.18", ", "," p < 0.025", sep ="" )), color="black", family ="serif", size =6)
dev.off()


#| PTEN and pAKT activity (for those with H-score noNA)
pdf("Results/Box_plots/PTEN_expression_and_pAKT_activitiry_147_PCa_samples_H_score_greater_than_0.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x=pAKT_positive, y =PTEN_Exp_log2, fill =pAKT_positive)) +
  geom_boxplot() +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("PTEN expression (Normalized log2 counts)")+ 
  ggtitle("Behavior of PTEN expression based on pAKT activity detected \n by IHC on PCa samples (147 with H-score >= 0)")+
  xlab("pAKT assessment by IHC") +
  scale_fill_manual(values =c("darkslateblue", "darkolivegreen3")) + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  labs(fill = "pAKT activity")
dev.off()

#| Mean expression and pAKT activity (for those with H-score noNA)
pdf("Results/Box_plots/Mean_expression_FOXO_and_pAKT_activitiry_147_PCa_samples_H_score_greater_than_0.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x=pAKT_positive, y =`Mean expression FOXO`, fill =pAKT_positive)) +
  geom_boxplot() +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) + 
  ylab("Mean expression FOXO signature (Normalized log2 counts)")+ 
  ggtitle("Behavior of Mean expression FOXO signature based on pAKT  \nactivity detected by IHC on PCa samples (147 with H-score >= 0)")+
  xlab("pAKT assessment by IHC") +
  scale_fill_manual(values =c("darkslateblue", "darkolivegreen3")) + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  labs(fill = "pAKT activity")
dev.off()


#| H-score and AKT1 expression
pdf("Results/Correlation_plots/AKT1_expression_and_H_score_147_PCa_samples.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= AKT1_Exp_log2)) + geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 Expression vs H-score PTEN Assessment \n 147 PCa samples") +
  ylab("AKT1 Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=50, y=9, label=paste("r = ", format(round(cor(sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$AKT1_Exp_log2[which(!is.na(sample_info_extracted$H_score))], method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6)
dev.off()

#| H-score and AKT1 expression (purity)
pdf("Results/Correlation_plots/AKT1_expression_and_H_score_147_PCa_samples_purity.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= AKT1_Exp_log2, colour =purity)) + geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 Expression vs H-score PTEN Assessment \n 147 PCa samples") +
  ylab("AKT1 Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=50, y=9, label=paste("r = ", format(round(cor(sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$AKT1_Exp_log2[which(!is.na(sample_info_extracted$H_score))], method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D")
dev.off()

#| H-score and AKT1 expression (only recurrence)
pdf("Results/Correlation_plots/AKT1_expression_and_H_score_Recurrence_36_PCa_samples.pdf")
ggplot(sample, aes(x= H_score, y= AKT1_Exp_log2)) + geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 Expression vs H-score PTEN Assessment \n Only Recurrence (36 PCa samples)") +
  ylab("AKT1 Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=50, y=9, label=paste("r = ", format(round(cor(sample$H_score, sample$AKT1_Exp_log2, method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6)
dev.off()

#| H-score and AKT1 expression (only recurrence) (purity)
pdf("Results/Correlation_plots/AKT1_expression_and_H_score_Recurrence_36_PCa_samples_purity.pdf")
ggplot(sample, aes(x= H_score, y= AKT1_Exp_log2, colour =purity)) + geom_point(size=2.5) + 
  geom_smooth(method = "lm", color ="black") + theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("AKT1 Expression vs H-score PTEN Assessment \n Only Recurrence (36 PCa samples)") +
  ylab("AKT1 Expression (Normalizede log2 counts)")+ 
  annotate(geom="text", x=50, y=9, label=paste("r = ", format(round(cor(sample$H_score, sample$AKT1_Exp_log2, method ="pearson"),2), nsmall =2),sep=""), color="black", family ="serif", size =6)+
  scale_color_viridis(discrete=F, option="D")
dev.off()
################################################################################


################ FINDING H-SCORE THRESHOLD BASED ON SIGNATURES #################
#| Is it possible to create a classification of H-Score based of the pAKT activity?
pdf("Results/Box_plots/H_score_and_pAKT_activity_to_determine_possible_threshold.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x=pAKT_positive, y =H_score, fill =pAKT_positive)) +
  geom_boxplot()+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("H-score")+ 
  ggtitle("Behavior of H-score based on pAKT activity detected by IHC \n on PCa samples (147)")+
  xlab("pAKT assessment by IHC") +
  scale_fill_manual(values =c("darkslateblue", "darkolivegreen3")) + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  labs(fill = "pAKT activity")
dev.off()

#| Can we set a threshold of the H-score based on the values of the pAKT obtained by IHC? Lets apply a logistic regression
pred <- prediction(sample_info_extracted$H_score[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$pAKT_positive[which(!is.na(sample_info_extracted$H_score))])
perf <- performance(pred,measure="tpr",x.measure="fpr")
par(family="serif")
plot(perf,
     colorize=FALSE,
     type="l",
     col ="blue",
     lwd=2,
     main = "AUC = 0.39")
abline(a=0,b=1)

# Área bajo la curva
AUC       <- performance(pred,measure="auc")
AUCaltura <- AUC@y.values

# Punto de corte óptimo
cost.perf <- performance(pred, measure ="cost")
opt.cut   <- pred@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
#coordenadas del punto de corte óptimo
x<-perf@x.values[[1]][which.min(cost.perf@y.values[[1]])]
y<-perf@y.values[[1]][which.min(cost.perf@y.values[[1]])]
points(x,y, pch=20, col="red")
sample_info_extracted <- sample_info
#| Can we set a threshold of the H-score based on the values of Mean expression FOXO?
threshold <- unique(sample_info_extracted$H_score)
threshold <- threshold[which(!is.na(threshold))]
threshold <- sort(threshold)
p_value_wilcox_FOXO <- c()
p_value_wilcox_FOXO_pathway <- c()
p_value_wilcox_PI3K_AKT_mTOR <- c()
p_value_wilcox_PTEN_loss <- c()
difference_wilcox_FOXO <- c()
difference_wilcox_FOXO_pathway <- c()
difference_wilcox_PI3K_AKT_mTOR <- c()
difference_wilcox_PTEN_loss <- c()

for (i in 1:(length(threshold)-1)){
  
  #| Creating a new column in the sample_info_extracted
  sample_info_extracted$H_score_threshold <- sample_info_extracted$H_score
  
  #| If values are below the threshold set as "PTEN loss", otherwise set "PTEN intact"
  sample_info_extracted$H_score_threshold[which(sample_info_extracted$H_score_threshold <= threshold[i] & !is.na(sample_info_extracted$H_score_threshold))] <- "PTEN loss"
  sample_info_extracted$H_score_threshold[which(sample_info_extracted$H_score_threshold != "PTEN loss" & !is.na(sample_info_extracted$H_score_threshold))] <- "PTEN intact"
  
  if (i == 1){
    cut_0 <- sample_info_extracted$H_score_threshold
  } else {
    if( i==2){
      cut_10<- sample_info_extracted$H_score_threshold
    } else {
      if(i== 3){
        cut_20 <- sample_info_extracted$H_score_threshold
      }else {
        if(i == 4){
          cut_30 <- sample_info_extracted$H_score_threshold
        }
      }
    }
  }
  
  #| Testing the mean differences for each signature
  #| 1) FOXO signature (+1200 genes)
  #test_wilcox_FOXO <- wilcox.test(sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$H_score_threshold == "PTEN loss")],sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )
  #p_value_wilcox_FOXO <- c(p_value_wilcox_FOXO, test_wilcox_FOXO$p.value) 
  #difference_wilcox_FOXO <- c(difference_wilcox_FOXO, mean(sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$H_score_threshold == "PTEN loss")]) - mean(sample_info_extracted$`Mean expression FOXO`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )) 
  #
  ##| 2) FOXO pathway signature (49 genes)
  #test_wilcox_FOXO_pathway <- wilcox.test(sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$H_score_threshold == "PTEN loss")],sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )
  #p_value_wilcox_FOXO_pathway <- c(p_value_wilcox_FOXO_pathway, test_wilcox_FOXO_pathway$p.value) 
  #difference_wilcox_FOXO_pathway <- c(difference_wilcox_FOXO_pathway, mean(sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$H_score_threshold == "PTEN loss")]) - mean(sample_info_extracted$`Mean expression FOXO pathway`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )) 
  
  #| 3) PI3K-AKT-mTOR signature (105 genes)
  test_wilcox_PI3K_AKT_mTOR <- wilcox.test(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "PTEN loss")],sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )
  p_value_wilcox_PI3K_AKT_mTOR <- c(p_value_wilcox_PI3K_AKT_mTOR, test_wilcox_PI3K_AKT_mTOR$p.value) 
  difference_wilcox_PI3K_AKT_mTOR <- c(difference_wilcox_PI3K_AKT_mTOR, mean(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "PTEN loss")]) - mean(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )) 
  
  #| 4) PTEN loss signature (812 genes)
  #test_wilcox_PTEN_loss <- wilcox.test(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$H_score_threshold == "PTEN loss")],sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )
  #p_value_wilcox_PTEN_loss <- c(p_value_wilcox_PTEN_loss, test_wilcox_PTEN_loss$p.value) 
  #difference_wilcox_PTEN_loss <- c(difference_wilcox_PTEN_loss, mean(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$H_score_threshold == "PTEN loss")]) - mean(sample_info_extracted$`Mean expression PTEN_loss`[which(sample_info_extracted$H_score_threshold == "PTEN intact")] )) 
  
}

difference_wilcox_FOXO[which(difference_wilcox_FOXO <0)] <- "Negative"
difference_wilcox_FOXO[which(difference_wilcox_FOXO >=0 & difference_wilcox_FOXO != "Negative")] <- "Positive"

difference_wilcox_FOXO_pathway[which(difference_wilcox_FOXO_pathway <0)] <- "Negative"
difference_wilcox_FOXO_pathway[which(difference_wilcox_FOXO_pathway >=0 & difference_wilcox_FOXO_pathway != "Negative")] <- "Positive"

difference_wilcox_PI3K_AKT_mTOR[which(difference_wilcox_PI3K_AKT_mTOR <0)] <- "Negative"
difference_wilcox_PI3K_AKT_mTOR[which(difference_wilcox_PI3K_AKT_mTOR >=0 & difference_wilcox_PI3K_AKT_mTOR != "Negative")] <- "Positive"

difference_wilcox_PTEN_loss[which(difference_wilcox_PTEN_loss <0)] <- "Negative"
difference_wilcox_PTEN_loss[which(difference_wilcox_PTEN_loss >=0 & difference_wilcox_PTEN_loss != "Negative")] <- "Positive"
length(p_value_wilcox_FOXO)

#| DATA FOR PLOTS 
data_frame_p_value_wilcoxon_threshold <- data.frame( threshold =threshold[-length(threshold)],
                                                     #p_value_wilcox_FOXO = p_value_wilcox_FOXO,
                                                     #p_value_wilcox_FOXO_pathway = p_value_wilcox_FOXO_pathway,
                                                     p_value_wilcox_PI3K_AKT_mTOR= p_value_wilcox_PI3K_AKT_mTOR,
                                                     #p_value_wilcox_PTEN_loss =p_value_wilcox_PTEN_loss,
                                                     #difference_wilcox_FOXO =difference_wilcox_FOXO,
                                                     #difference_wilcox_FOXO_pathway=difference_wilcox_FOXO_pathway,
                                                     difference_wilcox_PI3K_AKT_mTOR=difference_wilcox_PI3K_AKT_mTOR)
                                                     #difference_wilcox_PTEN_loss = difference_wilcox_PTEN_loss)

#| 1) FOXO signature (+1200 genes)
pdf("Results/Correlation_plots/FINDING_H_SCORE_THRESHOLD_p_value_FOXO.pdf")
ggplot(data_frame_p_value_wilcoxon_threshold, aes(x=threshold, y =p_value_wilcox_FOXO, color =difference_wilcox_FOXO)) +
  geom_point(size =2.5)+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("p-values of wilcoxon test FOXO signature")+ 
  ggtitle("Finding a threshold value of the H-score \n based on the FOXO signature")+
  xlab("threshold of H-score") +
  scale_color_manual(values = c("midnightblue", "green")) +
  labs(color = "Difference between means\n PTEN loss - PTEN intact")
dev.off()

#| 2) FOXO pathway signature (49 genes)
pdf("Results/Correlation_plots/FINDING_H_SCORE_THRESHOLD_p_value_FOXO_pathway.pdf")
ggplot(data_frame_p_value_wilcoxon_threshold, aes(x=threshold, y =p_value_wilcox_FOXO_pathway, color =difference_wilcox_FOXO_pathway)) +
  geom_point(size =2.5)+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("p-values of wilcoxon test FOXO pathway signature")+ 
  ggtitle("Finding a threshold value of the H-score \n based on the FOXO pathway signature")+
  xlab("threshold of H-score") +
  scale_color_manual(values = c("midnightblue", "green"))+
  labs(color = "Difference between means\n PTEN loss - PTEN intact")
dev.off()

#| 3) PI3K-AKT-mTOR signature (105 genes)
ggplot(data_frame_p_value_wilcoxon_threshold, aes(x=threshold, y =p_value_wilcox_PI3K_AKT_mTOR, color = difference_wilcox_PI3K_AKT_mTOR)) +
  geom_point(size =2.5)+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("p-values of wilcoxon test \n PI3K-AKT-mTOR pathway signature")+ 
  #ggtitle("Finding a threshold value of the H-score \n based on the PI3K-AKT-mTOR signature")+
  xlab("threshold of H-score") +
  scale_color_manual(values = c("black", "azure3")) +
  labs(color = "Difference between means\n PTEN loss - PTEN intact")
ggsave("Results/Correlation_plots/FINDING_H_SCORE_THRESHOLD_p_value_PI3K_AKT_mTOR.pdf", height = 4.5, width = 5.5)


#| 4) PTEN loss signature (812 genes)
pdf("Results/Correlation_plots/FINDING_H_SCORE_THRESHOLD_p_value_PTEN_loss.pdf")
ggplot(data_frame_p_value_wilcoxon_threshold, aes(x=threshold, y =p_value_wilcox_PTEN_loss, color =difference_wilcox_PTEN_loss)) +
  geom_point(size =2.5)+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("p-values of wilcoxon test PTEN loss signature")+ 
  ggtitle("Finding a threshold value of the H-score \n based on the PTEN loss signature")+
  xlab("threshold of H-score") +
  scale_color_manual(values = c("midnightblue", "green")) +
  labs(color = "Difference between means \n PTEN loss - PTEN intact")
dev.off()


sample_info_extracted$H_score_cut_0 <- cut_0
sample_info_extracted$H_score_cut_10 <- cut_10
sample_info_extracted$H_score_cut_20 <- cut_20
sample_info_extracted$H_score_cut_30 <- cut_30

sample_info_extracted$Mean.expression.PI3K.AKT.mTOR

ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = Mean.expression.PI3K.AKT.mTOR, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Mean expression PI3K-AKT-mTOR\n (Z score of Normalized log2 counts)")+ 
  scale_fill_manual(values =c("midnightblue", "orange")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_0_PI3K-AKT-mTOR_PTEN_intact_loss.pdf", height = 4.5, width = 5.5)



#| Purity as a function of PTEN loss and PTEN intact with H-score at 0
pdf("Results/Box_plots/Purity_vs_PTEN_Protein_Levels_H_score_thresholded_at_0.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x =H_score_cut_0, y =purity, fill =H_score_cut_0))+
  geom_boxplot() +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  theme(text=element_text(size=15,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("PTEN protein levels") +
  scale_fill_manual(values = c("midnightblue", "green")) +
  labs(fill = "PTEN status") +
  ggtitle("Behavior of PTEN status at protein level with H-score\n thresholded at 0 (147 PCa samples)")
dev.off()



cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$Diagnostico == "PCa")], sample_info_extracted$purity[which(sample_info_extracted$Diagnostico == "PCa")], method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 3)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Purity_vs_Mean_expression_PI3K-AKT-mTOR_PCa_198_samples.pdf")
ggplot(sample_info_extracted[which(sample_info_extracted$Diagnostico == "PCa"),], aes(x = `Mean expression PI3K-AKT-mTOR`, y =purity))+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("Mean expression PI3K-AKT-mTOR (Z score of normalized log2 counts)") +
  ggtitle("Purity as a function of PI3K-AKT-mTOR signature \n 198 PCa samples") + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") +
  annotate(geom="text", x=-0.25, y=0.7, label=paste( "r = ",r, ", ", "p < ", p, sep ="" ), color="black", family ="serif", size =6)
dev.off()


cor <- cor.test(sample_info_extracted$`Mean expression FOXO`[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$purity[which(!is.na(sample_info_extracted$H_score))], method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 3)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Purity_vs_Mean_expression_FOXO_PCa_198_samples.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x = `Mean expression FOXO`, y =purity))+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("Mean expression FOXO (Z score of normalized log2 counts)") +
  ggtitle("Purity as a function of FOXO signature \n 198 PCa samples") + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") +
  annotate(geom="text", x=-0.25, y=0.7, label=paste( "r = ",r, ", ", "p < ", p, sep ="" ), color="black", family ="serif", size =6)
dev.off()


cor <- cor.test(sample_info_extracted$`Mean expression FOXO pathway`[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$purity[which(!is.na(sample_info_extracted$H_score))], method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 2)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Purity_vs_Mean_expression_FOXO_pathway_PCa_198_samples.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x = `Mean expression FOXO pathway`, y =purity))+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("Mean expression FOXO pathway (Z score of normalized log2 counts)") +
  ggtitle("Purity as a function of FOXO pathway signature \n 198 PCa samples") + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") +
  annotate(geom="text", x=-0.25, y=0.7, label=paste( "r = ",r, ", ", "p < ", p, sep ="" ), color="black", family ="serif", size =6)
dev.off()


cor <- cor.test(sample_info_extracted$`Mean expression PTEN_loss`[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$purity[which(!is.na(sample_info_extracted$H_score))], method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 4)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Purity_vs_Mean_expression_PTEN_loss_PCa_198_samples.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x = `Mean expression PTEN_loss`, y =purity))+
  theme(text=element_text(size=16, family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("Mean expression PTEN loss (Z score of normalized log2 counts)") +
  ggtitle("Purity as a function of PTEN los signature \n 147 PCa samples") + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") +
  annotate(geom="text", x=-0.1, y=0.72, label=paste( "r = ",r, ", ", "p < ", p, sep ="" ), color="black", family ="serif", size =6)
dev.off()


cor <- cor.test(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(!is.na(sample_info_extracted$H_score))], sample_info_extracted$purity[which(!is.na(sample_info_extracted$H_score))], method ="pearson", exact=FALSE)
p  <- round(cor$p.value, 3)
r <- round(as.numeric(cor$estimate[[1]]),2)
pdf("Results/Correlation_plots/Purity_vs_Mean_expression_PI3K-AKT-mTOR_PCa_198_samples.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x = `Mean expression PTEN_loss`, y =purity))+
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Purity") +
  xlab("Mean expression PI3K-AKT-mTOR (Z score of normalized log2 counts)") +
  ggtitle("Purity as a function of PI3K-AKT-mTOR signature \n 147 PCa samples") + 
  geom_point(size=2.5, color ="midnightblue") + 
  geom_smooth(method = "lm", color ="green") +
  annotate(geom="text", x=-0.25, y=0.7, label=paste( "r = ",r, ", ", "p = ", p, sep ="" ), color="black", family ="serif", size =6)
dev.off()




pdf("Results/Correlation_plots/PTEN_Assessment_Correlation_PTEN_expression_H_score_pAKT_Activity.pdf")
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score)),], aes(x= H_score, y= PTEN_Exp_log2, color = pAKT_positive)) + geom_point(size=2.5) + 
  geom_smooth(color ="firebrick4") + 
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ggtitle("PTEN Expression vs H-score PTEN Assessment\n 147 PCa samples") +  
  ylab("PTEN Expression (Normalizede log2 counts)")+ 
  xlab("H-score") +
  labs(color ="pAKT activity") +
  scale_color_manual(values = c("midnightblue", "darkgoldenrod3")) +
  annotate(geom="text", x=150, y=10.2, label=expression(paste(rho, " = ", "-0.03", ", "," p < 0.75", sep ="" )), color="black", family ="serif", size =6)
dev.off()



#| What cut it is more suitable to use?
fraction_samples_PTEN_loss_intact <- data.frame(group =c("PTEN loss", "PTEN intact", "PTEN loss", "PTEN intact", "PTEN loss", "PTEN intact", "PTEN loss", "PTEN intact"),
                                                Samples = c(length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$H_score_cut_0== "PTEN loss")])/length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$H_score_cut_0== "PTEN intact")])/length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0))]), 
                                                            length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_10) & sample_info_extracted$H_score_cut_10== "PTEN loss")])/length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_10) & sample_info_extracted$H_score_cut_10== "PTEN intact")])/length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_20) & sample_info_extracted$H_score_cut_20== "PTEN loss")])/length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_20) & sample_info_extracted$H_score_cut_20== "PTEN intact")])/length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_30) & sample_info_extracted$H_score_cut_30== "PTEN loss")])/length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_30) & sample_info_extracted$H_score_cut_30== "PTEN intact")])/length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_0))])),
                                                cut = c("0", "0", "10", "10", "20", "20", "30", "30"))

fraction_samples_PTEN_loss_intact$round <- round(fraction_samples_PTEN_loss_intact$Samples, 2)
ggplot(fraction_samples_PTEN_loss_intact, aes(x =cut, y = Samples, fill =group)) +
  geom_bar(stat = "identity") + scale_fill_manual(values =c("midnightblue", "orange")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Fraction of samples") +
  xlab("H-score cutoff")+
  labs(fill ="PTEN protein status") +
  geom_text(aes(label=round),position="stack",vjust=1, colour = "white")#+ 
  #ggtitle("Fraction of samples for PTEN loss and PTEN intact changes depending on the \n H-score cut-off (147 PCa samples)")
ggsave("Results/Bar_plots/Fraction_samples_PTEN_loss_PTEN_intact_Different_H-score_cutoff_147_PCa_samples.pdf", heigh= 4.5, width = 5.5)
library(viridis)
################################################################################




########### CORRELATING CLINOPATHOLOGICAL VARIABLES WITH PTEN LOSS #############

#| Age
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = Edad, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Age (years)")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_0_PTEN_intact_loss_Age.pdf", height = 6, width = 6)

#| BMI
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = BMI, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("BMI (kg/m)")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_0_PTEN_intact_loss_BMI.pdf", height = 6, width = 6)

#| PSA
sample_info_extracted$PSA <- as.numeric(sample_info_extracted$PSA)
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0) & !is.na(sample_info_extracted$PSA)),], aes(x = H_score_cut_0, y = PSA, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("PSA")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  xlab("PTEN protein status")


#| Gleason score
data_fraction_gleason <- data.frame(group = c(length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 6 & sample_info_extracted$H_score_cut_0 == "PTEN loss")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 6 & sample_info_extracted$H_score_cut_0 == "PTEN intact")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 7 & sample_info_extracted$H_score_cut_0 == "PTEN loss")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 7 & sample_info_extracted$H_score_cut_0 == "PTEN intact")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 8 & sample_info_extracted$H_score_cut_0 == "PTEN loss")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 8 & sample_info_extracted$H_score_cut_0 == "PTEN intact")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 9 & sample_info_extracted$H_score_cut_0 == "PTEN loss")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                              length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$Gleason_score_pieza == 9 & sample_info_extracted$H_score_cut_0 == "PTEN intact")])/length(sample_info_extracted$Gleason_score_pieza[which(!is.na(sample_info_extracted$H_score_cut_0))])),
                                    value = c("Gleason == 6", "Gleason == 6","Gleason == 7", "Gleason == 7","Gleason == 8", "Gleason == 8","Gleason == 9", "Gleason == 9"),
                                    value2 = c("PTEN loss", "PTEN intact","PTEN loss", "PTEN intact","PTEN loss", "PTEN intact","PTEN loss", "PTEN intact"))
data_fraction_gleason$round <- round(data_fraction_gleason$group, 2)

ggplot(data_fraction_gleason, aes(x =value, y = group, fill =value2)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = c("aquamarine3", "gray17")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Fraction of samples") +
  xlab("Gleason score")+
  labs(fill ="PTEN protein status") +
  geom_text(aes(label=round),position="stack",vjust=1, colour = "white", size =2)+ 
  ggtitle("") +
  ylim(0,1)
ggsave("Results/Bar_plots/H_score_cut_0_PTEN_intact_loss_Gleason_score.pdf", height = 6, width = 7.3)

sample_info_extracted$Gleason_score_pieza <- as.numeric(sample_info_extracted$Gleason_score_pieza)
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0) & !is.na(sample_info_extracted$Gleason_score_pieza)),], aes(x =H_score_cut_0, y = Gleason_score_pieza, color =H_score_cut_0)) +
  geom_boxplot(outlier.shape = NA) + #scale_fill_manual(values = c("aquamarine3", "gray17")) +
  theme(text=element_text(size=16,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Gleason score") +
  xlab("PTEN protein") +
  geom_point(aes(fill =H_score_cut_0), size = 4, shape = 21, position = position_jitterdodge(0.3)) +
  scale_color_manual(values = c("midnightblue", "aquamarine3")) +
  scale_fill_manual(values = c("midnightblue", "aquamarine3"))+
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)+
  labs(fill ="PTEN protein status") +
  guides(color="none")
#geom_text(aes(label=round),position="stack",vjust=1, colour = "white", size =2)


#| Volume
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = Volume..æl., fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Volumen")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  xlab("PTEN protein status")


#| Stroma
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = stromal, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Stromal")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-0.5,family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/STROMA_PTEN_loss_vs_intact.pdf", height = 6, width = 6.3)

#| Immune
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = immune, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Immune")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-0.5,family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/IMMUNE_PTEN_loss_vs_intact.pdf", height = 6, width = 6.3)

#| Purity
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = purity, fill = H_score_cut_0)) +
  geom_boxplot()+
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("purity")+ 
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5) +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/Purity_PTEN_loss_vs_intact.pdf", height = 6, width = 6.3)


############| Can PTEN loss at protein level predicts recurrence?
#| DFS.TIME
sample_info_extracted$DFS.STATUS <- as.numeric(sample_info_extracted$DFS.STATUS)
sample_info_extracted$DFS <- sample_info_extracted$DFS.STATUS
sample_info_extracted$DFS[which(sample_info_extracted$DFS.STATUS ==0)] <- "No recurrence"
sample_info_extracted$DFS[which(sample_info_extracted$DFS.STATUS ==1)] <- "Recurrence"


colour <- factor(sample_info_extracted$H_score_cut_0, labels = c("midnightblue", "aquamarine3"))

ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0) & !is.na(sample_info_extracted$DFS)),], aes(x = DFS, y = DFS.TIME, color = H_score_cut_0)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill =H_score_cut_0), size = 5, shape = 21, position = position_jitterdodge(0.3)) +
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("DFS")+ 
  scale_color_manual(values = c("midnightblue", "aquamarine3")) +
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  labs(fill ="PTEN status") +
  stat_compare_means(aes(group = H_score_cut_0),vjust=-0.98,family="serif", size =5,label = "p.format" ) +
  xlab("PTEN protein status") +
  guides(color="none")
ggsave("Results/Box_plots/H_score_cut_0_PTEN_intact_loss_DFS_TIME_STATUS.pdf", height = 7, width = 8)



ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$DFS.STATUS==1),], aes(x = DFS.TIME, y = mTOR_Exp_log2, color = H_score_cut_0)) +
  geom_point() +
  theme(text=element_text(size=13,  family="serif"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) 


################################################################################


################################################################################
#| Saving sample information
write.table(sample_info_extracted, "Results/Sample_info_table/sample_info_extracted.txt", sep ="\t",row.names=T)
################################################################################




###################### STUDYING H-SCORE WITH ONLY ZEROS ########################
# what is the idea of studying just the ones with zeros? In this case, we can 
hist(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$H_score == 0)])
hist(sample_info_extracted$PTEN_Exp_log2[which(sample_info_extracted$H_score != 0 & !is.na(sample_info_extracted$H_score) )])

hist(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$H_score == 0)])
hist(sample_info_extracted$`Mean expression PI3K-AKT-mTOR`[which(sample_info_extracted$H_score != 0 & !is.na(sample_info_extracted$H_score) )])

ggplot(sample_info_extracted[which(sample_info_extracted$H_score == 0),], aes(x=pAKT_positive, y =`Mean expression PI3K-AKT-mTOR`)) +
  geom_boxplot() 


ggplot(sample_info_extracted[which(sample_info_extracted$H_score != 0 & !is.na(sample_info_extracted$H_score)),], aes(x=pAKT_positive, y =`Mean expression PI3K-AKT-mTOR`)) +
  geom_boxplot() +
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)


ggplot(sample_info_extracted[which(sample_info_extracted$H_score == 0),], aes(x=pAKT_positive, y =PTEN_Exp_log2)) +
  geom_boxplot() + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)


ggplot(sample_info_extracted[which(sample_info_extracted$H_score == 0),], aes(x=pAKT_positive, y =DFS.TIME)) +
  geom_boxplot() + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)

ggplot(sample_info_extracted[which(sample_info_extracted$H_score != 0  & !is.na(sample_info_extracted$H_score)),], aes(x=pAKT_positive, y =DFS.TIME)) +
  geom_boxplot() + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)


ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x=H_score_cut_0, y =purity)) +
  geom_boxplot() + 
  stat_compare_means(label ="p.format", vjust=-0.1, hjust=-1,family="serif", size =5)


ggplot(sample_info_extracted[which(sample_info_extracted$H_score == 0 ),], aes(x= `Mean expression PI3K-AKT-mTOR` , y =purity)) +
  geom_point() +ylim(0.75,1)

ggplot(sample_info_extracted[which(sample_info_extracted$H_score != 0  & !is.na(sample_info_extracted$H_score) ),], aes(x= `Mean expression PI3K-AKT-mTOR` , y =purity)) +
  geom_point() +ylim(0.75,1)


sample_info_extracted$DFS.TIME <- as.numeric(sample_info_extracted$DFS.TIME)
################################################################################