################################################################################
#################### SURVIVAL ANALYSIS OF THE BASURTO DATA #####################
################################################################################

#| This script is written to analyse a survival analysis on the samples for the
#| data AC-45_RNAseq-FFPE.

#| Survival analysis are classical statistical methods to explore the effect of 
#| parameters on a study across time. This analysis is used in several ways:

#| *) To describe the survival times of members of a group
#|          - Life tables
#|          - Kaplan-Meier curves
#|          - Survival function
#|          - Hazard function
#| *) To compare the survival times of two or more groups
#|          - Log-rank test
#| *) To describe the effect of categorical or quantitative variables on survival
#|          - Cox proportional hazards regression
#|          -Parametric survival models
#|          -Survival trees
#|          -Survival random forests

#| Proportional hazards: These are a class of survival models in statistics that 
#| related the time that passes, before some event ocurrs, to one or more covariates
#| that may be associated with that quantity of time. The models assumes that the
#| rate of change of the variable of interest across time is linear-exponentially
#| proportional. 

#| Cox analysis: The idea behind the cox proportional hazard model is to fit a 
#| linear model following the proportional hazard assumption, and to test if the 
#| covariates involved in the model have son effects on the differents between 
#| groups or quantitative variables. 

#| Survival function S(t): The probability that an individual survives longer than
#| time t.

#| Censored variables: These are variable that have been stopped followed in the
#| studied at a given time.

#| INTERPRETATION OF THE RESULTS FROM THE COX ANALYSIS:

#|      1) coef = regression coefficient who was estimated from the model. A negative
#|      value means more survival, a positive coeff means worse survival. This is 
#|      because of the ratio-nature of the harzard model.

#|      2) exp(coef) = exponential value of the estimated coefficient

#|      3) se(coef) = standard deviation of the coefficient


#| Note on statistical tests: The likelihood ratio test is almost always preferable 
#| to the Wald test, unless computational demands make it impractical to refit the model.          
################################################################################



################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))  
suppressMessages(library(survminer)) 
suppressMessages(library(ggfortify))
################################################################################



################################################################################
info.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_COUNTS_phenotype_Normalized_Log2_PTEN_loss_intact.txt"
dir.proj <- "W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/4_Survival_Analysis_Basurto/"
setwd(dir.proj)
#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("W:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| For plots
theme_set(theme_classic())
################################################################################




################################################################################
################################## DATA ########################################
################################################################################

#| Readtable Full counts. NOTE: it is not necessary to do any normalization 
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Readtable Sample information
sample_info <- read.table(info.file, sep ="\t")

#| Filtering BPH
sample_info_S <- sample_info[which(sample_info$Diagnostico == "PCa" & !is.na(sample_info$DFS.STATUS)),]
counts_data_S <- counts_data[,which(colnames(counts_data) %in% rownames(sample_info_S))]

#| Separating the groups for gleason score (two groups <=7 and >7 )
sample_info_S$Gleason_score_group <- sample_info_S$Gleason_score_pieza
sample_info_S$Gleason_score_group[sample_info_S$Gleason_score_group <= 7] <- "   < eq 7"
sample_info_S$Gleason_score_group[sample_info_S$Gleason_score_group > 7] <- "   > 7"

#| Setting DFS.S as factors
sample_info_S$Status <- sample_info_S$DFS.STATUS
sample_info_S$Status[sample_info_S$Status == "0"] <- "No recurrence"
sample_info_S$Status[sample_info_S$Status == "1"] <- "Recurrence"


################################################################################
#| Survival analysis
################################################################################

#| Gleason score
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ Gleason_score_group, data=sample_info_S)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 6), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point() +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + 
  ylab("Non-recurrence probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Gleason Score") +
  scale_color_manual(values =c("midnightblue", "orange"))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for Gleason")
ggsave("Results/Kaplan-Meier_Curve_Basurto_Gleason_score.pdf", height = 5, width = 6)

#| Gleason score ALL
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ Gleason_score_pieza, data=sample_info_S[which(!is.na(sample_info_S$Gleason_score_pieza)),])
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 8), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point() +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + 
  ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Gleason Score") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=17, y=0.05, label=paste("p = ","0.0001",sep=""), color="black", family ="serif", size =3.8) +
  ggtitle("Kaplan-Meier Curve for Gleason")
ggsave("Results/Kaplan-Meier_Curve_Basurto_Gleason_score_all.pdf", height = 5, width = 6)

ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           title="Kaplan-Meier Curve for Prostate Cancer Survival: Gleason score",# add title to plot
           ggtheme = theme_survminer( base_family="serif"), font.family = "serif"
           ) 
dev.off()

#| Estadio_pN_nodulospos_pieza
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ Estadio_pN_nodulospos_pieza, data=sample_info_S)
pdf("Results/SurvivalAnalysis/Kaplan-Meier_Curve_PCa_Estadio_pN_nodulospos_pieza_AC-45_RNAseq-FFPE_data.pdf")
ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           title="Kaplan-Meier Curve for Prostate Cancer Survival: Estadio pN nodulospos pieza",# add title to plot
           ggtheme = theme_survminer( base_family="serif"), font.family = "serif"
) 
dev.off()

#| H-score PTEN
sample_info$H_score_cut_0
#sample_info_S$quantile <- cut(as.numeric(sample_info_S$H_score),breaks=quantile(as.numeric(sample_info_S$H_score),probs =  c(0.5, 0.75,1), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS)~H_score_cut_0, data=sample_info_S)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.numeric(sfit$strata)
sfit$strata[which(sfit$strata==1)] <- "PTEN loss"
sfit$strata[which(sfit$strata==2)] <- "PTEN intact"

ggplot(sfit, aes(x=time, y= surv, color =strata)) +geom_line() + geom_point() +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Non-recurrence probability") + xlab("Time of recurrence (months)") +
  labs(color="PTEN status") +
  scale_color_manual(values =c("midnightblue", "orange"))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) 
  #ggtitle("Kaplan-Meier Curve for H-score PTEN Basurto")
ggsave("Results/Kaplan-Meier_Curve_H-score_Basurto.pdf", heigh=4.5, width = 6)

quantile(as.numeric(data_phenotype_basurto$PTEN_mrna), probs=c(0.25, 0.75),na.rm = TRUE)
#| PTEN expression
data_phenotype_basurto$quantile <- cut(as.numeric(data_phenotype_basurto$PTEN_mrna),breaks=quantile(as.numeric(data_phenotype_basurto$PTEN_mrna),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
#data_phenotype_basurto$quantile[which(data_phenotype_basurto$quantile ==2)] <-1
#data_phenotype_basurto$quantile[which(data_phenotype_basurto$quantile ==3)] <-4
sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ quantile, data=data_phenotype_basurto)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 6), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point() +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + 
  ylab("Non-recurrence probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles of \n PTEN mRNA") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) #+
  #ggtitle("Kaplan-Meier Curve for Gleason")
ggsave("Results/Kaplan-Meier_Curve_Basurto_PTEN_mRNA.pdf", height = 4.5, width = 6)

data_phenotype_basurto$H_score_cut_0
ggplot(data_phenotype_basurto[which(!is.na(data_phenotype_basurto$H_score_cut_0)),], aes(x = H_score_cut_0, y = PTEN_mrna, fill =H_score_cut_0)) +
  geom_boxplot(outlier.shape = NA) +
  theme(text=element_text(size=16,  family="serif"), legend.position = "none",plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  xlab("PTEN status") + 
  ylab("PTEN expression \n(Normalized log2 counts)") +
  scale_fill_manual(values =c("midnightblue", "orange"))+
  stat_compare_means(label ="p.format",family="serif", size =5)
ggsave("Results/PTEN_mRNA_PTEN_protein_loss_vs_intact.pdf", height = 4.2, width = 4.7)


################################################################################
library(viridis)
ggplot(sample_info_S, aes(x=H_score_cut_10, y =PTEN_Exp)) +
  geom_boxplot() +
  stat_compare_means(label ="p.format", family="serif", size =5)


data <- data.frame(counts = c(length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 6 & sample_info_S$H_score_cut_0 == "PTEN loss")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 7 & sample_info_S$H_score_cut_0 == "PTEN loss")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 8 & sample_info_S$H_score_cut_0 == "PTEN loss")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 9 & sample_info_S$H_score_cut_0 == "PTEN loss")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 6 & sample_info_S$H_score_cut_0 == "PTEN intact")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 7 & sample_info_S$H_score_cut_0 == "PTEN intact")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 8 & sample_info_S$H_score_cut_0 == "PTEN intact")]),
                               length(sample_info_S$Gleason_score_pieza[which(sample_info_S$Gleason_score_pieza == 9 & sample_info_S$H_score_cut_0 == "PTEN intact")])),
                   status =c("PTEN loss", "PTEN loss", "PTEN loss", "PTEN loss", "PTEN intact", "PTEN intact", "PTEN intact", "PTEN intact"),
                   score = c("6", "7", "8","9","6", "7", "8","9" ))

ggplot(data, aes(x =score, y =counts, fill=status)) +
  geom_bar(stat="identity") +
  xlab("Gleason score") +
  ylab("Counts") +
  labs(fill ="PTEN status") +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
  scale_fill_manual(values =c("midnightblue", "orange"))
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Bar_plots/Gleason_score_PTEN_Loss_vs_intact.pdf", height = 4.5, width = 5.5)

i <- sample_info[which(!is.na(sample_info$H_score_cut_0)),]
ggplot(i, aes(x =H_score_cut_0, y =PTEN_Exp, fill=H_score_cut_0)) +
  geom_boxplot() +
  xlab("PTEN status") +
  ylab("PTEN mRNA expression (Normalized log2 counts)") +
  labs(fill ="PTEN status") +
  theme(text=element_text(size=16,  family="serif"),legend.position = "none", plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
  scale_fill_manual(values =c("midnightblue", "orange")) +
  stat_compare_means(label ="p.format", family="serif", size =5)
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/PTEN_mRNA_PTEN_Loss_vs_intact.pdf", height = 4.5, width = 5.3)

################################################################################
#|  1)  Cox regression analysis for continuous variables
################################################################################

#| Fixing the cox model with all variables of interest to see the results we get and which of them are significant in nature.
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ Edad_Zscore + PTEN_Exp + DV200.value + BMI, data=sample_info_S)
summary(fit)
library(WGCNA)
packageVersion("gprofiler") 

#| Fixing the cox model with PTEN and Edad
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ Edad_Zscore + PTEN_Exp, data=sample_info_S)
summary(fit) #| From the results: As PTEN decreases, it reduces the chance of death by 0.04%

#| Fixing the cox model with Edad, DV200, H-score
sample_info_S_H_score_filtered <- sample_info_S[which(!is.na(sample_info_S$H_score)),]
sample_info_S_H_score_filtered$H_score_Zscore <- (sample_info_S_H_score_filtered$H_score - mean(sample_info_S_H_score_filtered$H_score))/sd(sample_info_S_H_score_filtered$H_score)
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ DV200_Zscore + Edad_Zscore + H_score_Zscore, data=sample_info_S_H_score_filtered)
summary(fit)

#| Fixing the cox model with only H-score
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~  H_score_Zscore, data=sample_info_S_H_score_filtered)
summary(fit)
################################################################################


