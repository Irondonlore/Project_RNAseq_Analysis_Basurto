#################### SURVIVAL ANALYSIS AT DIFFERENT LEVELS #####################
############################ PTEN LOSS AND RECURRENCE ##########################
sample_info$DFS_MONTHS <- as.numeric(sample_info$DFS_MONTHS)
sample_info$DFS_STATUS <- gsub(":DiseaseFree", "",sample_info$DFS_STATUS)
sample_info$DFS_STATUS <- gsub(":Recurred/Progressed", "",sample_info$DFS_STATUS)
sample_info$DFS_STATUS[which(sample_info$DFS_STATUS =="[Not Available]")] <- NA
sample_info$DFS_STATUS <- as.numeric(sample_info$DFS_STATUS)

#| Protein
sample_info$quantile <- cut(as.numeric(sample_info$PTEN_protein),breaks=quantile(as.numeric(sample_info$PTEN_protein), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_info)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN protein TCGA")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = TRUE,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_Protein.pdf", heigh = 6, width = 6)

#| Protein Quantile 1 and 4
sample_info$quantile <- cut(as.numeric(sample_info$PTEN_protein),breaks=quantile(as.numeric(sample_info$PTEN_protein), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sample_info$quantile[which(sample_info$quantile ==2)] <- 1
sample_info$quantile[which(sample_info$quantile ==3)] <- 4
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_info)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(2))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN protein TCGA")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = TRUE,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_Protein_Quantile_1_4.pdf", heigh = 6, width = 6)


#| Protein only CNA 0
sample <- sample_info[which(sample_info$PTEN_cna =="0"),]
sample$quantile <- cut(as.numeric(sample$PTEN_protein),breaks=quantile(as.numeric(sample$PTEN_protein), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  #annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN protein TCGA CNA 0")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = TRUE,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_Protein_CNA_0.pdf", heigh = 6, width = 6)


#| Protein only CNA 0 1 and 4 quantile
sample <- sample_info[which(sample_info$PTEN_cna =="0"),]
sample$quantile <- cut(as.numeric(sample$PTEN_protein),breaks=quantile(as.numeric(sample$PTEN_protein), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sample$quantile[which(sample$quantile ==2)] <- 1
sample$quantile[which(sample$quantile ==3)] <- 4
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(2))+
  #annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN protein TCGA CNA 0")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = TRUE,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_Protein_CNA_0_Quantile_1_4.pdf", heigh = 6, width = 6)


#| mRNA
sample_info$quantile <- cut(as.numeric(sample_info$PTEN_mrna),breaks=quantile(as.numeric(sample_info$PTEN_mrna), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_info)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]],3), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN mRNA TCGA")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = p_value,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_mRNA.pdf", heigh = 6, width = 6)

#| mRNA 1 and 4 quantile
sample_info$quantile <- cut(as.numeric(sample_info$PTEN_mrna),breaks=quantile(as.numeric(sample_info$PTEN_mrna), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sample_info$quantile[which(sample_info$quantile ==2)] <- 1
sample_info$quantile[which(sample_info$quantile ==3)] <- 4
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_info)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]],3), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(2))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN mRNA TCGA")
#ggsurvplot(sfit, 
#           data = sample_info,
#           ggtheme = theme_classic(),
#           size=0.8,
#           font.family="serif",
#           font.x = c(14, "black"),
#           font.y = c(14, "black"),
#           font.tickslab = c(12),
#           palette = viridis(4),
#           pval = p_value,
#           legend.title = "Quantiles",
#           xlab = "Time (months)",
#           censor.shape="|", censor.size = 3.5)
ggsave("Results/Kaplan-Meier_Curve_PTEN_mRNA_Quantile_1_4.pdf", heigh = 6, width = 6)


#| mRNA only CNA 0
sample <- sample_info[which(sample_info$PTEN_cna =="0"),]
sample$quantile <- cut(as.numeric(sample$PTEN_mrna),breaks=quantile(as.numeric(sample$PTEN_mrna), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  #annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for mRNA protein TCGA")
ggsave("Results/Kaplan-Meier_Curve_PTEN_mRNA_CNA_0.pdf", heigh = 6, width = 6)

#| mRNA only CNA 0 with 1 and 4 quantile
sample <- sample_info[which(sample_info$PTEN_cna =="0"),]
sample$quantile <- cut(as.numeric(sample$PTEN_mrna),breaks=quantile(as.numeric(sample$PTEN_mrna), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sample$quantile[which(sample$quantile ==2)] <- 1
sample$quantile[which(sample$quantile ==3)] <- 4
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(2))+
  #annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for mRNA protein TCGA")
ggsave("Results/Kaplan-Meier_Curve_PTEN_mRNA_CNA_0_Quantile_1_2.pdf", heigh = 6, width = 6)

#| CNA
#|sample_info$quantile <- cut(as.numeric(sample_info$PT),breaks=quantile(as.numeric(sample_info$PTEN_mrna), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~PTEN_cna, data=sample_info[which(sample_info$PTEN_cna == "-1" |sample_info$PTEN_cna == "-2" |sample_info$PTEN_cna == "0" ),])
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]],3), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point() +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="CNA") +
  scale_color_manual(values =viridis(6))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN CNA TCGA")
ggsave("Results/Kaplan-Meier_Curve_PTEN_CNA.pdf", heigh = 6, width = 6)









################### SURVIVAL ANALYSIS ON THE OVERLAPPING GENES #################

sample_info$DFS_MONTHS <- as.numeric(sample_info$DFS_MONTHS)
sample_info$DFS_STATUS <- gsub(":DiseaseFree", "",sample_info$DFS_STATUS)
sample_info$DFS_STATUS <- gsub(":Recurred/Progressed", "",sample_info$DFS_STATUS)
sample_info$DFS_STATUS[which(sample_info$DFS_STATUS =="[Not Available]")] <- NA
sample_info$DFS_STATUS <- as.numeric(sample_info$DFS_STATUS)

sample_info$CEACAM22P <- as.numeric(counts_data[which(rownames(counts_data) == "NRIP3"),])
sample_info$quantile <- cut(as.numeric(sample_info$CEACAM22P),breaks=quantile(as.numeric(sample_info$CEACAM22P), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
sample_info$quantile[which(sample_info$quantile == 2)] <- 1
sample_info$quantile[which(sample_info$quantile == 3)] <- 4

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_info)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 7), nsmall = 7)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for CCDC78")
ggsave("Results/SurvivalPlots/Kapla-Meier_CCDC78.pdf", heigh =6, width=7.5)

sample_PTEN_loss <- sample_info[which(sample_info$PTEN_cna == "-2"),]
sample_PTEN_loss$DFS_MONTHS <- as.numeric(sample_PTEN_loss$DFS_MONTHS)
sample_PTEN_loss$DFS_STATUS <- gsub(":DiseaseFree", "",sample_PTEN_loss$DFS_STATUS)
sample_PTEN_loss$DFS_STATUS <- gsub(":Recurred/Progressed", "",sample_PTEN_loss$DFS_STATUS)
sample_PTEN_loss$DFS_STATUS[which(sample_PTEN_loss$DFS_STATUS =="[Not Available]")] <- NA
sample_PTEN_loss$DFS_STATUS <- as.numeric(sample_PTEN_loss$DFS_STATUS)

sample_PTEN_loss$CEACAM22P <- as.numeric(counts_data[which(rownames(counts_data) == "IL21R-AS1"),sample_PTEN_loss$PATIENT_ID])
sample_PTEN_loss$quantile <- cut(as.numeric(sample_PTEN_loss$CEACAM22P),breaks=quantile(as.numeric(sample_PTEN_loss$CEACAM22P), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS)~quantile, data=sample_PTEN_loss)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for PTEN protein TCGA")


d <- data.frame(CCDC78 = as.numeric(counts_data[which(rownames(counts_data) == "CCDC78"),]),
                info =sample_info$PTEN_cna)

ggplot(d[which(d$info == "-2" |d$info == "0" ),], aes(x= info, y =  CCDC78)) +
  geom_boxplot() +
  ylim("0, 1000")



exprMatrix
data_trait


data_trait$CCDC78 <- as.numeric(exprMatrix[which(rownames(exprMatrix) == "CCDC78"),])
data_trait$quantile <- cut(as.numeric(data_trait$CCDC78),breaks=quantile(as.numeric(data_trait$CCDC78), na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_trait$quantile[which(data_trait$quantile == 2)] <- 1
data_trait$quantile[which(data_trait$quantile == 3)] <- 4

sfit <- survfit(Surv(DFS.TIME, DFS.STATUS)~quantile, data=data_trait)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 2), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)
ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line(size =0.9) + 
  geom_point(shape ="|", size= 2.4) +
  theme(text=element_text(size=16,  family="serif"), plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+ 
  ylim(c(0,1)) + ylab("Survival probability") + 
  xlab("Time of recurrence (months)") +
  labs(color="Quantiles") +
  scale_color_manual(values =viridis(4))+
  annotate(geom="text", x=20, y=0.2, label=paste("p = ",p_value,sep=""), color="black", family ="serif", size =6) +
  ggtitle("Kaplan-Meier Curve for CCDC78")

ggsave("Results/SurvivalPlots/Kapla-Meier_CCDC78_BASURTO.pdf", heigh =6, width=7.5)
################################################################################

