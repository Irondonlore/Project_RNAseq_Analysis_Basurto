################################################################################
#   STATISTICS OF THE TRIMMING
################################################################################

workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/0_STEP_TRIMMING/"
setwd(workingDir)

################################################################################
### LIBRARIES
################################################################################
library(ggplot2)
library(ggthemes) # Load

################################################################################
### DATA
################################################################################
samples <- list.files(pattern = "*_TRIM.out")
samples <- gsub("_TRIM.out", "", samples)
samples

value <- c()
for(s in 1:length(samples)){
  
  #Reading the file
  file <- readLines(paste(samples[s], "_TRIM.out", sep = "")) 
  
  #Obtaining the initial value of the reads
  initial <- grep('Total basepairs processed', file, value = TRUE)
  initial <- substring(initial[1], first = 28, last = 37)
  initial <- gsub(",","",initial)
  initial <- as.numeric(initial)

  #Obtaining the final value of the reads
  final <- grep('Total written', file, value = TRUE)
  final <- substring(final[2], first = 28, last = 37)
  final <- gsub(",","",final)
  final <- as.numeric(final)
  
  #Comparing the values (%)
  percentage <- (initial - final)*100/initial
  
  #Saving the values (%)
  value <- c(value, percentage)
}


################################################################################
### ANALYSIS
################################################################################
pdf(file="Percentage_trimmed_basepairs.pdf")
ggplot(as.data.frame(value), aes(x=value)) + 
  ggtitle("Percentage of trimmed basepairs") +
  xlab("Percentage (%)") + geom_histogram(color="black", fill="lightblue", bins = 60) + 
  theme_classic(base_family = "serif",base_size = 22)
dev.off()
#Mean value (%)
mean(value)

#Median value (%)
median(value)

#Standard deviation 
sd(value)

#Quartiles
quantile(value)     
