##########################################################################################################################################################
############### utilizzo questo script per caricare i maf file filtrati al 5% - utilizzando un ciclo FOR #################################################
##########################################################################################################################################################
# 17/09/21

library(R.utils)
library(maftools) # l'ho installato usando il file .tar.gz scaricato da BioConductur
library(tidyverse)
update.packages("maftools")

# carico i maf TOTALI
setwd(dir='/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/')
name_maf <- list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/maf/", pattern = "*.maf", full.names = T) 
merged_maf = merge_mafs(maf = name_maf)

# carico i maf relativi all'EX19
name_e19 <-list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/ex19_maf/", pattern = "*.maf", full.names = T)
merged_19 = merge_mafs(maf = name_e19) 

# carico i maf relativi all'EX21
name_e21 <- list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/ex21_maf/", pattern = "*.maf", full.names = T)
merged_21 = merge_mafs(maf = name_e21) 

# ciclo FOR per tenere solo un valore di AF
# per i MAF totali
vaf_old<- data.frame(merged_maf@data$t_AF)
length(merged_maf@data$t_AF)
class(vaf_old)
for (i in 1:length(merged_maf@data$t_AF)){
  temp <- unlist(strsplit(merged_maf@data$t_AF[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  # vaf_new$X1.3230[i] <- data.frame(temp_num[pos])
  vaf_old$y[i] <- data.frame(temp_num[pos])
} 
tail(vaf_old, n = 300)
class(as.character(vaf_old$y))
merged_maf@data$t_AF <- as.character(vaf_old$y)

# per i MAF EX19
vaf_old_19<- data.frame(merged_19@data$t_AF)
length(merged_19@data$t_AF)
class(vaf_old)
for (i in 1:length(merged_19@data$t_AF)){
  temp <- unlist(strsplit(merged_19@data$t_AF[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  vaf_old_19$y[i] <- data.frame(temp_num[pos])
}
class(as.character(vaf_old_19$y))
merged_19@data$t_AF <- as.character(vaf_old_19$y)

# per i MAF EX21
vaf_old_21<- data.frame(merged_21@data$t_AF)
length(merged_21@data$t_AF)
class(vaf_old)
for (i in 1:length(merged_21@data$t_AF)){
  temp <- unlist(strsplit(merged_21@data$t_AF[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  # vaf_new$X1.3230[i] <- data.frame(temp_num[pos])
  vaf_old_21$y[i] <- data.frame(temp_num[pos])
} 
class(as.character(vaf_old$y))
merged_21@data$t_AF <- as.character(vaf_old_21$y)

# carico i cnv-fusion

cnv_fusion <- read.csv('/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/cnv/tot_cnv.csv', sep = ",",  header = F) 
cnv_fusion <- cnv_fusion %>% 
  rename(
    Gene = V1,
    Sample_name = V2,
    CN = V3
  )

cnv_fusion19 <- read.csv('/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/cnv/CNV_EX19.csv', header = F, sep = ",")
cnv_fusion21 <- read.csv('/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/cnv/CNV_EX21.csv', header = F, sep = ",") 
cnv_fusion21 <- cnv_fusion21 %>% 
  rename(
    Gene = V1,
    Sample_name = V2,
    CN = V3
  )
cnv_fusion19 <- cnv_fusion19 %>% 
  rename(
    Gene = V1,
    Sample_name = V2,
    CN = V3
  )

# unisco i MAF con le CNV
merged_maf.plus.cnfus = merge_mafs(maf = name_maf, cnTable = cnv_fusion) 
merged_19.plus.cnfus = merge_mafs(maf = name_e19, cnTable = cnv_fusion19) 
merged_21.plus.cnfus = merge_mafs(maf = name_e21, cnTable = cnv_fusion21) 

# Plotting MAF summary, per vedere se tutto è andato bene
plotmafSummary(maf = merged_maf.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_19.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_21.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)

getClinicalData(merged_19)
getClinicalData(merged_21)

is.null(getClinicalData(merged_maf)) # controllo se ho dei valori nulli nei nomi dei campioni
