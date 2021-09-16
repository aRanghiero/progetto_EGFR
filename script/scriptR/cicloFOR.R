## per il file merged_maf
library(maftools)
setwd(dir='/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/')
name_maf <- list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/maf/", pattern = "*.maf", full.names = T) 
merged_maf = merge_mafs(maf = name_maf)

# name2_maf <- dir(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/maf/", pattern = "*.maf", full.names = T)
# merged_maf = merge_mafs(maf = name2_maf)

x <- merged_maf@data$t_AF
tail(x, n = 300)
vaf_old<- data.frame(x)
length(merged_maf@data$t_AF)
class(vaf_old)
# vaf_new <- data.frame(1:length(x))
for (i in 1:length(x)){
  temp <- unlist(strsplit(x[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  # vaf_new$X1.3230[i] <- data.frame(temp_num[pos])
  vaf_old$y[i] <- data.frame(temp_num[pos])
} 
tail(vaf_old, n = 300)
class(as.character(vaf_old$y))
merged_maf@data$t_AF <- as.character(vaf_old$y)
plotVaf(maf = merged_maf, vafCol = 't_AF', top = 10, keepGeneOrder = T, showN = T)

## per il file merged_19
name_e19 <-list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/ex19_maf/", pattern = "*.maf", full.names = T)
merged_19 = merge_mafs(maf = name_e19) 
x19 <- merged_19@data$t_AF
vaf_old_19<- data.frame(x19)
length(merged_19@data$t_AF)
class(vaf_old)

for (i in 1:length(x19)){
  temp <- unlist(strsplit(x19[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  vaf_old_19$y[i] <- data.frame(temp_num[pos])
}

class(as.character(vaf_old_19$y))
merged_19@data$t_AF <- as.character(vaf_old_19$y)
plotVaf(maf = merged_19, vafCol = 't_AF', top = 12, keepGeneOrder = T, showN = T)

## per il file merged_21
name_e21 <- list.files(path = "/media/albe/Volume/master/tirocinio/EGFRnew/progetto_EGFR/VCF/ex21_maf/", pattern = "*.maf", full.names = T)
merged_21 = merge_mafs(maf = name_e21) 
x21 <- merged_21@data$t_AF
vaf_old_21<- data.frame(x21)
length(merged_21@data$t_AF)
class(vaf_old)
# vaf_new <- data.frame(1:length(x21))
for (i in 1:length(x21)){
  temp <- unlist(strsplit(x21[i], ","))
  temp_num <- as.numeric(temp)
  pos <- which(temp_num==max(temp_num))
  # vaf_new$X1.3230[i] <- data.frame(temp_num[pos])
  vaf_old_21$y[i] <- data.frame(temp_num[pos])
} 

class(as.character(vaf_old$y))
merged_21@data$t_AF <- as.character(vaf_old_21$y)
tail(merged_21@data$t_AF)
plotVaf(maf = merged_21, vafCol = 't_AF', top = 12, keepGeneOrder = T, showN = T)

# dopo aver creato gli oggetti MAF
## creo un oggetto per settare il colore dei primi 12 geni per mutazione
bcol = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"))
names(bcol) = c(
  'PMS2',
  'BRCA2',
  'ATM',
  'EGFR',
  'CHEK1',
  'RNF43',
  'NOTCH3',
  'MYCL',
  'POLE', 
  'TP53',
  'ATR',
  'FANCI'
  )
top25Gene <- getGeneSummary(merged_maf)[1:25,1]
gene3 <- c("EGFR", "PDGFRB", "MDM2", "CCND3", "CDKN1B")
plotVaf(maf = merged_maf, vafCol = 't_AF', genes = top25Gene$Hugo_Symbol[1:12], keepGeneOrder = T, showN = T, color = bcol)
plotVaf(maf = merged_19, vafCol = 't_AF', genes = top25Gene$Hugo_Symbol[1:12], keepGeneOrder = T, showN = T)
plotVaf(maf = merged_21, vafCol = 't_AF', genes = top25Gene$Hugo_Symbol[1:12], keepGeneOrder = T, showN = T)
plotVaf(maf = merged_maf, vafCol = 't_AF', genes = c("EGFR", "PDGFRB", "MDM2", "CCND3", "TP53"), keepGeneOrder = T, showN = T, flip = F, )

plotmafSummary(maf = merged_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_19, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_21, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)

mafbarplot(maf = merged_maf)

#Detecting cancer driver genes based on positional clustering
gene_to_ing <- c("PMS2", "CHEK1", "ATM", "BRCA2", "RNF43", 
                 "NOTCH3", "RICTOR")
egfr.sig = oncodrive(maf = merged_maf, AACol = 'HGVSp', minMut = 10, ignoreGenes = gene_to_ing, pvalMethod = "zscore")
egfr.sig = oncodrive(maf = merged_maf, AACol = 'HGVSp', minMut = 10, ignoreGenes = gene_to_ing, pvalMethod = "poisson")
head(egfr.sig)
plotOncodrive(res = egfr.sig, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 0.6, colCode = c("#F27348", "#26474E"), )

# Lollipop plots for amino acid changes
## volevo evidenziare nel doppio lolliplot le mutazioni correlate a resistenza farmacologica
### purtroppo non lo si può fare
# install.packages("berryFunctions")
library(berryFunctions)
lolli_cols = RColorBrewer::brewer.pal(n = 3, name = 'Dark2')
names(lolli_cols) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp'
)

col_pastel = c("#218B82", "#EB96AA", "#98D4BB")
names(col_pastel) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp'
)

col_summer = c("#44A7B6", "#DC1E43", "#B07F54")
names(col_summer) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp'
)
lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, 
              gene = "EGFR", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              refSeqID = "NM_005228", 
              showDomainLabel = F,
              #colors = lolli_cols,
              colors = col_summer,
              roundedRect =T, 
              m1_label = "790",
              m2_label = "790",
              labPosAngle = 0,
              domainBorderCol = "white"
)

lollipopPlot(
  maf = merged_19,
  gene = 'EGFR',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  labelPos = 790,
  showDomainLabel = F,
  colors = col_summer,
  refSeqID = "NM_005228"
)

# somatic interation
genes_somaticInteractions <- getGeneSummary(x = merged_19)[7:35, Hugo_Symbol]
somaticInteractions(maf = merged_maf, top = 30, fontSize = 0.5, countType = "all", returnAll = T, pvalue = c(0.05, 0.01), genes = getGeneSummary(x = merged_maf)[6:35, Hugo_Symbol])# posso anche fornigli un elenco di geni
somaticInteractions(maf = merged_19, fontSize = 0.5, countType = "all", returnAll = T, top = 35, pvalue = c(0.05, 0.01), colPal = "RdBu", genes = getGeneSummary(x = merged_19)[7:35, Hugo_Symbol]) #genes = genes_somaticInteractions$Hugo_Symbol)# posso anche fornigli un elenco di geni
somaticInteractions(maf = merged_21, fontSize = 0.5, countType = "all", returnAll = T, pvalue = c(0.05, 0.01), top = 35, colPal = "RdBu", genes = getGeneSummary(x = merged_21)[8:35, Hugo_Symbol])

