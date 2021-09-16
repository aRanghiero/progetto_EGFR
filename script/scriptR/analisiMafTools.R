if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# carico le library

library(RColorBrewer)
library(R.utils)
library(maftools) # l'ho installato usando il file .tar.gz scaricato da BioConductur
library(tidyverse)
update.packages("maftools")

setwd(dir='/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/maf/maf_filtrati5/')
name_maf <- list.files(path = "/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/maf/maf_filtrati5/", pattern = "*.maf")
head(name_maf)
#merged_maf <- merge_mafs(mafs = name_maf) # sto usando dei file MAF filtrati al 5%
#merged_maf
name_e19 <-list.files(path = "/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/ex19/maf_ex19/", pattern = "*.maf")
name_e21 <- list.files(path = "/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/ex21/maf_ex21/", pattern = "*.maf")
#merged_e19 <- merge_mafs(mafs = name_e19)
#merged_e21 <- merge_mafs(mafs = name_e21)

# carico i cnv-fusion

cnv_fusion <- read.csv('/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/vcf_tot/cnv/tot_cnv.csv',sep = ";",  header = F) 
cnv_fusion <- cnv_fusion %>% 
  rename(
    Gene = V1,
    Sample_name = V2,
    CN = V3
  )

cnv_fusion19 <- read.csv('/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/vcf_tot/cnv/ex19_cnv.csv', header = F) # qui ho la fusione su ALK
cnv_fusion21 <- read.csv('/media/albe/Volume/master/tirocinio/tirocinioAlberto/VCF/vcf_tot/cnv/ex21_cnv.csv', header = F) # qui ho la fusione su ALK
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

# mergio i file maf
merged_maf.plus.cnfus = merge_mafs(maf = name_maf, cnTable = cnv_fusion) 
merged_19.plus.cnfus = merge_mafs(maf = name_e19, cnTable = cnv_fusion19) 
merged_21.plus.cnfus = merge_mafs(maf = name_e21, cnTable = cnv_fusion21) 
head(merged_maf.plus.cnfus@data$Hugo_Symbol)
# creo delle liste di geni
gene_to_ing <- c("PMS2", "BRCA1", "CHEK1", "ATM", "BRCA2", "RNF43", "TP53", "NOTCH3", "MYCL", "POLE", "ATR", "FANCA", "FANCI", "SETD2", "FGFR4", "CCND3", "RICTOR")
genes <- c("PMS2", "BRCA2", "ATM", "EGFR", "CHEK1", "RNF43", "TP53", "NOTCH3", "MYCL", "POLE") # 10 geni più mutati
gene2 <- c("PDGFRB", "MDM2", "CCND3", "MLH1", "NOTCH1", "NOTCH2", "CSF1R", "PIK3R1", "POLE","FANCA", "FGFR4", "EGFR", "TP53") # i geni che ottengo dall'analisi di mafCompare
gene3 <- c("EGFR", "PDGFRB", "MDM2", "CCND3", "CDKN1B", "TP53")
#Shows sample summry.
getSampleSummary(merged_maf.plus.cnfus)
#Shows gene summary.
getGeneSummary(merged_maf.plus.cnfus)[1:10]
#shows clinical data associated with samples
getClinicalData(merged_maf.plus.cnfus)
getClinicalData(merged_19.plus.cnfus)
getClinicalData(merged_21.plus.cnfus)
#Shows all fields in MAF
getFields(merged_maf.plus.cnfus)
#Shows gene summary.
top50Gene <- getGeneSummary(merged_maf.plus.cnfus)[1:50,1]
top25Gene <- getGeneSummary(merged_maf.plus.cnfus)[1:25,1]
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = merged_maf, basename = 'merged_maf')

# Plotting MAF summary.
plotmafSummary(maf = merged_maf.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_19.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)
plotmafSummary(maf = merged_21.plus.cnfus, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F)

# Drawing oncoplots
# Changing colors for variant classifications name = 'Paired', 'Spectral', 'Set3'
vc_cols = RColorBrewer::brewer.pal(n = 11, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Complex_Event',
  'Amp',
  'Fusion'
)

col_oncoplot <- c("#44A7B6", "#1487C8", "#123B49", "#DC1E43", "#677BB1", "#28B158", "#B8DBDA", "#CAA10F", "#337B8D", "#8F3A37", "#1E5579")
names(col_oncoplot) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Complex_Event',
  'Amp',
  'Fusion'
)

sample19 <- getClinicalData(merged_19.plus.cnfus)
sample19 <- sample19[-c(1, 11, 12, 21, 23, 30)]
Tumor_Sample_Barcode <- factor(c("19-C-001559_v1", "19-I-010290_v1", "19-I-010852_v1", 
                                 "19-I-022128_v1", "19-I-025402_v2", "20-I-000029_v1")) # metto da parte i campioni con simultaneamente una delezione del 19 ed una mutazione del 21
sample19_21 <- data.frame(Tumor_Sample_Barcode)
sample21 <- getClinicalData(merged_21.plus.cnfus)
sample19 <- rbind(sample19_21 ,sample19, fill = TRUE)
sample_tot <- rbind(sample19 ,sample21, fill = TRUE)
head(sample_tot$Tumor_Sample_Barcode)
class(sample_tot)
oncoplot(maf = merged_maf.plus.cnfus, 
         gene_mar = 5, 
         fontSize = 0.5, 
         legend_height = 6,
         legendFontSize = 1.4,
         showTitle = T,
         drawBox = F,
         cBioPortal = T,
         top = 10,
         #top = 85, # evidenzio i primi 85 geni per mutazioni
         sampleOrder = sample_tot$Tumor_Sample_Barcode,  # uso il dataframe sample_tot$Tumor_Sample_Barcode per dare un ordine di rappresentazione, prima l'E19 poi l'E21
         sortByMutation = TRUE, # Force sort matrix according mutations. Helpful in case of MAF was read along with copy number data. Default FALSE.
         altered = TRUE, # Default FALSE. Chooses top genes based on muatation status. If TRUE chooses top genes based alterations (CNV or mutation).
         #genesToIgnore = gene_to_ing, # scarto i geni mutati nella maggioranza dei campioni
          #draw_titv = TRUE, 
         #colors = vc_cols, # carico l'elenco del colore creato precedentemente
         colors = col_oncoplot,
         bgCol = "#F3EFF2"
         # genes = c("ALK", "TP53", "EGFR"),
         # genes = gene3
         ) # vedo anche la fusione!!!
oncoplot(maf = merged_19.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 10)
oncoplot(maf = merged_21.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 10)

# Co-onco plots con le cnv e fusioni
coOncoplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = 'EX19', m2Name = 'EX21', genes = genes, removeNonMutated = TRUE, colors = vc_cols)
coOncoplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = 'EX19', m2Name = 'EX21', genes = rev(gene2), removeNonMutated = TRUE, colors = vc_cols)
coOncoplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = 'EX19', m2Name = 'EX21', genes = gene3, removeNonMutated = TRUE, colors = vc_cols) # geni significativi
coOncoplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = 'EX19', m2Name = 'EX21', colors = vc_cols, genes = top50Gene$Hugo_Symbol)

# Transition and Transversions
merged.titv = titv(maf = merged_maf.plus.cnfus, plot = FALSE, useSyn = TRUE)
merged19.titv = titv(maf = merged_19.plus.cnfus, plot = FALSE, useSyn = TRUE)
merged21.titv = titv(maf = merged_21.plus.cnfus, plot = FALSE, useSyn = TRUE)
plotTiTv(res = merged.titv, sampleOrder = sample_tot$Tumor_Sample_Barcode, )
plotTiTv(res = merged19.titv, sampleOrder = sample19$Tumor_Sample_Barcode)
plotTiTv(res = merged21.titv, sampleOrder = sample21$Tumor_Sample_Barcode)

# Plotting VAF
## questo pezzo l'ho fatto nel R script che si chiama cicloFOR.R

# Somatic Interactions
#exclusive/co-occurance event analysis on top 10 mutated genes. 
# i file merged senza le cnv e le fusion si trovano dentro lo script che si chiama cicloFOR.R
?somaticInteractions
somaticInteractions(maf = merged_maf, top = 20, fontSize = 0.5, countType = "all", returnAll = T, pvalue = c(0.05, 0.01))# posso anche fornigli un elenco di geni
somaticInteractions(maf = merged_19, top = 20, pvalue = c(0.05, 0.01), countType = "all", fontSize = 0.5,)
somaticInteractions(maf = merged_21, top = 35, pvalue = c(0.05, 0.01), fontSize = 0.5) # impostando a 50 ho una esclusività con un gene
somaticInteractions(maf = merged_maf, genes = c("ATM", "BRCA2", "CHEK1", "EGFR", "PMS2", "NOTCH3", "ALK"),  pvalue = c(0.05, 0.01), fontSize = 0.5, countType = "all" )

# Detecting cancer driver genes based on positional clustering
# non ho ben capito a cosa possa servire
merged.sig = oncodrive(maf = merged_maf.plus.cnfus, AACol = 'HGVSp', minMut = 5, pvalMethod = 'zscore', bgEstimate = T) # aggiungere i geni da polimorfismi?
head(merged.sig)
plotOncodrive(res = merged.sig, fdrCutOff = 0.01,  useFraction = TRUE, labelSize = 1, bubbleSize = 0.8)
plotOncodrive(res = merged.sig, fdrCutOff = 0.05) # fdrCutOff = false discovery rate
merged.sig = oncodrive(maf = merged_maf, AACol = 'HGVSp', minMut = 5, pvalMethod = 'zscore', ignoreGenes = gene_to_ing) # aggiungere i geni da polimorfismi?
head(merged.sig)
plotOncodrive(res = merged.sig, fdrCutOff = 0.01,  useFraction = TRUE, labelSize = 0.2, )
plotOncodrive(res = merged.sig, fdrCutOff = 0.05) # fdrCutOff = false discovery rate
# size of the points proportional to the number of clusters found in the gene
# X-axis shows number of mutations (or fraction of mutations) observed in these clusters.


# Adding and summarizing pfam domains
#  Proteins are generally composed of one or more functional regions, commonly termed domains. Different combinations of domains give rise to the diverse range of proteins found in nature. 
# The identification of domains that occur within proteins can therefore provide insights into their function.
# Pfam also generates higher-level groupings of related entries, known as clans. A clan is a collection of Pfam entries which are related by similarity of sequence, structure or profile-HMM. 
merged.pfam = pfamDomains(maf = merged_maf.plus.cnfus, AACol = 'HGVSp', top = 10)

#  Comparing two cohorts (MAFs)
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.

e19_e21_cnv_fus <- mafCompare(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = 'E19', m2Name = 'E21', minMut = 5, useCNV = F)
e19VSe21 <- mafCompare(m1 = merged_19, m2 = merged_21, m1Name = 'E19', m2Name = 'E21', minMut = 5)
print(e19VSe21)

# Forest plots
forestPlot(mafCompareRes = e19_e21_cnv_fus, pVal = 0.3, geneFontSize = 0.8, ) 
forestPlot(mafCompareRes = e19VSe21, pVal = 0.3, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
# ottengo 2 geni (pdgfrb e ccnd3) significativi se uso il file con anche le cnv e le fusioni
# mentre ottengo tre geni (pdgfrb, mdm2 e ccnd3) se uso il file senza le cnv e le fusioni

#Lollipop plot-2
install.packages("berryFunctions")
library(berryFunctions)
lolli_cols = RColorBrewer::brewer.pal(n = 3, name = 'Dark2')
names(lolli_cols) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp'
)

lollipopPlot2(m1 = merged_19.plus.cnfus, 
              m2 = merged_21.plus.cnfus, 
              gene = "EGFR", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              refSeqID = "NM_005228", 
              showDomainLabel = F,
              colors = lolli_cols
              )

lollipopPlot2(m1 = merged_19.plus.cnfus, 
              m2 = merged_21.plus.cnfus, 
              gene = "TP53", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              refSeqID = "NM_005228", 
              showDomainLabel = F,
              colors = lolli_cols
)
lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, 
              gene = "CDKN1B", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              showDomainLabel = F,
              colors = lolli_cols
              #refSeqID = "NM_005228", 
)
lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, 
              gene = "PDGFRB", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              colors = lolli_cols,
              showDomainLabel = F
              #refSeqID = "NM_005228", 
)
lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, 
              gene = "MDM2", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              colors = lolli_cols,
              showDomainLabel = F
              #refSeqID = "NM_005228", 
)

lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, 
              gene = "CCND3", 
              AACol1 = "HGVSp_Short", 
              AACol2 = "HGVSp_Short", 
              m1_name = "ex19", 
              m2_name = "ex21", 
              refSeqID = "NM_001760", 
              colors = lolli_cols,
              showDomainLabel = F
)

# Drug-Gene Interactions
dgi = drugInteractions(maf = merged_maf.plus.cnfus, fontSize = 0.75, plotType = "bar", drugs = T)
egfr.dgi = drugInteractions(genes = "EGFR", drugs = TRUE)
egfr.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

# Oncogenic Signaling Pathways
OncogenicPathways(maf = merged_maf, panelWidths = c(2,4,0))
PlotOncogenicPathways(maf = merged_maf, tsgCol = "#AB1F34", ogCol = "#1487C8",  
                      #pathways = c("RTK-RAS", "PI3K", "Cell_Cycle"),
                      pathways = c("RTK-RAS", "PI3K", "Cell_Cycle"), 
                      showTumorSampleBarcodes = F,
                      fontSize = 0.5,
                      #pathways = c("RTK-RAS", "PI3K", "Cell_Cycle", "NOTCH"), 
                      #sampleOrder = sample21$Tumor_Sample_Barcode,
                      #sampleOrder = sample_tot$Tumor_Sample_Barcode
                      ) # vedo a vedere i geni mutati in che pathways si trovano

OncogenicPathways(maf = merged_19)
PlotOncogenicPathways(maf = merged_19, 
                      pathways = c("RTK-RAS", "PI3K"),
                      #pathways = "RTK-RAS", 
                      showTumorSampleBarcodes = F,
                      #pathways = c("RTK-RAS", "PI3K", "Cell_Cycle", "NOTCH"), 
                      #sampleOrder = sample21$Tumor_Sample_Barcode,
                      # sampleOrder = sample_tot$Tumor_Sample_Barcode
) # vedo a vedere i geni mutati in che pathways si trovano

OncogenicPathways(maf = merged_21)
PlotOncogenicPathways(maf = merged_maf.plus.cnfus, 
                      pathways = c("RTK-RAS", "PI3K"),
                      #pathways = "RTK-RAS", 
                      showTumorSampleBarcodes = F,
                      #pathways = c("RTK-RAS", "PI3K", "Cell_Cycle", "NOTCH"), 
                      #sampleOrder = sample21$Tumor_Sample_Barcode,
                      # sampleOrder = sample_tot$Tumor_Sample_Barcode
) # vedo a vedere i geni mutati in che pathways si trovano



# Heterogeneity in tumor samples
# This heterogeneity can be inferred by clustering variant allele frequencies.
library("mclust")
I19_015657.het = inferHeterogeneity(maf = merged_maf, 
                                      tsb = '19-I-015657',
                                      vafCol =  "t_AF") # secondo me il problema è la colonna della VAF
I19_010852.het = inferHeterogeneity(maf = merged_maf, 
                                    tsb = '19-I-010852_v1',
                                    vafCol =  "t_AF") 
print(I19_015657.het$clusterMeans)
print(I19_010852.het$clusterMeans)
#Visualizing results
plotClusters(clusters = I19_015657.het) 
plotClusters(clusters = I19_010852.het) 

barcode_19 <- merged_19@clinical.data
for (i in barcode_19) {
  i.het = inferHeterogeneity(maf = merged_19, 
                                tsb = i,
                                vafCol =  "t_AF") 
  het_tot <- data.frame(i.het)
}
plotClusters(clusters = i.het) 
library(ggplot2)
library(ggridges)
theme_set(theme_minimal())



# Compare mutation load against TCGA cohorts
laml.mutload = tcgaCompare(maf = c(merged_maf, merged_19, merged_21), cohortName = c('Studio EGFR', 'EGFR E19', 'EGFR E21'), logscale = TRUE, primarySite = F, rm_hyper = F, medianCol = "red")

# Draw two barplots side by side for cohort comparision.
maftools::coBarplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, m1Name = "ex19", m2Name = "ex21", colors = col_oncoplot, genes = rev(top25Gene$Hugo_Symbol), geneMar = 4, showPct = F)
maftools::coBarplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus, borderCol = "white", 
                    m1Name = "E19", m2Name = "E21", colors = col_oncoplot, genes =c("EGFR", "PDGFRB", "MDM2", "CCND3", "TP53"), geneMar = 4, showPct = F)

# Creates a Genotype Matrix for every variant
genotypeMatrix(merged_maf, genes = top25Gene$Hugo_Symbol)
matrix_e19 <- maftools::mutCountMatrix(merged_19, )
matrix_e21 <- maftools::mutCountMatrix(merged_21, )
matrix_tot <- mutCountMatrix(merged_maf)
head(matrix_e19)
class(matrix_e19)
write.table(matrix_tot, file = "matrix_tot.csv", sep = ",", row.names=T, quote=FALSE)
write.table(matrix_e19, file = "matrix_19.csv", sep = ",", row.names=T, quote=FALSE)
# calculates MATH (Mutant-Allele Tumor Heterogeneity) score.
math.score(merged_maf, vafCol = "t_AF")

# con queste funzioni posso scaricare delle coorti da tcga!
tcgaAvailable(repo = "gitee")
tcgaAvailable(repo = c("github", "gitee"))
tcgaDriverBP(merged_maf, )
luad.maf = tcgaLoad(study = "LUAD")
lusc.maf = tcgaLoad(study = "LUSC")
lung.maf = merge_mafs(mafs = c(luad.maf, lusc.maf)) # provo a caricare due cohorti per vedere se ottengo un p-value decente
luad.maf2 = TCGAmutations::tcga_load(study = "LUAD")
plotmafSummary(maf = luad.maf, rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE)
oncoplot(maf = luad.maf, top = 10)
somaticInteractions(lung.maf,  genes = c("EGFR", "TP53", "MDM2", "PDGFRB", "KRAS"))
lollipopPlot(lung.maf, gene = "EGFR", cBioPortal = F, refSeqID = "NM_005228", showDomainLabel = F, repel = F, labelPos = 790)
mafSurvival(luad.maf, isTCGA = F, genes = c("EGFR", "MDM2"), time = "days_to_death", Status = "surv_event", ) # il p-value è ancora non significativo!
head(luad.maf@clinical.data$days_to_death)
head(luad.maf@clinical.data$days_to_last_followup)
