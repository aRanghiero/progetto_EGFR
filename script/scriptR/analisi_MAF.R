##########################################################################################################################################################
############### uso questo script per analizzare i MAF dopo averli caricati  #############################################################################
##########################################################################################################################################################
# 17/09/21

library(RColorBrewer)
library(R.utils)
library(maftools) # l'ho installato usando il file .tar.gz scaricato da BioConductur
library(tidyverse)

# creo una palette di colori
# 1 palette
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
# 2 palette
col_oncoplot <- c("#0caeb8", "#1487C8", "#123B49", "#DC1E43", "#677BB1", "#c1d5e3", "#B8DBDA", "#CAA10F", "#337B8D", "#8F3A37", "#1E5579")
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


# creo una lista di nomi per suddividere i campioni nell'oncoplot tra EX19 e Ex21
# sample19 <- getClinicalData(merged_19.plus.cnfus)
#sample19_MH <- sample19[c(1,2,4,11,33,35,37,38,40,41,52,53,45)]
#sample19 <- sample19[-c(1,2,4,11,33,35,37,38,40,41,52,53,45)]
#sample21 <- getClinicalData(merged_21.plus.cnfus)
#sample21_MH <- sample21[c(4,5,13,22,28,31,39,51,52)]
#sample21 <- sample21[-c(4,5,13,22,28,31,39,51,52)]

#sample19 <- rbind(data.frame(sample19_MH) ,sample19, fill = T)
#sample21 <- rbind(data.frame(sample21_MH) ,sample21, fill = TRUE)
#sample_tot <- rbind(sample19 ,sample21, fill = TRUE)
#head(sample_tot$Tumor_Sample_Barcode)
#class(sample_tot)
# esporto la lista di nomi in un file csv in modo da non dover rifare questo pezzo ogni volta
#write_csv(sample_tot, "NameFolder/sample_tot.csv", col_names = F, ) 
#write_csv(sample19, "NameFolder/sample_E19.csv", col_names = F, )
#write_csv(sample21, "NameFolder/sample_E21.csv", col_names = F, )
#dora in poi dovrò eseguire solo queste tre stringhe ;)
sample_tot <- read_csv(file = "NameFolder/sample_tot.csv", col_names = "sample_name", )
sample_e19 <- read_csv(file = "NameFolder/sample_E19.csv", col_names = "sample_name", )
sample_e21 <- read_csv(file = "NameFolder/sample_E21.csv", col_names = "sample_name", )
sample_tot <- sample_tot$sample_name[-c(65, 122, 123)]
sample_e19 <- sample_e19$sample_name[-c(65)]
sample_e21 <- sample_e21$sample_name[-c(57)]
#genero un oncoplot MAF TOT
oncoplot(maf = merged_maf.plus.cnfus, 
         gene_mar = 5, 
         fontSize = 0.5, 
         legend_height = 6,
         legendFontSize = 1.4,
         showTitle = T,
         drawBox = F,
         cBioPortal = T,
         # top = 10,
         top = 50, # evidenzio i primi 50 geni per mutazioni
         sampleOrder = sample_tot,  # uso il dataframe sample_tot$Tumor_Sample_Barcode per dare un ordine di rappresentazione, prima l'E19 poi l'E21
         sortByMutation = TRUE, # Force sort matrix according mutations. Helpful in case of MAF was read along with copy number data. Default FALSE.
         altered = TRUE, # Default FALSE. Chooses top genes based on muatation status. If TRUE chooses top genes based alterations (CNV or mutation).
         #genesToIgnore = gene_to_ing, # scarto i geni mutati nella maggioranza dei campioni
         #draw_titv = TRUE, 
         #colors = vc_cols, # carico l'elenco del colore creato precedentemente, predominanza azzurro e verde
         colors = col_oncoplot, # carico l'elenco del colore creato precedentemente, predominanza fuxia e verde petrolio
         bgCol = "#F3EFF2", # uso il bianco come background
         # genes = c("ALK", "TP53", "EGFR"),
         # genes = gene3
) 

oncoplot(maf = merged_19.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 10, sampleOrder = sample_e19, showTumorSampleBarcodes = T, barcode_mar = 5,)
# vedo che il campione 64_EGFR 21-S-000015 presenta un numero elevato di mutazioni. Probabilmente da scartare
oncoplot(maf = merged_21.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 10, sampleOrder = sample_e21, showTumorSampleBarcodes = T, barcode_mar = 5,)
# vedo che il campione 6_EGFR 15-I-000801 presenta un numero elevato di mutazioni. Probabilmente da scartare

########### posso decidere di eliminare questi due file ##############
sample_tot2 <- sample_tot[-c(9, 71)]
subsetMAF_tot <- subsetMaf(merged_maf, tsb = sample_tot2) # in questo modo ho eliminato i due campioni con un elevato numerdo di mutazioni

# Transition and Transversions.
merged_maf.titv = titv(maf = merged_maf, plot = FALSE, useSyn = TRUE)
E19_maf.titv = titv(maf = merged_19, plot = FALSE, useSyn = TRUE)
E21_maf.titv = titv(maf = merged_21, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = merged_maf.titv)
plotTiTv(res = E19_maf.titv)
plotTiTv(res = E21_maf.titv)
# provo a fare la stessa cosa usando uno studio di TCGA
LUAD.titv = titv(maf = LUAD, plot = FALSE, useSyn = TRUE)
LUAD_fh.titv = titv(maf = LUAD_fh, plot = FALSE, useSyn = TRUE)
plotTiTv(res = LUAD.titv)
plotTiTv(res = LUAD_fh.titv)
### confrontando i mie campioni rispetto a quelli del TCGA si nota una differenza tra il numero di Transition e Transversions.
### ed anche una differenza nel tipo di alterazione, più legata al fumo

# Lollipop plots for amino acid changes

install.packages("berryFunctions") # carico un pacchetto per selezionare i colori
library(berryFunctions)
lolli_cols = RColorBrewer::brewer.pal(n = 3, name = 'Dark2')
names(lolli_cols) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp'
)

lollipopPlot( # plotto tutte le mutazioni su EGFR della casistica totale
  maf = merged_maf, 
  roundedRect = T,
  gene = 'EGFR',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  refSeqID = "NM_005228", 
  showDomainLabel = F,
  colors = lolli_cols,
  labelPos = 790,
  
)

# plotProtein(gene = "EGFR", refSeqID = "NM_005228") # posso anche plottare soltanto EGFR

# Compare mutation load against TCGA cohorts
# si vede chiaramente che le mie coorti presentano un numero di mutazioni inferiori se comparate ad altre coorti del polmone
EGFR_tot.mutload = tcgaCompare(maf = merged_maf.plus.cnfus, cohortName = 'EGFR_study', logscale = TRUE, capture_size = 50, col = c("#A3D2CA", "#5EAAA8"), medianCol = "#EB5E0B")
EGFR_e19.mutload = tcgaCompare(maf = merged_19.plus.cnfus, cohortName = 'EGFR_study', logscale = TRUE, capture_size = 50, col = c("#A3D2CA", "#5EAAA8"), medianCol = "#EB5E0B" )
EGFR_e21.mutload = tcgaCompare(maf = merged_21.plus.cnfus, cohortName = 'EGFR_study', logscale = TRUE, capture_size = 50, col = c("#A3D2CA", "#5EAAA8"), medianCol = "#EB5E0B" )

# Plotting VAF

bcol = c(RColorBrewer::brewer.pal(n = 12, name = "Paired")) #creo una palette
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
top25Gene <- getGeneSummary(merged_maf)[1:25,1] # ed una lista dei primi 12 geni mutati nel merged_maf
plotVaf(maf = merged_maf, vafCol = 't_AF', showN = T, color = bcol, top = 12, 
        #genes = top25Gene$Hugo_Symbol[1:12]
        )
plotVaf(maf = merged_19, vafCol = 't_AF', showN = T, color = bcol, top = 12, 
        #genes = top25Gene$Hugo_Symbol[1:12]
)
plotVaf(maf = merged_21, vafCol = 't_AF', showN = T, color = bcol, top = 12, 
        #genes = top25Gene$Hugo_Symbol[1:12]
)

# Somatic Interactions - presenza di co-occurance
somaticInteractions(maf = merged_maf, top = 40, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T)
somaticInteractions(maf = merged_19, top = 40, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T)
somaticInteractions(maf = merged_21, top = 40, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T)

# Comparing two cohorts (MAFs)
E19VSE21 <- mafCompare(m1 = merged_19, m2 = merged_21, m1Name = 'E19', m2Name = 'E21', minMut = 5, useCNV = F)
print(E19VSE21)
forestPlot(mafCompareRes = E19VSE21, pVal = 0.12, ) # ci sono 3 geni mutati in maniera significativa: CCND3, MDM2, SMARCA4

# Co-onco plots
# vado a vedere le differenze tra le due coorti
coOncoplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus,
           m1Name = 'E19', m2Name = 'E21',  
           removeNonMutated = TRUE, colors = col_oncoplot, 
           sampleOrder1 = sample_e19, sampleOrder2 = sample_e21,
           genes = c("EGFR", "TP53", "CCND3", "MDM2", "SMARCA4"), 
           geneNamefont = 0.7, gene_mar = 1.45, bgCol = "#eeeeee")

# Co-bar plots
# oppure posso vederla in questo modo
coBarplot(m1 = merged_19.plus.cnfus, m2 = merged_21.plus.cnfus,
          m1Name = 'E19', m2Name = 'E21',
          colors = col_oncoplot, genes = rev(c("EGFR", "TP53", "CCND3", "MDM2", "SMARCA4")), borderCol = "white", 
          geneMar = 5, showPct = F, legendTxtSize = 0.9, axisSize = 0.9 
            )

# Lollipop plot-2
lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, gene = "EGFR", labPosSize = 0.6, 
              m1_label = 790, m2_label = 790, 
              m1_name = 'E19', m2_name = 'E21', 
              domainBorderCol = "black", 
              roundedRect = T, 
              AACol1 = 'HGVSp_Short', legendTxtSize = 0.8,
              AACol2 = 'HGVSp_Short', 
              refSeqID = "NM_005228", 
              showDomainLabel = F,
              colors = lolli_cols,
              )

# Oncogenic Signaling Pathways
OncogenicPathways(maf = merged_maf, )
OncogenicPathways(maf = merged_19, )
OncogenicPathways(maf = merged_21, )

PlotOncogenicPathways(maf = merged_maf, pathways = "RTK-RAS", 
                      #sampleOrder = sample_tot
                      )
PlotOncogenicPathways(maf = merged_19, pathways = "RTK-RAS", 
                      #sampleOrder = sample_e19
)
PlotOncogenicPathways(maf = merged_21, pathways = "RTK-RAS", 
                      #sampleOrder = sample_e21
)
