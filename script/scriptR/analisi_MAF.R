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
# ora si trova dentro il file loading_maf.R
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
#sample_tot <- read_csv(file = "NameFolder/sample_tot.csv", col_names = "sample_name", )
#sample_e19 <- read_csv(file = "NameFolder/sample_E19.csv", col_names = "sample_name", )
#sample_e21 <- read_csv(file = "NameFolder/sample_E21.csv", col_names = "sample_name", )
#sample_tot <- sample_tot$sample_name[-c(65, 122, 123)]
#sample_e19 <- sample_e19$sample_name[-c(65)]
#sample_e21 <- sample_e21$sample_name[-c(57)]
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

oncoplot(maf = merged_19.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 50, sampleOrder = sample_e19, showTumorSampleBarcodes = F, barcode_mar = 5)

oncoplot(maf = merged_21.plus.cnfus, colors = col_oncoplot, altered = TRUE, top = 10, sampleOrder = sample_e21, showTumorSampleBarcodes = F, barcode_mar = 5,)


########### posso decidere di eliminare questi due MAF ##############
# i MAF 64_EGFR 21-S-000015, 6_EGFR 15-I-000801 sono stati eliminati dall'analisi 
# questo processo viene effetuato nello script loading_maf.R
# sample_tot2 <- sample_tot[-c(9, 71)]
# subsetMAF_tot <- subsetMaf(merged_maf, tsb = sample_tot2) # in questo modo ho eliminato i due campioni con un elevato numerdo di mutazioni

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

# install.packages("berryFunctions") # carico un pacchetto per selezionare i colori
library(berryFunctions)
lolli_cols = RColorBrewer::brewer.pal(n = 5, name = 'Dark2')
names(lolli_cols) = c(
  'Missense_Mutation',
  'In_Frame_Del',
  'Amp', 
  'Nonsense_Mutation',
  'Splice_Site'
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
EGFR_tot.mutload = tcgaCompare(maf = merged_maf.plus.cnfus, 
                               cohortName = 'EGFR_study', 
                               logscale = TRUE, capture_size = 50, 
                               col = c("#4E89AE", "#43658B"), 
                               medianCol = "#ED6663", 
                               bg_col = c("#FAF3F3", "#A7BBC7"), 
                               # tcga_cohorts = c("LUSC", "LUAD")
                               )
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
top50Gene <- getGeneSummary(merged_maf)[1:50,1]
somatic_gene <- top50Gene[(-c(1, 2, 3, 4))] # in questo modo elimino i 4 primi geni più mutati
top50Gene_e19 <- getGeneSummary(merged_19)[1:50,1]
somatic_gene_e19 <- top50Gene_e19[(-c(1, 2, 3, 4))] # in questo modo elimino i 4 primi geni più mutati
top50Gene_e21 <- getGeneSummary(merged_21)[1:50,1]
somatic_gene_e21 <- top50Gene_e21[(-c(1, 2, 3, 4, 5, 6))] # in questo modo elimino i 6 primi geni più mutati
somaticInteractions(maf = merged_maf, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T, genes = somatic_gene$Hugo_Symbol)
somaticInteractions(maf = merged_19, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T, genes = somatic_gene_e19$Hugo_Symbol)
somaticInteractions(maf = merged_21, pvalue = c(0.05, 0.1), fontSize = 0.5, countType = "all", returnAll = T, genes = somatic_gene_e21$Hugo_Symbol)

# Comparing two cohorts (MAFs)
E19VSE21 <- mafCompare(m1 = merged_19, m2 = merged_21, m1Name = 'E19', m2Name = 'E21', minMut = 5, useCNV = F, )
print(E19VSE21)
forestPlot(mafCompareRes = E19VSE21, pVal = 0.1) # ci sono 2 geni mutati in maniera significativa: CCND3, SMARCA4

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

lollipopPlot2(m1 = merged_19, 
              m2 = merged_21, gene = "SMARCA4", labPosSize = 0.6, 
              #m1_label = 790, m2_label = 790, 
              m1_name = 'E19', m2_name = 'E21', 
              domainBorderCol = "black", 
              roundedRect = T, 
              AACol1 = 'HGVSp_Short', legendTxtSize = 0.8,
              AACol2 = 'HGVSp_Short', 
              refSeqID = "NM_001128849", 
              showDomainLabel = F,
              colors = lolli_cols,
)

# Oncogenic Signaling Pathways
OncogenicPathways(maf = merged_maf, )
OncogenicPathways(maf = merged_19, )
OncogenicPathways(maf = merged_21, )

PlotOncogenicPathways(maf = merged_maf, pathways = c("RTK-RAS", "Cell_Cycle"),  
                      #sampleOrder = sample_tot
                      )
PlotOncogenicPathways(maf = merged_19, pathways = "RTK-RAS", 
                      #sampleOrder = sample_e19
)
PlotOncogenicPathways(maf = merged_21, pathways = "RTK-RAS", 
                      #sampleOrder = sample_e21
)

# Detecting cancer driver genes based on positional clustering
merged.sig = oncodrive(maf = merged_maf, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(merged.sig)
plotOncodrive(res = merged.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5, colCode = c('#F5A962', '#3C8DAD'))

# Drug-Gene Interactions
dgi = drugInteractions(maf = merged_maf, fontSize = 0.75)

# Heterogeneity in tumor samples
library("mclust")
EGFR105.het = inferHeterogeneity(maf = merged_maf, tsb = sample_e19[1:4], vafCol = 't_AF', minVaf = 0) # provo a analizzare soltanto 4 campioni presi dal gruppo e19
print(EGFR105.het)
plotClusters(clusters = EGFR105.het, )

# Mutational Signatures
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
egfr.tnm = trinucleotideMatrix(maf = merged_maf, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = egfr.tnm, maf = merged_maf, pVal = 0.2) # None of the samples are enriched for APOBEC. Nothing to plot

# Signature analysis
library('NMF')
egfr.sign = estimateSignatures(mat = egfr.tnm, nTry = 6) # Draw elbow plot to visualize and decide optimal number of signatures from above results.
# Best possible value is the one at which the correlation value on the y-axis drops significantly. In this case it appears to be at n = 3
egfr.sig = extractSignatures(mat = egfr.tnm, n = 3)
# Compare detected signatures to COSMIC Legacy or SBS signature database.

egfr.og30.cosm = compareSignatures(nmfRes = egfr.sig , sig_db = "legacy") # Compate against original 30 signatures 
egfr.v3.cosm = compareSignatures(nmfRes = egfr.sig, sig_db = "SBS") #Compate against updated version3 60 signatures 
# compareSignatures returns full table of cosine similarities against COSMIC signatures, 
# which can be further analysed. Below plot shows comparison of similarities of detected signatures against validated signatures.

install.packages('pheatmap')
library('pheatmap')
pheatmap::pheatmap(mat = egfr.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures") # non funziona perchè si blocca la sessione
# Finally plot signatures
maftools::plotSignatures(nmfRes = egfr.sig, title_size = 1.2, sig_db = "legacy", show_title = T)
maftools::plotSignatures(nmfRes = egfr.sig, title_size = 1.2, sig_db = "SBS")
# If you fancy 3D barpots, you can install barplot3d package and visualize the results with legoplot3d function.
