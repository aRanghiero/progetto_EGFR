# creo uno script per caricare uno studio di TCGA
# 21/09/21

library(R.utils)
library(maftools) # l'ho installato usando il file .tar.gz scaricato da BioConductur
library(tidyverse)
update.packages("maftools")

tcga_avail = tcgaAvailable()
head(tcga_avail, 50)

# carico lo studio LUAD che comprende tumori Lung_adenocarcinoma
# By default MAF from MC3 project will be loaded
LUAD = tcgaLoad(study = "LUAD", source = "M") # Samples 517
# Change the source to Firehose
LUAD_fh = tcgaLoad(study = "LUAD", source = "Firehose") # Samples 533
