###### cBioPortalData

# Bioconductor version: Release (3.13)
# The cBioPortalData package takes compressed resources from repositories such as cBioPortal and assembles a MultiAssayExperiment object with Bioconductor classes.
# Author: Levi Waldron [aut], Marcel Ramos [aut, cre]
# Maintainer: Marcel Ramos <marcel.ramos at roswellpark.org>
# Citation (from within R, enter citation("cBioPortalData")):
# Ramos M, Geistlinger L, Oh S, Schiffer L, Azhar R, Kodali H, de Bruijn I, Gao J, Carey VJ, Morgan M, Waldron L (2020). “Multiomic Integration of Public Oncology Databases in Bioconductor.” JCO Clinical Cancer Informatics, 958-971. doi: 10.1200/CCI.19.00119, PMID: 33119407, https://doi.org/10.1200/CCI.19.00119, https://doi.org/10.1200/CCI.19.00119.

BiocManager::install("cBioPortalData")
BiocManager::install("AnVIL") # The AnVIL is a cloud computing resource developed in part by the National Human Genome Research Institute
library(cBioPortalData)
library(AnVIL)

# To see what datasets are currently not building, we can look at the studiesTable dataset that comes with the package. We can invoke it using:
data("studiesTable", package = "cBioPortalData")
studiesTable[136:145,] 
studiesTable$cancer_study_id
# The last two columns will show the availability of each cancer_study_id for either download method 
# (pack_build for cBioDataPack and api_build for cBioPortalData).

# Choosing download method: cBioDataPack \\ cBioPortalData
# cBioDataPack: Obtain Study Data as Zipped Tarballs
## Use ask=FALSE for non-interactive use
luad <- cBioDataPack("luad_tcga", ask = FALSE)
luad

# cBioPortalData: Obtain data from the cBioPortal API

cbio <- cBioPortal()
tags(cbio)
getStudies(cbio) # Get the list of studies available
clinicalData(cbio, "luad_tcga")
mols <- molecularProfiles(cbio, "luad_tcga")
mols[["molecularProfileId"]]
genePanels(cbio)
getGenePanel(cbio, "ARCHER-SOLID-CV1")

luad_api <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "luad_tcga", genePanelId = "ARCHER-SOLID-CV1", molecularProfileIds = c("luad_tcga_rppa", "luad_tcga_mutations"))

## Example Analysis: Kaplan-Meier Plot
library(survival)
library(survminer)
table(colData(luad)$OS_STATUS)
class(colData(luad)$OS_MONTHS)

# Now, we clean the data a bit to ensure that our variables are of the right type for the subsequent survival model fit.
colluad <- colData(luad)
colluad[colluad$OS_MONTHS == "[Not Available]", "OS_MONTHS"] <- NA
colluad$OS_MONTHS <- as.numeric(colluad$OS_MONTHS)
colData(luad) <- colluad

# We specify a simple survival model using SEX as a covariate and we draw the K-M plot.
fit <- survfit(
  Surv(OS_MONTHS, as.numeric(substr(OS_STATUS, 1, 1))) ~ SEX,
  data = colData(luad)
)
ggsurvplot(fit, data = colData(luad), risk.table = TRUE) # si blocca

laml <- cBioDataPack("laml_tcga", ask = FALSE)
laml
table(colData(laml)$OS_STATUS)
class(colData(laml)$OS_MONTHS)
collaml <- colData(laml)
collaml[collaml$OS_MONTHS == "[Not Available]", "OS_MONTHS"] <- NA
collaml$OS_MONTHS <- as.numeric(collaml$OS_MONTHS)
colData(laml) <- collaml
fit <- survfit(
  Surv(OS_MONTHS, as.numeric(substr(OS_STATUS, 1, 1))) ~ SEX,
  data = colData(laml)
)
ggsurvplot(fit, data = colData(laml), risk.table = TRUE)