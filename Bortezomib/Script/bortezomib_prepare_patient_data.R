# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")
library(hgu133a.db) # unload "GEOquery" if this fails to load
library(jetset)

source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load Data ---------------------------------------------------------------
# uncomment this line to download the data directly from GEO.
# library(GEOquery)
# bortezomib_mas5 <- getGEO("GSE9782")
load("Bortezomib/WS/bortGeo.RData")

exprDataU133a <- cbind(exprs(bortezomib_mas5[[1]]), exprs(bortezomib_mas5[[2]]))
p <- rbind(pData(phenoData(bortezomib_mas5[[1]])), pData(phenoData(bortezomib_mas5[[2]])))

# get the index for bortezomib --------------------------------------------
bortIndex <- which(p$"characteristics_ch1.1" == "treatment = PS341")

# get the gene symbol -----------------------------------------------------
mapped_probes <- mappedkeys(hgu133aENSEMBL)
affy2ensg <- as.list(hgu133aENSEMBL[mapped_probes])

ensembl_vector <- unlist(affy2ensg)
names(ensembl_vector) <- gsub("_at[0-9]*", "_at", names(ensembl_vector))

best_probesets <- jmap("hgu133a", ensembl = ensembl_vector)
affyBest2ensg <- names(best_probesets)
names(affyBest2ensg) <- best_probesets

rownames(exprDataU133a) <- affyBest2ensg[rownames(exprDataU133a)]

# get the expression level ------------------------------------------------
NARows <- which(is.na(rownames(exprDataU133a)))
bortezomibExpr <- exprDataU133a[-NARows, bortIndex]
bortezomibExpr = t(log2(bortezomibExpr))
stopifnot(sum(is.na(colnames(bortezomibExpr))) == 0)
dim(bortezomibExpr)

# double check that there are no duplicate gene symbols ---------------------
temp.data = comGENE(bortezomibExpr, bortezomibExpr)
dim(temp.data[[1]])
stopifnot(dim(bortezomibExpr) == temp.data[[1]])

# get the binary response -------------------------------------------------
response = as.vector(p$"characteristics_ch1.8")[bortIndex]
binaryResponse = rep(0, length(response))
binaryResponse[which(response == "PGx_Responder = NR")] = 0
binaryResponse[which(response == "PGx_Responder = R")] = 1
binaryResponse[which(response == "PGx_Responder = IE")] = NA

# exclude NAs -------------------------------------------------------------
nonIEind = which(!is.na(binaryResponse))
binaryResponse = binaryResponse[nonIEind]
bortezomibExpr = bortezomibExpr[nonIEind, ]

# get the batch labels ----------------------------------------------------
studyCode = as.vector(p$"characteristics_ch1")[bortIndex]
codeLabel = rep(0, length(studyCode))

# for reference
codeLabel[which(studyCode == "studyCode = 40")] = 1
codeLabel[which(studyCode == "studyCode = 39")] = 2 # hybridized differently
codeLabel[which(studyCode == "studyCode = 25")] = 3
codeLabel[which(studyCode == "studyCode = 24")] = 4

codeLabel = codeLabel[nonIEind]
stopifnot(sum(codeLabel == 0) == 0)

# apply ComBat ------------------------------------------------------------
stopifnot(nrow(bortezomibExpr) == length(codeLabel))
show_pca(input_data = bortezomibExpr, label = codeLabel)
show_pca(input_data = bortezomibExpr, label = binaryResponse)

bortezomib.patient_ComBat = t(ComBat(dat=t(bortezomibExpr), batch = codeLabel, 
                                     mod=NULL, par.prior=TRUE, prior.plots=FALSE))
show_pca(input_data = bortezomib.patient_ComBat, label = codeLabel)
show_pca(input_data = bortezomib.patient_ComBat, label = binaryResponse)

# save the workspace file -------------------------------------------------
save(bortezomib.patient_ComBat, binaryResponse, file = "Bortezomib/WS/bortezomib.patient.RData")
