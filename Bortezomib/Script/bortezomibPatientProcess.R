# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")
library(hgu133a.db) # unload "GEOquery" if this fails to load

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
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
names(mapped_probes) <- as.character(x[mapped_probes])
affy2sym <- as.character(x[mapped_probes])
sym2affy <- names(affy2sym)
rownames(exprDataU133a) <- affy2sym[rownames(exprDataU133a)]

# get the expression level ------------------------------------------------
nonNARows <- which(!is.na(rownames(exprDataU133a)))
bortezomibExpr <- exprDataU133a[nonNARows, bortIndex]
bortezomibExpr = t(log2(bortezomibExpr))

# deal with duplicate gene symbols by taking the mean ---------------------
temp.data = comGENE(bortezomibExpr, bortezomibExpr)
dim(temp.data[[1]])
bortezomibExpr = temp.data[[1]]

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
