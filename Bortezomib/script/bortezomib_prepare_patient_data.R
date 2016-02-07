# Load Libraries ----------------------------------------------------------
library("sva")

setwd("~/Git/CP2P/")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')

# Load Data ---------------------------------------------------------------
# uncomment these 2 line to download the data directly from GEO.
# library(GEOquery)
# bortezomib_mas5 <- getGEO("GSE9782")
load("Bortezomib/WS/bortGeo.RData")

exprDataU133a <- cbind(exprs(bortezomib_mas5[[1]]), exprs(bortezomib_mas5[[2]]))
p <- rbind(pData(phenoData(bortezomib_mas5[[1]])), pData(phenoData(bortezomib_mas5[[2]])))

# get the index for bortezomib --------------------------------------------
bortIndex <- which(p$"characteristics_ch1.1" == "treatment = PS341")

# get the gene symbol -----------------------------------------------------
HGU133A_Hs_ENSG_mapping <- read.delim("GatherData/BrainArray/HGU133A_Hs_ENSG_20.0.0/HGU133A_Hs_ENSG_mapping.txt")

affy2ensembl <- gsub("_at", "", HGU133A_Hs_ENSG_mapping$Probe.Set.Name)
names(affy2ensembl) <- HGU133A_Hs_ENSG_mapping$Affy.Probe.Set.Name
affy2ensembl <- affy2ensembl[which(!duplicated(affy2ensembl))]
rownames(exprDataU133a) <- affy2ensembl[rownames(exprDataU133a)]

# get the expression level ------------------------------------------------
NARows <- which(is.na(rownames(exprDataU133a)))
print(length(NARows)) # 10463
bortezomibExpr <- exprDataU133a[-NARows, bortIndex]
bortezomibExpr = t(log2(bortezomibExpr))
stopifnot(sum(is.na(rownames(bortezomibExpr))) == 0)
dim(bortezomibExpr) # 188 11820

# double check that there are no duplicate gene symbols ---------------------
temp.data = comGENE(bortezomibExpr, bortezomibExpr)
dim(temp.data[[1]]) # 188 11820
stopifnot(dim(bortezomibExpr) == dim(temp.data[[1]]))

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
