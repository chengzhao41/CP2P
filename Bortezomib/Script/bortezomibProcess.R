rm(list = ls()); # clears all workspace
gc()

library(GEOquery)
#load("E:/CELLINE_project/bortezomibAnalysis/bortezomibData/bortGeo.RData")
#load("/home/cheng/Dropbox/ccle_cgp_bortezomib/bortezomibData/bortGeo.RData")
load("~/Dropbox/SNF_DRUG_PROJECT/CLINICAL/bortezomibData/bortGeo.RData")

stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Bortezomib/Script")

exprDataU133a <- cbind(exprs(bortezomib_mas5[[1]]), exprs(bortezomib_mas5[[2]]))
# View(exprDataU133a)
p <- rbind(pData(phenoData(bortezomib_mas5[[1]])), pData(phenoData(bortezomib_mas5[[2]])))

bortIndex <- which(p$"characteristics_ch1.1" == "treatment = PS341")
dexIndex <- which(p$"characteristics_ch1.1" == "treatment = Dex")
library("hgu133a.db")
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
names(mapped_probes) <- as.character(x[mapped_probes])
affy2sym <- as.character(x[mapped_probes])
sym2affy <- names(affy2sym)
rownames(exprDataU133a) <- affy2sym[rownames(exprDataU133a)]

nonNARows <- which(!is.na(rownames(exprDataU133a)))
bortezomibExpr <- exprDataU133a[nonNARows, bortIndex]
bortezomibExpr = t(log2(bortezomibExpr))

response = as.vector(p$"characteristics_ch1.8")[bortIndex]
binaryResponse = rep(0, length(response))
binaryResponse[which(response == "PGx_Responder = NR")] = 0
binaryResponse[which(response == "PGx_Responder = R")] = 1
binaryResponse[which(response == "PGx_Responder = IE")] = NA

nonIEind = which(!is.na(binaryResponse))
binaryResponse = binaryResponse[nonIEind]



studyCode = as.vector(p$"characteristics_ch1")[bortIndex]
codeLabel = rep(0, length(studyCode))

codeLabel[which(studyCode == "studyCode = 40")] = 1
codeLabel[which(studyCode == "studyCode = 39")] = 2
codeLabel[which(studyCode == "studyCode = 25")] = 3
codeLabel[which(studyCode == "studyCode = 24")] = 4

bortezomibExpr = bortezomibExpr[nonIEind, ]
codeLabel = codeLabel[nonIEind]
library(sva)

bortezomib.patient_ComBat = ComBat(dat=t(bortezomibExpr), batch = codeLabel, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

save(bortezomib.patient_ComBat, binaryResponse, file = "../WS/bortezomib.patient.RData")
