# Script to Process GEX using BrainArray CDF
library(gdata)
library(genefu)
library(parallel)
library(limma)
library(affy)

library(Biobase)
library(GEOquery)

# [HG-U133_Plus_2]
library(hgu133plus2hsensgcdf) #for EnsemblGeneId

process_CEL <- function (celFiles) {
  affyRaw <- ReadAffy(cdfname='hgu133plus2hsensgcdf', filenames=celFiles)
  print(affyRaw)
  
  eset <- expresso(affyRaw,
                   bgcorrect.method="rma",
                   normalize.method="quantiles",
                   pmcorrect.method="pmonly",
                   summary.method="medianpolish"
  )
  # Full Expression Matrix (without being part of an object)
  datamatrix <- exprs(eset)
  
  ### for Ensembl IDs 
  # filter for probes corresponding to genes
  datamatrix <- datamatrix[grep("ENSG", rownames(datamatrix)), ]
  # The CDF always represents the Ensembl IDs at ###_at, so you need to remove the "_at" from those IDs!
  # Polish the rownames (remove the _at from the Ensembl IDs)
  rownames(datamatrix) <- gsub(rownames(datamatrix), pattern="_at", replacement="")
  #print(rownames(datamatrix))  
  
  #   ### for Entrez IDs 
  #   # The CDF always represents the entrez IDs at ###_at, so you need to remove the "_at" from those IDs!
  #   # Polish the rownames (remove the _at from the Entrez IDs)
  #   rownames(datamatrix) <- gsub(rownames(datamatrix), pattern="_at", replacement="")
  #   #print(rownames(datamatrix))
  
  #   # Translate EntrezIDs to EnsemblIDs using biomaRt
  #   # !!! NOTE: general gene ID mapping is not unique, thus for different gene IDs use
  #   # appropriate customCDF like hgu133plus2hsensgcdf !!!
  #   library(biomaRt)
  #   genesEntrez <- rownames(datamatrix)
  #   ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  #   idtable <- getBM(filters="entrezgene",
  #                    attributes=c("ensembl_gene_id", "entrezgene"),
  #                    values=genesEntrez,
  #                    mart=ensembl)
  
  # Polish the data names
  colnames(datamatrix) <- sapply(strsplit(colnames(datamatrix), split=".CEL"),`[`,1)
  # Depending on dataset, sometimes need to remove suffix after '_'
  #colnames(datamatrix) <- sapply(strsplit(colnames(datamatrix), split="_"),`[`,1)
  
  return(datamatrix)
}

# take GEOid as an argument of the Rscript
#args <- commandArgs(TRUE)
#GEOid <- args[1]

# hardcode the GEOid
GEOid <- "Epirubicin/WS/GSE16446_RAW/"

# Get the CEL files list (of just downloaded files or those already saved)
celFiles <- list.celfiles(GEOid, full.names = TRUE)

# process with BrainArray
exprs_data <- process_CEL(celFiles)

desmedt2009_demo <- read.csv("Epirubicin/WS/desmedt2009_demo.csv", header=TRUE)
temp.order <- match(colnames(exprs_data), desmedt2009_demo$geo_accn)
stopifnot(desmedt2009_demo$geo_accn[temp.order] == colnames(exprs_data))

temp.labels <- desmedt2009_demo$pCR[temp.order]
names(temp.labels) <- desmedt2009_demo$geo_accn[temp.order]
temp.data <- exprs_data[, which(!is.na(temp.labels))]
temp.labels <- temp.labels[which(!is.na(temp.labels))]
stopifnot(names(temp.labels) == colnames(temp.data))

epirubicin.patient <- t(exprs_data)
# the response labels
binaryResponse = rep(FALSE, length(temp.labels))
binaryResponse[temp.labels == "YES"] = TRUE
names(binaryResponse) = names(temp.labels)
stopifnot(rownames(epirubicin.patient) == names(binaryResponse))
table(binaryResponse)
# FALSE  TRUE 
# 101    17 

# Save the expression matrix
save(epirubicin.patient, binaryResponse, file="Epirubicin/WS/epirubicin.patient.RData")

