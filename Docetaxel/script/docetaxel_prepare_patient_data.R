# Script to Process GEX using BrainArray CDF
library(gdata)
library(genefu)
library(parallel)
library(limma)
library(affy)

library(Biobase)
library(GEOquery)

# HG_U95Av2]
library(hgu95av2hsensgcdf) #for EnsemblGeneId

process_CEL <- function (celFiles) {
  affyRaw <- ReadAffy(cdfname='hgu95av2hsensgcdf', filenames=celFiles)
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
GEOid_res <- "Docetaxel/WS/GSE349_RAW/"
GEOid_sens <- "Docetaxel/WS/GSE350_RAW/"

# Get the CEL files list (of just downloaded files or those already saved)
celFiles <- list.celfiles(c(GEOid_res, GEOid_sens), full.names = TRUE)

# process with BrainArray
exprs_data <- process_CEL(celFiles)
docetaxel.patient <- t(exprs_data)

# the response labels
celFiles[1:14]
celFiles[15:24]
binaryResponse <- c(rep(FALSE, 14), rep(TRUE, 10))

stopifnot(nrow(docetaxel.patient) == length(binaryResponse))

# Save the expression matrix
save(exprs_data, binaryResponse, file="Docetaxel/WS/docetaxel.patient.RData")

