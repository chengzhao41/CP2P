# Script to Process GEX using BrainArray CDF
library(gdata)
library(genefu)
library(parallel)
library(limma)
library(affy)

library(Biobase)
library(GEOquery)

# HuGene-1_0-st
library(hugene10sthsensgcdf) #for EnsemblGeneId

process_CEL <- function (celFiles) {
  affyRaw <- ReadAffy(cdfname='hugene10sthsensgcdf', filenames=celFiles)
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
GEOid <- "Erlotinib/WS/GSE33072_RAW/"

# Get the CEL files list (of just downloaded files or those already saved)
celFiles <- list.celfiles(GEOid, full.names = TRUE)

# process with BrainArray
exprs_data <- process_CEL(celFiles)
dim(t(exprs_data))
# [1]   131 21637

# the response labels
# G2 indicates treatment using erlotinib only
# group names for all samples
# X marks samples that have not been treated with erlotinib
# G1 marks samples that have progress-free survival time < 2 months
# G2 marks samples that have progress-free survival time >= 2 months
temp.labels <- c("X","G2","X","X","G2","X","X","X","X","G2","G1","X","X","X","X",
         "X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X",
         "X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X",
         "X","X","X","X","X","G1","G1","X","X","G1","X","G2","X","G2","G1",
         "X","X","X","X","X","X","X","X","G1","X","X","X","X","G2","G1","G2",
         "X","X","X","X","X","X","G1","X","X","X","X","X","X","X","X","G1",
         "X","X","X","X","G2","G1","X","X","X","X","G2","X","X","X","G1",
         "G2","G2","X","X","X","X","X","X","G1","G1","X","X","X","G1","X",
         "X","X","X","X","X","X","X")
rownames(erlotinib.patient)[temp.labels == "G2"]

erlotinib.patient <- t(exprs_data[, temp.labels != "X"])
dim(erlotinib.patient)
#[1]    25 21637

binaryResponse <- temp.labels[which(temp.labels != "X")]
binaryResponse[which(binaryResponse == "G1")] = FALSE
binaryResponse[which(binaryResponse == "G2")] = TRUE
binaryResponse <- as.logical(binaryResponse)
names(binaryResponse) <- rownames(erlotinib.patient)
table(binaryResponse)
# FALSE  TRUE 
# 14    11 

stopifnot(nrow(erlotinib.patient) == length(binaryResponse))

# Save the expression matrix
save(erlotinib.patient, binaryResponse, file="Erlotinib/WS/erlotinib.patient.RData")

