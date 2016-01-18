# based on script by Deena M.A. Gendoo
# Script to to Process GEX using BrainArray CDF

library(gdata)
library(genefu)
library(parallel)
library(limma)
library(affy)

library(Biobase)
library(GEOquery)
library(hgu133plus2hsentrezgcdf) # Hg U133 Plus 2.0 CDF

process_CEL_hgu133plus2 <- function (celFiles) {
  affyRaw <- ReadAffy(cdfname='hgu133plus2hsentrezgcdf',filenames=celFiles)
  print(affyRaw)
  
  eset <- expresso(affyRaw,
  				bgcorrect.method="rma",
  				normalize.method="quantiles",
  				pmcorrect.method="pmonly",
  				summary.method="medianpolish"
  				)
  
  # Full Expression Matrix (without being part of an object)
  datamatrix <- exprs(eset)
  
  # The CDF always represents the entrez IDs at ###_at, so you need to remove the "_at" from those IDs!
  # Polish the rownames (remove the _at from the Entrez IDs)
  rownames(datamatrix) <- gsub(rownames(datamatrix), pattern="_at", replacement="")
  print(rownames(datamatrix))
  
  colnames(datamatrix) <- sapply(strsplit(colnames(datamatrix), split=".CEL"),`[`,1)
  # Depending on dataset, sometimes need to remove suffix after '_'
  #colnames(datamatrix) <- sapply(strsplit(colnames(datamatrix), split="_"),`[`,1)
  
  return(datamatrix)
}

# take GEOid as an argument of the Rscript
#args <- commandArgs(TRUE)
#GEOid <- args[1]

# hardcode the GEOid
GEOid <- "GSE20711"

# Download the raw CEL files from supplement MAY TAKE LONG TIME!
getGEOSuppFiles(GEOid)
# Unpack the downloaded file and get the list of files
untar(paste0(GEOid, "/", GEOid, "_RAW.tar"), exdir=GEOid)
# Get the CEL files list (of just downloaded files or those already saved)
celFiles <- list.celfiles(GEOid, full.names = TRUE)

# process with BrainArray
exprs_data <- process_CEL_hgu133plus2(celFiles)

# Save the expression matrix
save(exprs_data, file=paste0(GEOid, "exprs.RData"))

