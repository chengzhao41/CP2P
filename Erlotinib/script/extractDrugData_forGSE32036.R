extractDrugData_forGSE32036 <- function (drugID, dataPset, mDataType='rna') {
  ## custom function just for GSE32036

  
  ## get IC50 for the drug
  label.IC50 <- summarizeSensitivityProfiles(dataPset,
                                             drugs=drugID,
                                             sensitivity.measure="ic50_published",
                                             summary.stat="median")
  label.IC50 <- label.IC50[1, !is.na(label.IC50)]
  
  ## get expression data of cell lines with sensitivity measurement for the drug
  data.expression <- summarizeMolecularProfiles(GSE32036,
                                                # cell.lines=cc,
                                                mDataType='rna',
                                                verbose=FALSE)
  ## rename genes to EnsemblGeneId
  gInfo <- featureInfo(dataPset, mDataType)
  rna.all <- exprs(data.expression)
  class(rna.all) <- "numeric"
  dim(rna.all)

  # get the mapping of Entrez_Gene_ID -> ensembl
  require(org.Hs.eg.db)
  x <- org.Hs.egENSEMBL
  # Get the entrez gene IDs that are mapped to an Ensembl ID
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  
  # remove NA genes
  temp.entrez_gene_id <- gInfo$Entrez_Gene_ID
  temp.rna <- rna.all[!is.na(gInfo$Entrez_Gene_ID), ]
  temp.entrez_gene_id <-gInfo$Entrez_Gene_ID[!is.na(gInfo$Entrez_Gene_ID)]
  
  temp.ensembl <- xx[temp.entrez_gene_id]
  temp.ensembl <- unlist(temp.ensembl)
  stopifnot(sum(is.null(temp.ensembl)) == 0)
  
  temp.ensembl_assignment <- temp.ensembl[temp.entrez_gene_id]
  
  temp.rna <- temp.rna[!is.na(temp.ensembl_assignment), ]
  rownames(temp.rna) <- temp.ensembl_assignment[!is.na(temp.ensembl_assignment)]

  # remove samples for which no valid gene measurement exists
  # and remove any genes which has any NA values
  temp.rna2 <- temp.rna[, -which(colSums(is.na(temp.rna)) == dim(temp.rna)[1])]
  rna.all <- temp.rna2[rowSums(is.na(temp.rna2)) == 0, ]
   
  ## get list of cell lines
  cl.all <- colnames(rna.all)
  cl.IC50 <- intersect(names(label.IC50), cl.all)
 
  ## filter for cell lines with both expression and measurement data
  rna.IC50 <- rna.all[, cl.IC50]
  rna.IC50 <- rna.IC50[, colSums(is.na(rna.IC50))==0]
  cl.IC50 <- colnames(rna.IC50)
  label.IC50 <- label.IC50[cl.IC50]
 
  stopifnot(all(colnames(rna.IC50) == colnames(label.IC50)))
  
  return(list('rna.IC50'=rna.IC50, 'IC50.cont'=label.IC50))
}
