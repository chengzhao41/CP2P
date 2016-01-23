require(PharmacoGx)

## function to compute Slope measure for a PSet using computeSlope() from PharmacoGx
computeSlopeMeasure <- function (dataPset){
  slopes <- list()
  len <- dim(dataPset@sensitivity$raw)[[1]]
  for (i in 1:len) {
    slopes[i] <- computeSlope(dataPset@sensitivity$raw[i, , "Dose"], dataPset@sensitivity$raw[i, , "Viability"])
  }
  dataPset@sensitivity$profiles$slope_recomputed <- as.numeric(slopes)
  return(dataPset)
}