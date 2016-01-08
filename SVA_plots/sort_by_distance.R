rootdir <- "~/Dropbox/CP2P"
setwd(rootdir)

load("Bortezomib/WS/bortezomib_data.RData");

#### following functions are from SNF codebase
dist2 <- function(X,C) {
  ndata = nrow(X)
  ncentres = nrow(C)
  
  sumsqX = rowSums(X^2)
  sumsqC = rowSums(C^2)
  
  XC = 2 * (X %*% t(C))
  
  res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
  res[res < 0] = 0
  return(res)
}

affinityMatrix <- function(Diff,K=20,sigma=0.5) {
  ###This function constructs similarity networks.
  N = nrow(Diff)
  
  Diff = (Diff + t(Diff)) / 2
  diag(Diff) = 0;
  sortedColumns = as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }
  means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
  
  avg <- function(x,y) ((x+y)/2)
  Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
  densities = dnorm(Diff,0,sigma*Sig,log = FALSE)
  
  W = (densities + t(densities)) / 2
  return(W)
}


# bortezomib$combined_slope.sva
# bortezomib.labels$slope_combined
# bortezomib.labels$slope_combined.source

## Compute the affinities 
K <- 20
alpha <- 0.5
data <- scale(bortezomib$combined_slope.sva)
Dist1 <- dist2(data, data)
A <- affinityMatrix(Dist1, K, alpha)

# mean affinity to patients
W <- rowMeans(A[, bortezomib.labels$slope_combined.source == 'patient'])

# subset cell lines only and order them by mean affinity to patients
W <- W[bortezomib.labels$slope_combined.source == 'cgp']
CL.data <- bortezomib$combined_slope.sva[bortezomib.labels$slope_combined.source == 'cgp', ]
CL.labels <- bortezomib.labels$slope_combined[bortezomib.labels$slope_combined.source == 'cgp']

similarity.order <- rev(order(W))
CL.data <- CL.data[similarity.order, ]
CL.labels <- CL.labels[similarity.order]

# select top 100 similar cell lines
CL.data <- CL.data[1:100, ]
CL.labels <- CL.labels[1:100]

# join with patients
new.source <- c(rep('cgp', length(CL.labels)), rep('patient', sum(bortezomib.labels$slope_combined.source == 'patient')))
new.data <- rbind(CL.data, bortezomib$combined_slope.sva[bortezomib.labels$slope_combined.source == 'patient', ])
new.labels <- c(CL.labels, bortezomib.labels$slope_combined[bortezomib.labels$slope_combined.source == 'patient'])



