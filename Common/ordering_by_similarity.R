# generates the random partitions with the test set common to all approaches
order_by_similarity <- function(
  data,
  labels,
  source_labels,
  K = 20,
  alpha = 0.5) 
{
  stopifnot(length(which(source_labels == "patient")) > 0)
  stopifnot(length(which(source_labels == "cgp")) > 0)
  stopifnot(dim(data)[1] == length(labels))
  stopifnot(length(labels) == length(source_labels))
  
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
  
  ## Compute the affinities 
  Dist1 <- dist2(data, data)
  A <- affinityMatrix(Dist1, K, alpha)
  
  # mean affinity to patients
  W <- rowMeans(A[, source_labels == 'patient'])
  
  # subset cell lines only and order them by mean affinity to patients
  W <- W[source_labels == 'cgp']
  
  similarity.order <- rev(order(W))
  return (similarity.order)
}
