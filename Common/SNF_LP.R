SNF_LP <- function(train, test, groups, K=0, alpha=0.5, t=20, method=1) {  
  
  require(SNFtool)
  
  ### This function is used to predict the subtype of new patients.
  # train and test have the same number of view and the same number of columns
  # group is the label for the train data
  # K, alpha, t are the prameters for SNF. 
  # K is the number of neighbors
  # alpha is the hyperparameter used in constructing similarity network
  # t is the number of iterations
  # method is a indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
  
  stopifnot(K > 1)
  
  # Normalize each column of x to have mean 0 and standard deviation 1.
  Standard_Normalization = function(x) {
    x = as.matrix(x);
    mean = apply(x, 2, mean)
    sd = apply(x, 2, sd)
    sd[sd==0] = 1
    xNorm = t((t(x) - mean) / sd)
    return(xNorm)
  }  
  
  csPrediction <- function(W,Y0,method,NLabel){
    ###This function implements the label propagation to predict the label(subtype) for new patients.	
    ### note method is an indicator of which semi-supervised method to use
    # method == 0 indicates to use the local and global consistency method
    # method >0 indicates to use label propagation method.
    
    alpha=0.9;
    P= W/rowSums(W)
    if(method==0){
      Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
    } else {
      Y=Y0;
      for (i in 1:1000){
        Y=P%*%Y;
        Y[1:NLabel,]=Y0[1:NLabel,];
      }
    }
    return(Y);
  }
  
  view= (rbind(train,test));
  view = Standard_Normalization(view);
  Dist1 = SNFtool::dist2(view, view);
  W = SNFtool::affinityMatrix(Dist1, K);
  
  U = sort(unique(groups));
  groupU = groups;
  
  for (i in 1:length(U)) {
    groupU[which(groups == U[i])] = i;
  }
  
  Y0=matrix(0,nrow(view), max(groupU));
  for (i in 1:length(groups)) Y0[i,groupU[i]]=1;
  NLabel = nrow(train);
  Y=csPrediction(W,Y0,method, NLabel);
  if (length(U)>2){
    newgroups=rep(0,nrow(view));
    for (i in 1:nrow(Y)) newgroups[i]=which(Y[i,]==max(Y[i,]));
  }
  else{
    newgroups = (Y[,2]+0.00000001)/(Y[,1]+0.00000001)/2
  }
  
  return (list(newgroups, Y[, 2]))
}
