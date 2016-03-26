comGENE <- function(data1, data2) {

  sym1 <- colnames(data1)
  sym2 <- colnames(data2)

ind1 = match(sym1, sym2, nomatch = NA_integer_, incomparables = NA) 

com = !is.na(ind1)
comSym = c()
comInd = c()
for (i in 1 : length(sym1))
{
  if (com[i])
  {
    comSym = c(comSym, sym1[i])   #Take the symbols with duplicates here
    comInd = c(comInd, i)           #Record the current index
  }
}

uniqueComSym = unique(comSym) ## Take the unique symbols

numUnique = length(uniqueComSym) ## Number of unique symbols

comdata1 = matrix(0, dim(data1)[1], numUnique)
comdata2 = matrix(0, dim(data2)[1], numUnique)
for( symInd in 1 : numUnique)
{
  curSym = uniqueComSym[symInd]
  
  data1ind = match(sym1, curSym)
  nonNadata1ind = which(!is.na(data1ind))
  
  if(length(nonNadata1ind) > 1)
  {
    comdata1[, symInd] = rowMeans(data1[, nonNadata1ind]) #Average over different probes for the same gene for each cell line
  }else
  {
    comdata1[, symInd] = data1[, nonNadata1ind]
  }
  
  
  
  data2ind = match(sym2, curSym)
  nonNadata2ind = which(!is.na(data2ind))
  if(length(nonNadata2ind) > 1)
  {
    comdata2[, symInd] = rowMeans(data2[, nonNadata2ind])
  }else
  {
    comdata2[, symInd] = data2[, nonNadata2ind]
  }
}
rownames(comdata1) = rownames(data1)
colnames(comdata1) = uniqueComSym

rownames(comdata2) = rownames(data2)
colnames(comdata2) = uniqueComSym

output = list(comdata1, comdata2)
return(output)
}