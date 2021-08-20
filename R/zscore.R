#' @title zscore functions
#'
#'
#'
#' @description zscore function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return Z score matrix
#' @import data.table
#' @export
#'
calZscore = function(data,
                  gSets) {
  cell_names = colnames(data)
  data = data.table::data.table(t(data))
  data[,(colnames(data)):=lapply(.SD,function(x) (x-mean(x))/sd(x))]
  data = t(data)
  colnames(data) = cell_names

  es = sapply(names(gSets), function(gSetName)
    combinez(gSet=gSets[[gSetName]], Z=data))

  if(length(gSets) == 1){
    es = matrix(es, nrow=1)
  }
  #rownames(es) = names(gSets)
  #colnames(es) = colnames(X)

  t(es)
}



combinez = function(gSet, Z){
  gSet = unique(gSet)
  nGenes = length(gSet)
  gSet = intersect(gSet, rownames(Z))

  if(nGenes == 1){
    Z[gSet,]
  }else{
    colSums(Z[gSet,]) / sqrt(length(nGenes))
  }
}

