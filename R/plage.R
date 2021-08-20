#' @title plage functions
#'
#'
#'
#' @description plage function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return plage score matrix
#' @import stats
#' @import data.table
#' @export
#'
calPlage = function(data,
                    gSets){
  cell_names = colnames(data)
  if(is(data,'sparseMatrix')){
    data = as.matrix(data)
  }
  data = data.table::data.table(t(data))
  data[,(colnames(data)):=lapply(.SD,function(x) (x-mean(x))/sd(x))]
  data = t(data)
  colnames(data) = cell_names


  es = sapply(names(gSets), function(gSetName)
    rightsingularsvdvectorgset(gSet=gSets[[gSetName]],Z=data))
  if (length(gSets) == 1)
    es = matrix(es, ncol=1)
  es = t(es)
  colnames(es) = colnames(data)
  es
}


rightsingularsvdvectorgset = function(gSet, Z){
  gSet = unique(gSet)
  gSet = intersect(gSet, rownames(Z))
  nGenes = length(gSet)
  if(nGenes==1){
    return(rep(NA,ncol(Z)))
  }else{
    s = svd(Z[gSet, ])
    return(s$v[, 1])
  }
}

