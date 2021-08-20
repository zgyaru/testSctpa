#' @title ssGSEA functions
#'
#'
#'
#' @description ssGSEA function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return ssgsea score matrix
#' @import data.table
#' @export
#'
calSSgsea = function(data,
                     gSets,
                     alpha=0.25,
                     normalization=T){
  if(is(data,'sparseMatrix')){
    data = as.matrix(data)
  }
  data = convertData(data)
  gset.idx.list = lapply(gSets,
                         function(x,y) na.omit(match(x, y)),
                         data$rn)


  R = rankGenes(data)
  O = orderGenes(data)

  m = t(sapply(gset.idx.list,
               ks_matrix_ssgsea,
               as.matrix(R[,-"rn"]),
               as.matrix(O[,-"rn"]),
               tau=alpha))
  colnames(m) = colnames(data)[-1]
  if (normalization) {
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    score_range = range(m)[2] - range(m)[1]
    m = m/score_range
  }
  m
}
