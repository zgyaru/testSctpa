#' @title GSVA functions
#'
#'
#'
#' @description GSVA function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return gsva score matrix
#' @import data.table
#' @importFrom stats na.omit
#' @export
#'
calGsva = function(data,
                   gSets,
                   alpha=0.25,
                   rnaseq=F){
  if(is(data,'sparseMatrix')){
    data = as.matrix(data)
  }
  data = convertData(data)
  num_genes = nrow(data)
  gset.idx.list = lapply(gSets,
                         function(x,y) na.omit(match(x, y)),
                         data$rn)

  gene.density = compute.gene.density(data[,-"rn"])
  rank.scores = rep(0, num_genes)


  O = orderGenes(data)
  R = apply(O[,-"rn"], 2, compute_rank_score, num_genes)

  m = t(sapply(gset.idx.list,
               ks_matrix_gsva,
               as.matrix(R),
               as.matrix(O[,-"rn"]),
               tau=1,
               mx_diff=1,
               abs_rnk=0))

  colnames(m) = colnames(data)[-1]

  m
}



compute.gene.density = function(expr,
                                rnaseq=FALSE,
                                kernel=TRUE){
  gene.density = NA
  gene.density = matrix_d(as.matrix(expr), as.integer(rnaseq))

  rownames(gene.density) = rownames(expr)
  colnames(gene.density) = colnames(expr)

  return(gene.density)
}

compute_rank_score = function(sort_idx_vec,
                              num_genes){
  tmp = rep(0, num_genes)
  tmp[sort_idx_vec] = abs(seq(from=num_genes,to=1) - num_genes/2)
  return (tmp)
}
