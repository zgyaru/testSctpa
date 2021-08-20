#' @title AUCell functions
#'
#'
#'
#' @description calculate PAS by AUCell function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return AUC score matrix
#' @import data.table
#' @export
#'
calAUC = function(data,
                  gSets){
  # step1 rank genes in each cell
  data = convertData(data)
  data.ranking = rankGenes(data)
  aucMaxRank=ceiling(0.05*nrow(data.ranking))
  # step2 calculate AUC score
  gSetName = NULL
  # calculate each gene set
  aucMatrix = sapply(names(gSets), function(gSetName)
    .AUC.geneSet(gSet=gSets[[gSetName]], rankings=data.ranking,
                 aucMaxRank=aucMaxRank))
  aucMatrix = t(aucMatrix)
}



#' @description calculate a gene set auc score for all cells
#' @param gSet a gene set in character vector format
#' @param rankings ranked gene matrix in data.table format
#' @param aucMaxRank max auc score for normalization
#' @return auc score vector for a gene set

.AUC.geneSet = function(gSet, rankings, aucMaxRank){
  gSet = unique(gSet)
  nGenes = length(gSet)
  gSet = intersect(gSet, rankings$rn)
  #missing = nGenes-length(gSet)

  # filtering ranked gene matrix by genes in gene set
  gSetRanks = rankings[rn %in% gSet,]
  rm(rankings)

  ########### NEW version:  #######################
  x_th = 1:nrow(gSetRanks)
  x_th = sort(x_th[x_th<aucMaxRank])
  y_th = seq_along(x_th)
  maxAUC = sum(diff(c(x_th, aucMaxRank)) * y_th)
  ############################################

  # Apply by columns (i.e. to each ranking)
  auc = apply(gSetRanks[,-1], 2, .auc, aucMaxRank, maxAUC)

  return(auc)
}


#' @description calculate a gene set auc score for a cell
#' @param oneRanking ranked gene vector for a cell
#' @param aucMaxRank max auc score for normalization
#' @param maxAUC sum auc score
#' @return auc score vector for a gene set
.auc <- function(oneRanking, aucMaxRank, maxAUC)
{
  x = unlist(oneRanking)
  x = sort(x[x<aucMaxRank])
  y = seq_along(x)
  sum(diff(c(x, aucMaxRank)) * y)/maxAUC
}
