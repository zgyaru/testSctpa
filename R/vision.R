#' @title VISION functions
#'
#'
#' @description calculate PAS by Vision function
#' @param data expression matrix or data.frame with row in genes and column in cells
#' @param gSets gene sets in list format
#' @return vision score matrix
#' @import data.table
#' @export
#'
calVision = function(data,
                    gSets){
  
  data = as.matrix(data)
  #latenSpace = computeLatentSpace(data)
  #clusters = clusterCells(object)
  normExpr_list = getNormalizedCopySparse(data)
  sigScores = innerEvalSignatureBatchNorm(gSets, normExpr_list)
  sigScores
}




getNormalizedCopySparse = function(data, func='znorm_columns') {
  ### just it....
  data <- matLog2(data)
  
  rowOffsets <- NULL
  colOffsets <- NULL
  rowScaleFactors <- NULL
  colScaleFactors <- NULL
  
  if (func == "znorm_columns") {
    colOffsets <- colMeans(data) * -1
    colScaleFactors <- colVarsSp(data) ** -0.5
    colScaleFactors[is.infinite(colScaleFactors)] <- 1
  }
  
  rowOffsets = numeric(nrow(data))
  rowScaleFactors = numeric(nrow(data)) + 1
  
  return(list(norm_data = data,
              rowOffsets = rowOffsets,
              colOffsets = colOffsets,
              rowScaleFactors = rowScaleFactors, 
              colScaleFactors = colScaleFactors))
}



#' Used in inner loop of batchSigEvalNorm
#'
#' Computes signature scores without inflating the genes/cells matrix
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#'
#' @param normData NormData row/column normalization factors
#' @param sigs List of Signature to be evalauting
#' @return matrix containing signature values (sigs x cells)
innerEvalSignatureBatchNorm = function(gSets,
                                       normExpr_list){
  norm_data = normExpr_list$norm_data
  rowScaleFactors = normExpr_list$rowScaleFactors
  colScaleFactors = normExpr_list$colScaleFactors
  rowOffsets = normExpr_list$rowOffsets
  colOffsets = normExpr_list$colOffsets
  
  sigSparseMatrix = sigsToSparseMatrix(gSets, norm_data)
  NCells = ncol(norm_data)
  NGenes = nrow(norm_data)
  Rs = Diagonal(x = rowScaleFactors)
  Cs = Diagonal(x = colScaleFactors)
  Rog = Matrix(rowOffsets, ncol = 1)
  Roc = Matrix(1, nrow = 1, ncol = NCells)
  
  SRs <- (sigSparseMatrix %*% Rs)
  SRsE <- SRs %*% norm_data
  SRsRo <- (SRs %*% Rog) %*% Roc
  
  # Note: this requires sparse=TRUE so OMP/MKL won't use many threads
  # for the next multiply (e.g., avoid multiplying dense x dense).  This
  # is important because this runs inside a parallel loop already
  Cog <- Matrix(1, ncol = 1, nrow = NGenes, sparse = TRUE)
  Coc <- Matrix(colOffsets, nrow = 1, sparse = TRUE)
  
  SCo <- (sigSparseMatrix %*% Cog) %*% Coc
  
  C <- (SRsE + SRsRo + SCo) %*% Cs
  
  sigScores <- as.matrix(C)
  
  denom <- rowSums(abs(sigSparseMatrix)) # denom is vector of length N_sigs
  
  sigScores <- sigScores / denom
  
  return(sigScores)
}

