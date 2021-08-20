#' @title functions to convert data types
#' @import data.table
#'
#'
#' @description Reading a .gmt file containing pathways or gene sets
#' @details Input a .gmt file path and return a list which including gene sets
#' @param gmt_file_path .gmt file path
#' @param n_gene_thre minimum number of genes in gene sets
#' @export
#' @example
#' gsets = getGMT('/path/kegg.gmt')
#'
getGMT = function(gmt_file_path,
                  n_gene_thre = 0){
  if (!file.exists(gmt_file_path)) {
    stop(paste0("Cannot find file: ", gmt_file_path))
  }
  paths = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(paths)){
    t = strsplit(paths[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    genes = unique(genes)
    if(length(genes)>n_gene_thre){
      gsets[[t[1]]] = genes
    }
  }
  return (gsets)
}


#' @description Writing a .gmt file containing pathways or gene sets
#' @details Input a gsets list and .gmt file path, write these gsets into gmt format
#' @importFrom utils write.table
#' @param gSets a list which containing gene sets
#' @param gmt_file_path .gmt file path
#' importFrom("stats", "na.omit", "sd")
#' importFrom("utils", "write.table")
#' @export
#' @example
#' toGMT(kegg.gsets,'/path/kegg.gmt')
#'
toGMT = function(gSets,
                 gmt_file_path){
  if(file.exists(gmt_file_path)){
    file.remove(gmt_file_path)
    warning(paste0("There are already exit the ",
                   gmt_file_path,
                   ". Automatically removing it"))
  }
  # replace "-", ",", "/", " ", "(" or ")" with "_"
  aa = gsub('-|,|/| |\\(|\\)','_',names(gSets))
  aa = gsub('_+','_',aa)
  aa = gsub('_\\b','',aa)
  aa = gsub(',','',aa)
  names(gSets) = aa
  for(i in 1:length(gSets)){
    genes = paste(gSets[[i]],collapse ='\t')
    aa = paste(names(gSets)[i],names(gSets)[i],genes,sep='\t')
    write.table(aa, file=gmt_file_path, row.names=F,col.names = F,quote=F,append = T)
  }
}



len_gSets = function(gSets){
  res = sapply(gSets,length)
  unlist(res)
}


# convert data to data.table format
convertData = function(data){
  data = data.table::data.table(data, keep.rownames = T)
  data = data.table::setkey(data, "rn")
  return(data)
}


#' @description rank genes (decreasing) with in cells according it's expression
#' @import data.table
#' @param data data.table matrix with row in genes and column in cells
#' @return ranked data in data.table format
#'
rankGenes = function(data){
  #print(class(data))
  colsNam = colnames(data)[-1]
  # Similar to rank but much faster
  # pay attention to the random......
  data[, (colsNam):=lapply(-.SD, data.table::frank, ties.method="random", na.last="keep"),
       .SDcols=colsNam]
  return(data)
}

#' @description order genes (decreasing) with in cells according it's expression
#' @import data.table
#' @param data data.table matrix with row in genes and column in cells
#' @return ordered data in data.table format
#'
orderGenes = function(data){
  #print(class(data))
  colsNam = colnames(data)[-1]
  # Similar to rank but much faster
  # pay attention to the random......
  data[, (colsNam):=lapply(-.SD, order),
       .SDcols=colsNam]
  return(data)
}


#' log2-scale transform a dense OR sparse matrix
#'
#' This avoids the creation of a dense intermediate matrix when
#' operating on sparse matrices
#'
#' Either performs result <- log2(spmat+1) or if scale = TRUE
#' returns result <- log2(spmat/colSums(spmat)*scaleFactor + 1)
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @importFrom methods is
#' @param spmat sparse Matrix
#' @param scale boolean - whether or not to scale the columns to sum to `scale_factor`
#' @param scaleFactor if scale = TRUE, columns are scaled to sum to this number
#' @return logmat sparse Matrix
matLog2 <- function(spmat, scale = FALSE, scaleFactor = 1e6) {
  if (scale == TRUE) {
    spmat <- t( t(spmat) / colSums(spmat)) * scaleFactor
  }
  if (is(spmat, "sparseMatrix")) {
    matsum = Matrix::summary(spmat)
    
    logx <- log2(matsum$x + 1)
    
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j,
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  } else {
    logmat <- log2(spmat + 1)
  }
  return(logmat)
}

#' Compute col-wise variance on matrix without densifying
#'
#' @importFrom Matrix colMeans
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom matrixStats colVars
#'
#' @param x expression matrix
#' @return numeric vector col-wise variance
colVarsSp <- function(x) {
  if (is(x, "matrix")) {
    out <- matrixStats::colVars(x)
    names(out) <- colnames(x)
  } else {
    rm <- Matrix::colMeans(x)
    out <- Matrix::colSums(x ^ 2)
    out <- out - 2 * Matrix::rowSums(t(x) * rm)
    out <- out + nrow(x) * rm ^ 2
    out <- out / (nrow(x) - 1)
  }
  return(out)
}


#' Utility method to load signatures into a sparse matrix
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats na.omit
#' @param gSets List of gene sets
#' @param expression numeric Matrix Genes x Cells
#' @return sparseMatrix containing signature matched values
sigsToSparseMatrix = function(gSets, expression) {
  
  sigMatches <- lapply(seq(length(gSets)), function(i) {
    sig = gSets[[i]]
    indices = na.omit(match(sig, rownames(expression)))
    values = rep(1,length(indices))
    ii = rep(i, length(indices))
    return(list(indices, ii, values))
  })
  
  j <- unlist(lapply(sigMatches, function(x) x[[1]]))
  i <- unlist(lapply(sigMatches, function(x) x[[2]]))
  x <- unlist(lapply(sigMatches, function(x) x[[3]]))
  
  dn = list( names(gSets),
              rownames(expression))
  
  sigSparseMatrix <- sparseMatrix(i = i, j = j, x = x,
                                  dims = c(length(gSets),
                                           nrow(expression)),
                                  dimnames = dn)
  return(sigSparseMatrix)
}


setAs("data.frame", "Matrix", function(from) {
  mat = do.call(cbind, lapply(from, as, "Matrix"))
  colnames(mat) <- colnames(from)
  rownames(mat) <- rownames(from)
  mat
})

#roxygen2::roxygenize(package.dir = ".")
