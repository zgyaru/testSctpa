#' examples for test
#' 
#' 
#' loading test data 
#' @return a sparse matrix
#' @export
#' 
load_counts = function(){
  folder = system.file("data",package = "testSctpa")
  if(folder == ""){
    stop("could not find test data directory, try-re-installing 'testSctpa'")
  }
  readRDS(file.path(folder,'test.rds'))
}

#' loading pathways or gene sets 
#' @param species species; human or mouse
#' @param pathway abbreviation for pathway database
#'                detials see https://github.com/zgyaru/PASBench
#' @return a list containing pathways or gene sets
#' @export
#' 
getPathways = function(species, pathway){
  gSets_path = system.file(file.path("gmtFiles",species,
                                     paste0(pathway,'.gmt')),
                           package = "testSctpa")
  if(gSets_path == ""){
    stop("could not find pathway file, try-re-installing 'testSctpa' or read a patway gmt file through 'getGMT' function.")
  }
  getGMT(gSets_path)
}