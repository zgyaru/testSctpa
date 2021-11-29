
#' calculation Pathway Activity Score
#'
#' parameters `counts`,  `gSets_path` and `gSets` are consistent across all tools
#'
#' other parameters see details of tools
#'
#'
#'
#' @param counts gene expresion matrix, rows are genes and cols are cells
#' @param tool tool name
#' @param gmt_file pathways in GMT format
#' @param species species; human or mouse
#' @param pathway abbreviation for pathway database
#' @param filter whether filtering for genes expressed in less than 5 percent cells
#' @param normalize normalization method
#'

#' @export
cal_PAS = function(seurat_object,
                   tool = c('AUCell',
                            'Vision',
                            'GSVA','ssGSEA',
                            'plage','zscore'),
                   species = 'none',
                   pathway = 'none',
                   gmt_file = 'none',
                   normalize = 'sctransform',
                   n_cores = 3,
                   rand_seed = 123){
  
  if(gmt_file == 'none'){
    if(species == 'none' || pathway == 'none'){
      stop("please define species and pathway when do not assign gmt_file")
    }
    
    gSets_path = system.file(file.path("gmtFiles",species,
                                       paste0(pathway,'.gmt')),
                             package = "testSctpa")
  }else{
    if(species != 'none' || pathway != 'none'){
      warning(" 'gmt_file' is already present.
              Ignoring 'species' and 'pathway'.")
    }
    gSets_path = gmt_file
  }
  
  if(tool %in% c('pagoda2','vision')){
    warning(" 'log' transform will counter error for 'pagoda2' or 'vision'.
            force the 'normalize' to 'none' or you could modify your code
            to call other normalization function.")
    normalize = 'none'
  }
  
  score = cal_all_tools(seurat_object,
                        gSets_path,
                        tools = tool,
                        normalize = normalize,
                        n_cores = n_cores,
                        rand_seed = rand_seed
  )
  
  PAS = CreateAssayObject(data = score[[tool]])
  seurat_object[['PAS']] = PAS
  DefaultAssay(seurat_object) = 'PAS'
  warning(class(seurat_object))
  return(seurat_object)
}







#' calculation seven tools
#'
#' parameters `counts`,  `gSets_path` and `gSets` are consistent across all tools
#'
#' other parameters see details of tools
#'
#'
#'
#' @param counts gene expresion matrix, rows are genes and cols are cells
#' @param gSets_path patwhays/gene sets in clasical GMT format
#' @param tool select PAS tools
#' @param filter whether filtering for genes expressed in less than 5% cells
#' @param normalize normalization method

cal_all_tools = function(seurat_object,
                         gSets_path,
                         #cells_label,
                         tools = c('AUCell',
                                   'Vision',
                                   'GSVA','ssGSEA',
                                   'plage','zscore'),
                         normalize = c('log','CLR','RC','scran','none'),
                         mat_type = 'counts',
                         n_cores = 3,
                         rand_seed = 123){
  
  
  
  ## normalization
  tryCatch({
    if(normalize == 'log'){
      seurat_object = Seurat::NormalizeData(seurat_object,normalization.method = 'LogNormalize',verbose=0)
      #seurat_object = Seurat::ScaleData(seurat_object)
    }else if(normalize == 'CLR'){
      seurat_object = Seurat::NormalizeData(seurat_object,normalization.method = 'CLR',verbose=0)
      #seurat_object = Seurat::ScaleData(seurat_object)
    }else if(normalize == 'RC'){
      seurat_object = Seurat::NormalizeData(seurat_object,normalization.method = 'RC',verbose=0)
      #seurat_object = Seurat::ScaleData(seurat_object)
    }else if(normalize == 'scran'){
      counts = GetAssayData(object = seurat_object, assay = "RNA", slot = "counts")
      if(mat_type == 'counts'){
        sc = SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = counts)
        )
        if(ncol(counts)>300){
          clusters = scran::quickCluster(sc, assay.type = "counts")
          sc = scran::computeSumFactors(sc, clusters=clusters, assay.type = "counts")
        }else{
          sc = scran::computeSumFactors(sc,assay.type = "counts")
        }
        
        sc = scater::normalize(sc,exprs_values  ='counts')
        seurat_object@assays$data = sc@assays$data$logcounts
      }else{
        counts = GetAssayData(object = seurat_object, assay = "RNA", slot = "counts")
        sc = SingleCellExperiment::SingleCellExperiment(
          assays = list(tpm = counts)
        )
        if(ncol(counts)>300){
          clusters = scran::quickCluster(sc, assay.type = "tpm")
          sc = scran::computeSumFactors(sc, clusters=clusters, assay.type = "tpm")
        }else{
          sc = scran::computeSumFactors(sc,assay.type = "tpm")
        }
        
        sc = scater::normalize(sc,exprs_values  ='tpm')
        seurat_object@assays$data = sc@assays$data$logcounts
      }
      
      
    }else if(normalize == 'sctransform'){
      seurat_object = Seurat::SCTransform(seurat_object,verbose=FALSE)
    }else if(normalize == 'scnorm_9'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=9,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
      },error=function(e){
        print("scnorm error")
        print(e)
        return("error")
      })
    }else if(normalize == 'scnorm_5'){
      tryCatch({
        counts = na.omit(counts)
        DataNorm = SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K=5,NCores=n_cores)
        counts = DataNorm@assays$data$normcounts
        rm(DataNorm)
      },error=function(e){
        print("scnorm error")
        print(e)
        return("error")
      })
      
    }},error=function(e){
      print("normalize error")
      return("error")
    })
  cat("normalize success\n")
  
  
  
  eval_tools = vector(mode="list")
  
  n_cores = n_cores
  gSets = getGMT(gSets_path)
  
  counts = as.matrix(counts)
  for(i in 1:length(tools)){
    tool = tools[i]
    t_start = Sys.time()
    gc()
    score = switch(tool,
                   AUCell = cal_AUCell(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                       gSets,
                                       n_cores),  #0,0.8
                   pagoda2 = cal_pagoda2(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                         gSets,
                                         n_cores=n_cores),  #-13,65
                   fscLVM = cal_fscLVM(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                       gSets_path,
                                       type = mat_type),
                   Vision = cal_vision(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                       gSets_path,
                                       n_cores),  #-0.6895973, 3.32
                   ROMA = cal_ROMA(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                   gSets,
                                   n_cores),   #-0.07, 0.4
                   GSVA = cal_gsva(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                   gSets,
                                   n_cores), #-0.98,0.94
                   ssGSEA = cal_ssgsea(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                       gSets,
                                       n_cores),#-0.5,0.5
                   plage = cal_plage(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                     gSets,
                                     n_cores),  #-1,1
                   zscore = cal_zscore(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                       gSets,
                                       n_cores),   #-4,29
                   seurat = cal_SeuratScore(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"),
                                            gSets),  #-1753,6958
                   scSigScore = cal_scSigScore(),
                   gene = counts)
    
    
    eval_tools[[i]] = score
    
  }
  names(eval_tools) = tools
  eval_tools
}

