#' calculate PAS by nine tools
#'
#' parameters `counts`,  `gSets_path` and `gSets` are consistent across all tools
#'
#' other parameters see details of tools
#' 




cal_fscLVM = function(counts,
                      # counts matrix; recommend log-transform
                      gSets_path,
                      # *.gmt file path
                      type = 'counts',
                      n_hidden=1,
                      # number of hidden factors
                      ifLog=T,
                      # wether using log-transfrom
                      rand_seed = 123){
  
  ## matrix could be numeric or integer
  if(ncol(counts)>5000){
    if(ncol(counts)>10000){
      n_hidden = 3
    }else{
      n_hidden = 2
    }
  }else{
    n_hidden = 1
  }
  print(n_hidden)
  print(type)
  # no parallel
  if(type == 'counts'){
    if(ifLog){
      sc_mat = SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = log2(counts[,]+1))
      )
      slot = 'logcounts'
    }else{
      sc_mat = SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts[,])
      )
      slot = 'counts'
    }
  }else if(type == 'tpm'){
    if(ifLog){
      sc_mat = SingleCellExperiment::SingleCellExperiment(
        assays = list(logtpm = log2(counts[,]+1))
      )
      slot = 'logtpm'
    }else{
      sc_mat = SingleCellExperiment::SingleCellExperiment(
        assays = list(tpm = counts[,])
      )
      slot = 'tpm'
    }
  }
  print(slot)
  tryCatch({
    genesets = GSEABase::getGmt(gSets_path)
    # @must step1: construct model
    model = slalom::newSlalomModel(sc_mat,
                                   # a SingleCellExperiment object
                                   genesets,
                                   # a GeneSetCollection object
                                   n_hidden = n_hidden,
                                   # number of hidden factors, recommend 2-5
                                   min_genes = 1,
                                   # minum genes of every gene Sets cotained matrix genes
                                   assay_name = slot,
                                   # slot name of SingleCellExperiment
                                   prune_genes = T,
                                   # whether filtter out genes not annotated in any gene sets
                                   design = NULL,
                                   # numeric matrix for covariates
                                   anno_fpr = 0.01,
                                   # false positive rate (FPR) for assigning genes to factors (pathways)
                                   anno_fnr = 0.001,
                                   # false negative rate (FNR) for assigning genes to factors (pathways)
                                   verbose = F
    )
    print('model construction')
    # @must step2: initialize model
    model = slalom::initSlalom(model,
                               noise_model = 'gauss',
                               # just gauss now, future: support "hurdle", "poisson"
                               alpha_priors = NULL,
                               epsilon_priors = NULL,
                               design = NULL,
                               pi_prior = NULL,
                               n_hidden = NULL,
                               # number of hidden factors, required if pi_prior is not null
                               seed=rand_seed,
                               verbose=F
    )
    print('initial successful')
    # @must step3: train model
    model = slalom::trainSlalom(model,
                                nIterations = 1000,
                                # maximum number of iterations to use in training the model
                                minIterations = 300,
                                # minimum number of iterations to perform
                                shuffle=T,
                                # whether should the order in which factors are updated be shuffled between iterations
                                # generally helps speed up covergence
                                seed=rand_seed,
                                tolerance = 1e-05,
                                forceIterations = FALSE,
                                # whether should the model be forced to update nIterations times
                                pretrain = TRUE,
                                # should the model be "pre-trained" to achieve faster convergence
                                verbose = F,
                                drop_factors = TRUE
                                # should factors be dropped from the model if designed be be not relevant
    )
    if(model@.xData$converged){
      score = model$X_E1
      colnames(score) = model$termNames
      rownames(score) = colnames(counts)
      score = score[,(n_hidden+1):ncol(score)]
      score = t(score)
      return(score)
    }else{
      return("not converged")
    }
    
  },error = function(e){
    print(e)
    return("error")
  })
}



cal_AUCell = function(counts,
                      gSets,
                      n_cores
                      # number of cores used for parallel
){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    print(gc())
    # @must step1: rank genes
    ac_rankings = AUCell::AUCell_buildRankings(counts[,],
                                               # matrix, dgCMatrix, SummarizedExperiment, ExpressionSet
                                               nCores=n_cores,
                                               plotStats=FALSE,
                                               # plot the expression boxplots or histograms
                                               #assayName = NULL,
                                               # slot name of assay containing expression matrix
                                               verbose = F)
    # @must step2: calculate AUC scores
    print(gc())
    sc_AUC = AUCell::AUCell_calcAUC(gSets,
                                    ac_rankings,
                                    normAUC=T,
                                    # Wether to normalize the maximum possible AUC to 1
                                    aucMaxRank=ceiling(0.05 * nrow(ac_rankings)),
                                    # the number of genes (maximum ranking) that is used to computation
                                    # default: 5% of pathways; recommend: 1%-20%
                                    verbose = F
    )
    print(gc())
    score = AUCell::getAUC(sc_AUC)
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


cal_ROMA = function(counts,
                    gSets,
                    n_cores){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    gSets2 = del_geneSets_roma(gSets)
    roma = rRoma.R(ExpressionMatrix = data.matrix(counts),
                   ModuleList = gSets2,
                   ClusType = "FORK",
                   MinGenes = 1,   ## min: 4
                   #OutGeneSpace = NULL,
                   UseParallel = TRUE,
                   GeneOutThr = 1,
                   OutGeneNumber = 1,
                   nCores = n_cores,
                   ### do not change follow parameters
                   centerData=F,
                   # should the gene expression values be centered over the samples
                   ExpFilter=F,
                   # logical, should the samples be filtered
                   UseWeigths=F,
                   # logical, should the weigths be used for PCA calculation
                   DefaultWeight=1,
                   OutGeneSpace=NULL,
                   ApproxSamples=5,
                   nSamples=100,
                   Ncomp=10,
                   FixedCenter= TRUE,
                   GeneOutDetection = "L1OutExpOut",
                   GeneSelMode = "All",
                   SampleFilter = FALSE,
                   MoreInfo = FALSE,
                   PlotData = FALSE,
                   PCADims = 2,
                   PCSignMode = "none",
                   PCSignThr = NULL,
                   SamplingGeneWeights = NULL,
                   Grouping = NULL,
                   FullSampleInfo = FALSE,
                   GroupPCSign = FALSE,
                   CorMethod = "pearson",
                   PCAType = "DimensionsAreGenes")
    score = roma$SampleMatrix
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


cal_vision = function(counts,
                      gSets_path,
                      n_cores){
  
  
  ## matrix could be integer or numeric
  
  ### other parameters for Vsion:
  # proteinData = NULL, unnormalizedData = NULL, meta = NULL,
  # projection_genes = c("fano"), min_signature_genes = 5,
  # sig_gene_threshold = 0.001, threshold = 0.05, perm_wPCA = FALSE,
  # projection_methods = c("tSNE30"),
  # sig_norm_method = c("znorm_columns", "none", "znorm_rows",
  #                     "znorm_rows_then_columns", "rank_norm_columns"),
  # pool = "auto",
  # cellsPerPartition = 10, name = NULL, num_neighbors = NULL,
  # latentSpace = NULL, latentSpaceName = NULL,
  # latentTrajectory = NULL, pools = list()
  
  
  ## counts recommend for scale or normalized, but not transformed
  ## The expression data should not be log-transformed prior to loading into VISION.
  
  
  tryCatch({
    #if(ncol(counts)<100){
    #  projection_method = 'tSNE10'
    #}else{
    #  projection_method = 'tSNE30'
    #}
    print(gc())
    vis = VISION::Vision(counts,            ## Gene X Cell
                         # data.frame; sparseMatrix; dgeMatrix; ExpressionSet; SummarizedExperiment; Seurat
                         signatures = gSets_path,
                         projection_method = 'UMAP',
                         sig_gene_threshold=0)
    print(gc())
    options(mc.cores=n_cores)
    vis = VISION::analyze(vis)
    print(gc())
    score = t(vis@SigScores)    ## pathway X cell
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


cal_pagoda2 = function(counts,
                       gSets,
                       trim = 5,
                       n_cores){
  
  
  ### must be counts matrix !!!!!
  
  ### other parameters for knn.error.models
  # min.nonfailed = 5, min.count.threshold = 1,
  # max.model.plots = 50,
  # min.size.entries = 2000, min.fpm = 0, cor.method = "pearson",
  # verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE,
  # local.theta.fit = linear.fit, theta.fit.range = c(0.01, 100),
  # alpha.weight.power = 1/2
  
  ### other parameters for pagoda.varnorm
  # batch = NULL, prior = NULL,
  # fit.genes = NULL, minimize.underdispersion = FALSE,
  # n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9,
  # verbose = 0, weight.df.power = 1, smooth.df = -1,
  # theta.range = c(0.01, 100), gene.length = NULL
  print(gc())
  nPcs = min(round(ncol(counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  tryCatch({
    p2 = Pagoda2$new(counts, n.cores = n_cores,log.scale=F)
    print(gc())
    p2$adjustVariance(plot=F)
    print(gc())
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    print(gc())
    path_names = c()
    env = new.env(parent=globalenv())
    invisible(lapply(1:length(gSets),function(i) {
      genes = intersect(gSets[[i]],rownames(counts))
      name = paste0(names(gSets[i]),i)
      if(length(genes)>3){
        assign(name, genes, envir = env)
        path_names = c(path_names, name)
      }
    }))
    print(gc())
    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1)
    print(gc())
    path_names = names(p2@.xData$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
    rownames(score) = path_names
    colnames(score) = colnames(counts)
    for(i in 1:length(p2@.xData$misc$pwpca)){
      if(!is.null(p2@.xData$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(p2@.xData$misc$pwpca[[i]]$xp$scores)
      }
    }
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


cal_gsva = function(counts,
                    gSets,
                    n_cores){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    print(gc())
    score = GSVA::gsva(counts,
                       gSets,
                       method='gsva',
                       parallel.sz=n_cores,
                       verbose=T)
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


cal_ssgsea = function(counts,
                      gSets,
                      n_cores){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    print(gc())
    score = GSVA::gsva(counts,
                       gSets,
                       method='ssgsea',
                       parallel.sz=n_cores,
                       verbose=T)
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

cal_plage = function(counts,
                     gSets,
                     n_cores){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    print(gc())
    score = GSVA::gsva(counts,
                       gSets,
                       method='plage',
                       parallel.sz=n_cores,
                       verbose=T)
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

cal_zscore = function(counts,
                      gSets,
                      n_cores){
  
  ## matrix could be integer or numeric
  
  tryCatch({
    print(gc())
    score = GSVA::gsva(counts,
                       gSets,
                       method='zscore',
                       parallel.sz=n_cores,
                       verbose=T)
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

