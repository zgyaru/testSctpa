#data = read.csv('D:/0research/scRNA-seq/scTPA2/scTPA_local-v7/example/e1_scRNA_UMI_count.csv',
#                row.names = 1)
#data = as.matrix(data)
#kegg = getGMT('D:/0research/scRNA-seq/scTPA2/scTPA_local-v7/data/pathway/homo/c2.kegg.gmt')


###########################################
########### test calVision         ########
###########################################
#kegg_sig = readSignaturesInput('D:/0research/scRNA-seq/scTPA2/scTPA_local-v7/data/pathway/homo/c2.kegg.gmt')
#normExpr = getNormalizedCopySparse(as.matrix(data))
#sigsSpars = sigsToSparseMatrix(kegg_sig, normExpr$norm_data)
#sigsSpars = as.matrix(sigsSpars)
#sigMatch = sapply(kegg,function(x,y) na.omit(match(x, y)),ownames(data))
##sigMatch = do.call('rbind',sigMatch)
#sigsSpars2 = sigsToSparseMatrix(kegg, normExpr$norm_data)
##identical(as(sigsSpars,'sparseMatrix'), sigsSpars2)
#sigScores <- innerEvalSignatureBatchNorm(kegg, normExpr)

#vision_score = calVsion(data,kegg)


###########################################
########### test calZscore         ########
###########################################
#zscore_score = calZscore(data,kegg)
#zscore_score[1:5,1:5]


###########################################
########### test calPlage         ########
###########################################
#lage_score = calPlage(data,kegg)
#plage_score[1:5,1:5]

###########################################
########### test calSSgsea         ########
###########################################
#ssgsea_score = calSSgsea(data,kegg)
#ssgsea_score[1:5,1:5]

###########################################
########### test calGsva         ########
###########################################
#gsva_score = calGSVA(data,kegg)
#plage_score[1:5,1:5]
#aucell_score = calAUC(data,kegg)
#aucell_score[1:5,1:5]


#Rcpp.package.skeleton("sctpaBeta",
#                      example_code=FALSE,
#                      attributes=TRUE,
#                      cpp_files=c("ker.cpp","ks_test.cpp"))
