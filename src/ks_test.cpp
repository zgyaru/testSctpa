#include <Rcpp.h>
#include <iostream>

#include <typeinfo>
using namespace std;
using namespace Rcpp;



double ks_sample_gsva(IntegerVector x, 
                      IntegerVector x_sort_indxs,
                      IntegerVector geneset_idxs,
                      IntegerVector geneset_mask,
                      int n_genes,
                      int n_geneset,
                      double dec,
                      double tau, 
                      int mx_diff, 
                      int abs_rnk){
  

  //NumericVector walk(n_genes);

  double sum_gset = 0.0;
  //Rprintf("%f\n", dec);
  for(int i = 0; i < n_geneset; i++){
    //Rprintf("%f\n",pow(x[geneset_idxs[i]-1], tau));
    sum_gset += pow(x[geneset_idxs[i]-1], tau);
  }
  //Rprintf("%f\n", sum_gset);
  
  //double mx_value = 0.0;
  double mx_value_sign = 0.0;
  double cum_sum = 0.0;
  
  double mx_pos = 0.0;
  double mx_neg = 0.0;
  
  int idx;
  for(int i = 0; i < n_genes; i++){
    idx = x_sort_indxs[i]-1;
    
    if(geneset_mask[idx] == 1){
      //Rprintf("%d\t", i);
      cum_sum += pow(x[idx], tau) / sum_gset;
    }else{
      cum_sum -= dec;
    }
    //walk[i] = cum_sum;
    
    if(cum_sum > mx_pos){ mx_pos = cum_sum; }
    if(cum_sum < mx_neg){ mx_neg = cum_sum; }
  }
  
  if (mx_diff != 0) {
    mx_value_sign = mx_pos + mx_neg;
    if (abs_rnk != 0)
      mx_value_sign = mx_pos - mx_neg;
  } else {
    mx_value_sign = (mx_pos > fabs(mx_neg)) ? mx_pos : mx_neg;
  }
  return mx_value_sign;
}

// [[Rcpp::export]]
NumericVector ks_matrix_gsva(IntegerVector geneset_idxs,
                             IntegerMatrix expr,
                             IntegerMatrix sort_idxs,
                             double tau, 
                             int mx_diff, 
                             int abs_rnk){
  
  int n_genes = expr.nrow();
  int n_samples = expr.ncol();
  int n_geneset = geneset_idxs.size();
  IntegerVector geneset_mask(n_genes);
  NumericVector es_geneset(n_samples);
  for(int i = 0; i < n_genes; i++){
    geneset_mask[i] = 0;
  }
  
  for(int i = 0; i < n_geneset; i++){
    geneset_mask[geneset_idxs[i]-1] = 1;
  }
  
  double dec = 1.0 / (n_genes - n_geneset);
  for(int i = 0; i < n_samples; i++){
    //Rprintf("%s\t", typeid(expr(_,i)).name());
    //int i=0;
    //IntegerVector x(n_genes);
    //x = expr(_,i);
    //IntegerVector x_sort_idxs(n_genes);
    //x_sort_idxs = sort_idxs(_,i);
    //Rprintf("%f\n",x[2]);
    es_geneset[i] = ks_sample_gsva(expr(_,i), sort_idxs(_,i), geneset_idxs,geneset_mask, 
                                   n_genes,n_geneset,dec,
                                   tau, mx_diff, abs_rnk);
  }
  return es_geneset;
}

// [[Rcpp::export]]
double ks_sample_ssgsea(IntegerVector x, 
                        IntegerVector x_sort_indxs, 
                        IntegerVector geneset_idxs, 
                        IntegerVector geneset_mask,
                        int n_genes,
                        int n_geneset,
                        double dec,
                        double tau){
  

  double walk_sum=0;

  double sum_gset = 0.0;

  for(int i = 0; i < n_geneset; i++){
    //Rprintf("%f\n",pow(x[geneset_idxs[i]-1], tau));
    sum_gset += pow(x[geneset_idxs[i]-1], tau);
  }

  double cum_sum = 0.0;

  
  int idx;
  for(int i = 0; i < n_genes; i++){
    idx = x_sort_indxs[i]-1;
    if(geneset_mask[idx] == 1){
      cum_sum += pow(x[idx], tau) / sum_gset;
    }else{
      cum_sum -= dec;
    }
    walk_sum += cum_sum;
  }

  return walk_sum;
}

// [[Rcpp::export]]
NumericVector ks_matrix_ssgsea(IntegerVector geneset_idxs,
                               IntegerMatrix expr,
                               IntegerMatrix sort_idxs,
                               double tau){
  
  int n_genes = expr.nrow();
  int n_samples = expr.ncol();
  int n_geneset = geneset_idxs.size();
  IntegerVector geneset_mask(n_genes);
  NumericVector es_geneset(n_samples);
  for(int i = 0; i < n_genes; i++){
    geneset_mask[i] = 0;
  }
  
  for(int i = 0; i < n_geneset; i++){
    geneset_mask[geneset_idxs[i]-1] = 1;
  }
  
  double dec = 1.0 / (n_genes - n_geneset);
  for(int i = 0; i < n_samples; i++){
    es_geneset[i] = ks_sample_ssgsea(expr(_,i), sort_idxs(_,i), 
                                     geneset_idxs,geneset_mask, 
                                     n_genes,n_geneset,dec,tau);
  }
  return es_geneset;
}

