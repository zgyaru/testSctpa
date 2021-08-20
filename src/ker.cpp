#include <Rcpp.h>
using namespace Rcpp;

void initCdfs(void);
double precomputedCdf(double x, double sigma);

#define SIGMA_FACTOR 4.0
#define PRECOMPUTE_RESOLUTION 10000
#define MAX_PRECOMPUTE 10.0

double precomputed_cdf[PRECOMPUTE_RESOLUTION+1];
int is_precomputed = 0;

double sd(NumericVector x, int n) {
  int         i, n1;
  double      mean, sd;
  long double sum = 0.0;
  long double tmp;

  for (i=0; i < n; i++)
    sum += x[i];
  tmp = sum / n;
  if (R_FINITE((double) tmp)) {
    sum = 0.0;
    for (i=0; i < n; i++)
      sum += x[i] - tmp;
    tmp = tmp + sum / n;
  }
  mean = tmp;
  n1 = n - 1;

  sum = 0.0;
  for (i=0; i < n; i++)
    sum += (x[i] - mean) * (x[i] - mean);
  sd = sqrt((double) (sum / ((long double) n1)));

  return(sd);
}

void initCdfs(void){
  double divisor = PRECOMPUTE_RESOLUTION * 1.0;
  for(int i = 0; i <= PRECOMPUTE_RESOLUTION; ++i)
    precomputed_cdf[i] = R::pnorm(MAX_PRECOMPUTE * ((double) i) / divisor, 0.0, 1.0, TRUE, FALSE);
  /* standard normal distribution function, lower.tail=TRUE, log.p=FALSE */
}

inline double precomputedCdf(double x, double sigma){
  double v = x / sigma;
  if(v < (-1 * MAX_PRECOMPUTE)){
    return 0;
  }else if(v > MAX_PRECOMPUTE){
    return 1;
  }else{
    double cdf = precomputed_cdf[(int)(fabs(v) / MAX_PRECOMPUTE * PRECOMPUTE_RESOLUTION)];
    if(v < 0){
      return 1.0 - cdf;
    }else{
      return cdf;
    }
  }
}


// [[Rcpp::export]]
NumericVector row_d(NumericVector x, int rnaseq){
  int size_density_n = x.size();
  //Rprintf("%d\n",size_density_n);
  int size_test_n = size_density_n;
  NumericVector r(size_density_n);
  double bw = rnaseq ? 0.5 : (sd(x, size_density_n) / SIGMA_FACTOR);
  //Rprintf("%d\n",is_precomputed);

  if (!rnaseq && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
    //Rprintf("%d\n",is_precomputed);
  }

  for(int j = 0; j < size_test_n; ++j){
    double left_tail = 0.0;

    for(int i = 0; i < size_density_n; ++i){
      left_tail += rnaseq ? R::ppois(x[j], x[i]+bw, TRUE, FALSE) : precomputedCdf(x[j]-x[i], bw);
    }
    left_tail = left_tail / size_density_n;
    r[j] = -1.0 * log((1.0-left_tail)/left_tail);
  }
  return(r);
}

// [[Rcpp::export]]
NumericMatrix matrix_d(NumericMatrix expr, int rnaseq){
  int n_genes = expr.nrow();
  NumericMatrix m(n_genes,expr.ncol());
  for(int i = 0; i < n_genes; i++){
    m(i,_) = row_d(expr(i,_), rnaseq);
  }
  return(m);
}
