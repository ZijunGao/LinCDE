#include <Rcpp.h>
using namespace Rcpp;

//' LinCDECdf
//'
//' This function computes the cdf. for a LinCDE boosting model in C++.
//'
//' @param cellProb matrix of cell probabilities, of dimension nobs x number of bins.
//' @param y values to compute cdf. at (length nobs).
//' @param splitPointY split points of responses (length number of bins + 1).
//'
//' @return This function return a vector of cdf. (length nobs).
//'
//[[Rcpp::export]]
NumericVector LinCDECdf(const NumericMatrix& cellProb,
                              const NumericVector& y,
                              const NumericVector& splitPointY){
  // preprocessing
  int n = cellProb.nrow();
  int numberBin = cellProb.ncol();
  // initialization
  NumericVector cdfs(n);

  for (int i = 0; i < n; i++){
    double probTemp = 0;
    for(int j = 0; j < numberBin; j++){
      if(y[i] >= splitPointY[j+1]){
        probTemp += cellProb(i, j);
      } else {
        cdfs[i] += probTemp + cellProb(i, j) * (y[i] - splitPointY[j])/(splitPointY[j+1] - splitPointY[j]);
        break;
      }
      if(j == (numberBin-1)){
        cdfs[i] = 1;
      }
    }
  }
  return cdfs;
}
