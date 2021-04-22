#include <Rcpp.h>
using namespace Rcpp;

// function: compute cdf for LinCDE (C++)
// input:
// cellProb: matrix of cell Probabilities, each row sum to one (n by numberBin)
// y: values to compute cdf at (length n)
// splitPointY: split points of responses (length numberBin+1)

// output:
//   cdf: vector of cdf (length n)

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
