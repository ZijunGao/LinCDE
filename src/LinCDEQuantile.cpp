#include <Rcpp.h>
using namespace Rcpp;

//' LinCDEQuantiles
//'
//' This function computes the quantiles for a LinCDE boosting model in C++.
//'
//' @param cellProb matrix of cell probabilities, of dimension nobs x number of bins.
//' @param probs levels of quantiles.
//' @param splitPointY split points of responses (length number of bins + 1).
//'
//' @return This function return a matrix of quantiles (dimension nobs x  number of quantile levels).
//' @export
//[[Rcpp::export]]
NumericMatrix LinCDEQuantiles(const NumericMatrix& cellProb,
                              const NumericVector& probs,
                              const NumericVector& splitPointY){
  // preprocessing
  int n = cellProb.nrow();
  int numberBin = cellProb.ncol();
  int numQuantiles = probs.length();
  // initialization
  NumericMatrix quantiles(n, numQuantiles);

  for (int i = 0; i < n; i++){
    for(int j = 0; j < numQuantiles; j++){
      double probTemp = 0; double probTempPrev = 0; int counter = 0;
      while(1 > 0){
        if(probTemp < probs[j]){probTempPrev = probTemp;
          probTemp += cellProb(i, counter); counter += 1;
        } else {
          if(counter == 0){quantiles(i, j) = splitPointY[0]; break;}
          else if(counter == numberBin + 1){quantiles(i, j) = splitPointY[numberBin]; break;}
          else {
            quantiles(i, j) = splitPointY[counter] * (probs[j] - probTempPrev)
            /(probTemp - probTempPrev) + splitPointY[counter - 1] * (probTemp - probs[j])
            /(probTemp - probTempPrev); break;
            }
          }
        }
      }
    }
    return quantiles;
  }
