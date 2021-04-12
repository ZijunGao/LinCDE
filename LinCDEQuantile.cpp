#include <Rcpp.h>
using namespace Rcpp;

// function: compute quantiles for LinCDE (C++)
// input:
  // cellProb: matrix of cell Probabilities, each row sum to one (n by numberBin)
  // probs: probabilities to evaluate quantiles
  // splitPointY: split points of responses (length numberBin+1)
    
// output:
//   quantiles: matrix of quantiles (n by length(probs))

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
