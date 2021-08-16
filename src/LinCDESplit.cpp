#include <Rcpp.h>
using namespace Rcpp;

IntegerVector orderIndex(NumericVector arr) {
  // if (is_true(any(duplicated(x)))) {
  //   Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  // }
  NumericVector sorted = clone(arr).sort();
  return match(sorted, arr);
}

//' LinCDESplit
//'
//' This function computes the approximately optimal split in C++. A split is invalid if the sample size on the left or right of the split falls below 10
//'
//' @param X input matrix, of dimension nobs x nvars; each row represents an observation vector.
//' @param yIndex discretized response vector, of length nobs.
//' @param cellProb cell probability matrix, of dimension nobs x \code{numberBin}.
//' @param z sufficient statistics matrix, of dimension \code{numberBin} x number of sufficient statistics.
//' @param covMatrixInv inverse of psi'' (number of sufficient statistics x number of sufficient statistics).
//' @param splitPoint a candidate split list (length nvars). Each element is a vector corresponding to a certain variable's candidate splits (including the left and right end points).
//' @param numberBin the number of bins for response discretization.
//'
//' @return The function returns \code{splitVar}: the index of the variable to split at; \code{splitVal}: the cut-point of the split; \code{improvement}: the contribution of the split to the objective.
//'
//[[Rcpp::export]]
NumericVector LinCDESplit(const NumericMatrix& X, const IntegerVector& yIndex,
                          const NumericMatrix& cellProb, const NumericMatrix& z,
                          const NumericMatrix& covMatrixInv,
                          const List& splitPoint, const int& numberBin){
  // preprocessing
  int n = X.nrow();
  int k = z.ncol();
  int d = X.ncol();
  // initialization
  List improvement(d);
  double maxImprovement = -1; int splitVar = 0; double splitVal = 0;

  // loop over potential splits
  for (int m = 0; m < d;  m++){
    NumericVector currX(n);
    NumericVector currSplitPoint = splitPoint[m];
    int s = currSplitPoint.length();
    currSplitPoint[0] -= 100; currSplitPoint[s-1] += 100;
    NumericVector currImprovement(s-1);
    for(int i = 0; i < n; i++){currX[i] = X(i,m);}
    NumericVector sortX = clone(currX);
    NumericVector indexXTemp = clone(currX);
    std::sort(sortX.begin(), sortX.begin() + n);
    IntegerVector indexX = orderIndex(indexXTemp) - 1;
    // compute cell proportions, cell mean of z, cell mean of r
    NumericVector propCell(s);
    IntegerMatrix yIndexPropCell(numberBin, s);
    NumericMatrix cellProbCell(numberBin, s);
    int counter = 0;
    for (int i = 0; i < n; i++){
      if(sortX[i] > currSplitPoint[counter]){
        counter += 1;
        if(sortX[i] > currSplitPoint[counter]){i -= 1; continue;}
      }
      propCell[counter] += 1;
      yIndexPropCell(yIndex[indexX[i]], counter) += 1;
      cellProbCell(_, counter) = cellProbCell(_, counter) + cellProb(indexX[i], _);
    }

    NumericMatrix zCellMean(k, s);
    NumericMatrix rCellMean(k, s);
    for (int i = 0; i < s; i++){
      for (int j = 0; j < k;  j++){
        double zTemp = 0; double rTemp = 0;
        for (int l = 0; l < numberBin; l++){
          zTemp += yIndexPropCell(l, i) * z(l, j);
          rTemp += cellProbCell(l, i) * z(l, j);
        }
        if(propCell[i] > 0){
          zCellMean(j, i) = zTemp / propCell[i];
          rCellMean(j, i) = rTemp / propCell[i];
        }
      }
    }

    // find the minimum and maximum valid splits
    int minSplitIndex = 0; int maxSplitIndex = 0;
    int countL = 0; int countR = 0;
    for(int i = 0; i < s; i++){
      countL += propCell[i];
      if (countL > 10){ // minimum samples on either side of the split
        minSplitIndex = i; break;
      }
    }
    for(int i = s-1; i >= 0; i--){
      countR += propCell[i];
      if (countR > 10){
        maxSplitIndex = i-1; break;
      }
    }
    propCell = propCell/n;

    if(minSplitIndex > maxSplitIndex){continue;}

    double pL = 0; double pR = 1;
    NumericVector zL(k); NumericVector zR(k);
    NumericVector rL(k); NumericVector rR(k);
    for(int i = 0; i < minSplitIndex; i++){
      pL += propCell[i];
      zL += propCell[i] * zCellMean(_, i);
      rL += propCell[i] * rCellMean(_, i);
    }
    if(pL > 0){
      zL = zL / pL; rL = rL / pL;
    }
    pR = 1 - pL;
    for(int i = s-1; i >= minSplitIndex; i--){
      zR += propCell[i] * zCellMean(_, i);
      rR += propCell[i] * rCellMean(_, i);
    }
    if(pR > 0){
      zR = zR / pR; rR = rR / pR;
    }

    // compute the improvement in the objective
    for (int i = minSplitIndex; i <= maxSplitIndex; i++){
      zL = zL * pL / (pL + propCell[i]) + zCellMean(_, i) * propCell[i] / (pL + propCell[i]);
      zR = zR * pR / (pR - propCell[i]) - zCellMean(_, i) * propCell[i] / (pR - propCell[i]);
      rL = rL * pL / (pL + propCell[i]) + rCellMean(_, i) * propCell[i] / (pL + propCell[i]);
      rR = rR * pR / (pR - propCell[i]) - rCellMean(_, i) * propCell[i] / (pR - propCell[i]);
      pL = pL + propCell[i]; pR = pR - propCell[i];
      NumericVector temp = (zL - rL) - (zR - rR);
      for (int j = 0; j < k; j++){
        for (int l = 0; l < k; l++){
          currImprovement[i] += temp[j] * temp[l] * covMatrixInv(j, l);
        }
      }
      currImprovement[i] = pL * pR * currImprovement[i];
      if (currImprovement[i] > maxImprovement){
        splitVar = m; splitVal = currSplitPoint[i]; maxImprovement = currImprovement[i];
      }
    }
    improvement[m] = currImprovement;
  }
  NumericVector result(3);
  result[0] = splitVar + 1; result[1] = splitVal; result[2] = maxImprovement * n / 2;
  return result;
}


