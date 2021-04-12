# LinCDE
LinCDE: Conditional Density Estimation via Lindsey's Method

Description:

Conditional density estimation is a fundamental problem in statistics, with various scientific and practical applications such as genetics and economics. We propose a conditional density estimator based on tree boosting and Lindsey's method (LinCDE). LinCDE admits flexible modeling of the density family and captures distributional characteristics like modality and shape. In particular, LinCDE always produces smooth and nonnegative density estimates. Furthermore, in the presence of nuisances, LinCDE identifies the influential covariates to the response distribution. 

Code:

To start with, run the helper.R, LinCDEBoosting.R, and LinCDEPredict.R, and source LinCDESplit.cpp, LinCDEQuantile.cpp, and LinCDECdf.cpp. Next, run the one-line code in exampleWrapper.R with the default parameters. You will be able to see the conditional density estimates produced by LinCDE boosting at several representative subregions of the covariate space. Enjoy :)
