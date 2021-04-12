# wrapper of LinCDEExample

# load helper.R, LinCDEBoosting.R, LinCDEPredict.R
# setwd("...") # fill in the address of the LinCDE
# sourceCpp("./LinCDESplit.cpp")
# sourceCpp("./LinCDEQuantile.cpp")
# sourceCpp("./LinCDECdf.cpp")

LinCDEExample(setting = "GM", m = 1) #GM, GLM; m: number of simulations


