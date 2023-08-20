# intRlearner
Integrative R learner by Wu and Yang (2022)


library(glmnet)
library(mvtnorm)

library(stringr)
library(caret)
library(ncvreg)
source("./R/rlasso.R")
seed <- 0
set.seed(seed)
# If KPS_data = False: the data generating process in Simulation I in the paper;
#   otherwise Simulation II
KPS_data <- FALSE
if (KPS_data) {
  m <- 2000  # the sample size of RWE
  n <- 500  # the sample size of RCT
  N <- 1e5
  b <- NA
  n0 <- NA
  cf <- 0.75
  multidim <- FALSE
  p <- 2
  p0 <- NA
} else {
  m <- 500
  n <- NA
  N <- 1e5
  b <- 2.5
  n0 <- 7
  cf <- NA
  multidim <- NA
  p <- 10
  p0 <- 1# p-2
}

dataset <- gen_data(m, N, b=b, p=p, lambda="lambda.min", n0=n0, n=n, KPS_data=KPS_data, multidim=multidim, cf=cf, p0=p0)
x_rct <- dataset$x_rct
x_rwe <- dataset$x_rwe
a_rct <- dataset$a_rct
a_rwe <- dataset$a_rwe
y_rct <- dataset$y_rct
y_rwe <- dataset$y_rwe
newx_rlearner <- dataset$newx_rlearner
tau <- dataset$tau
# R-learner method with RCT
rlasso_fit_rct <- rlasso(x=x_rct, w=a_rct, y=y_rct, p=p)
rlasso_est_rct <- predict(rlasso_fit_rct, newx_rlearner)
# R-learner method with RWE
rlasso_fit_rwe <- rlasso(x=x_rwe, w=a_rwe, y=y_rwe, p=p)
rlasso_est_rwe <- predict(rlasso_fit_rwe, newx_rlearner)
# Naive data combination
rlasso_fit_naive <- rlasso(x=rbind(x_rct, x_rwe), w=c(a_rct, a_rwe), y=c(y_rct, y_rwe), p=p)
rlasso_est_naive <- predict(rlasso_fit_naive, newx_rlearner)
# R-learner for the combined data
rlasso_fit_cmb <- intrlearner(x=rbind(x_rct, x_rwe), w=c(a_rct, a_rwe),
                              y=c(y_rct, y_rwe), s=c(rep(1, length(a_rct)), rep(0, m)), p=p)
rlasso_est_cmb <- predict(rlasso_fit_cmb, newx_rlearner)
cmb_rlasso <- sqrt(mean((rlasso_est_cmb - tau)**2))

rct_rlasso <- sqrt(mean((rlasso_est_rct - tau)**2))
rwe_rlasso <- sqrt(mean((rlasso_est_rwe - tau)**2))
naive_rlasso <- sqrt(mean((rlasso_est_naive - tau)**2))
print("RMSE of int_rlearner, rct_rlearner, os_rlearner, naive_rlearner:")
print(c(cmb_rlasso, rct_rlasso, rwe_rlasso, naive_rlasso))




