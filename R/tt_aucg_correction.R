# Simulation 1

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse)
library(splines) 
library(ROCR)

# source code
path <- "/work/bj91/BMC/"
source(paste0(path, "source/simulate_data.R"))

########
# Data #
########
m <- 30
J <- 150

set.seed(123)
seedsave <- sample(10000, 50, replace=FALSE)

tt_aucg = tt_aucg0 = tt_aucgi = rep(NA, 50)
save <- 1000

for (iter in 1:50) {
  seed = seedsave[iter]
  gendata = generate_data(m, J, d=3, seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  misdata$Y <- misdata$orgY
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij[missing_idx])

  # BMC
  out = readRDS(paste0(path, "data/sim1_BMC_res_", iter, ".rds"))$out
  gamma_ij.save = out$gamma_ij.save 
  z.save = lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
  gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
  gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
  
  ## 3. AUC for predicted gamma 
  tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), tt_gamma)
  tt_rocs = performance(tt_pred, measure = "auc")
  tt_aucg[iter] = tt_rocs@y.values[[1]]
  
  # BMC0
  out0 = readRDS(paste0(path, "data/sim1_BMC0_res_", iter, ".rds"))$out
  gamma_ij.save0 = out0$gamma_ij.save
  gamma_ij.save0[is.na(gamma_ij.save0)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, out0$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postm0 = rowMeans(gamma_ij.save0, dim=2, na.rm=TRUE)
  
  tt_pred0 = prediction(as.numeric(gamma_ij.postm0[missing_idx]), tt_gamma)
  tt_rocs0 = performance(tt_pred0, measure = "auc")
  tt_aucg0[iter] = tt_rocs0@y.values[[1]]
}

BMC_new = readRDS(paste0(path, "data/sim1_BMC_new.rds"))
BMC_new$tt_aucg = tt_aucg
saveRDS(BMC_new, paste0(path, "data/sim1_BMC_new2.rds"))

BMC0_new = readRDS(paste0(path, "data/sim1_BMC0_new.rds"))
BMC0_new$tt_aucg = tt_aucg0
saveRDS(BMC0_new, paste0(path, "data/sim1_BMC0_new2.rds"))
