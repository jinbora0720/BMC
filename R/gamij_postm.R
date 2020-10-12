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

save <- 1000

iter <- 50
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
out = readRDS(paste0(path, "data/sim1_BMC_res_50.rds"))$out
gamma_ij.save = out$gamma_ij.save 
z.save = lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)

## 3. AUC for predicted gamma 
tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), tt_gamma)
tt_rocs = performance(tt_pred, measure = "auc")
tt_aucg = tt_rocs@y.values[[1]]

saveRDS(list(tt_aucg = tt_aucg, gamma_ij.postm = gamma_ij.postm), paste0(path, "data/sim1_BMC_res_50_re.rds"))
