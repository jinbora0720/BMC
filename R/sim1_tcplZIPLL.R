# Simulation 1

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 
library(ROCR)
library(tcpl)
# library(devtools)
# devtools::install_github("AnderWilson/ZIPLL")
library(ZIPLL)

# source code
path <- "~/Documents/BoraJin2018~/Research/DoseResponse/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler2.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simulate_data.R"))

########
# Data #
########
m <- 30
J <- 150

set.seed(123)
seedsave <- sample(10000, 50, replace=FALSE)

#################
# Store results #
#################
rmse_tcpl = tr_aucg_tcpl <- rep(NA, 50)
rmse_ZIPLL = tr_aucg_ZIPLL <- rep(NA, 50)

###################
# MCMC parameters #
###################
thin <- 10
burnin <- 10000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

###############
# Simulations #
###############
pb <- txtProgressBar(style=3,width=50)
for (iter in 1:50) {
  seed = seedsave[iter]
  gendata = generate_data(m, J, d=3, seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  names(misdata)[5] <- "Y"
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij)[missing_idx_col]
  tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
  tr_t = as.numeric(truth$t_ij)[-missing_idx_col]
  
  X = misdata$X
  orgX = misdata$orgX
  Y = misdata$Y
  m_j = misdata$m_j
  idx_j = misdata$idx_j
  Start = misdata$Start
  End = misdata$End
  
  ########
  # tcpl #
  ########
  rmse_TCPL = hitc = matrix(NA, m, J)
  for (j in 1:J) {
    for (i in 1:m_j[j]) {
      params <- tcplFit(logc = orgX[[j]][Start[i,j]:End[i,j]],
                        resp = Y[[j]][Start[i,j]:End[i,j]],
                        bmad = (median(Y[[j]]) + sd(Y[[j]]))/3, bidirectional=TRUE)
      rmse_TCPL[idx_j[[j]][i],j] = c(params$cnst_rmse, params$hill_rmse, params$gnls_rmse)[which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))]
      hitc[idx_j[[j]][i],j] = ifelse(which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))==1, 0, 1)
    }
  }
  
  ## 1. overall RMSE 
  rmse_tcpl[iter] = mean(rmse_TCPL, na.rm=TRUE)
  
  ## 2. AUC for gamma
  tr_pred_tcpl = prediction(as.numeric(hitc)[-missing_idx_col], tr_gamma)
  tr_rocs_tcpl = performance(tr_pred_tcpl, measure = "auc")
  tr_aucg_tcpl[iter] = tr_rocs_tcpl@y.values[[1]]
  
  #########
  # ZIPLL #
  #########
  # formulate data
  dat <- NULL
  for(i in 1:m){
    for (j in 1:J) {
      dat <- rbind(dat,cbind(i,j,simdata$orgX[[j]][simdata$Start[i,j]:simdata$End[i,j]],
                             simdata$orgY[[j]][simdata$Start[i,j]:simdata$End[i,j]]))
    }
  }
  
  # missing index
  ZIPLL_missing = NULL
  for (s in 1:nrow(missing_idx)) {
    i = missing_idx[s,1]; j = missing_idx[s,2]
    ZIPLL_missing = c(ZIPLL_missing, which(dat[,1]==i & dat[,2]==j))
  }
  
  # fit ZIPLL
  fit <- ZIPLL(dat, nitter=burnin+thin*save, burnin=burnin)
  
  ## 1. overall RMSE 
  rmse_ZIPLL[iter] = sqrt(mean((fit$dat[-ZIPLL_missing,4] - fit$dat[-ZIPLL_missing,5])^2, na.rm=TRUE))
  
  ## 2. AUC for train gamma 
  tr_pred_ZIPLL = prediction(fit$parms[,7][-test_idx], as.numeric(t(truth$gamma_ij))[-test_idx])
  tr_rocs_ZIPLL = performance(tr_pred_ZIPLL, measure = "auc")
  tr_aucg_ZIPLL[iter] = tr_rocs_ZIPLL@y.values[[1]]
  
  setTxtProgressBar(pb, iter/50) 
}
close(pb)

################
# Save results #
################
res_tcpl <- list(rmse = rmse_tcpl, tr_aucg = tr_aucg_tcpl, seedsave = seedsave)
res_ZIPLL <- list(rmse = rmse_ZIPLL, tr_aucg = tr_aucg_ZIPLL, seedsave = seedsave, fit_ZIPLL = fit)

saveRDS(res_tcpl, file.path(path, "data/sim1_tcpl.rds"))
saveRDS(res_ZIPLL, file.path(path, "data/sim1_ZIPLL.rds"))

lapply(res_tcpl[1:2], mean)
lapply(res_tcpl[1:2], sd)
lapply(res_ZIPLL[1:2], mean)
lapply(res_ZIPLL[1:2], sd)
