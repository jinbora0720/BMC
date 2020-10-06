# Process simulation 1 results 

# dependencies 
library(tidyverse)
library(splines)

# source code 
source("metadata.R")
source("simulate_data.R")

# simulation 1 data (run either this) ############################################################

res <- readRDS("sim1_BMC.rds")
out <- readRDS("sim1_BMC_out1.rds")
out$beta_ij.save <- c(readRDS("sim1_BMC_out2.rds"), readRDS("sim1_BMC_out3.rds"))

m <- 30; J <- 150
gendata <- generate_data(m, J, d=3, seed=res$seedsave[50])
simdata <- gendata$simdata
misdata <- data_missing(simdata, prob_missing=0.03, missing_idx=NULL, seed=res$seedsave[50])
missing_idx <- misdata$missing_idx 
X <- misdata$X
orgX <- misdata$orgX
Y <- misdata$orgY
m_j <- misdata$m_j
idx_j <- misdata$idx_j
Start <- misdata$Start
End <- misdata$End

##################################################################################################

# obese data (or this) ###########################################################################

# datalist <- readRDS("~/Documents/GitHub/BMC/data/datalist.rds")
# obese_data <- datalist$obese_data
# hit_vec <- datalist$hit_vec
# obese_out <- readRDS("obese_out1.rds")
# out <- obese_out$out
# missing_idx <- obese_out$missing_idx
# out$beta_ij.save <- c(readRDS("obese_out2.rds"), readRDS("obese_out3.rds"))
# neuro_active_aenm <- readRDS("neuro_active_aenm.rds")
# 
# # arrage data
# meta <- metadata(obese_data)
# uniq_chnm <- meta$uniq_chnm
# m <- length(uniq_chnm)
# uniq_aenm <- meta$uniq_aenm
# J <- length(uniq_aenm)
# simdata <- list(X = meta$X, orgX = lapply(meta$orgData, function(x) x[,1]), 
#                 orgY = meta$Y, K_ij = meta$K_ij, m_j = meta$m_j, 
#                 idx_j = meta$idx_j, Start = meta$Start, End = meta$End)
# 
# misdata <- data_missing(simdata = simdata, missing_idx = missing_idx, seed = 330)
# 
# X <- misdata$X
# orgX <- misdata$orgX
# Y <- misdata$orgY
# m_j <- misdata$m_j
# idx_j <- misdata$idx_j
# Start <- misdata$Start
# End <- misdata$End

##################################################################################################

# MCMC parameter
save <- 1000

# posterior predictive Y
predY.save <- lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]) + 
        rnorm(End[i,j]-Start[i,j]+1, 0, 1)*exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)*sqrt(out$sigj_sq.save[j,s])
    }) %>% unlist() 
  }) 
})
predY.postL <- lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
predY.postU <- lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))

# full model residual 
Yhat.save <- lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s])/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
Yhat.postm <- lapply(Yhat.save, rowMeans)

Ytil.save <- lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      Y[[j]][Start[i,j]:End[i,j]]/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
res.save <- mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE)
res.postm <- lapply(res.save, rowMeans)

# dose-response function estimate
Fit.save <- lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]
    }) %>% unlist() 
  }) 
})
Fit.postm <- lapply(Fit.save, rowMeans)
Fit.postL <- lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
Fit.postU <- lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))
naive_res <- mapply('-', Y, Fit.postm, SIMPLIFY = FALSE)

# save simulation 1 data (run either this) #######################################################

saveRDS(list(predY.postL = predY.postL, predY.postU = predY.postU, 
             Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU, 
             Yhat.postm = Yhat.postm, res.postm = res.postm, naive_res = naive_res), 
        "sim1_processed.rds")

##################################################################################################

# obese data (or this) ###########################################################################

# saveRDS(list(predY.postL = predY.postL, predY.postU = predY.postU, 
#              Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU, 
#              Yhat.postm = Yhat.postm, res.postm = res.postm, naive_res = naive_res), 
#         "obese_processed.rds")

##################################################################################################
