# Process simulation 1 results 

# dependencies 
library(tidyverse)

source("simulate_data.R")
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

# MCMC parameter
save  = 1000

# posterior predictive Y
predY.save <- lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]) + 
        rnorm(End[i,j]-Start[i,j]+1, 0, 1)*exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)*sqrt(out$sigj_sq.save[j,s])
    }) %>% unlist() 
  }) 
})
predY.postL = lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
predY.postU = lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))

# full model residual 
Yhat.save = lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s])/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
Yhat.postm = lapply(Yhat.save, rowMeans)

Ytil.save = lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      Y[[j]][Start[i,j]:End[i,j]]/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
res.save = mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE)
res.postm = lapply(res.save, rowMeans)

# dose-response function estimate
Fit.save = lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]
    }) %>% unlist() 
  }) 
})
Fit.postm = lapply(Fit.save, rowMeans)
Fit.postL = lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
Fit.postU = lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))
naive_res = mapply('-', Y, Fit.postm, SIMPLIFY = FALSE)

saveRDS(list(predY.postL = predY.postL, predY.postU = predY.postU, 
             Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU, 
             Yhat.postm = Yhat.postm, res.postm = res.postm, naive_res = naive_res), 
        "sim1_processed.rds")
