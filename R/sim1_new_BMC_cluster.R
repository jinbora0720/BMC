# Simulation 1
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse)
library(splines) 
library(ROCR)

# source code
path <- "~/BMC/BMC/"
source(paste0(path, "source/bmc_new.R"))
source(paste0(path, "source/bmc_sampler_new.R"))
sourceCpp(paste0(path, "source/bmc_sampler2_new.cpp"))
sourceCpp(paste0(path, "source/samplerBits_new.cpp"))
source(paste0(path, "source/simulate_data_new.R"))

########
# Data #
########
m <- 30
J <- 150

set.seed(123)
seedsave <- sample(10000, 50, replace=FALSE)

###################
# MCMC parameters #
###################
thin <- 10
burnin <- 10000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

# parallelising
iter <- taskID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

  seed = seedsave[iter]
  gendata = generate_data(m, J, d=3, q=2, xi=0, alpha=c(-.5, 1.3), seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  misdata$Y <- misdata$orgY
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij[missing_idx])
  tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
  tr_t = as.numeric(truth$t_ij)[-missing_idx_col]
  
  X = misdata$X
  orgX = misdata$orgX
  Y = misdata$Y
  m_j = misdata$m_j
  idx_j = misdata$idx_j
  Start = misdata$Start
  End = misdata$End
  
  #########
  # hyper #
  #########
  p = 7
  a = p+2
  Sigj = sapply(1:J, function(j) {
    covtemp = sapply(1:m_j[j], function(i) {
      coef(lm(Y[[j]][Start[i,j]:End[i,j]]~
                bs(orgX[[j]][Start[i,j]:End[i,j]],df=p,intercept=TRUE)-1))    
    }) %>% t(.) 
    if (nrow(na.omit(covtemp))>1) { cov(covtemp, use="complete.obs") }
    else {rep(0, p*p)}
  }) 
  Sigj[,colSums(Sigj)==0] = rowSums(Sigj)/(J-sum(colSums(Sigj)==0))
  invSigj = apply(Sigj, 2, function(x) c(ginv(matrix(x, p, p))))
  R = matrix((a-p-1)*rowMeans(Sigj, na.rm=TRUE), p, p)
  v_d = var(unlist(Y))
  hyper = list(a = a, R = R, v_d = v_d)
  
  ##################
  # initial values #
  ##################
  q = 5 
  ols = lm(scale(unlist(Y))~bs(unlist(orgX), df=p)-1)
  beta_ij = list()
  for (j in 1:J) {
    beta_ij[[j]] = matrix(coef(ols), nrow=p, ncol=m_j[j])
  }
  init = list(q = q, beta_ij = beta_ij)
  
  #######
  # bmc # 
  #######
  out <- bmc_new(Data = misdata, MCMC = MCMC, hyper = hyper, init = init,
                 hetero = TRUE, hetero_simpler = TRUE,
                 gamma_simpler = FALSE, adapt = FALSE, 
                 update_xi = FALSE, apply_cutoff = FALSE, verbose = TRUE)
  
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
  res.postm = lapply(mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE), rowMeans)
  
  ## 1. overall RMSE
  rmse = sqrt(mean(unlist(res.postm)^2))
  
  ## 2. AUC for gamma 
  gamma_ij.save = out$gamma_ij.save 
  le.save = mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE)
  z.save = lapply(le.save, function(x) matrix(rnorm(m*J, x, 1), m, J))
  gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
  gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
  
  tr_pred = prediction(as.numeric(gamma_ij.postm)[-missing_idx_col], tr_gamma)
  tr_rocs = performance(tr_pred, measure = "auc")
  tr_aucg = tr_rocs@y.values[[1]]
  
  ## 3. AUC for predicted gamma 
  tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), tt_gamma)
  tt_rocs = performance(tt_pred, measure = "auc")
  tt_aucg = tt_rocs@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.postm = rowMeans(out$t_ij.save, dims=2)
  tr_pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], tr_t)
  tr_rocs_t = performance(tr_pred_t, measure = "auc")
  tr_auct = tr_rocs_t@y.values[[1]]
  
################
# Save results #
################
  out_name <- paste0(path, "data/sim1_new_BMC_res_", iter, ".RDS")
  saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
               tr_auct = tr_auct, seed = seed), file.path(out_name))



