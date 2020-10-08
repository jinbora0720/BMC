# Simulation 1

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse)
library(splines) 
library(ROCR)

# source code
path <- "/work/bj91/BMC/"
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
  gendata = generate_data(m, J, d=3, seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  misdata$Y <- misdata$orgY
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij)[missing_idx]
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
  
  ########
  # bmci # 
  ########
  outi <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = TRUE, adapt = FALSE, simpler = "bmci", verbose = FALSE)
  
  Yhat.savei = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        (X[[j]][Start[i,j]:End[i,j],]%*%outi$beta_ij.save[[j]][,i,s])/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*outi$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  Yhat.postmi = lapply(Yhat.savei, rowMeans)
  
  Ytil.savei = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        Y[[j]][Start[i,j]:End[i,j]]/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*outi$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  res.postmi = lapply(mapply('-', Ytil.savei, Yhat.savei, SIMPLIFY = FALSE), rowMeans)
  
  ## 1. overall RMSE
  rmsei = sqrt(mean(unlist(res.postmi)^2))
  
  ## 2. AUC for gamma
  gamma_ij.savei = outi$gamma_ij.save
  gamma_ij.savei[is.na(gamma_ij.savei)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outi$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postmi = rowMeans(gamma_ij.savei, dim=2, na.rm=TRUE)
  
  tr_predi = prediction(as.numeric(gamma_ij.postmi)[-missing_idx_col], tr_gamma)
  tr_rocsi = performance(tr_predi, measure = "auc")
  tr_aucgi = tr_rocsi@y.values[[1]]
  
  ## 3. AUC for predicted gamma
  tt_predi = prediction(as.numeric(gamma_ij.postmi[missing_idx]), tt_gamma)
  tt_rocsi = performance(tt_predi, measure = "auc")
  tt_aucgi = tt_rocsi@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.postmi = rowMeans(outi$t_ij.save, dims=2)
  pred_ti = prediction(as.numeric(t_ij.postmi)[-missing_idx_col], tr_t)
  rocs_ti = performance(pred_ti, measure = "auc")
  tr_aucti = rocs_ti@y.values[[1]]
  
  ################
  # Save results #
  ################
  out_name <- paste0(path, "data/sim1_BMCi_res_", iter, ".rds")
  saveRDS(list(rmse = rmsei, tr_aucg = tr_aucgi, tt_aucg = tt_aucgi, tr_auct = tr_aucti,
               seed = seed, out = outi), file.path(out_name))
  






