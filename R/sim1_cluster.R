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
  out <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
             hetero = TRUE, adapt = FALSE, simpler = FALSE, verbose = FALSE)
  
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
  z.save = lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
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
  pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], tr_t)
  rocs_t = performance(pred_t, measure = "auc")
  tr_auct = rocs_t@y.values[[1]]
  
  ########
  # bmc0 # 
  ########
  out0 <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = TRUE, adapt = FALSE, simpler = "bmc0", verbose = FALSE)
  
  Yhat.save0 = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        (X[[j]][Start[i,j]:End[i,j],]%*%out0$beta_ij.save[[j]][,i,s])/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*out0$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  Yhat.postm0 = lapply(Yhat.save0, rowMeans)
  
  Ytil.save0 = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        Y[[j]][Start[i,j]:End[i,j]]/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*out0$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  res.postm0 = lapply(mapply('-', Ytil.save0, Yhat.save0, SIMPLIFY = FALSE), rowMeans)
  
  ## 1. overall RMSE
  rmse0 = sqrt(mean(unlist(res.postm0)^2))
  
  ## 2. AUC for gamma
  gamma_ij.save0 = out0$gamma_ij.save
  gamma_ij.save0[is.na(gamma_ij.save0)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, out0$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postm0 = rowMeans(gamma_ij.save0, dim=2, na.rm=TRUE)
  
  tr_pred0 = prediction(as.numeric(gamma_ij.postm0)[-missing_idx_col], tr_gamma)
  tr_rocs0 = performance(tr_pred0, measure = "auc")
  tr_aucg0 = tr_rocs0@y.values[[1]]
  
  ## 3. AUC for predicted gamma
  tt_pred0 = prediction(as.numeric(gamma_ij.postm0[missing_idx]), tt_gamma)
  tt_rocs0 = performance(tt_pred0, measure = "auc")
  tt_aucg0 = tt_rocs0@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.postm0 = rowMeans(out0$t_ij.save, dims=2)
  pred_t0 = prediction(as.numeric(t_ij.postm0)[-missing_idx_col], tr_t)
  rocs_t0 = performance(pred_t0, measure = "auc")
  tr_auct0 = rocs_t0@y.values[[1]]
  
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
  
  ########
  # bmcj # 
  ########
  outj <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = TRUE, adapt = FALSE, simpler = "bmcj", verbose = FALSE)
  
  Yhat.savej = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        (X[[j]][Start[i,j]:End[i,j],]%*%outj$beta_ij.save[[j]][,i,s])/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*outj$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  Yhat.postmj = lapply(Yhat.savej, rowMeans)
  
  Ytil.savej = lapply(1:J, function(j) {
    sapply(1:save, function(s) {
      lapply(1:m_j[j], function(i) {
        Y[[j]][Start[i,j]:End[i,j]]/
          exp(orgX[[j]][Start[i,j]:End[i,j]]*outj$d_ij.save[idx_j[[j]][i],j,s]/2)
      }) %>% unlist()
    })
  })
  res.postmj = lapply(mapply('-', Ytil.savej, Yhat.savej, SIMPLIFY = FALSE), rowMeans)
  
  ## 1. overall RMSE
  rmsej = sqrt(mean(unlist(res.postmj)^2))

  ## 2. AUC for gamma
  gamma_ij.savej = outj$gamma_ij.save
  gamma_ij.savej[is.na(gamma_ij.savej)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outj$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postmj = rowMeans(gamma_ij.savej, dim=2, na.rm=TRUE)
  
  tr_predj = prediction(as.numeric(gamma_ij.postmj)[-missing_idx_col], tr_gamma)
  tr_rocsj = performance(tr_predj, measure = "auc")
  tr_aucgj = tr_rocsj@y.values[[1]]
  
  ## 3. AUC for predicted gamma
  tt_predj = prediction(as.numeric(gamma_ij.postmj[missing_idx]), tt_gamma)
  tt_rocsj = performance(tt_predj, measure = "auc")
  tt_aucgj = tt_rocsj@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.postmj = rowMeans(outj$t_ij.save, dims=2)
  pred_tj = prediction(as.numeric(t_ij.postmj)[-missing_idx_col], tr_t)
  rocs_tj = performance(pred_tj, measure = "auc")
  tr_auctj = rocs_tj@y.values[[1]]

################
# Save results #
################
  out_name <- paste0(path, "data/sim1_res_", iter, ".rds")
  saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, tr_auct = tr_auct,
               rmse0 = rmse0, tr_aucg0 = tr_aucg0, tt_aucg0 = tt_aucg0, tr_auct0 = tr_auct0,
               rmsei = rmsei, tr_aucgi = tr_aucgi, tt_aucgi = tt_aucgi, tr_aucti = tr_aucti,
               rmsej = rmsej, tr_aucgj = tr_aucgj, tt_aucgj = tt_aucgj, tr_auctj = tr_auctj,
               seed = seed, out = out), file.path(out_name))



