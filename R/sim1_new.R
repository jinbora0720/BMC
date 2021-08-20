# Simulation 1

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 
library(ROCR)
library(tcpl)

# source code
path <- "~/Documents/BoraJin2018~/Research/DoseResponse/BMC/"
source(paste0(path, "source/bmc_new.R"))
source(paste0(path, "source/bmc_sampler_new.R"))
sourceCpp(paste0(path, "source/bmc_sampler2.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simulate_data_new.R"))
source(paste0(path, "source/dosres_plot.R"))

########
# Data #
########
m <- 30
J <- 150

# set.seed(123)
# seedsave <- sample(10000, 50, replace=FALSE)
# 
# #################
# # Store results #
# #################
# rmse = tr_aucg = tt_aucg = tr_auct <- rep(NA, 50)
# rmse0 = tr_aucg0 = tt_aucg0 = tr_auct0 <- rep(NA, 50)
# rmsei = tr_aucgi = tt_aucgi = tr_aucti <- rep(NA, 50)
# rmsej = tr_aucgj = tt_aucgj = tr_auctj <- rep(NA, 50)
# rmse_tcpl = tr_aucg_tcpl <- rep(NA, 50)
# rmse_ZIPLL = tr_aucg_ZIPLL <- rep(NA, 50)

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
# pb <- txtProgressBar(style=3,width=50)
# for (iter in 1:50) {
#   seed = seedsave[iter]
  seed = 123
  gendata = generate_data(m, J, d=3, seed)
  simdata = gendata$simdata
  truth = gendata$truth
  simdata$rxf = truth$rxf
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  misdata$Y <- misdata$orgY
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij[missing_idx])
  tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
  tt_t = as.numeric(truth$t_ij[missing_idx])
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
             hetero = TRUE, adapt = FALSE, simpler = FALSE, verbose = TRUE)
  
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
  # rmse[iter] = sqrt(mean(unlist(res.postm)^2))
  rmse = sqrt(mean(unlist(res.postm)^2))
  
  ## 2. AUC for gamma 
  gamma_ij.save = out$gamma_ij.save 
  le.save = mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE)
  z.save = lapply(le.save, function(x) matrix(rnorm(m*J, x, 1), m, J))
  gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
  gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
  
  tr_pred = prediction(as.numeric(gamma_ij.postm)[-missing_idx_col], tr_gamma)
  tr_rocs = performance(tr_pred, measure = "auc")
  # tr_aucg[iter] = tr_rocs@y.values[[1]]
  tr_aucg = tr_rocs@y.values[[1]]
  
  ## 3. AUC for predicted gamma 
  tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), tt_gamma)
  tt_rocs = performance(tt_pred, measure = "auc")
  # tt_aucg[iter] = tt_rocs@y.values[[1]]
  tt_aucg = tt_rocs@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.save = out$t_ij.save 
  lea = lapply(1:save, function(s) out$alpha.save[1,s] + out$alpha.save[2,s]*le.save[[s]])
  u.save = lapply(lea, function(x) matrix(rnorm(m*J, x, 1), m, J))
  t_ij.save[is.na(t_ij.save)] = sapply(u.save, function(x) x[missing_idx]>0)
  t_ij.postm = rowMeans(t_ij.save, dim=2, na.rm=TRUE)
  
  tr_pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], tr_t)
  tr_rocs_t = performance(tr_pred_t, measure = "auc")
  # tr_auct[iter] = tr_rocs_t@y.values[[1]]
  tr_auct = tr_rocs_t@y.values[[1]]
  
  ## 5. AUC for predicted t_ij 
  tt_pred_t = prediction(as.numeric(t_ij.postm[missing_idx]), tt_t)
  tt_rocs_t = performance(tt_pred_t, measure = "auc")
  # tt_auct[iter] = tt_rocs_t@y.values[[1]]
  tt_auct = tt_rocs_t@y.values[[1]]

  saveRDS(list(out = out, rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
               tr_auct = tr_auct, tt_auct = tt_auct), paste0(path, "data/sim1_new.RDS"))
  
#   ########
#   # bmc0 # 
#   ########
#   out0 <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
#               hetero = TRUE, adapt = FALSE, simpler = "bmc0", verbose = FALSE)
#   
#   Yhat.save0 = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         (X[[j]][Start[i,j]:End[i,j],]%*%out0$beta_ij.save[[j]][,i,s])/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*out0$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   Yhat.postm0 = lapply(Yhat.save0, rowMeans)
#   
#   Ytil.save0 = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         Y[[j]][Start[i,j]:End[i,j]]/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*out0$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   res.postm0 = lapply(mapply('-', Ytil.save0, Yhat.save0, SIMPLIFY = FALSE), rowMeans)
#   
#   ## 1. overall RMSE
#   rmse0[iter] = sqrt(mean(unlist(res.postm0)^2))
#   
#   ## 2. AUC for gamma
#   gamma_ij.save0 = out0$gamma_ij.save
#   gamma_ij.save0[is.na(gamma_ij.save0)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, out0$pi_ij.save[,,s][missing_idx]))
#   gamma_ij.postm0 = rowMeans(gamma_ij.save0, dim=2, na.rm=TRUE)
#   
#   tr_pred0 = prediction(as.numeric(gamma_ij.postm0)[-missing_idx_col], tr_gamma)
#   tr_rocs0 = performance(tr_pred0, measure = "auc")
#   tr_aucg0[iter] = tr_rocs0@y.values[[1]]
#   
#   ## 3. AUC for predicted gamma
#   tt_pred0 = prediction(as.numeric(gamma_ij.postm0[missing_idx]), tt_gamma)
#   tt_rocs0 = performance(tt_pred0, measure = "auc")
#   tt_aucg0[iter] = tt_rocs0@y.values[[1]]
#   
#   ## 4. AUC for t_ij
#   t_ij.postm0 = rowMeans(out0$t_ij.save, dims=2)
#   pred_t0 = prediction(as.numeric(t_ij.postm0)[-missing_idx_col], tr_t)
#   rocs_t0 = performance(pred_t0, measure = "auc")
#   tr_auct0[iter] = rocs_t0@y.values[[1]]
#   
#   ########
#   # bmci # 
#   ########
#   outi <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
#               hetero = TRUE, adapt = FALSE, simpler = "bmci", verbose = FALSE)
#   
#   Yhat.savei = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         (X[[j]][Start[i,j]:End[i,j],]%*%outi$beta_ij.save[[j]][,i,s])/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*outi$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   Yhat.postmi = lapply(Yhat.savei, rowMeans)
#   
#   Ytil.savei = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         Y[[j]][Start[i,j]:End[i,j]]/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*outi$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   res.postmi = lapply(mapply('-', Ytil.savei, Yhat.savei, SIMPLIFY = FALSE), rowMeans)
#   
#   ## 1. overall RMSE
#   rmsei[iter] = sqrt(mean(unlist(res.postmi)^2))
#   
#   ## 2. AUC for gamma
#   gamma_ij.savei = outi$gamma_ij.save
#   gamma_ij.savei[is.na(gamma_ij.savei)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outi$pi_ij.save[,,s][missing_idx]))
#   gamma_ij.postmi = rowMeans(gamma_ij.savei, dim=2, na.rm=TRUE)
#   
#   tr_predi = prediction(as.numeric(gamma_ij.postmi)[-missing_idx_col], tr_gamma)
#   tr_rocsi = performance(tr_predi, measure = "auc")
#   tr_aucgi[iter] = tr_rocsi@y.values[[1]]
#   
#   ## 3. AUC for predicted gamma
#   tt_predi = prediction(as.numeric(gamma_ij.postmi[missing_idx]), tt_gamma)
#   tt_rocsi = performance(tt_predi, measure = "auc")
#   tt_aucgi[iter] = tt_rocsi@y.values[[1]]
#   
#   ## 4. AUC for t_ij
#   t_ij.postmi = rowMeans(outi$t_ij.save, dims=2)
#   pred_ti = prediction(as.numeric(t_ij.postmi)[-missing_idx_col], tr_t)
#   rocs_ti = performance(pred_ti, measure = "auc")
#   tr_aucti[iter] = rocs_ti@y.values[[1]]
#   
#   ########
#   # bmcj # 
#   ########
#   outj <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
#               hetero = TRUE, adapt = FALSE, simpler = "bmcj", verbose = FALSE)
#   
#   Yhat.savej = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         (X[[j]][Start[i,j]:End[i,j],]%*%outj$beta_ij.save[[j]][,i,s])/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*outj$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   Yhat.postmj = lapply(Yhat.savej, rowMeans)
#   
#   Ytil.savej = lapply(1:J, function(j) {
#     sapply(1:save, function(s) {
#       lapply(1:m_j[j], function(i) {
#         Y[[j]][Start[i,j]:End[i,j]]/
#           exp(orgX[[j]][Start[i,j]:End[i,j]]*outj$d_ij.save[idx_j[[j]][i],j,s]/2)
#       }) %>% unlist()
#     })
#   })
#   res.postmj = lapply(mapply('-', Ytil.savej, Yhat.savej, SIMPLIFY = FALSE), rowMeans)
#   
#   ## 1. overall RMSE
#   rmsej[iter] = sqrt(mean(unlist(res.postmj)^2))
#   
#   ## 2. AUC for gamma
#   gamma_ij.savej = outj$gamma_ij.save
#   gamma_ij.savej[is.na(gamma_ij.savej)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outj$pi_ij.save[,,s][missing_idx]))
#   gamma_ij.postmj = rowMeans(gamma_ij.savej, dim=2, na.rm=TRUE)
#   
#   tr_predj = prediction(as.numeric(gamma_ij.postmj)[-missing_idx_col], tr_gamma)
#   tr_rocsj = performance(tr_predj, measure = "auc")
#   tr_aucgj[iter] = tr_rocsj@y.values[[1]]
#   
#   ## 3. AUC for predicted gamma
#   tt_predj = prediction(as.numeric(gamma_ij.postmj[missing_idx]), tt_gamma)
#   tt_rocsj = performance(tt_predj, measure = "auc")
#   tt_aucgj[iter] = tt_rocsj@y.values[[1]]
#   
#   ## 4. AUC for t_ij
#   t_ij.postmj = rowMeans(outj$t_ij.save, dims=2)
#   pred_tj = prediction(as.numeric(t_ij.postmj)[-missing_idx_col], tr_t)
#   rocs_tj = performance(pred_tj, measure = "auc")
#   tr_auctj[iter] = rocs_tj@y.values[[1]]
#   
#   ########
#   # tcpl #
#   ########
#   rmse_TCPL = hitc = matrix(NA, m, J)
#   for (j in 1:J) {
#     for (i in 1:m_j[j]) {
#       params <- tcplFit(logc = orgX[[j]][Start[i,j]:End[i,j]],
#                         resp = Y[[j]][Start[i,j]:End[i,j]],
#                         bmad = (median(Y[[j]]) + sd(Y[[j]]))/3, bidirectional=TRUE)
#       rmse_TCPL[idx_j[[j]][i],j] = c(params$cnst_rmse, params$hill_rmse, params$gnls_rmse)[which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))]
#       hitc[idx_j[[j]][i],j] = ifelse(which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))==1, 0, 1)
#     }
#   }
#   
#   ## 1. overall RMSE 
#   rmse_tcpl[iter] = mean(rmse_TCPL, na.rm=TRUE)
#   
#   ## 2. AUC for gamma
#   tr_pred_tcpl = prediction(as.numeric(hitc)[-missing_idx_col], tr_gamma)
#   tr_rocs_tcpl = performance(tr_pred_tcpl, measure = "auc")
#   tr_aucg_tcpl[iter] = tr_rocs_tcpl@y.values[[1]]
#   
#   #########
#   # ZIPLL #
#   #########
#   # formulate data
#   dat <- NULL
#   for(i in 1:m){
#     for (j in 1:J) {
#       dat <- rbind(dat,cbind(i,j,simdata$orgX[[j]][simdata$Start[i,j]:simdata$End[i,j]],
#                              simdata$orgY[[j]][simdata$Start[i,j]:simdata$End[i,j]]))
#     }
#   }
#   
#   # missing index
#   ZIPLL_missing = NULL
#   for (s in 1:nrow(missing_idx)) {
#     i = missing_idx[s,1]; j = missing_idx[s,2]
#     ZIPLL_missing = c(ZIPLL_missing, which(dat[,1]==i & dat[,2]==j))
#   }
#   
#   # fit ZIPLL
#   fit <- ZIPLL(dat, nitter=burnin+thin*save, burnin=burnin)
#   
#   ## 1. overall RMSE 
#   rmse_ZIPLL[iter] = sqrt(mean((fit$dat[-ZIPLL_missing,4] - fit$dat[-ZIPLL_missing,5])^2, na.rm=TRUE))
#   
#   ## 2. AUC for train gamma 
#   tr_pred_ZIPLL = prediction(fit$parms[,7][-test_idx], as.numeric(t(truth$gamma_ij))[-test_idx])
#   tr_rocs_ZIPLL = performance(tr_pred_ZIPLL, measure = "auc")
#   tr_aucg_ZIPLL[iter] = tr_rocs_ZIPLL@y.values[[1]]
#   
#   setTxtProgressBar(pb, iter/50) 
# }
# close(pb)

# ################
# # Save results #
# ################
# res <- list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, tr_auct = tr_auct, 
#             seedsave = seedsave)
# res0 <- list(rmse = rmse0, tr_aucg = tr_aucg0, tt_aucg = tt_aucg0, tr_auct = tr_auct0,  
#              seedsave = seedsave)
# resi <- list(rmse = rmsei, tr_aucg = tr_aucgi, tt_aucg = tt_aucgi, tr_auct = tr_aucti,  
#              seedsave = seedsave)
# resj <- list(rmse = rmsej, tr_aucg = tr_aucgj, tt_aucg = tt_aucgj, tr_auct = tr_auctj,  
#              seedsave = seedsave)
# res_tcpl <- list(rmse = rmse_tcpl, tr_aucg = tr_aucg_tcpl, seedsave = seedsave)
# res_ZIPLL <- list(rmse = rmse_ZIPLL, tr_aucg = tr_aucg_ZIPLL, seedsave = seedsave)
# 
# # save result
# # saveRDS(res, file.path(path, "data/sim1_BMC.rds"))
# # saveRDS(list(gamma_ij.save = out$gamma_ij.save, 
# #              gamma_ij.postm = gamma_ij.postm, 
# #              Lambda.save = out$Lambda.save, 
# #              eta.save = out$eta.save, 
# #              covMean = out$covMean, 
# #              t_ij.save = out$t_ij.save, 
# #              d_ij.save = out$d_ij.save), file.path(path, "data/sim1_BMC_out1.rds"))
# # saveRDS(out$beta_ij.save[c(1:75)], file.path(path, "data/sim1_BMC_out2.rds"))
# # saveRDS(out$beta_ij.save[c(76:150)], file.path(path, "data/sim1_BMC_out3.rds"))
# # saveRDS(res0, file.path(path, "data/sim1_BMC0.rds"))
# # saveRDS(resi, file.path(path, "data/sim1_BMCi.rds"))
# # saveRDS(resj, file.path(path, "data/sim1_BMCj.rds"))
# # saveRDS(res_tcpl, file.path(path, "data/sim1_tcpl.rds"))
# # saveRDS(res_ZIPLL, file.path(path, "data/sim1_ZIPLL.rds"))
# 
# # call result
# # res <- readRDS(paste0(path, "data/sim1_BMC.rds")) # BMC
# # out <- readRDS(paste0(path, "data/sim1_BMC_out1.rds")) # BMC results at iter = 50
# # gamma_ij.postm <- out$gamma_ij.postm
# # out$beta_ij.save <- c(readRDS(paste0(path, "data/sim1_BMC_out2.rds")), readRDS(paste0(path, "data/sim1_BMC_out3.rds")))
# # res0 <- readRDS(paste0(path, "data/sim1_BMC0.rds")) # BMC_0
# # resi <- readRDS(paste0(path, "data/sim1_BMCi.rds")) # BMC_i
# # resj <- readRDS(paste0(path, "data/sim1_BMCj.rds")) # BMC_j
# # res_tcpl <- readRDS(paste0(path, "data/sim1_tcpl.rds")) # tcpl
# # res_ZIPLL <- readRDS(paste0(path, "data/sim1_ZIPLL.rds")) # ZIPLL

###################
# Process results #
###################
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
prc_res <- list(predY.postL = predY.postL, predY.postU = predY.postU, 
                Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU, 
                Yhat.postm = Yhat.postm, res.postm = res.postm, naive_res = naive_res)

# save result
# saveRDS(prc_res, paste0(path, "data/sim1_processed.rds"))

# call result
# prc_res <- readRDS(paste0(path, "data/sim1_processed.rds")) 

#----------------------------------------------------------------------------------------------------------------------------------
# Table 1: summary of results (make mean and sd tables separately)
mres <- matrix(c(round(unlist(lapply(res[1:4], mean)), 3), 
                 round(unlist(lapply(res0[1:4], mean)), 3),
                 round(unlist(lapply(resi[1:4], mean)), 3),
                 round(unlist(lapply(resj[1:4], mean)), 3), 
                 c(round(unlist(lapply(res_ZIPLL[1:2], mean)), 3), NA, NA), 
                 c(round(unlist(lapply(res_tcpl[1:2], mean)), 3), NA, NA)), 4, 6)
rownames(mres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 'In-sample AUC for tij')
colnames(mres) <- c('BMC', 'BMC0', 'BMCi', 'BMCj', 'ZIPLL', 'tcpl')

sdres <- matrix(c(round(unlist(lapply(res[1:4], sd)), 3), 
                  round(unlist(lapply(res0[1:4], sd)), 3), 
                  round(unlist(lapply(resi[1:4], sd)), 3), 
                  round(unlist(lapply(resj[1:4], sd)), 3), 
                  c(round(unlist(lapply(res_ZIPLL[1:2], sd)), 3), NA, NA), 
                  c(round(unlist(lapply(res_tcpl[1:2], sd)), 3), NA, NA)), 4, 6)
rownames(sdres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 'In-sample AUC for tij')
colnames(sdres) <- c('BMC', 'BMC0', 'BMCi', 'BMCj', 'ZIPLL', 'tcpl')

cat('mean table\n'); mres
cat('sd table\n'); sdres

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S5: heat map of the estimated and true correlation matrix 
z_cor_T <- cov2cor(tcrossprod(truth$Lambda) + diag(1, m))
z_cor_data <- reshape2::melt(z_cor_T)
cor.postm <- cov2cor(out$covMean)
cor.postm_data <- data.frame(reshape2::melt(cor.postm), truevalue = z_cor_data$value)

rng <- range(cor.postm, z_cor_T)
t1 <- ggplot(cor.postm_data)+
  geom_tile(aes(reorder(Var2, truevalue), reorder(Var1, truevalue), fill=value)) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white", 
                       limits=c(floor(rng[1]), ceiling(rng[2]))) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  labs(x="Chemical", y="Chemical", title="Estimate") 

t2 <- ggplot(z_cor_T %>% reshape2::melt())+
  geom_tile(aes(reorder(Var2, value), reorder(Var1, value), fill=value)) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white",
                       space = "Lab", na.value = "grey50", guide = "colourbar",
                       aesthetics = "fill", name="Correlation", 
                       limits=c(floor(rng[1]), ceiling(rng[2]))) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  labs(y="Chemical", x="Chemical", title="Truth") 

gridExtra::grid.arrange(t1, t2, nrow=1, widths=c(1.275/3,1.725/3))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S6: the estimated and true entries of loading matrix 
lambda.postm <- Reduce("+", out$Lambda.save)/length(out$Lambda.save)
lambda.postm[,1:2] <- lambda.postm[,c(2,1)] # for better visuatlization
lambda.postm[,1] <- (-1)*lambda.postm[,1] # for better visuatlization

t3 <- ggplot(lambda.postm %>% reshape2::melt()) +
  geom_tile(aes(Var2, Var1, fill=value), color=1) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white", midpoint = 0, 
                       name="Loading") + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20), 
        axis.title = element_text(size=20), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15)) +
  scale_y_reverse(breaks=c(1,10,20,30)) +
  labs(x="Factor", y="Chemical", title=expression(paste("Estimated ", Lambda)))

t4 <- ggplot(cbind(truth$Lambda,matrix(0,30,3)) %>% reshape2::melt()) + # for better visuatlization
  geom_tile(aes(Var2, Var1, fill=value), color=1) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white", midpoint = 0, 
                       name="Loading") + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20), 
        axis.title = element_text(size=20), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15)) +
  scale_y_reverse(breaks=c(1,10,20,30)) +
  scale_x_continuous(breaks=c(1,2)) + 
  labs(x="Factor", y="Chemical", title=expression(paste("True ", Lambda))) 

gridExtra::grid.arrange(t3, t4, nrow=1)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 3: heat map of estimate and true profiles of the mean effect 
pred <- which(is.na(tr_gamma_ij.postm[6:10,6:10]), arr.ind = TRUE)
frames <- data.frame(Var1=pred[,1], Var2=pred[,2])
t5 <- reshape2::melt(gamma_ij.postm[6:10,6:10]) %>% 
  rename(Chem = Var1, Assay = Var2) %>% 
  ggplot() + geom_tile(aes(Assay, Chem, fill=value), color="grey50") + 
  geom_text(aes(Assay, Chem, label=round(value,3))) + 
  scale_fill_gradient2(high = "red3", low = "white", mid = "#E6F598", midpoint = 0.5, limits=c(0,1), na.value = "black") + 
  scale_y_reverse() + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  geom_rect(data=frames, size=1, fill=NA, colour="black",
            aes(xmin=Var2 - 0.5, xmax=Var2 + 0.5, ymin=Var1 - 0.5, ymax=Var1 + 0.5)) + 
  labs(y="Chemical", x="Assay endpoint", title="Estimate", fill="Prob.Active") 

t6 <- reshape2::melt(truth$gamma_ij[6:10,6:10]) %>% 
  rename(Chem = Var1, Assay = Var2) %>% 
  ggplot() + geom_tile(aes(Assay, Chem, fill=factor(value)), color="grey50") + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  scale_y_reverse() +
  scale_fill_manual(values=c("1"="red3", "0"="white")) + 
  labs(y="Chemical", x="Assay endpoint", title="Truth", fill=expression(gamma[ij])) 

gridExtra::grid.arrange(t5, t6, nrow=1, widths=c(1.345/3,1.655/3))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 2: dose-response curves
data_f7 <- list(misdata = misdata, truth = truth)
res_f7 <- list(res_ZIPLL = res_ZIPLL, 
               out = list(gamma_ij.postm = gamma_ij.postm, t_ij.postm = t_ij.postm), 
               prc_res = prc_res)

t7 <- dosres_plot(i = 26, j = 67, data = data_f7, result = res_f7)
t8 <- dosres_plot(i = 5, j = 2, data = data_f7, result = res_f7)
t9 <- dosres_plot(i = 24, j = 104, data = data_f7, result = res_f7)

ggpubr::ggarrange(t7, t8, t9, nrow=1, ncol=3, common.legend=TRUE, legend="bottom", labels="AUTO")

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S7: residuals versus fitted values 
data_f8 <- misdata
res_f8 <- list(res_ZIPLL = res_ZIPLL, prc_res = prc_res)

t10 <- resid_plot(i = 26, j = 67, data = data_f8, result = res_f8)
t11 <- resid_plot(i = 5, j = 2, data = data_f8, result = res_f8)

ggpubr::ggarrange(t10, t11, nrow=1, common.legend=TRUE, legend="bottom", labels="AUTO")
