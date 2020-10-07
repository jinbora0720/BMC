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
path <- "~/Documents/GitHub/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler.cpp"))
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
rmse = tr_aucg = tt_aucg = tr_auct = tt_auct <- rep(NA, 50)
rmse0 = tr_aucg0 = tt_aucg0 = tr_auct0 = tt_auct0 <- rep(NA, 50)
rmsei = tr_aucgi = tt_aucgi = tr_aucti = tt_aucti <- rep(NA, 50)
rmsej = tr_aucgj = tt_aucgj = tr_auctj = tt_auctj <- rep(NA, 50)
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
  gendata = generate_data(nchem=m, nassay=J, d=3, seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.03
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  
  X = misdata$X
  orgX = misdata$orgX
  Y = misdata$orgY
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
  
  ### check from here
  tr_pred = prediction(as.numeric(gamma_ij.postm)[-missing_idx_col], as.numeric(truth$gamma_ij)[-missing_idx_col])
  tr_rocs = performance(tr_pred, measure = "auc")
  tr_aucg = tr_rocs@y.values[[1]]
  
  ## 3. AUC for predicted gamma 
  tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), as.numeric(truth$gamma_ij[missing_idx]))
  tt_rocs = performance(tt_pred, measure = "auc")
  tt_aucg = tt_rocs@y.values[[1]]
  
  ## 4. AUC for t_ij
  t_ij.postm = rowMeans(out$t_ij.save, dims=2)
  pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], as.numeric(truth$t_ij)[-missing_idx_col])
  rocs_t = performance(pred_t, measure = "auc")
  auct = rocs_t@y.values[[1]]
  
  ########
  # bmc0 # 
  ########
  out0 <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = FALSE, adapt = FALSE, simpler = "bmc0", verbose = FALSE)
  
  beta_ij.postm0 = lapply(out0$beta_ij.save, function(x) rowMeans(x,dim=2))
  Fit.postm0 = lapply(1:J, function(j) {
    sapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%beta_ij.postm0[[j]][,i]
    })
  })
  
  ## 1. overall RMSE
  rmse0[iter] = sqrt(mean(unlist(mapply('-', Y, Fit.postm0))^2))
  
  ## 2. AUC for gamma 
  gamma_ij.save0 = out0$gamma_ij.save 
  gamma_ij.save0[is.na(gamma_ij.save0)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, out0$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postm0 = rowMeans(gamma_ij.save0, dim=2, na.rm=TRUE)
  
  tr_pred0 = prediction(as.numeric(t(gamma_ij.postm0))[-test_idx], tr_gamma)
  tr_rocs0 = performance(tr_pred0, measure = "auc")
  tr_aucg0[iter] = tr_rocs0@y.values[[1]]
  
  ## 4. AUC for predicted gamma 
  tt_pred0 = prediction(as.numeric(t(gamma_ij.postm0))[test_idx], tt_gamma)
  tt_rocs0 = performance(tt_pred0, measure = "auc")
  tt_aucg0[iter] = tt_rocs0@y.values[[1]]
  
  ########
  # bmci # 
  ########
  outi <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = FALSE, adapt = FALSE, simpler = "bmci", verbose = FALSE)
  
  beta_ij.postmi = lapply(outi$beta_ij.save, function(x) rowMeans(x,dim=2))
  Fit.postmi = lapply(1:J, function(j) {
    sapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%beta_ij.postmi[[j]][,i]
    })
  })
  
  ## 1. overall RMSE
  rmsei[iter] = sqrt(mean(unlist(mapply('-', Y, Fit.postmi))^2))
  
  ## 2. AUC for gamma 
  gamma_ij.savei = outi$gamma_ij.save 
  gamma_ij.savei[is.na(gamma_ij.savei)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outi$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postmi = rowMeans(gamma_ij.savei, dim=2, na.rm=TRUE)
  
  tr_predi = prediction(as.numeric(t(gamma_ij.postmi))[-test_idx], tr_gamma)
  tr_rocsi = performance(tr_predi, measure = "auc")
  tr_aucgi[iter] = tr_rocsi@y.values[[1]]
  
  ## 4. AUC for predicted gamma 
  tt_predi = prediction(as.numeric(t(gamma_ij.postmi))[test_idx], tt_gamma)
  tt_rocsi = performance(tt_predi, measure = "auc")
  tt_aucgi[iter] = tt_rocsi@y.values[[1]]
  
  ########
  # bmcj # 
  ########
  outj <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
              hetero = FALSE, adapt = FALSE, simpler = "bmcj", verbose = FALSE)
  
  beta_ij.postmj = lapply(outj$beta_ij.save, function(x) rowMeans(x,dim=2))
  Fit.postmj = lapply(1:J, function(j) {
    sapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%beta_ij.postmj[[j]][,i]
    })
  })
  
  ## 1. overall RMSE
  rmsej[iter] = sqrt(mean(unlist(mapply('-', Y, Fit.postmj))^2))
  
  ## 2. AUC for gamma 
  gamma_ij.savej = outj$gamma_ij.save 
  gamma_ij.savej[is.na(gamma_ij.savej)] = sapply(1:save, function(s) rbinom(nrow(missing_idx), 1, outj$pi_ij.save[,,s][missing_idx]))
  gamma_ij.postmj = rowMeans(gamma_ij.savej, dim=2, na.rm=TRUE)
  
  tr_predj = prediction(as.numeric(t(gamma_ij.postmj))[-test_idx], tr_gamma)
  tr_rocsj = performance(tr_predj, measure = "auc")
  tr_aucgj[iter] = tr_rocsj@y.values[[1]]
  
  ## 4. AUC for predicted gamma 
  tt_predj = prediction(as.numeric(t(gamma_ij.postmj))[test_idx], tt_gamma)
  tt_rocsj = performance(tt_predj, measure = "auc")
  tt_aucgj[iter] = tt_rocsj@y.values[[1]]
  
  ########
  # tcpl #
  ########
  rmse_TCPL = hitc = matrix(NA, m, J)
  for (j in 1:J) {
    for (i in 1:m) {
      these <- (dat[,1]==i & dat[,2]==j)
      params <- tcplFit(logc = dat[these,3],
                        resp = dat[these,4],
                        bmad = (median(dat[which(dat[,2]==j),4]) + sd(dat[which(dat[,2]==j),4]))/3, bidirectional=TRUE)
      rmse_TCPL[i,j] = c(params$cnst_rmse, params$hill_rmse, params$gnls_rmse)[which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))]
      hitc[i,j] = ifelse(which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))==1, 0, 1)
    }
  }
  
  ## 1. overall RMSE
  rmse_tcpl[iter] = mean(rmse_TCPL[-missing_idx_ij], na.rm=TRUE)
  
  ## 2. AUC for gamma
  tr_pred_tcpl = prediction(as.numeric(t(hitc))[-test_idx], tr_gamma)
  tr_rocs_tcpl = performance(tr_pred_tcpl, measure = "auc")
  tr_aucg_tcpl[iter] = tr_rocs_tcpl@y.values[[1]]
  
  #########
  # ZIPLL #
  #########
  fit <- ZIPLL(dat, nitter=burnin+thin*save, burnin=burnin)
  
  ## 1. overall RMSE 
  rmse_ZIPLL[iter] = sqrt(mean((fit$dat[-ZIPLL_missing,4] - fit$dat[-ZIPLL_missing,5])^2, na.rm=TRUE))
  
  ## 2. AUC for train gamma 
  tr_pred_ZIPLL = prediction(fit$parms[,7][-test_idx], tr_gamma)
  tr_rocs_ZIPLL = performance(tr_pred_ZIPLL, measure = "auc")
  tr_aucg_ZIPLL[iter] = tr_rocs_ZIPLL@y.values[[1]]
  
  setTxtProgressBar(pb, iter/50) 
}
close(pb)

################
# Save results #
################
res <- list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
            seedsave = seedsave[iterpool], out = out)
res0 <- list(rmse = rmse0, tr_aucg = tr_aucg0, tt_aucg = tt_aucg0, 
             seedsave = seedsave[iterpool], out = out0)
resi <- list(rmse = rmsei, tr_aucg = tr_aucgi, tt_aucg = tt_aucgi, 
             seedsave = seedsave[iterpool], out = outi)
resj <- list(rmse = rmsej, tr_aucg = tr_aucgj, tt_aucg = tt_aucgj, 
             seedsave = seedsave[iterpool], out = outj)
res_tcpl <- list(rmse = rmse_tcpl, tr_aucg = tr_aucg_tcpl, seedsave = seedsave[iterpool])
res_ZIPLL <- list(rmse = rmse_ZIPLL, tr_aucg = tr_aucg_ZIPLL, seedsave = seedsave[iterpool])

saveRDS(list(res = res, res0 = res0, resi = resi, resj = resj, 
             res_tcpl = res_tcpl, res_ZIPLL = res_ZIPLL),
        file.path(path, "data/sim2_results.rds"))
