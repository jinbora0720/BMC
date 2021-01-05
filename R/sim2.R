# Simulation 2

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
sourceCpp(paste0(path, "source/bmc_sampler.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simdata.R"))
source(paste0(path, "source/simulate_data.R"))

########
# Data #
########
m <- 15
J <- 15

set.seed(123)
seedsave <- sample(10000, 100, replace=FALSE)
iterpool <- c(1,2,3,4,5,6,9,10,
              11,13,15,16,17,18,19,20, 
              21,23,24,25,26,27,28,29,30,
              31,32,33,34,35,36,37,38, 
              41,43,44,45,46,47,48,49,50,
              51,52,53,54,58,59,60,
              61)

# hold out bottom right 3 x 3 cells 
missing_idx <- rbind(c(13,13), c(13,14), c(13,15),
                     c(14,13), c(14,14), c(14,15),
                     c(15,13), c(15,14), c(15,15))
test_idx <- apply(missing_idx, 1, function(x) x[2]+J*(x[1]-1))
missing_idx_ij <- apply(missing_idx, 1, function(x) x[1] + m*(x[2]-1))
simout <- simdata(nchem = m, nassay = J, seed = 123)
ZIPLL_missing <- NULL
for (s in 1:nrow(missing_idx)) {
  i <- missing_idx[s,1]
  j <- missing_idx[s,2]
  ZIPLL_missing <- c(ZIPLL_missing, which(simout$dat[,1]==i & simout$dat[,2]==j))
}

#################
# Store results #
#################
rmse = tr_aucg = tt_aucg <- rep(NA, 50)
rmse0 = tr_aucg0 = tt_aucg0 <- rep(NA, 50)
rmsei = tr_aucgi = tt_aucgi <- rep(NA, 50)
rmsej = tr_aucgj = tt_aucgj <- rep(NA, 50)
rmse_tcpl = tr_aucg_tcpl <- rep(NA, 50)
rmse_ZIPLL = tr_aucg_ZIPLL <- rep(NA, 50)

###################
# MCMC parameters #
###################
thin <- 1
burnin <- 10000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

###############
# Simulations #
###############
pb <- txtProgressBar(style=3,width=50)
for (iter in 1:50) {
  seed = seedsave[iterpool[iter]]
  simout = simdata(nchem=m, nassay=J, seed)
  dat = simout$dat
  truth = simout$truth
  
  tt_gamma = truth$gamma_ij[test_idx]
  tr_gamma = truth$gamma_ij[-test_idx]
  
  Y = orgX = idx_j = list()
  m_j = rep(NA, J)
  for (j in 1:J) {
    these = (dat[,2]==j)
    Y[[j]] = dat[these,4]
    orgX[[j]] = dat[these,3]
    idx_j[[j]] = unique(dat[these,1])
    m_j[j] = length(idx_j[[j]])
  }
  
  K_ij = matrix(8, m, J)   
  for (j in 1:J) {
    K_ij[-idx_j[[j]], j] = 0
  }
  
  Start = End = matrix(0, m, J)
  for (j in 1:J) {
    Start[1:m_j[j],j] = cumsum(c(1,K_ij[,j]))[idx_j[[j]]]
    End[1:m_j[j],j] = cumsum(K_ij[,j])[idx_j[[j]]]
  }
  
  p = 7
  X = list() 
  for (j in 1:J) {
    X[[j]] = matrix(NA, max(End[,j]), p)
    for (i in 1:m_j[j]) {
      logc = orgX[[j]][Start[i,j]:End[i,j]] 
      X[[j]][Start[i,j]:End[i,j],] = apply(bs(logc, df=p, intercept=TRUE), 2, 
                                           function(x) scale(x, center=TRUE, scale=FALSE))
    }
  }
  
  simuldata = list(X = X, orgX = orgX, orgY = Y, K_ij = K_ij, 
                   m_j = m_j, idx_j = idx_j, Start = Start, End = End)
  
  misdata = data_missing(simdata = simuldata, missing_idx = missing_idx, seed = seed)
  names(misdata)[5] <- "Y"
  
  X = misdata$X
  Y = misdata$Y
  orgX = misdata$orgX
  Start = misdata$Start
  End = misdata$End 
  idx_j = misdata$idx_j
  m_j = misdata$m_j
  
  #########
  # hyper #
  #########
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
  hyper = list(a = a, R = R)
  
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
             hetero = FALSE, adapt = FALSE, simpler = FALSE, verbose = FALSE)
  
  beta_ij.postm = lapply(out$beta_ij.save, function(x) rowMeans(x,dim=2))
  Fit.postm = lapply(1:J, function(j) {
    sapply(1:m_j[j], function(i) {
      X[[j]][Start[i,j]:End[i,j],]%*%beta_ij.postm[[j]][,i]
    })
  })
  
  ## 1. overall RMSE
  rmse[iter] = sqrt(mean(unlist(mapply('-', Y, Fit.postm))^2))
  
  ## 2. AUC for gamma 
  gamma_ij.save = out$gamma_ij.save 
  z.save = lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
  gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
  gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
  
  tr_pred = prediction(as.numeric(t(gamma_ij.postm))[-test_idx], tr_gamma)
  tr_rocs = performance(tr_pred, measure = "auc")
  tr_aucg[iter] = tr_rocs@y.values[[1]]
  
  ## 3. AUC for predicted gamma 
  tt_pred = prediction(as.numeric(t(gamma_ij.postm))[test_idx], tt_gamma)
  tt_rocs = performance(tt_pred, measure = "auc")
  tt_aucg[iter] = tt_rocs@y.values[[1]]
  
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
  
  ## 3. AUC for predicted gamma 
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
  
  ## 3. AUC for predicted gamma 
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
  
  ## 3. AUC for predicted gamma 
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

# save result
# saveRDS(list(res = res, res0 = res0, resi = resi, resj = resj, 
#              res_tcpl = res_tcpl, res_ZIPLL = res_ZIPLL),
#         paste0(path, "data/sim2_results.rds"))

# call result 
# all <- readRDS(paste0(path, "data/sim2_results.rds"))
# res <- all$res
# res0 <- all$res0
# resi <- all$resi
# resj <- all$resj
# res_tcpl <- all$res_tcpl
# res_ZIPLL <- all$res_ZIPLL

#----------------------------------------------------------------------------------------------------------------------------------
# Table S1: Summary of results from Simulation 2

mres <- matrix(c(round(unlist(lapply(res[1:3], mean)), 3), 
                  round(unlist(lapply(res0[1:3], mean)), 3),
                  round(unlist(lapply(resi[1:3], mean)), 3), 
                  round(unlist(lapply(resj[1:3], mean)), 3),
                  c(round(unlist(lapply(res_ZIPLL[1:2], mean)), 3), NA), 
                  c(round(unlist(lapply(res_tcpl[1:2], mean)), 3), NA)), 3, 6)
rownames(mres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij')
colnames(mres) <- c('BMC', 'BMC0', 'BMCi', 'BMCj', 'ZIPLL', 'tcpl')

sdres <- matrix(c(round(unlist(lapply(res[1:3], sd)), 3), 
                   round(unlist(lapply(res0[1:3], sd)), 3), 
                   round(unlist(lapply(resi[1:3], sd)), 3), 
                   round(unlist(lapply(resj[1:3], sd)), 3), 
                   c(round(unlist(lapply(res_ZIPLL[1:2], sd)), 3), NA), 
                   c(round(unlist(lapply(res_tcpl[1:2], sd)), 3), NA)), 3, 6)
rownames(sdres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij')
colnames(sdres) <- c('BMC', 'BMC0', 'BMCi', 'BMCj', 'ZIPLL', 'tcpl')

cat('mean table\n'); mres
cat('sd table\n'); sdres
