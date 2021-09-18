# Simulation 1
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse)
library(splines) 
library(ROCR)

# source code
# path <- "~/BMC/BMC/"
path <- "~/BMC/"
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
q <- 2
alpha <- c(-0.1,1.2)
xi <- 0

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
  gendata = generate_data(m=m, J=J, d=3, q=q, 
                          xi=xi, alpha=alpha, delta_mean=1.5, seed=seed)
  simdata = gendata$simdata
  truth = gendata$truth
  
  prob_missing = 0.05
  misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
  misdata$Y <- misdata$orgY
  missing_idx = misdata$missing_idx 
  missing_idx_col = misdata$missing_idx_col 
  test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
  
  tt_gamma = as.numeric(truth$gamma_ij[missing_idx])
  tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
  tt_t = as.numeric(truth$t_ij[missing_idx])
  tr_t = as.numeric(truth$t_ij)[-missing_idx_col]
  actprob = 1-(1-truth$gamma_ij)*(1-truth$t_ij)
  tt_actprob = as.numeric(actprob[missing_idx])
  tr_actprob = as.numeric(actprob)[-missing_idx_col]
  
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
  
  ## informative prior
  mu_xi = qnorm(sum(truth$gamma_ij)/(m*J)) # xi 
  sigsq_xi = 1
  
  mu_alpha = c(0,1)
  sigsq_alpha = c(5,5)
  
  hyper = list(a = a, R = R, v_d = v_d,
               mu_xi = mu_xi, sigsq_xi = sigsq_xi,
               mu_alpha = mu_alpha, sigsq_alpha = sigsq_alpha)
  
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
                 hetero = TRUE, hetero_simpler = FALSE,
                 gamma_simpler = FALSE, mgsp_adapt = FALSE, 
                 update_xi = TRUE, apply_cutoff = FALSE, verbose = TRUE)
  
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
  z.save = list()
  for (i in 1:MCMC$save) {
    z.save[[i]] = matrix(rnorm(m*J, le.save[[i]] + out$xi.save[i], 1), m, J)
  }
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
  t_ij.save = out$t_ij.save 
  lea = lapply(1:save, function(s) out$alpha.save[1,s] + out$alpha.save[2,s]*le.save[[s]])
  u.save = lapply(lea, function(x) matrix(rnorm(m*J, x, 1), m, J))
  t_ij.save[is.na(t_ij.save)] = sapply(u.save, function(x) x[missing_idx]>0)
  t_ij.postm = rowMeans(t_ij.save, dim=2, na.rm=TRUE)
  
  tr_pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], tr_t)
  tr_rocs_t = performance(tr_pred_t, measure = "auc")
  tr_auct = tr_rocs_t@y.values[[1]]
  
  ## 5. AUC for predicted t_ij 
  tt_pred_t = prediction(as.numeric(t_ij.postm[missing_idx]), tt_t)
  tt_rocs_t = performance(tt_pred_t, measure = "auc")
  tt_auct = tt_rocs_t@y.values[[1]]

  ## 6. AUC for gamma_ij = 1 or t_ij = 1
  actprob.postm = 1-rowMeans((1-gamma_ij.save)*(1-t_ij.save), dim=2)
  tr_pred_ap = prediction(as.numeric(actprob.postm)[-missing_idx_col], tr_actprob)
  tr_rocs_ap = performance(tr_pred_ap, measure = "auc")
  tr_aucap = tr_rocs_ap@y.values[[1]]
  
  ## 7. AUC for predicted gamma_ij = 1 or t_ij = 1
  tt_pred_ap = prediction(as.numeric(actprob.postm[missing_idx]), tt_actprob)
  tt_rocs_ap = performance(tt_pred_ap, measure = "auc")
  tt_aucap = tt_rocs_ap@y.values[[1]]
  
  ## 8. alphas 
  alpha.postm <- rowMeans(out$alpha.save)
  
  ## 9. xi
  xi.postm <- mean(out$xi.save)
  
################
# Save results #
################
  out_name <- paste0(path, "data/sim1_new2_BMCnew_res_", iter, ".RDS")
  saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg,
               tr_auct = tr_auct, tt_auct = tt_auct,
               tr_aucap = tr_aucap, tt_aucap = tt_aucap,
               alpha = alpha.postm, xi = xi.postm, seed = seed), file.path(out_name))


