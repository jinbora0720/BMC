# Simulation 4: behaviour of Pr(gamma_ij = 1 U t_ij = 1)
# highly correlated
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 
library(ROCR)

# source code
path <- "~/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simulate_data.R"))

set.seed(1)
seedsave <- sample(10000, 30, replace=FALSE)

########
# Data #
########
m <- 30
J <- 10
q <- 1
alpha <- c(0.3,1)
xi <- 0.8

# ################
# # Save results #
# ################
# rmse = tr_aucg = tt_aucg = tr_auct = tt_auct = tr_aucap = tt_aucap = xi.postm = matrix(0, nrow = 4, ncol = 30)
# alpha.postm = array(0, dim = c(2,4,30))
# 
# ###################
# # MCMC parameters #
# ###################
# thin <- 5
# burnin <- 5000
# save <- 1000
# MCMC <- list(thin = thin, burnin = burnin, save = save)
# 
# probs = c(0.1, 0.2, 0.3, 0.5)
# 
# pb <- txtProgressBar(style=3,width=50)
# for (sim in 1:30) {
#   seed = seedsave[sim]
#   
#   ########
#   # Data #
#   ########
#   gendata = generate_data(m, J, d=5, q, xi=xi, alpha=alpha, 
#                           delta_mean=1, df=0.1, seed=seed) # df = 10 (new2)
#   simdata = gendata$simdata
#   truth = gendata$truth
#   simdata$rxf = truth$rxf
#   
#   for (iter in 1:4) {
#     prob_missing = probs[iter]
#     misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
#     misdata$Y = misdata$orgY
#     missing_idx = misdata$missing_idx 
#     missing_idx_col = misdata$missing_idx_col 
#     test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))
#     
#     tt_gamma = as.numeric(truth$gamma_ij[missing_idx])
#     tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
#     tt_t = as.numeric(truth$t_ij[missing_idx])
#     tr_t = as.numeric(truth$t_ij)[-missing_idx_col]
#     actprob = 1-(1-truth$gamma_ij)*(1-truth$t_ij)
#     tt_actprob = as.numeric(actprob[missing_idx])
#     tr_actprob = as.numeric(actprob)[-missing_idx_col]
#     
#     X = misdata$X
#     orgX = misdata$orgX
#     Y = misdata$Y
#     m_j = misdata$m_j
#     idx_j = misdata$idx_j
#     Start = misdata$Start
#     End = misdata$End
#     
#     #########
#     # hyper #
#     #########
#     p = 7
#     a = p+2
#     Sigj = sapply(1:J, function(j) {
#       covtemp = sapply(1:m_j[j], function(i) {
#         coef(lm(Y[[j]][Start[i,j]:End[i,j]]~
#                   bs(orgX[[j]][Start[i,j]:End[i,j]],df=p,intercept=TRUE)-1))    
#       }) %>% t(.) 
#       if (nrow(na.omit(covtemp))>1) { cov(covtemp, use="complete.obs") }
#       else {rep(0, p*p)}
#     }) 
#     Sigj[,colSums(Sigj)==0] = rowSums(Sigj)/(J-sum(colSums(Sigj)==0))
#     invSigj = apply(Sigj, 2, function(x) c(ginv(matrix(x, p, p))))
#     R = matrix((a-p-1)*rowMeans(Sigj, na.rm=TRUE), p, p)
#     v_d = var(unlist(Y))
#     
#     ## informative prior
#     mu_xi = qnorm(sum(truth$gamma_ij)/(m*J)) # xi 
#     sigsq_xi = 1
#     
#     mu_alpha = c(qnorm(sum(truth$t_ij)/(m*J)), 
#                  sum(truth$gamma_ij == truth$t_ij)/(m*J)) # alpha
#     sigsq_alpha = c(1,1)
#     
#     hyper = list(a = a, R = R, v_d = v_d,
#                  mu_xi = mu_xi, sigsq_xi = sigsq_xi,
#                  mu_alpha = mu_alpha, sigsq_alpha = sigsq_alpha)
#     
#     ##################
#     # initial values #
#     ##################
#     ols = lm(scale(unlist(Y))~bs(unlist(orgX), df=p)-1)
#     beta_ij = list()
#     for (j in 1:J) {
#       beta_ij[[j]] = matrix(coef(ols), nrow=p, ncol=m_j[j])
#     }
#     init = list(q = q, beta_ij = beta_ij)
#     
#     #######
#     # BMC #
#     #######
#     out <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init,
#                hetero = TRUE, hetero_simpler = FALSE,
#                gamma_simpler = FALSE, mgsp_adapt = FALSE,
#                update_xi = TRUE, apply_cutoff = FALSE, verbose = FALSE)
#     
#     # ##########
#     # # checks #
#     # ##########
#     # rowMeans(out$alpha.save); alpha
#     # apply(out$alpha.save, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
#     # plot(as.numeric(rowMeans(out$d_ij.save, dims = 2) - truth$d_ij))
#     # plot(rowMeans(out$sigj_sq.save) - truth$sigj_sq)
#     # mean(out$xi.save); xi
#     # quantile(out$xi.save, probs = c(0.025, 0.975));
#     # 
#     # plot(out$xi.save, type="l")
#     # plot(out$alpha.save[1,], type="l")
#     # plot(out$alpha.save[2,], type="l")
#     # plot(out$sigj_sq.save[1,], type="l")
#     # plot(out$d_ij.save[sample.int(m,1),1,], type="l")
#     # for (i in 1:p) {
#     #   plot(out$beta_ij.save[[1]][i,sample.int(m,1),], type="l")
#     # }
#     
#     ###########
#     # results #
#     ###########
#     ### dose-response
#     Yhat.save = lapply(1:J, function(j) {
#       sapply(1:save, function(s) {
#         lapply(1:m_j[j], function(i) {
#           (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s])/
#             exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
#         }) %>% unlist() 
#       }) 
#     })
#     Yhat.postm = lapply(Yhat.save, rowMeans)
#     
#     Ytil.save = lapply(1:J, function(j) {
#       sapply(1:save, function(s) {
#         lapply(1:m_j[j], function(i) {
#           Y[[j]][Start[i,j]:End[i,j]]/
#             exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
#         }) %>% unlist() 
#       }) 
#     })
#     res.postm = lapply(mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE), rowMeans)
#     
#     ## 1. overall RMSE
#     rmse[iter,sim] = sqrt(mean(unlist(res.postm)^2))
#     
#     ## 2. AUC for gamma 
#     gamma_ij.save = out$gamma_ij.save 
#     le.save = mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE)
#     z.save = list()
#     for (i in 1:MCMC$save) {
#       z.save[[i]] = matrix(rnorm(m*J, le.save[[i]] + out$xi.save[i], 1), m, J)
#     }
#     gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
#     gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
#     
#     tr_pred = prediction(as.numeric(gamma_ij.postm)[-missing_idx_col], tr_gamma)
#     tr_rocs = performance(tr_pred, measure = "auc")
#     tr_aucg[iter,sim] = tr_rocs@y.values[[1]]
#     
#     ## 3. AUC for predicted gamma 
#     tt_pred = prediction(as.numeric(gamma_ij.postm[missing_idx]), tt_gamma)
#     tt_rocs = performance(tt_pred, measure = "auc")
#     tt_aucg[iter,sim] = tt_rocs@y.values[[1]]
#     
#     ## 4. AUC for t_ij
#     t_ij.save = out$t_ij.save 
#     lea = lapply(1:save, function(s) out$alpha.save[1,s] + out$alpha.save[2,s]*le.save[[s]])
#     u.save = lapply(lea, function(x) matrix(rnorm(m*J, x, 1), m, J))
#     t_ij.save[is.na(t_ij.save)] = sapply(u.save, function(x) x[missing_idx]>0)
#     t_ij.postm = rowMeans(t_ij.save, dim=2, na.rm=TRUE)
#     
#     tr_pred_t = prediction(as.numeric(t_ij.postm)[-missing_idx_col], tr_t)
#     tr_rocs_t = performance(tr_pred_t, measure = "auc")
#     tr_auct[iter,sim] = tr_rocs_t@y.values[[1]]
#     
#     ## 5. AUC for predicted t_ij 
#     tt_pred_t = prediction(as.numeric(t_ij.postm[missing_idx]), tt_t)
#     tt_rocs_t = performance(tt_pred_t, measure = "auc")
#     tt_auct[iter,sim] = tt_rocs_t@y.values[[1]]
#     
#     ## 6. alphas 
#     alpha.postm[,iter,sim] = rowMeans(out$alpha.save)
#     
#     ## 7. xi 
#     xi.postm[iter,sim] = mean(out$xi.save)
#     
#     ## 8. AUC for gamma_ij = 1 or t_ij = 1
#     actprob.postm = 1-rowMeans((1-gamma_ij.save)*(1-t_ij.save), dim=2)
#     tr_pred_ap = prediction(as.numeric(actprob.postm)[-missing_idx_col], tr_actprob)
#     tr_rocs_ap = performance(tr_pred_ap, measure = "auc")
#     tr_aucap[iter,sim] = tr_rocs_ap@y.values[[1]]
#     
#     tt_pred_ap = prediction(as.numeric(actprob.postm[missing_idx]), tt_actprob)
#     tt_rocs_ap = performance(tt_pred_ap, measure = "auc")
#     tt_aucap[iter,sim] = tt_rocs_ap@y.values[[1]]
#   }
#   setTxtProgressBar(pb, sim/30) 
# }
# close(pb)

################
# Save results #
################
out_name = paste0(path, "data/sim4_highcor.RDS")
# saveRDS(list(out = out,
#              rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
#              tr_auct = tr_auct, tt_auct = tt_auct, 
#              tr_aucap = tr_aucap, tt_aucap = tt_aucap, 
#              alpha = alpha.postm, xi = xi.postm, seed = seedsave), file.path(out_name))

###########
# results #
###########
# Figure S10: Heat map of correlation matrix when chemicals are highly correlated
seed = seedsave[1]
truth = generate_data(m, J, d=5, q, xi=xi, alpha=alpha, 
                        delta_mean=1, df=0.1, seed=seed)$truth
z_cor_T <- cov2cor(tcrossprod(truth$Lambda) + diag(1, m))
ggplot(z_cor_T %>% reshape2::melt())+
  geom_tile(aes(reorder(Var2, value), reorder(Var1, value), fill=value)) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white",
                       space = "Lab", na.value = "grey50", guide = "colourbar",
                       aesthetics = "fill", name="Correlation") +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  labs(y="Chemical", x="Chemical", title="Truth")
# ggsave("~/BMC/Figure/sim4_highcorrelation.eps")
  
#----------------------------------------------------------------------------------------------------------------------------------
# Table S2: Summary of results from Simulation 4 under highly correlated structure
res = readRDS(out_name)
tab = rbind(do.call("rbind", lapply(res[2:8], rowMeans)), 
            rowMeans(res$alpha, dim=2), rowMeans(res$xi)) 
rownames(tab)[8:10] = c("alpha0", "alpha1", "xi")
colnames(tab) = c("10%", "20%", "30%", "50%")

tab_sd = rbind(do.call("rbind", lapply(res[2:8], function(x) apply(x, 1, sd))), 
               apply(res$alpha, c(1,2), sd), apply(res$xi, 1, sd)) 
rownames(tab_sd)[8:10] = c("alpha0", "alpha1", "xi")
colnames(tab_sd) = c("10%", "20%", "30%", "50%")

cat('mean table\n'); round(tab[c(6,7,2,3,4,5,1,8,9,10),], 3)
cat('sd table\n'); round(tab_sd[c(6,7,2,3,4,5,1,8,9,10),], 3)
