# Reproduce obesity-relevant results 

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 

# source code
path <- "~/Documents/BoraJin2018~/Research/DoseResponse/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler2.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/metadata.R"))
source(paste0(path, "source/simulate_data.R"))
source(paste0(path, "source/dosres_plot.R"))

########
# Data #
########
datalist <- readRDS(paste0(path, "data/datalist.rds"))
obese_data <- datalist$obese_data
hit_vec <- datalist$hit_vec
meta <- metadata(obese_data)
uniq_chnm <- meta$uniq_chnm
m <- length(uniq_chnm)
uniq_aenm <- meta$uniq_aenm
J <- length(uniq_aenm)
simdata <- list(X = meta$X, orgX = lapply(meta$orgData, function(x) x[,1]), 
                orgY = meta$Y, K_ij = meta$K_ij, m_j = meta$m_j, 
                idx_j = meta$idx_j, Start = meta$Start, End = meta$End)

# pool for hold-out cells 
# one can replace this by any non-missing cells 
pred_idx <- readRDS(paste0(path, "data/obese_pred_idx.rds"))
prob_missing <- 0.03
set.seed(330); missing_idx <- pred_idx[sample(nrow(pred_idx), prob_missing*m*J),]
misdata <- data_missing(simdata = simdata, missing_idx = missing_idx, seed = 330)
misdata$Y <- misdata$orgY

###################
# MCMC parameters #
###################
thin <- 10
burnin <- 30000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

###################
# hyperparameters #
###################
dat <- data.frame(resp = unlist(misdata$Y), logc = unlist(misdata$orgX))   
ddat <- dat %>% filter(logc!=0) %>% group_by(logc) %>% 
  summarise(v = sd(resp)) %>% mutate(d = 2*log(v)/logc)
v_d <- ((max(ddat$d) - min(ddat$d))/4)^2
hyper <- list(v_d = v_d)

##################
# initial values #
##################
q <- 11 
init <- list(q = q)

###########
# run bmc # 
###########
out <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init, 
           hetero = TRUE, adapt = TRUE, simpler = FALSE, verbose = TRUE)

###########
# results #
###########
# split result due to file size
# saveRDS(list(out = list(gamma_ij.save = out$gamma_ij.save,
#                         Lambda.save = out$Lambda.save, eta.save = out$eta.save,
#                         Sigj.save = out$Sigj.save, sigj_sq.save = out$sigj_sq.save,
#                         pi_t.save = out$pi_t.save, t_ij.save = out$t_ij.save,
#                         d_ij.save = out$d_ij.save, covMean = out$covMean, accept = out$accept),
#              missing_idx = missing_idx),
#         paste0(path, "data/obese_out1.rds"))
# saveRDS(out$beta_ij.save[1:135], paste0(path, "data/obese_out2.rds"))
# saveRDS(out$beta_ij.save[136:271], paste0(path, "data/obese_out3.rds"))

# call result
# obese_out <- readRDS(paste0(path, "data/obese_out1.rds"))
# out <- obese_out$out
# out$beta_ij.save <- c(readRDS(paste0(path, "data/obese_out2.rds")),
#                       readRDS(paste0(path, "data/obese_out3.rds")))

X <- misdata$X
orgX <- misdata$orgX
Y <- misdata$Y
m_j <- misdata$m_j
idx_j <- misdata$idx_j
Start <- misdata$Start
End <- misdata$End

# arrange posterior samples
set.seed(330)
gamma_ij.save <- out$gamma_ij.save 
z.save <- lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
for (s in 1:save) {
  gamma_ij.save[,,s][missing_idx] <- (z.save[[s]][missing_idx]>0)
}
tr_gamma_ij.postm <- rowMeans(out$gamma_ij.save, dim=2, na.rm=TRUE)
gamma_ij.postm <- rowMeans(gamma_ij.save, dim=2)
t_ij.postm <- rowMeans(out$t_ij.save, dims=2)

uniq_code <- paste("C",gsub("-", "", meta$uniq_casn), sep="")
hit_mat <- hit_vec %>% filter(aenm %in% uniq_aenm, code %in% uniq_code) %>% 
  reshape2::dcast(code~aenm, value.var = "hitc") 
rownames(hit_mat) <- hit_mat[,1]
hit_mat <- as.matrix(hit_mat[uniq_code,uniq_aenm])
hit_mat[which(is.na(gamma_ij.postm))] <- NA

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

prc_res <- list(predY.postL = predY.postL, predY.postU = predY.postU, 
                Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU, 
                Yhat.postm = Yhat.postm, res.postm = res.postm)

# save result
# saveRDS(prc_res, paste0(path, "data/obese_processed.rds"))

# call result
# prc_res <- readRDS(paste0(path, "data/obese_processed.rds"))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S9: Chemical ranks by the average active probability 
# Prob.Active
actprob <- 1-rowMeans((1-out$gamma_ij.save)*(1-out$t_ij.save), dim=2, na.rm=TRUE)

colors <- c("BMC" = "black", "EPA" = "blue") 
shapes <- c("BMC" = 19, "EPA" = 17)

data.frame(act = rowMeans(actprob, na.rm=TRUE), 
           averhit = rowMeans(hit_mat, na.rm=TRUE), chnm=uniq_chnm) %>% 
  ggplot() + geom_point(aes(reorder(chnm, act), act, shape="BMC", color="BMC"), size=5) + 
  geom_point(aes(reorder(chnm,act), averhit, shape="EPA", color="EPA"), size=3) +
  labs(x="Chemical", y="Active probability", title="Obesity", color=NULL, shape=NULL) +
  theme_minimal() + theme_bw() + 
  scale_color_manual(values = colors) + 
  scale_shape_manual(values = shapes) + 
  theme(axis.text.x = element_text(angle=65, hjust=1, size=15), 
        legend.text = element_text(size=15), 
        axis.title = element_text(size=15), 
        plot.title = element_text(size=20)) 

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S12: Ranks of assay endpoints likely to be activated by the top 5 chemicals 
obese_ascend_chnm <- uniq_chnm[order(rowMeans(actprob, na.rm=TRUE))]
obese_active_aenm <- data.frame(act = colMeans(actprob[which(uniq_chnm %in% obese_ascend_chnm[c(26:30)]),], na.rm=TRUE), 
                                averhit = colMeans(hit_mat, na.rm=TRUE), aenm=uniq_aenm) %>% 
  filter(0.9 < act) %>% select(aenm)
# saveRDS(obese_active_aenm, paste0(path, "data/obese_active_aenm.rds"))
neuro_active_aenm <- readRDS(paste0(path, "data/neuro_active_aenm.rds"))

data.frame(act = colMeans(actprob[which(uniq_chnm %in% obese_ascend_chnm[c(26:30)]),], na.rm=TRUE), 
           averhit = colMeans(hit_mat, na.rm=TRUE), aenm=uniq_aenm, 
           common = uniq_aenm %in% neuro_active_aenm$aenm) %>% 
  filter(0.9 < act) %>% 
  ggplot() + geom_point(aes(reorder(aenm, act), act, color=factor(common), shape=factor(common)), size=5) + 
  scale_color_manual(values=c("TRUE"="red","FALSE"="blue")) + 
  labs(title="Obesity", x="Assay endpoint", y="Activated probability", 
       color="Involved in both diseases", 
       shape="Involved in both diseases") +
  theme_minimal() + theme_bw() + 
  theme(axis.text.x = element_text(angle=75, hjust=1, size=10), 
        legend.position = "bottom", 
        legend.title = element_text(size=15),
        legend.text = element_text(size=15), 
        axis.title = element_text(size=15), 
        plot.title = element_text(size=20))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (top row): Results for active pairs 
data_f <- list(misdata = misdata, meta = meta)
result_f <- list(out = list(gamma_ij.postm = gamma_ij.postm, t_ij.postm = t_ij.postm), 
                 hit_mat = hit_mat, prc_res = prc_res)

t1 <- dosres_plot(i = which(uniq_chnm == "Bisphenol A"), 
                  j = which(uniq_aenm == "NVS_TR_hDAT"),
                  data_f, result_f, realdata = TRUE)
t2 <- dosres_plot(i = which(uniq_chnm == "p,p'-DDE"), 
                  j = which(uniq_aenm == "TOX21_ERR_Antagonist"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t1, t2, nrow=1)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (middle row): Results for active pairs with mean effect
t3 <- dosres_plot(i = which(uniq_chnm == "Dibutyl phthalate"), 
                  j = which(uniq_aenm == "ATG_PPARg_TRANS_up"),
                  data_f, result_f, realdata = TRUE)
t4 <- dosres_plot(i = which(uniq_chnm == "Bisphenol A"), 
                  j = which(uniq_aenm == "ATG_PPARg_TRANS_up"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t3, t4, nrow=1)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (bottom row): Results for active pairs with variance effect
t5 <- dosres_plot(i = which(uniq_chnm == "Diisobutyl phthalate"), 
                  j = which(uniq_aenm == "TOX21_CAR_Antagonist"),
                  data_f, result_f, realdata = TRUE)
t6 <- dosres_plot(i = which(uniq_chnm == "Di(2-ethylhexyl) phthalate"), 
                  j = which(uniq_aenm == "TOX21_Aromatase_Inhibition"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t5, t6, nrow=1)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 5: Results for hold-out pairs
data_f12 <- list(simdata = simdata, meta = meta)
result_f12 <- gamma_ij.postm

t7 <- dosres_plot(i = which(uniq_chnm == "p,p'-DDE"), 
                  j = which(uniq_aenm == "BSK_SAg_CD38_down"),
                  data_f12, result_f12, realdata = TRUE, pred = TRUE)
t8 <- dosres_plot(i = which(uniq_chnm == "p,p'-DDE"), 
                  j = which(uniq_aenm == "BSK_3C_IL8_down"),
                  data_f12, result_f12, realdata = TRUE, pred = TRUE)

gridExtra::grid.arrange(t7, t8, nrow=1)

# "The chemical p,p'-DDE: average probability 0.808"
round(mean(gamma_ij.postm[which(uniq_chnm == "p,p'-DDE"),], na.rm = TRUE), 3)
# "BSK_SAg_CD38_down: the average posterior probability 0.874"  
round(mean(gamma_ij.postm[,which(uniq_aenm == "BSK_SAg_CD38_down")], na.rm = TRUE), 3)
# "BSK_3C_IL8_down: the average posterior probability 0.694".
round(mean(gamma_ij.postm[,which(uniq_aenm == "BSK_3C_IL8_down")], na.rm = TRUE), 3)

# 95% highest density intervals for the estimated and predicted $\gamma_{ij}$'s are 
# (0.205, 1) and (0.510, 0.863), respectively.
HDInterval::hdi(as.numeric(tr_gamma_ij.postm))
HDInterval::hdi(as.numeric(gamma_ij.postm[missing_idx]))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S8: Results for hold-out pairs
t9 <- dosres_plot(i = which(uniq_chnm == "Dichlorodiphenyltrichloroethane"), 
                  j = which(uniq_aenm == "ATG_IR1_CIS_dn"),
                  data_f12, result_f12, realdata = TRUE, pred = TRUE)
t10 <- dosres_plot(i = which(uniq_chnm == "Benzyl butyl phthalate"), 
                   j = which(uniq_aenm == "ATG_PPARg_TRANS_dn"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE)
t11 <- dosres_plot(i = which(uniq_chnm == "Bisphenol A"), 
                   j = which(uniq_aenm == "TOX21_ERb_BLA_Antagonist_ratio"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE)
t12 <- dosres_plot(i = which(uniq_chnm == "Di(2-ethylhexyl) phthalate"), 
                   j = which(uniq_aenm == "TOX21_FXR_BLA_antagonist_ratio"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE)

gridExtra::grid.arrange(t9, t10, t11, t12, nrow=2)
