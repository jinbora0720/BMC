# Reproduce obesity-relevant results 
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
theme_set(theme_bw())
library(splines) 
library(tcpl)

# source code
path <- "~/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simulate_data.R"))
source(paste0(path, "source/metadata.R"))
source(paste0(path, "source/dosres_plot.R"))
sourceCpp(paste0(path, "source/list_mean.cpp")) # From Michele Peruzzi: https://github.com/mkln/meshgp/blob/master/src/list_mean.cpp
sourceCpp(paste0(path, "source/msf.cpp"))

########
# Data #
########
datalist <- readRDS(paste0(path, "data/datalist.rds"))
obese_data <- datalist$obese_data
hit_vec <- datalist$hit_vec

# standardise cutoff by each pair 
assay_cutoff <- readRDS(paste0(path, "data/assay_cutoff.RDS"))
# meta <- metadata(obese_data, aeid_coff = assay_cutoff)
meta <- readRDS(paste0(path, "data/obese_meta.RDS"))
meta$coff_ij

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

## informative prior
uniq_code <- paste("C",gsub("-", "", meta$uniq_casn), sep="")
hit_mat <- hit_vec %>% filter(aenm %in% uniq_aenm, code %in% uniq_code) %>%
  reshape2::dcast(code~aenm, value.var = "hitc")
rownames(hit_mat) <- hit_mat[,1]
hit_mat <- as.matrix(hit_mat[uniq_code,uniq_aenm])

mu_xi <- qnorm(mean(hit_mat, na.rm = TRUE)) # xi 
sigsq_xi <- 1

mu_alpha <- c(0,1) # alpha
sigsq_alpha <- c(5,5)

hyper = list(v_d = v_d, cutoff = meta$coff_ij,
             mu_xi = mu_xi, sigsq_xi = sigsq_xi,
             mu_alpha = mu_alpha, sigsq_alpha = sigsq_alpha)

##################
# initial values #
##################
q <- 11 
init <- list(q = q)

###########
# run bmc # 
###########
# out <- bmc(Data = misdata, MCMC = MCMC, hyper = hyper, init = init,
#            hetero = TRUE, hetero_simpler = FALSE,
#            gamma_simpler = FALSE, mgsp_adapt = FALSE,
#            update_xi = TRUE, apply_cutoff = TRUE, verbose = TRUE)

###########
# results #
###########
# split result due to file size
# saveRDS(list(out = list(gamma_ij.save = out$gamma_ij.save,
#                         kappa_ij.save = out$kappa_ij.save,
#                         Sigj.save = out$Sigj.save, 
#                         sigj_sq.save = out$sigj_sq.save,
#                         xi.save = out$xi.save, 
#                         Lambda.save = out$Lambda.save, 
#                         eta.save = out$eta.save,
#                         covMean = out$covMean,
#                         alpha.save = out$alpha.save, 
#                         t_ij.save = out$t_ij.save,
#                         d_ij.save = out$d_ij.save,
#                         accept = out$accept),
#              missing_idx = missing_idx),
#         paste0(path, "data/obese_out1.RDS"))
# saveRDS(out$beta_ij.save[1:135], paste0(path, "data/obese_out2.RDS"))
# saveRDS(out$beta_ij.save[136:271], paste0(path, "data/obese_out3.RDS"))

# call result
obese_out <- readRDS(paste0(path, "data/obese_out1.RDS"))
out <- obese_out$out
out$beta_ij.save <- c(readRDS(paste0(path, "data/obese_out2.RDS")),
                      readRDS(paste0(path, "data/obese_out3.RDS")))

X <- misdata$X
orgX <- misdata$orgX
Y <- misdata$Y
m_j <- misdata$m_j
idx_j <- misdata$idx_j
Start <- misdata$Start
End <- misdata$End

# convergence check 
plot(out$xi.save, type="l")
plot(out$alpha.save[1,], type="l")
plot(out$alpha.save[2,], type="l")
plot(out$sigj_sq.save[sample.int(J,1),], type="l")
hist(out$gamma_ij.save)
hist(out$kappa_ij.save)
hist(out$t_ij.save)

#----------------------------------------------------------------------------------------------------------------------------------
# activity calls 
## mean activity
set.seed(330)
gamma_ij.save <- out$gamma_ij.save
le.save <- mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE)
z.save <- list()
for (s in 1:save) {
  z.save[[s]] <- matrix(rnorm(m*J, le.save[[s]] + out$xi.save[s], 1), m, J)
}
missing_all <- is.na(gamma_ij.save[,,1])
for (s in 1:save) {
  gamma_ij.save[,,s][missing_all] <- (z.save[[s]][missing_all]>0)
}
tr_gamma_ij.postm <- rowMeans(out$gamma_ij.save, dim=2, na.rm=TRUE)
gamma_ij.postm <- rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)

## mean activity
set.seed(330)
kappa_ij.save <- out$kappa_ij.save
missing_all <- is.na(kappa_ij.save[,,1])
for (s in 1:save) {
  kappa_ij.save[,,s][missing_all] <- gamma_ij.save[,,s][missing_all]
}
tr_kappa_ij.postm <- rowMeans(out$kappa_ij.save, dim=2, na.rm=TRUE)
kappa_ij.postm <- rowMeans(kappa_ij.save, dim=2, na.rm=TRUE)

## xi coefficients
xi.postm <- mean(out$xi.save)
xi.postm %>% round(3)
round(quantile(out$xi.save, probs = c(0.025, 0.975)), 3)

## variance activity
set.seed(330)
t_ij.save <- out$t_ij.save
lea <- lapply(1:save, function(s) out$alpha.save[1,s] + out$alpha.save[2,s]*le.save[[s]])
u.save <- lapply(lea, function(x) matrix(rnorm(m*J, x, 1), m, J))
for (s in 1:save) {
  t_ij.save[,,s][missing_all] <- (u.save[[s]][missing_all]>0)
}
tr_t_ij.postm <- rowMeans(out$t_ij.save, dim=2, na.rm=TRUE)
t_ij.postm <- rowMeans(t_ij.save, dim=2, na.rm=TRUE)

## Prob.Active
tr_actprob <- 1-rowMeans((1-out$gamma_ij.save)*(1-out$t_ij.save), dim=2, na.rm=TRUE)
tr_actprob_cutoff <- 1-rowMeans((1-out$kappa_ij.save)*(1-out$t_ij.save), dim=2, na.rm=TRUE)
actprob <- 1-rowMeans((1-kappa_ij.save)*(1-t_ij.save), dim=2, na.rm=TRUE)

## alpha coefficients
alpha.postm <- rowMeans(out$alpha.save)
alpha.postm %>% round(3)
out$alpha.save %>% apply(1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

#----------------------------------------------------------------------------------------------------------------------------------
# Post-processing 

# # posterior predictive Y
# predY.save <- lapply(1:J, function(j) {
#   sapply(1:save, function(s) {
#     lapply(1:m_j[j], function(i) {
#       (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]) +
#         rnorm(End[i,j]-Start[i,j]+1, 0, 1)*exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)*sqrt(out$sigj_sq.save[j,s])
#     }) %>% unlist()
#   })
# })
# predY.postL <- lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
# predY.postU <- lapply(predY.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))
# 
# # full model residual
# Yhat.save <- lapply(1:J, function(j) {
#   sapply(1:save, function(s) {
#     lapply(1:m_j[j], function(i) {
#       (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s])/
#         exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
#     }) %>% unlist()
#   })
# })
# Yhat.postm <- lapply(Yhat.save, rowMeans)
# 
# Ytil.save <- lapply(1:J, function(j) {
#   sapply(1:save, function(s) {
#     lapply(1:m_j[j], function(i) {
#       Y[[j]][Start[i,j]:End[i,j]]/
#         exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
#     }) %>% unlist()
#   })
# })
# res.save <- mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE)
# res.postm <- lapply(res.save, rowMeans)
# 
# # dose-response function estimate
# Fit.save <- lapply(1:J, function(j) {
#   sapply(1:save, function(s) {
#     lapply(1:m_j[j], function(i) {
#       X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s]
#     }) %>% unlist()
#   })
# })
# Fit.postm <- lapply(Fit.save, rowMeans)
# Fit.postL <- lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.025)))
# Fit.postU <- lapply(Fit.save, function(x) apply(x, 1, function(y) quantile(y, prob=0.975)))
# 
# prc_res <- list(predY.postL = predY.postL, predY.postU = predY.postU,
#                 Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU,
#                 Yhat.postm = Yhat.postm, res.postm = res.postm)
# 
# # save result
# saveRDS(prc_res, paste0(path, "data/obese_processed.RDS"))

# call result
prc_res <- readRDS(paste0(path, "data/obese_processed.RDS"))

# rmse
rmse <- sqrt(mean(unlist(prc_res$res.postm)^2))

#----------------------------------------------------------------------------------------------------------------------------------
# Section 5.2: Empirical coverage of 95% posterior predictive intervals for observations in heteroscedastic cells 
predY.postL <- prc_res$predY.postL
predY.postU <- prc_res$predY.postU

hetero_cells <- which(tr_t_ij.postm > 0.5, arr.ind = TRUE)
CIcov <- list()
for (k in 1:nrow(hetero_cells)) {
  i <- hetero_cells[k,1]
  j <- hetero_cells[k,2]
  inew <- which(idx_j[[j]] == i)
  resp = Y[[j]][Start[inew,j]:End[inew,j]]
  respLqt = predY.postL[[j]][Start[inew,j]:End[inew,j]]
  respUqt = predY.postU[[j]][Start[inew,j]:End[inew,j]]
  CIcov[[k]] <- sum(resp > respLqt & resp < respUqt)/length(resp)
}
CIcov %>% unlist() %>% mean
CIcov %>% unlist() %>% summary()

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (top row)
data_f <- list(misdata = misdata, meta = meta)
result_f <- list(out = list(gamma_ij.postm = gamma_ij.postm, 
                            kappa_ij.postm = kappa_ij.postm,
                            t_ij.postm = t_ij.postm),
                 hit_mat = hit_mat, prc_res = prc_res)

t1 <- dosres_plot(i = which(uniq_chnm == "Bisphenol A"),
                  j = which(uniq_aenm == "NVS_TR_hDAT"),
                  data_f, result_f, realdata = TRUE)
t2 <- dosres_plot(i = which(uniq_chnm == "p,p'-DDE"),
                  j = which(uniq_aenm == "TOX21_ERR_Antagonist"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t1, t2, nrow=1)

# g <- gridExtra::arrangeGrob(t1, t2, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/ToxCast_obese_hitcbad1.pdf"), g)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (middle row)
t3 <- dosres_plot(i = which(uniq_chnm == "Dibutyl phthalate"),
                  j = which(uniq_aenm == "ATG_PPARg_TRANS_up"),
                  data_f, result_f, realdata = TRUE)
t4 <- dosres_plot(i = which(uniq_chnm == "Bisphenol A"),
                  j = which(uniq_aenm == "ATG_PPARg_TRANS_up"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t3, t4, nrow=1)

# g <- gridExtra::arrangeGrob(t3, t4, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/ToxCast_obese_new_hitcbad3.pdf"), g)

# Why is hit-call 0 for (Dibutyl phthalate, ATG_PPARg_TRANS_up) and (Bisphenol A, ATG_PPARg_TRANS_up)? 
i = which(uniq_chnm == "Dibutyl phthalate")
# i = which(uniq_chnm == "Bisphenol A")
j = which(uniq_aenm == "ATG_PPARg_TRANS_up")

aeid_j <- meta$uniq_aeid[j]
aenm_j <- uniq_aenm[j]
bmad_j <- assay_cutoff[which(assay_cutoff$aeid == aeid_j),"bmad"]
coff_j <- assay_cutoff[which(assay_cutoff$aeid == aeid_j),"coff"]
chnm_i <- uniq_chnm[i]

datatmp <- obese_data %>% filter(chnm == chnm_i, aenm == aenm_j)
## negative response data need an inverse transformation
resp <- abs(datatmp$resp)
## fitting
params <- tcplFit(datatmp$logc, resp, bmad=bmad_j, 
                  force.fit = FALSE, bidirectional = FALSE, verbose = TRUE)
max(params$hill_modl) > coff_j

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4 (bottom row)
t5 <- dosres_plot(i = which(uniq_chnm == "Diisobutyl phthalate"),
                  j = which(uniq_aenm == "TOX21_CAR_Antagonist"),
                  data_f, result_f, realdata = TRUE)
t6 <- dosres_plot(i = which(uniq_chnm == "Di(2-ethylhexyl) phthalate"),
                  j = which(uniq_aenm == "TOX21_Aromatase_Inhibition"),
                  data_f, result_f, realdata = TRUE)

gridExtra::grid.arrange(t5, t6, nrow=1)

# g <- gridExtra::arrangeGrob(t5, t6, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/ToxCast_obese_hitcbad2.pdf"), g)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 6: Chemical ranks by the average active probability
tr_hit_mat <- hit_mat
tr_hit_mat[which(is.na(tr_gamma_ij.postm))] <- NA

colors <- c("BMC_cutoff" = "black", "EPA" = "blue")
shapes <- c("BMC_cutoff" = 19, "EPA" = 17)

data.frame(act_cutoff = rowMeans(tr_actprob_cutoff, na.rm=TRUE),
           averhit = rowMeans(tr_hit_mat, na.rm=TRUE), chnm=uniq_chnm) %>%
  ggplot() + geom_point(aes(reorder(chnm, act_cutoff), act_cutoff, 
                            shape="BMC_cutoff", color="BMC_cutoff"), size=5) +
  geom_point(aes(reorder(chnm,act_cutoff), averhit, shape="EPA", color="EPA"), size=3) +
  labs(x="Chemical", y="Active probability",
       title="Obesity", color=NULL, shape=NULL) +
  theme_minimal() + theme_bw() +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=15),
        legend.text = element_text(size=15),
        axis.title = element_text(size=15),
        plot.title = element_text(size=20))
# ggsave(paste0(path, "Figure/ToxCast_obese_disruptchem.pdf"))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S7: Detection technology is one of the latent factors
eta.postm <- out$eta.save %>% list_mean()

obese_assays <- readxl::read_xlsx(paste0(path, "data/Obese_Assays.xlsx")) %>%
  rename("stricter_obesity" = 'Matches to stricter obesity?')
obese_assays <- left_join(data.frame(assay_component_endpoint_name = uniq_aenm),
                          obese_assays)

scales_17 <- scales::hue_pal()(17)
detection_cols <- c("AlphaLISA immunoassay" = scales_17[1],               
                    "CALUX luciferase quantitation" = scales_17[2],                  
                    "Colorimetric" = scales_17[3],                                   
                    "ELISA" = scales_17[4],                                          
                    "Filter-based radiodetection" = scales_17[5],                    
                    "Fluorescence" = scales_17[6],                                   
                    "Fluorescence and electrophoretic mobility shift" = scales_17[7],
                    "Fluorescence Polarization" = scales_17[8],                      
                    "GAL4 b-lactamase reporter gene" = scales_17[9],                 
                    "Luciferase" = scales_17[10],                                     
                    "Luciferase-coupled ATP quantitation" = scales_17[11],            
                    "Lysate-based radiodetection" = scales_17[12],                    
                    "Microscopy" = scales_17[13],                                     
                    "NA" = "#8c8c8c",                                             
                    "Protein-fragment Complementation" = scales_17[14],               
                    "RT-CES" = scales_17[15],                                         
                    "RT-PCR and Capillary electrophoresis" = scales_17[16],           
                    "TR-FRET" = scales_17[17])

data.frame(obese_assays,
           V1 = eta.postm[,1], V2 = eta.postm[,2]) %>%
  filter(detection_technology %in% c("ELISA",
                                     "Luciferase",
                                     "Microscopy",
                                     "AlphaLISA immunoassay",
                                     "Protein-fragment Complementation",
                                     "RT-PCR and Capillary electrophoresis",
                                     "RT-CES", 
                                     "Fluorescence Polarization")) %>%
  ggplot() + geom_vline(xintercept = 0, color="red", linetype = "dashed") + 
  geom_point(aes(V1, V2, color = detection_technology)) +
  labs(x = "Factor1", y = "Factor2", color = "Detection technology") + 
  scale_color_manual(values = detection_cols) +
  xlim(-2,1) + ylim(-0.085, 0.075)
# ggsave(paste0(path, "Figure/ToxCast_obese_latfac_left.eps"))

data.frame(obese_assays,
           V1 = eta.postm[,1], V2 = eta.postm[,2]) %>%
  filter(detection_technology %in% c("CALUX luciferase quantitation", 
                                     "Colorimetric", 
                                     "Fluorescence and electrophoretic mobility shift",
                                     "Fluorescence",
                                     "Filter-based radiodetection",
                                     "TR-FRET",
                                     "GAL4 b-lactamase reporter gene",
                                     "Luciferase-coupled ATP quantitation",
                                     "Lysate-based radiodetection"
  )) %>%
  ggplot() + geom_vline(xintercept = 0, color="red", linetype = "dashed") + 
  geom_point(aes(V1, V2, color = detection_technology)) +
  labs(x = "Factor1", y = "Factor2", color = "Detection technology") +
  scale_color_manual(values = detection_cols) + 
  xlim(-2,1) + ylim(-0.085, 0.075) 
# ggsave(paste0(path, "Figure/ToxCast_obese_latfac_right.eps"))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S9: Ranks of assay endpoints likely to be activated by the top 5 chemicals 
obese_ascend_chnm <- uniq_chnm[order(rowMeans(tr_actprob_cutoff, na.rm=TRUE))]
obese_active_aenm <- data.frame(act = colMeans(tr_actprob_cutoff[which(uniq_chnm %in% obese_ascend_chnm[c(26:30)]),], na.rm=TRUE), 
                                aenm = uniq_aenm) %>% 
  filter(0.9 < act) %>% select(aenm)
# saveRDS(obese_active_aenm, paste0(path, "data/obese_active_aenm.RDS"))
neuro_active_aenm <- readRDS(paste0(path, "data/neuro_active_aenm.RDS"))

data.frame(act_cutoff = colMeans(tr_actprob_cutoff[which(uniq_chnm %in% obese_ascend_chnm[c(26:30)]),], na.rm=TRUE), 
           aenm = uniq_aenm, 
           common = uniq_aenm %in% neuro_active_aenm$aenm) %>% 
  filter(0.9 < act_cutoff) %>% 
  ggplot() + geom_point(aes(reorder(aenm, act_cutoff), act_cutoff, 
                            color=factor(common), shape=factor(common)), size=5) + 
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
# Figure S11: Results for hold-out pairs 
data_f12 <- list(simdata = simdata, meta = meta)
result_f12 <- list(gamma_ij.postm = gamma_ij.postm, 
                   t_ij.postm = t_ij.postm)

t9 <- dosres_plot(i = which(uniq_chnm == "Cyfluthrin"),
                   j = which(uniq_aenm == "NVS_GPCR_h5HT7"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE) +
  ylim(-1.9,2.9) + geom_hline(yintercept = 0, col="red", linetype="dashed")
t10 <- dosres_plot(i = which(uniq_chnm == "2-Hydroxy-4-methoxybenzophenone"),
                  j = which(uniq_aenm == "ACEA_AR_antagonist_80hr"),
                  data_f12, result_f12, realdata = TRUE, pred = TRUE) + 
  ylim(-1.9,2.9) + geom_hline(yintercept = 0, col="red", linetype="dashed")
t11 <- dosres_plot(i = which(uniq_chnm == "Dichlorodiphenyltrichloroethane"),
                   j = which(uniq_aenm == "TOX21_ARE_BLA_agonist_ratio"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE) +
  geom_hline(yintercept = 0, col="red", linetype="dashed")
t12 <- dosres_plot(i = which(uniq_chnm == "Benzyl butyl phthalate"),
                   j = which(uniq_aenm == "TOX21_ERa_BLA_Agonist_ratio"),
                   data_f12, result_f12, realdata = TRUE, pred = TRUE) +
  geom_hline(yintercept = 0, col="red", linetype="dashed")

gridExtra::grid.arrange(t9, t10, t11, t12, nrow=2)

# g <- gridExtra::arrangeGrob(t9, t10, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/ToxCast_obese_pred2.pdf"), g)
# g <- gridExtra::arrangeGrob(t11, t12, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/ToxCast_obese_pred3.pdf"), g)
