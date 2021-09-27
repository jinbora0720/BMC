# Reproduce results related to neurodevelopmental disorders 
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
theme_set(theme_bw())
library(splines) 

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
neuro_data <- datalist$neuro_data
hit_vec <- datalist$hit_vec

# standardise cutoff by each pair 
assay_cutoff <- readRDS(paste0(path, "data/assay_cutoff.RDS"))
# meta <- metadata(neuro_data, aeid_coff = assay_cutoff)
meta <- readRDS(paste0(path, "data/neuro_meta.RDS"))
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
pred_idx <- readRDS(paste0(path, "data/neuro_pred_idx.rds"))
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

mu_xi <- qnorm(mean(hit_mat, na.rm = TRUE))  
sigsq_xi <- 1

mu_alpha <- c(0,1) 
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
# save result
# saveRDS(list(out = out, missing_idx = missing_idx), paste0(path, "data/neuro_out.RDS"))

# call result
neuro_out <- readRDS(paste0(path, "data/neuro_out.RDS"))
out <- neuro_out$out

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
# Figure 5: Chemical ranks by the average active probability
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
       title="Neurodevelopmental Disorders", color=NULL, shape=NULL) +
  theme_minimal() + theme_bw() +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=15),
        legend.text = element_text(size=15),
        axis.title = element_text(size=15),
        plot.title = element_text(size=20))
# ggsave(paste0(path, "Figure/ToxCast_neuro_disruptchem.pdf"))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S8: Ranks of assay endpoints likely to be activated by the top 5 chemicals
# ("Bisphenol A", "p,p'-DDE", "2,4,5-Trichlorophenol", "Dichlorodiphenyltrichloroethane", "Triclosan")
neuro_ascend_chnm <- uniq_chnm[order(rowMeans(tr_actprob_cutoff, na.rm=TRUE))]
neuro_active_aenm <- data.frame(act_cutoff = colMeans(tr_actprob_cutoff[which(uniq_chnm %in% neuro_ascend_chnm[c(26:30)]),], na.rm=TRUE),
                                aenm = uniq_aenm) %>%
  filter(0.9 < act_cutoff) %>% select(aenm)
# saveRDS(neuro_active_aenm, paste0(path, "data/neuro_active_aenm.RDS"))
obese_active_aenm <- readRDS(paste0(path, "data/obese_active_aenm.RDS"))

data.frame(act_cutoff = colMeans(tr_actprob_cutoff[which(uniq_chnm %in% neuro_ascend_chnm[c(26:30)]),], na.rm=TRUE),
           aenm=uniq_aenm, common = uniq_aenm %in% obese_active_aenm$aenm) %>%
  filter(0.9 < act_cutoff) %>%
  ggplot() + geom_point(aes(reorder(aenm, act_cutoff), act_cutoff, 
                            color=factor(common), shape=factor(common)), size=5) + 
  scale_color_manual(values=c("TRUE"="red","FALSE"="blue")) +
  labs(title="Neurodevelopmental Disorders", x="Assay endpoint", 
       y="Activated probability", color="Involved in both diseases", 
       shape="Involved in both diseases") +
  theme_minimal() + theme_bw() +
  theme(axis.text.x = element_text(angle=75, hjust=1, size=15),
        legend.position = "bottom",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title = element_text(size=15),
        plot.title = element_text(size=20))
# ggsave(paste0(path, "Figure/ToxCast_neuro_disruptassay.pdf"))
