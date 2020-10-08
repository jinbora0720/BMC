# Reproduce results related to neurodevelopmental disorders 

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 

# source code
# path <- "~/Documents/GitHub/BMC/"
path <- "/home/grad/bj91/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler2.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/metadata.R"))
source(paste0(path, "source/simulate_data.R"))

########
# Data #
########
datalist <- readRDS(paste0(path, "data/datalist.rds"))
neuro_data <- datalist$neuro_data
hit_vec <- datalist$hit_vec
meta <- metadata(neuro_data)
uniq_chnm <- meta$uniq_chnm
m <- length(uniq_chnm)
uniq_aenm <- meta$uniq_aenm
J <- length(uniq_aenm)
simdata <- list(X = meta$X, orgX = lapply(meta$orgData, function(x) x[,1]), 
                orgY = meta$Y, K_ij = meta$K_ij, m_j = meta$m_j, 
                idx_j = meta$idx_j, Start = meta$Start, End = meta$End)

# pool for hold-out cells 
# one can replace this by any cells with measurements 
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
# saveRDS(list(out = out, missing_idx = missing_idx), file.path(paste0(path, "data/neuro_out.rds")))

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
gamma_ij.postm <- rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)

uniq_code <- paste("C",gsub("-", "", meta$uniq_casn), sep="")
hit_mat <- hit_vec %>% filter(aenm %in% uniq_aenm, code %in% uniq_code) %>% 
  reshape2::dcast(code~aenm, value.var = "hitc") 
rownames(hit_mat) <- hit_mat[,1]
hit_mat <- as.matrix(hit_mat[uniq_code,uniq_aenm])
hit_mat[which(is.na(gamma_ij.postm))] <- NA

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S4: Chemical ranks by the average active probability 
# Prob.Active
actprob <- 1-rowMeans((1-out$gamma_ij.save)*(1-out$t_ij.save), dim=2, na.rm=TRUE)

colors <- c("BMC" = "black", "EPA" = "blue") 
shapes <- c("BMC" = 19, "EPA" = 17)

data.frame(act = rowMeans(actprob, na.rm=TRUE), 
           averhit = rowMeans(hit_mat, na.rm=TRUE), chnm=uniq_chnm) %>% 
  ggplot() + geom_point(aes(reorder(chnm, act), act, shape="BMC", color="BMC"), size=5) + 
  geom_point(aes(reorder(chnm,act), averhit, shape="EPA", color="EPA"), size=3) +
  labs(x="Chemical", y="Active probability", title="Neurodevelopmental Disorders", color=NULL, shape=NULL) +
  theme_minimal() + theme_bw() + 
  scale_color_manual(values = colors) + 
  scale_shape_manual(values = shapes) + 
  theme(axis.text.x = element_text(angle=65, hjust=1, size=15), 
        legend.text = element_text(size=15), 
        axis.title = element_text(size=15), 
        plot.title = element_text(size=20)) 

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 14: Ranks of assay endpoints likely to be activated by the top 5 chemicals 
# ("Bisphenol A", "p,p'-DDE", "2,4,5-Trichlorophenol", "Dichlorodiphenyltrichloroethane", "Triclosan")
neuro_ascend_chnm <- uniq_chnm[order(rowMeans(actprob, na.rm=TRUE))]
neuro_active_aenm <- data.frame(act = colMeans(actprob[which(uniq_chnm %in% neuro_ascend_chnm[c(24,27:30)]),], na.rm=TRUE), 
                                averhit = colMeans(hit_mat, na.rm=TRUE), aenm=uniq_aenm) %>% 
  filter(0.9 < act) %>% select(aenm)
# saveRDS(neuro_active_aenm, paste0(path, "data/neuro_active_aenm.rds"))
obese_active_aenm <- readRDS(paste0(path, "data/obese_active_aenm.rds"))

data.frame(act = colMeans(actprob[which(uniq_chnm %in% neuro_ascend_chnm[c(24,27:30)]),], na.rm=TRUE), 
           averhit = colMeans(hit_mat, na.rm=TRUE), aenm=uniq_aenm, 
           common = uniq_aenm %in% obese_active_aenm$aenm) %>% 
  filter(0.9 < act) %>% 
  ggplot() + geom_point(aes(reorder(aenm, act), act, color=factor(common), shape=factor(common)), size=5) + 
  scale_color_manual(values=c("TRUE"="red","FALSE"="blue")) + 
  labs(title="Neurodevelopmental Disorders", x="Assay endpoint", y="Activated probability", 
       color="Involved in both diseases", 
       shape="Involved in both diseases") +
  theme_minimal() + theme_bw() + 
  theme(axis.text.x = element_text(angle=75, hjust=1, size=15), 
        legend.position = "bottom", 
        legend.title = element_text(size=15),
        legend.text = element_text(size=15), 
        axis.title = element_text(size=15), 
        plot.title = element_text(size=20)) 
 