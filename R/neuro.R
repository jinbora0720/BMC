# Reproduce results related to neurodevelopmental disorders 

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 

# source code
path <- "~/Documents/GitHub/BMC/"
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
