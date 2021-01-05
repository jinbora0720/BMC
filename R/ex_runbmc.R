# Example code to show how to run bmc 

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

########
# Data #
########
datalist <- readRDS(paste0(path, "data/datalist.rds"))
neuro_data <- datalist$neuro_data
hit_vec <- datalist$hit_vec
meta <- readRDS(paste0(path, "data/neuro_meta.rds"))
meta$orgX <- lapply(meta$orgData, function(x) x[,1])

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
dat <- data.frame(resp = unlist(meta$Y), logc = unlist(lapply(meta$orgData, function(x) x[,1])))  
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
out <- bmc(Data = meta, MCMC = MCMC, hyper = hyper, init = init, 
           hetero = TRUE, adapt = TRUE, simpler = FALSE, verbose = TRUE)
