# Simulation 2

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 

# source code
source("bmc.R")
source("bmc_sampler.R")
sourceCpp("bmc_sampler2.cpp")
sourceCpp("samplerBits.cpp")
source("simdata.R")

########
# Data #
########
m <- 15
J <- 15
rmse = tr_aucg = tt_aucg <- rep(NA, 50)
set.seed(123)
seedsave <- sample(10000, 100, replace=FALSE)
missing_idx <- rbind(c(13,13), c(13,14), c(13,15),
                    c(14,13), c(14,14), c(14,15),
                    c(15,13), c(15,14), c(15,15))
test_idx <- apply(missing_idx, 1, function(x) x[2]+J*(x[1]-1))
iterpool <- c(1,2,3,4,5,6,9,10,
             11,13,15,16,17,18,19,20, 
             21,23,24,25,26,27,28,29,30,
             31,32,33,34,35,36,37,38, 
             41,43,44,45,46,47,48,49,50,
             51,52,53,54,58,59,60,
             61)

seed = seedsave[iterpool[iter]]
simout = simdata(nchem=m, nassay=J, seed)
dat = simout$dat
truth = simout$truth

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
    Xtemp = apply(bs(logc, df=p, intercept=TRUE), 2, function(x) scale(x, center=TRUE, scale=FALSE)) # centred version
    X[[j]][Start[i,j]:End[i,j],] = Xtemp
  }
}

simuldata = list(X=X, orgX=orgX, orgY=Y, K_ij=K_ij, m_j=m_j, idx_j=idx_j, Start=Start, End=End)

misdata = data_missing(simdata=simuldata, missing_idx=missing_idx, seed=seed)

X = misdata$X
orgX = misdata$orgX
Y = misdata$orgY
m_j = misdata$m_j
idx_j = misdata$idx_j
Start = misdata$Start
End = misdata$End

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