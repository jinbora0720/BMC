# Simulation 3: multiplicity adjustment 
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

fp_cutoff <- function(pred, truth, cutoff) {
  k <- length(cutoff) 
  fp <- rep(0, k)
  for (i in 1:k) {
    pred_01 = ifelse(pred >= cutoff[i], 1, 0)
    fp[i] = sum(pred_01 == 1 & truth == 0)
  }
  return(fp)
}

########
# Data #
########
m <- 5
J <- 5
q <- 5
alpha <- c(0.3,1)
xi <- 0.8

seed = 1
gendata = generate_data(m, J, d=5, q, xi=xi, alpha=alpha, delta_mean=1, seed=seed)
simdata = gendata$simdata
truth = gendata$truth

Lambda = truth$Lambda
eta = truth$eta
z_ij = truth$z_ij
u_ij = truth$u_ij
gamma_ij = truth$gamma_ij
idx_g1 <- which(gamma_ij == 1)
t_ij = truth$t_ij
idx_t1 <- which(t_ij == 1)

X = simdata$X
orgX = simdata$orgX
Y = simdata$Y = simdata$orgY
m_j = simdata$m_j
idx_j = simdata$idx_j
Start = simdata$Start
End = simdata$End

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
mu_xi = qnorm(sum(gamma_ij)/(m*J)) # xi 
sigsq_xi = 1

mu_alpha = c(qnorm(sum(t_ij)/(m*J)), sum(gamma_ij == t_ij)/(m*J)) # alpha
sigsq_alpha = c(1,1)

hyper = list(a = a, R = R, v_d = v_d,
             mu_xi = mu_xi, sigsq_xi = sigsq_xi,
             mu_alpha = mu_alpha, sigsq_alpha = sigsq_alpha)

##################
# initial values #
##################
ols = lm(scale(unlist(Y))~bs(unlist(orgX), df=p)-1)
beta_ij = list()
for (j in 1:J) {
  beta_ij[[j]] = matrix(coef(ols), nrow=p, ncol=m_j[j])
}
init = list(q = q, beta_ij = beta_ij)

###################
# MCMC parameters #
###################
thin <- 5
burnin <- 10000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

#######
# BMC #
#######
out <- bmc(Data = simdata, MCMC = MCMC, hyper = hyper, init = init,
           hetero = TRUE, hetero_simpler = FALSE,
           gamma_simpler = FALSE, mgsp_adapt = FALSE,
           update_xi = TRUE, apply_cutoff = FALSE, verbose = TRUE)

###########
# results #
###########
### gamma
gamma_ij.postm = rowMeans(out$gamma_ij.save, dim=2, na.rm=TRUE)

### t
t_ij.postm = rowMeans(out$t_ij.save, dim=2, na.rm=TRUE)

### posterior probability 
ppg = data.frame(gamma_signal = z_ij[idx_g1], gamma_postprob = gamma_ij.postm[idx_g1])
ppt = data.frame(t_signal = u_ij[idx_t1], t_postprob = t_ij.postm[idx_t1])

### number of false positive
cutoff = seq(0, 1, by = 0.1)
fpg = fp_cutoff(gamma_ij.postm, gamma_ij, cutoff) 
fpt = fp_cutoff(t_ij.postm, t_ij, cutoff) 

saveRDS(list(out = out, 
             ppg = ppg, ppt = ppt, 
             fpg = fpg, fpt = fpt), 
        "~/BMC/data/sim3_m5J5.RDS")

