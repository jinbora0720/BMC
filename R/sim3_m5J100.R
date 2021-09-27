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
J <- 100
q <- 5

d <- 5
alpha <- c(0.3,1)
xi <- 0.8
delta_mean <- 1

seed = 1
gendata = generate_data(m, J=5, d, q, xi, alpha, delta_mean, seed)
truth = gendata$truth
gamma_ij = truth$gamma_ij
t_ij = truth$t_ij
z_ij = truth$z_ij
u_ij = truth$u_ij
idx_g1_orig <- which(gamma_ij == 1)
idx_t1_orig <- which(t_ij == 1)

# append 0 in gamma and t
set.seed(seed)

m_j = rep(m, J) 
idx_j = list()
for (j in 1:J) {
  idx_j[[j]] = sample(m, m_j[j], replace=FALSE) %>% sort()
}
K_ij = matrix(d*8, m, J)   
for (j in 1:J) {
  K_ij[-idx_j[[j]], j] = 0
}
Start = End = matrix(0, m, J)
for (j in 1:J) {
  Start[1:m_j[j],j] = cumsum(c(1,K_ij[,j]))[idx_j[[j]]]
  End[1:m_j[j],j] = cumsum(K_ij[,j])[idx_j[[j]]]
}

gamma_ij = cbind(gamma_ij, matrix(0, nrow=m, ncol=J-5))
for (j in 1:J) {
  gamma_ij[-idx_j[[j]], j] = NA
}

nu0 = 5; sig0_sq = 0.1
sigj_sq = 1/rgamma(J, nu0/2, nu0*sig0_sq/2)   
resid = lapply(1:J, function(j) {
  rnorm(max(End[,j]), 0, sqrt(sigj_sq[j]))
})

p = 9
beta_pool = list()
beta_pool[[1]] = c(0.5, 0.5, 1, 1.3, 1.5, 2, 2.5, 3, 3) # hill
beta_pool[[2]] = c(0, 0, 0, 0.3, 0.8, 1.3, 2.5, 1.5, 0.8) # gain-loss
beta_pool[[3]] = c(0, 0, 0, 0, 0, 0, -0.5, -0.7, -2) # decreasing

beta_ij = list()
if (J > 2) {
  beta_ij[1:floor(J/3)] = lapply(1:floor(J/3), function(j) {
    matrix(beta_pool[[1]], p, m_j[j])
  }) 
  beta_ij[(1+floor(J/3)):(2*floor(J/3))] = lapply((1+floor(J/3)):(2*floor(J/3)), function(j) {
    matrix(beta_pool[[2]], p, m_j[j])
  }) 
  beta_ij[(1+2*floor(J/3)):J] = lapply((1+2*floor(J/3)):J, function(j) {
    matrix(beta_pool[[3]], p, m_j[j])
  }) 
} else {
  beta_ij[1:J] = lapply(1:J, function(j) {
    matrix(beta_pool[[1]], p, m_j[j])
  }) 
}

for (j in 1:J) {
  for (i in 1:m_j[j]) {
    if (gamma_ij[idx_j[[j]][i],j]==0) beta_ij[[j]][,i] = rep(0,p)
  }
}

# heteroscedasticity 
t_ij = cbind(t_ij, matrix(0, nrow=m, ncol=J-5))
for (j in 1:J) {
  t_ij[-idx_j[[j]],j] = NA
}
d_ij = matrix(NA, m, J)
d_ij[which(t_ij==1, arr.ind=TRUE)] = rnorm(sum(t_ij==1, na.rm=TRUE), delta_mean, 0.1)
d_ij[which(t_ij==0, arr.ind=TRUE)] = 0

# Data
logc = rep(c(0.301, 0.477, 0.602, 0.845, 1.000, 1.301, 1.602, 2.000), d)
Xtemp = apply(bs(logc, df=p, intercept=TRUE), 2, function(x) scale(x, center=TRUE, scale=FALSE)) 

orgData = rxf = list()
for (j in 1:J) {
  orgData[[j]] = matrix(NA, max(End[,j]), 2)
  rxf[[j]] = rep(NA, max(End[,j]))
  for (i in 1:m_j[j]) {
    rxf[[j]][Start[i,j]:End[i,j]] = Xtemp%*%matrix(beta_ij[[j]][,i], ncol=1)
    resp = rxf[[j]][Start[i,j]:End[i,j]]+exp(logc*d_ij[idx_j[[j]][i], j]/2)*resid[[j]][Start[i,j]:End[i,j]]
    orgData[[j]][Start[i,j]:End[i,j],] = cbind(logc, resp)
  }
}
orgX = lapply(orgData, function(x) x[,1]) 
orgY = lapply(orgData, function(x) x[,2])

p1 = 7
CX = apply(bs(logc, df=p1, intercept=TRUE), 2, function(x) scale(x, center=TRUE, scale=FALSE))
X = replicate(J, 
              do.call("rbind", replicate(m, CX, simplify = FALSE)), 
              simplify = FALSE)

simdata = list(orgX=orgX, orgY=orgY, X=X, K_ij=K_ij, 
               m_j=m_j, idx_j=idx_j, Start=Start, End=End)
truth = list(gamma_ij=gamma_ij, rxf=rxf, 
             t_ij=t_ij, d_ij=d_ij, sigj_sq=sigj_sq)

gamma_ij = truth$gamma_ij
idx_g1 <- which(gamma_ij == 1)
sum(gamma_ij)
t_ij = truth$t_ij
idx_t1 <- which(t_ij == 1)
sum(t_ij)

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
ppg = data.frame(gamma_signal = z_ij[idx_g1_orig], 
                 gamma_postprob = gamma_ij.postm[idx_g1])
ppt = data.frame(t_signal = u_ij[idx_t1_orig], 
                 t_postprob = t_ij.postm[idx_t1])

### number of false positive
cutoff = seq(0, 1, by = 0.1)
fpg = fp_cutoff(gamma_ij.postm, gamma_ij, cutoff) 
fpt = fp_cutoff(t_ij.postm, t_ij, cutoff) 

saveRDS(list(out = out, 
             ppg = ppg, ppt = ppt, 
             fpg = fpg, fpt = fpt), 
        "~/BMC/data/sim3_m5J100.RDS")
