# Simulation 1
rm(list = ls())

# dependencies 
library(MASS)
library(Rcpp) 
library(tidyverse) 
library(splines) 
library(ROCR)
library(tcpl)

# source code
path <- "~/BMC/"
source(paste0(path, "source/bmc.R"))
source(paste0(path, "source/bmc_sampler.R"))
sourceCpp(paste0(path, "source/bmc_sampler.cpp"))
sourceCpp(paste0(path, "source/samplerBits.cpp"))
source(paste0(path, "source/simulate_data.R"))
source(paste0(path, "source/dosres_plot.R"))

########
# Data #
########
m <- 30
J <- 150
q <- 2
alpha <- c(-0.1,1.2)
xi <- 0

seed = 123
gendata = generate_data(m=m, J=J, d=3, q=q, 
                        xi=xi, alpha=alpha, delta_mean=1.5, seed=seed)
simdata = gendata$simdata
truth = gendata$truth
simdata$rxf = truth$rxf

prob_missing = 0.05
misdata = data_missing(simdata, prob_missing, missing_idx=NULL, seed=seed)
misdata$Y <- misdata$orgY
missing_idx = misdata$missing_idx 
missing_idx_col = misdata$missing_idx_col 
test_idx = apply(missing_idx,1,function(x) x[2]+J*(x[1]-1))

tt_gamma = as.numeric(truth$gamma_ij[missing_idx])
tr_gamma = as.numeric(truth$gamma_ij)[-missing_idx_col]
tt_t = as.numeric(truth$t_ij[missing_idx])
tr_t = as.numeric(truth$t_ij)[-missing_idx_col]
actprob = 1-(1-truth$gamma_ij)*(1-truth$t_ij)
tt_actprob = as.numeric(actprob[missing_idx])
tr_actprob = as.numeric(actprob)[-missing_idx_col]

X = misdata$X
orgX = misdata$orgX
Y = misdata$Y
m_j = misdata$m_j
idx_j = misdata$idx_j
Start = misdata$Start
End = misdata$End

###################
# MCMC parameters #
###################
thin <- 10
burnin <- 10000
save <- 1000
MCMC <- list(thin = thin, burnin = burnin, save = save)

################
# Call results #
################
# BMC
resnew <- readRDS(paste0(path, "data/sim1_BMC_all.RDS")) 
# result split due to file size
bmc123 <- readRDS(paste0(path, "data/sim1_BMC_seed123_out1.RDS"))
out <- bmc123$out
out$beta_ij.save <- c(readRDS(paste0(path, "data/sim1_BMC_seed123_out2.RDS")),
                      readRDS(paste0(path, "data/sim1_BMC_seed123_out3.RDS")))

Yhat.save = lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      (X[[j]][Start[i,j]:End[i,j],]%*%out$beta_ij.save[[j]][,i,s])/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
Yhat.postm = lapply(Yhat.save, rowMeans)

Ytil.save = lapply(1:J, function(j) {
  sapply(1:save, function(s) {
    lapply(1:m_j[j], function(i) {
      Y[[j]][Start[i,j]:End[i,j]]/
        exp(orgX[[j]][Start[i,j]:End[i,j]]*out$d_ij.save[idx_j[[j]][i],j,s]/2)
    }) %>% unlist() 
  }) 
})
res.postm = lapply(mapply('-', Ytil.save, Yhat.save, SIMPLIFY = FALSE), rowMeans)

## 1. overall RMSE
rmse = sqrt(mean(unlist(res.postm)^2))

## 2. gamma 
gamma_ij.save = out$gamma_ij.save 
le.save = mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE)
z.save = list()
for (i in 1:MCMC$save) {
  z.save[[i]] = matrix(rnorm(m*J, le.save[[i]] + out$xi.save[i], 1), m, J)
}
gamma_ij.save[is.na(gamma_ij.save)] = sapply(z.save, function(x) x[missing_idx]>0)
gamma_ij.postm = rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)

## 3. t_ij
t_ij.save = out$t_ij.save 
lea = lapply(1:save, function(s) out$alpha.save[1,s] + out$alpha.save[2,s]*le.save[[s]])
u.save = lapply(lea, function(x) matrix(rnorm(m*J, x, 1), m, J))
t_ij.save[is.na(t_ij.save)] = sapply(u.save, function(x) x[missing_idx]>0)
t_ij.postm = rowMeans(t_ij.save, dim=2, na.rm=TRUE)

## 4. gamma_ij = 1 or t_ij = 1
actprob.postm = 1-rowMeans((1-gamma_ij.save)*(1-t_ij.save), dim=2)

## 5. alphas 
alpha.postm <- rowMeans(out$alpha.save)

## 6. xi
xi.postm <- mean(out$xi.save)

# ### posterior predictive Y
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
# ### dose-response function estimate
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
# naive_res <- mapply('-', Y, Fit.postm, SIMPLIFY = FALSE)
# prc_res <- list(predY.postL = predY.postL, predY.postU = predY.postU,
#                 Fit.postm = Fit.postm, Fit.postL = Fit.postL, Fit.postU = Fit.postU,
#                 Yhat.postm = Yhat.postm, res.postm = res.postm, naive_res = naive_res)
# 
# ### save result
# saveRDS(prc_res, paste0(path, "data/sim1_BMC_processed.RDS"))

### call result
prcnew <- readRDS(paste0(path, "data/sim1_BMC_processed.RDS"))

# BMC_i
resi <- readRDS(paste0(path, "data/sim1_BMCi_all.RDS")) 

# tcpl
res_tcpl <- readRDS(paste0(path, "data/sim1_tcpl.RDS")) 

# ZIPLL
res_ZIPLL <- readRDS(paste0(path, "data/sim1_ZIPLL.RDS")) 
ZIPLL <- readRDS(paste0(path, "data/sim1_ZIPLL_seed123.RDS"))$fit_ZIPLL

#----------------------------------------------------------------------------------------------------------------------------------
# Table 1: summary of results (make mean and sd tables separately)
mres <- matrix(c(round(unlist(lapply(resnew[1:7], mean)), 3),
                 round(unlist(lapply(resi[1:7], mean)), 3),
                 c(round(unlist(lapply(res_ZIPLL[1:2], mean)), 3), rep(NA,5)),
                 c(round(unlist(lapply(res_tcpl[1:2], mean)), 3), rep(NA,5))), 7, 4)
rownames(mres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 
                    'In-sample AUC for tij', 'Out-of-sample AUC for tij', 
                    'In-sample AUC for rij = 1 U tij = 1', 
                    'Out-of-sample AUC for rij = 1 U tij = 1')
colnames(mres) <- c('BMC', 'BMCi', 'ZIPLL', 'tcpl')

sdres <- matrix(c(round(unlist(lapply(resnew[1:7], sd)), 3),
                  round(unlist(lapply(resi[1:7], sd)), 3),
                  c(round(unlist(lapply(res_ZIPLL[1:2], sd)), 3), rep(NA,5)),
                  c(round(unlist(lapply(res_tcpl[1:2], sd)), 3), rep(NA,5))), 7, 4)
rownames(sdres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 
                    'In-sample AUC for tij', 'Out-of-sample AUC for tij', 
                    'In-sample AUC for rij = 1 U tij = 1', 
                    'Out-of-sample AUC for rij = 1 U tij = 1')
colnames(sdres) <- c('BMC', 'BMCi', 'ZIPLL', 'tcpl')

cat('mean table\n'); mres
cat('sd table\n'); sdres

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 2: dose-response curves
data_f7 <- list(misdata = misdata, truth = truth)
res_f7 <- list(res_ZIPLL = ZIPLL,
               out = list(gamma_ij.postm = gamma_ij.postm, t_ij.postm = t_ij.postm),
               prc_res = prcnew)

t7 <- dosres_plot(i = 23, j = 69, data = data_f7, result = res_f7)
t8 <- dosres_plot(i = 18, j = 10, data = data_f7, result = res_f7)
t9 <- dosres_plot(i = 27, j = 105, data = data_f7, result = res_f7)

ggpubr::ggarrange(t7, t8, t9, nrow=1, ncol=3, common.legend=TRUE, legend="bottom", labels="AUTO")
# ggsave(paste0(path, "Figure/Simulation1_typicalcurves_lab.pdf"))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 3: heat map of estimate and true profiles of activity
Ind <- 1:5; Ind2 <- 38:42

tr_gamma_ij.postm <- rowMeans(out$gamma_ij.save, dim=2, na.rm=TRUE)
pred <- which(is.na(tr_gamma_ij.postm[Ind,Ind2]), arr.ind = TRUE)
frames <- data.frame(Var1=pred[,1], Var2=pred[,2])
t5g <- reshape2::melt(actprob.postm[Ind,Ind2]) %>%
  rename(Chem = Var1, Assay = Var2) %>%
  ggplot() + geom_tile(aes(Assay, Chem, fill=value), color="grey50") +
  geom_text(aes(Assay, Chem, label=round(value,3))) +
  scale_fill_gradient2(high = "red3", low = "white", mid = "#E6F598", midpoint = 0.5, limits=c(0,1), na.value = "black") +
  scale_y_reverse() +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  geom_rect(data=frames, size=1, fill=NA, colour="black",
            aes(xmin=Var2 - 0.5, xmax=Var2 + 0.5, ymin=Var1 - 0.5, ymax=Var1 + 0.5)) +
  labs(y="Chemical", x="Assay endpoint", title="Estimate", fill="Prob.Active")

t6g <- reshape2::melt(actprob[Ind,Ind2]) %>%
  rename(Chem = Var1, Assay = Var2) %>%
  ggplot() + geom_tile(aes(Assay, Chem, fill=factor(value)), color="grey50") +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  scale_y_reverse() +
  scale_fill_manual(values=c("1"="red3", "0"="white")) +
  labs(y="Chemical", x="Assay endpoint", title="Truth", 
       fill=expression(paste("1(",gamma[ij],"=1 U ",t[ij],"=1)")))

gridExtra::grid.arrange(t5g, t6g, nrow=1, widths=c(1.2/3,1.8/3))

# g <- gridExtra::arrangeGrob(t5g, t6g, nrow=1, widths=c(1.2/3,1.8/3))
# ggsave(paste0(path, "Figure/Simulation1_actprofile.eps"), g)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S4: heat map of the estimated and true correlation matrix
z_cor_T <- cov2cor(tcrossprod(truth$Lambda) + diag(1, m))
z_cor_data <- reshape2::melt(z_cor_T)
cor.postm <- cov2cor(out$covMean)
cor.postm_data <- data.frame(reshape2::melt(cor.postm), truevalue = z_cor_data$value)

rng <- range(cor.postm, z_cor_T)
t1 <- ggplot(cor.postm_data)+
  geom_tile(aes(reorder(Var2, truevalue), reorder(Var1, truevalue), fill=value)) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white",
                       limits=c(floor(rng[1]), ceiling(rng[2]))) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  labs(x="Chemical", y="Chemical", title="Estimate")

t2 <- ggplot(z_cor_T %>% reshape2::melt())+
  geom_tile(aes(reorder(Var2, value), reorder(Var1, value), fill=value)) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white",
                       space = "Lab", na.value = "grey50", guide = "colourbar",
                       aesthetics = "fill", name="Correlation",
                       limits=c(floor(rng[1]), ceiling(rng[2]))) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  labs(y="Chemical", x="Chemical", title="Truth")

gridExtra::grid.arrange(t1, t2, nrow=1, widths=c(1.265/3,1.735/3))

# g <- gridExtra::arrangeGrob(t1, t2, nrow=1, ncol=2, widths=c(1.265/3,1.735/3))
# ggsave(paste0(path, "Figure/Simulation1_correlation.pdf"), g)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S5: the estimated and true entries of loading matrix
lambda.postm <- Reduce("+", out$Lambda.save)/length(out$Lambda.save)
lambda.postm <- (-1)*lambda.postm # for better visualization

t3 <- ggplot(lambda.postm %>% reshape2::melt()) +
  geom_tile(aes(Var2, Var1, fill=value), color=1) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white", midpoint = 0,
                       name="Loading") +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20),
        axis.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  scale_y_reverse(breaks=c(1,10,20,30)) +
  labs(x="Factor", y="Chemical", title=expression(paste("Estimated ", Lambda)))

t4 <- ggplot(cbind(truth$Lambda,matrix(0,30,3)) %>% reshape2::melt()) + # for better visuatlization
  geom_tile(aes(Var2, Var1, fill=value), color=1) +
  scale_fill_gradient2(low = "#56B1F7", high = "red", mid = "white", midpoint = 0,
                       name="Loading") +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20),
        axis.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  scale_y_reverse(breaks=c(1,10,20,30)) +
  scale_x_continuous(breaks=c(1,2)) +
  labs(x="Factor", y="Chemical", title=expression(paste("True ", Lambda)))

gridExtra::grid.arrange(t3, t4, nrow=1)

# g <- gridExtra::arrangeGrob(t3, t4, nrow=1, ncol=2)
# ggsave(paste0(path, "Figure/Simulation1_loadings.pdf"), g)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S6: residuals versus fitted values
data_f8 <- misdata
res_f8 <- list(res_ZIPLL = ZIPLL, prc_res = prcnew)

t10 <- resid_plot(i = 23, j = 69, data = data_f8, result = res_f8)
t11 <- resid_plot(i = 18, j = 10, data = data_f8, result = res_f8)

ggpubr::ggarrange(t10, t11, nrow=1, common.legend=TRUE, legend="bottom", labels="AUTO")
# ggsave(paste0(path, "Figure/Simulation1_fitvsres_lab.pdf"))
