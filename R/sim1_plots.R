# Simulation 1: Table 1 and Figure 4, 5, 6, 7, 8

# dependencies
library(tidyverse)
library(reshape2)
library(gridExtra)
library(Rcpp)
library(ggpubr)
library(tcpl)

# source code
sourceCpp("msf.cpp")
source("simulate_data.R")
source("dosres_plot.R")

# call data
res <- readRDS("sim1_BMC.rds") # BMC
out <- readRDS("sim1_BMC_out1.rds") # BMC results 
res0 <- readRDS("sim1_BMC0.rds") # BMC_0
resi <- readRDS("sim1_BMCi.rds") # BMC_i
resj <- readRDS("sim1_BMCj.rds") # BMC_j
res_tcpl <- readRDS("sim1_tcpl.rds") # tcpl
res_ZIPLL <- readRDS("sim1_ZIPLL.rds") # ZIPLL
prc_res <- readRDS("sim1_processed.rds") # pre-processed posterior predictive results from BMC

#----------------------------------------------------------------------------------------------------------------------------------
# Table 1: summary of results (make mean and sd tables separately)
mres <- matrix(c(round(unlist(lapply(res[1:4], mean)), 3), 
                 round(unlist(lapply(resj[1:4], mean)), 3), 
                 c(round(unlist(lapply(res_ZIPLL[1:2], mean)), 3), NA, NA), 
                 c(round(unlist(lapply(res_tcpl[1:2], mean)), 3), NA, NA)), 4, 4)
rownames(mres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 'In-sample AUC for tij')
colnames(mres) <- c('BMC', 'BMCj', 'ZIPLL', 'tcpl')

sdres <- matrix(c(round(unlist(lapply(res[1:4], sd)), 3), 
                  round(unlist(lapply(resj[1:4], sd)), 3), 
                  c(round(unlist(lapply(res_ZIPLL[1:2], sd)), 3), NA, NA), 
                  c(round(unlist(lapply(res_tcpl[1:2], sd)), 3), NA, NA)), 4, 4)
rownames(sdres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij', 'In-sample AUC for tij')
colnames(sdres) <- c('BMC', 'BMCj', 'ZIPLL', 'tcpl')

cat('mean table\n'); mres
cat('sd table\n'); sdres

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 4: heat map of the estimated and true correlation matrix 
m <- 30; J <- 150
gendata <- generate_data(m, J, d=3, seed=res$seedsave[50])
simdata <- gendata$simdata
truth <- gendata$truth
simdata$rxf <- truth$rxf
misdata <- data_missing(simdata, prob_missing=0.03, missing_idx=NULL, seed=res$seedsave[50])
missing_idx <- misdata$missing_idx 

z_cor_T <- cov2cor(tcrossprod(truth$Lambda) + diag(1, m))
z_cor_data <- melt(z_cor_T)
cor.postm <- cov2cor(out$covMean)
cor.postm_data <- data.frame(melt(cor.postm), truevalue = z_cor_data$value)

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

gridExtra::grid.arrange(t1, t2, nrow=1, widths=c(1.275/3,1.725/3))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 5: the estimated and tru entries of loading matrix 
Rotout <- jointRot(out$Lambda.save, out$eta.save)
lamRot.postm <- Reduce("+", Rotout$lambda)/length(Rotout$lambda)
lamRot.postm[,2] <- (-1)*lamRot.postm[,2] # for better visuatlization

t3 <- ggplot(lamRot.postm %>% reshape2::melt()) +
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

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 6: heat map of estimate and tru profiles of the mean effect 
gamma_ij.save <- out$gamma_ij.save 
z.save <- lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
gamma_ij.save[is.na(gamma_ij.save)] <- sapply(z.save, function(x) x[missing_idx]>0)
gamma_ij.postm <- rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
tr_gamma_ij.postm <- rowMeans(out$gamma_ij.save, dim=2, na.rm=TRUE)
pred <- which(is.na(tr_gamma_ij.postm[6:10,6:10]), arr.ind = TRUE)

frames <- data.frame(Var1=pred[,1], Var2=pred[,2])
t5 <- reshape2::melt(gamma_ij.postm[6:10,6:10]) %>% 
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

t6 <- reshape2::melt(truth$gamma_ij[6:10,6:10]) %>% 
  rename(Chem = Var1, Assay = Var2) %>% 
  ggplot() + geom_tile(aes(Assay, Chem, fill=factor(value)), color="grey50") + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5), aspect.ratio = 1) +
  scale_y_reverse() +
  scale_fill_manual(values=c("1"="red3", "0"="white")) + 
  labs(y="Chemical", x="Assay endpoint", title="Truth", fill=expression(gamma[ij])) 

gridExtra::grid.arrange(t5, t6, nrow=1, widths=c(1.345/3,1.655/3))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 7: dose-response curves
data_f7 <- list(simdata = simdata, truth = truth)
res_f7 <- list(res_ZIPLL = res_ZIPLL, out = out, prc_res = prc_res)

t7 <- dosres_plot(i = 26, j = 67, data = data_f7, result = res_f7)
t8 <- dosres_plot(i = 5, j = 2, data = data_f7, result = res_f7)
t9 <- dosres_plot(i = 24, j = 104, data = data_f7, result = res_f7)

ggarrange(t7, t8, t9, nrow=1, ncol=3, common.legend=TRUE, legend="bottom", labels="AUTO")

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 8: residuals versus fitted values 
data_f8 <- simdata
res_f8 <- list(res_ZIPLL = res_ZIPLL, prc_res = prc_res)

t10 <- resid_plot(i = 26, j = 67, data = data_f8, result = res_f8)
t11 <- resid_plot(i = 5, j = 2, data = data_f8, result = res_f8)

ggarrange(t10, t11, nrow=1, common.legend=TRUE, legend="bottom", labels="AUTO")
