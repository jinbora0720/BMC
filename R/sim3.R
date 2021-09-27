# Simulation 3: multiplicity adjustment 
rm(list = ls())

## m = 5, J = 5
res5 <- readRDS("~/BMC/data/sim3_m5J5.RDS")
out5 <- res5$out
le5.save <- mapply("tcrossprod", out5$Lambda.save, out5$eta.save, SIMPLIFY = FALSE)
le5.postm <- meshed::summary_list_mean(le5.save)
xi5 <- mean(out5$xi.save)
alpha5 <- rowMeans(out5$alpha.save)

## m = 5, J = 6
resJ6 <- readRDS("~/BMC/data/sim3_m5J6.RDS")
outJ6 <- resJ6$out
leJ6.save <- mapply("tcrossprod", outJ6$Lambda.save, outJ6$eta.save, SIMPLIFY = FALSE)
leJ6.postm <- meshed::summary_list_mean(leJ6.save)
xiJ6 <- mean(outJ6$xi.save)
alphaJ6 <- rowMeans(outJ6$alpha.save)

## m = 5, J = 20
resJ20 <- readRDS("~/BMC/data/sim3_m5J20.RDS")
outJ20 <- resJ20$out
leJ20.save <- mapply("tcrossprod", outJ20$Lambda.save, outJ20$eta.save, SIMPLIFY = FALSE)
leJ20.postm <- meshed::summary_list_mean(leJ20.save)
xiJ20 <- mean(outJ20$xi.save)
alphaJ20 <- rowMeans(outJ20$alpha.save)

## m = 5, J = 50
resJ50 <- readRDS("~/BMC/data/sim3_m5J50.RDS")
outJ50 <- resJ50$out
leJ50.save <- mapply("tcrossprod", outJ50$Lambda.save, outJ50$eta.save, SIMPLIFY = FALSE)
leJ50.postm <- meshed::summary_list_mean(leJ50.save)
xiJ50 <- mean(outJ50$xi.save)
alphaJ50 <- rowMeans(outJ50$alpha.save)

## m = 5, J = 100
resJ100 <- readRDS("~/BMC/data/sim3_m5J100.RDS")
outJ100 <- resJ100$out
leJ100.save <- mapply("tcrossprod", outJ100$Lambda.save, outJ100$eta.save, SIMPLIFY = FALSE)
leJ100.postm <- meshed::summary_list_mean(leJ100.save)
xiJ100 <- mean(outJ100$xi.save)
alphaJ100 <- rowMeans(outJ100$alpha.save)

##########
# tables #
##########
## Table2: Multiplicity adjustment
tab_J <- matrix(0, nrow = 5, ncol = 4)
rownames(tab_J) <- c("J=5", "J=6", "J=20", "J=50", "J=100")
colnames(tab_J) <- c("FP_gamma", "FP_t", "alpha0", "xi")
tab_J[,1] <- c(res5$fpg[6], resJ6$fpg[6], resJ20$fpg[6], resJ50$fpg[6], resJ100$fpg[6])
tab_J[,2] <- c(res5$fpt[6], resJ6$fpt[6], resJ20$fpt[6], resJ50$fpt[6], resJ100$fpt[6])
tab_J[,3] <- c(alpha5[1], alphaJ6[1], alphaJ20[1], alphaJ50[1], alphaJ100[1])
tab_J[,4] <- c(xi5, xiJ6, xiJ20, xiJ50, xiJ100)
round(tab_J, 3)
