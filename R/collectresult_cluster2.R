# BMCnew
rmse = rep(NA, 30)
tr_aucg = rep(NA, 30)
tt_aucg = rep(NA, 30)
tr_auct = rep(NA, 30)
tt_auct = rep(NA, 30)
tr_aucap = rep(NA, 30)
tt_aucap = rep(NA, 30)
xi = rep(NA, 30)
alpha = matrix(NA, nrow = 2, ncol = 30)
seedsave = rep(NA, 30)
for (i in 1:30) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new2_BMCnew_res_", i, ".RDS"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  tt_auct[i] <- res$tt_auct
  tr_aucap[i] <- res$tr_aucap
  tt_aucap[i] <- res$tt_aucap
  alpha[,i] <- res$alpha
  xi[i] <- res$xi
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
             tr_auct = tr_auct, tt_auct = tt_auct, 
             tr_aucap = tr_aucap, tt_aucap = tt_aucap, 
             alpha = alpha, xi = xi, seedsave = seedsave), 
        "~/BMC/BMC/data/sim1_new2_BMCnew_all.RDS")

# BMC0
rmse = rep(NA, 30)
tr_aucg = rep(NA, 30)
tt_aucg = rep(NA, 30)
tr_auct = rep(NA, 30)
tt_auct = rep(NA, 30)
tr_aucap = rep(NA, 30)
tt_aucap = rep(NA, 30)
seedsave = rep(NA, 30)
for (i in 1:30) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new2_BMC0_res_", i, ".RDS"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  tt_auct[i] <- res$tt_auct
  tr_aucap[i] <- res$tr_aucap
  tt_aucap[i] <- res$tt_aucap
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
             tr_auct = tr_auct, tt_auct = tt_auct, 
             tr_aucap = tr_aucap, tt_aucap = tt_aucap,
             seedsave = seedsave), 
        "~/BMC/BMC/data/sim1_new2_BMC0_all.RDS")
