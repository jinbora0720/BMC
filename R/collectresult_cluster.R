# BMCnew
rmse = rep(NA, 50)
tr_aucg = rep(NA, 50)
tt_aucg = rep(NA, 50)
tr_auct = rep(NA, 50)
tt_auct = rep(NA, 50)
tr_aucap = rep(NA, 50)
tt_aucap = rep(NA, 50)
xi = rep(NA, 50)
alpha = matrix(NA, nrow = 2, ncol = 50)
seedsave = rep(NA, 50)
for (i in 1:50) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new_BMCnew_res_", i, ".RDS"))
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
        "~/BMC/BMC/data/sim1_new_BMCnew_all.RDS")

# # BMC
# rmse = rep(NA, 50)
# tr_aucg = rep(NA, 50)
# tt_aucg = rep(NA, 50)
# tr_auct = rep(NA, 50)
# seedsave = rep(NA, 50)
# for (i in 1:50) {
#   res <- readRDS(paste0("~/BMC/BMC/data/sim1_new_BMC_res_", i, ".RDS"))
#   rmse[i] <- res$rmse
#   tr_aucg[i] <- res$tr_aucg
#   tt_aucg[i] <- res$tt_aucg
#   tr_auct[i] <- res$tr_auct
#   seedsave[i] <- res$seed
# }
# saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, tr_auct = tr_auct, 
#              seedsave = seedsave), 
#         "~/BMC/BMC/data/sim1_new_BMC_all.RDS")

# BMC0
rmse = rep(NA, 50)
tr_aucg = rep(NA, 50)
tt_aucg = rep(NA, 50)
tr_auct = rep(NA, 50)
tr_aucap = rep(NA, 50)
seedsave = rep(NA, 50)
for (i in 1:50) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new_BMC0_res_", i, ".RDS"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  tr_aucap[i] <- res$tr_aucap
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
             tr_auct = tr_auct, tr_aucap = tr_aucap, 
             seedsave = seedsave), 
        "~/BMC/BMC/data/sim1_new_BMC0_all.RDS")

# BMCi
rmse = rep(NA, 50)
tr_aucg = rep(NA, 50)
tt_aucg = rep(NA, 50)
tr_auct = rep(NA, 50)
tr_aucap = rep(NA, 50)
seedsave = rep(NA, 50)
for (i in 1:50) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new_BMCi_res_", i, ".RDS"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  tr_aucap[i] <- res$tr_aucap
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
             tr_auct = tr_auct, tr_aucap = tr_aucap, 
             seedsave = seedsave), 
        "~/BMC/BMC/data/sim1_new_BMCi_all.RDS")

# BMCj
rmse = rep(NA, 50)
tr_aucg = rep(NA, 50)
tt_aucg = rep(NA, 50)
tr_auct = rep(NA, 50)
tr_aucap = rep(NA, 50)
seedsave = rep(NA, 50)
for (i in 1:50) {
  res <- readRDS(paste0("~/BMC/BMC/data/sim1_new_BMCj_res_", i, ".RDS"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  tr_aucap[i] <- res$tr_aucap
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, 
             tr_auct = tr_auct, tr_aucap = tr_aucap, 
             seedsave = seedsave), 
        "~/BMC/BMC/data/sim1_new_BMCj_all.RDS")

