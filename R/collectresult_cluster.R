rmse = rep(NA, 50)
tr_aucg = rep(NA, 50)
tt_aucg = rep(NA, 50)
tr_auct = rep(NA, 50)
seedsave = rep(NA, 50)
for (i in 1:50) {
  res <- readRDS(paste0("sim1_BMC0_res_", i, ".rds"))
  rmse[i] <- res$rmse
  tr_aucg[i] <- res$tr_aucg
  tt_aucg[i] <- res$tt_aucg
  tr_auct[i] <- res$tr_auct
  seedsave[i] <- res$seed
}
saveRDS(list(rmse = rmse, tr_aucg = tr_aucg, tt_aucg = tt_aucg, tr_auct = tr_auct, seedsave = seedsave), 
        "sim1_BMC0_new.rds")
