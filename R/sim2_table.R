# Simulation 2
# Table S1: summary of results from simulation 2

# call data
path <- "~/Documents/GitHub/BMC/"
AllRes <- readRDS(paste0(path, "data/sim2_results.rds"))

# BMC
res <- AllRes$res
# BMC_0
res0 <- AllRes$res0
# BMC_i
resi <- AllRes$resi
# BMC_j
resj <- AllRes$resj
# tcpl
res_tcpl <- AllRes$res_tcpl
# ZIPLL
res_ZIPLL <- AllRes$res_ZIPLL

# make mean and sd tables separately
mres <- matrix(c(round(unlist(lapply(res[1:3], mean)), 3), 
                 round(unlist(lapply(resi[1:3], mean)), 3), 
                 c(round(unlist(lapply(res_ZIPLL[1:2], mean)), 3), NA), 
                 c(round(unlist(lapply(res_tcpl[1:2], mean)), 3), NA)), 3, 4)
rownames(mres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij')
colnames(mres) <- c('BMC', 'BMCi', 'ZIPLL', 'tcpl')

sdres <- matrix(c(round(unlist(lapply(res[1:3], sd)), 3), 
                  round(unlist(lapply(resi[1:3], sd)), 3), 
                  c(round(unlist(lapply(res_ZIPLL[1:2], sd)), 3), NA), 
                  c(round(unlist(lapply(res_tcpl[1:2], sd)), 3), NA)), 3, 4)
rownames(sdres) <- c('RMSE', 'In-sample AUC for rij', 'Out-of-sample AUC for rij')
colnames(sdres) <- c('BMC', 'BMCi', 'ZIPLL', 'tcpl')

cat('mean table\n'); mres
cat('sd table\n'); sdres
