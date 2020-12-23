# EDA plots 
# from obesity: Figure S1
# from neurodevelopmental disorders: Figure 1, S2, S3, S4

# dependencies 
library(tidyverse)

# call data 
path <- "~/Documents/GitHub/BMC/"
datalist <- readRDS(paste0(path, "data/datalist.rds"))
neuro_data <- datalist$neuro_data
obese_data <- datalist$obese_data

# arrage data
obese_meta <- readRDS(paste0(path, "data/obese_meta.rds"))
K_ij <- obese_meta$K_ij

neuro_meta <- readRDS(paste0(path, "data/neuro_meta.rds"))
uniq_chnm <- neuro_meta$uniq_chnm
uniq_aenm <- neuro_meta$uniq_aenm
J <- length(uniq_aenm)
idx_j <- neuro_meta$idx_j
Start <- neuro_meta$Start
End <- neuro_meta$End
Y <- neuro_meta$Y
orgX <- lapply(neuro_meta$orgData, function(x) x[,1])

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S1: Heat map of the number of observations in ToxCast/Tox21 data for obesity
reshape2::melt(K_ij) %>% 
  mutate(countfactor=cut(value, breaks=c(-1,0,1,10,100,500,max(value)), 
                         labels=c("0","1","2-10","11-100","101-500","501-870"))) %>% 
  mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor)))) %>% 
  rename(Chem = Var1, Assay = Var2) %>% 
  ggplot() + geom_tile(aes(Assay, Chem, fill=countfactor)) +
  guides(fill=guide_legend(title="Number of \nObservations")) +
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","white")) +
  scale_y_reverse() + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size=20), 
        legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  labs(y="Chemical", x="Assay endpoint") + scale_x_continuous(position = 'top')

rm(obese_meta, K_ij)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure 1: Detailed illustration of the ToxCast/Tox21 data structure. 
chnm_sub <- uniq_chnm[c(9, 11, 16, 22, 23, 25, 27)]
aenm_sub <- uniq_aenm[c(2, 30, 101, 120, 125)]
neuro_sub <- neuro_data %>% filter(chnm %in% chnm_sub, aenm %in% aenm_sub) %>% select(logc, resp, chnm, aenm) 
neuro_sub$aenm <- plyr::revalue(neuro_sub$aenm, c("NVS_ADME_hCYP19A1_Activator"="NVS_ADME_\nhCYP19A1_Activator", 
                                                  "CEETOX_H295R_OHPROG_up"="CEETOX_H295R_\nOHPROG_up", 
                                                  "TOX21_AR_LUC_MDAKB2_Antagonist"="TOX21_AR_LUC_\nMDAKB2_Antagonist", 
                                                  "OT_AR_ARSRC1_0960"="OT_AR_\nARSRC1_0960")) 
neuro_sub$chnm <- plyr::revalue(neuro_sub$chnm, c("Diethyl phthalate"="Diethyl\nphthalate",
                                                  "Dibutyl phthalate"="Dibutyl\nphthalate", 
                                                  "Cypermethrin"="Cyper-\nmethrin"))
neuro_sub %>% 
  ggplot() + geom_point(aes(logc, resp)) + facet_grid(chnm~aenm, scales="free") +
  theme_bw() + 
  labs(x=expression(paste("Dose (", log[10], " ", mu, "M)", sep="")), y="Response") +
  theme(axis.title= element_text(size=20))

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S3: Scatter plots of two chemicals on TOX21_ERa_LUC_BG1_Agonist
i <- which(uniq_chnm == "Di-n-octyl phthalate")
j <- which(uniq_aenm == "TOX21_ERa_LUC_BG1_Agonist")

inew <- which(idx_j[[j]] == i)
t1 <- data.frame(resp = Y[[j]][Start[inew,j]:End[inew,j]], 
                 logc = orgX[[j]][Start[inew,j]:End[inew,j]]) %>% 
  ggplot(aes(logc, resp)) + geom_point() + 
  geom_smooth(formula = y~x, method = loess) + 
  labs(x=expression(paste("Dose (", log[10], " ", mu, "M)", sep="")), y="Response") +
  theme_minimal() + theme_bw() + theme(aspect.ratio = 0.7) + 
  labs(title=uniq_chnm[i])

i <- which(uniq_chnm == "2,4,5-Trichlorophenol")
j <- which(uniq_aenm == "TOX21_ERa_LUC_BG1_Agonist")

inew <- which(idx_j[[j]] == i)
t2 <- data.frame(resp = Y[[j]][Start[inew,j]:End[inew,j]], 
                 logc = orgX[[j]][Start[inew,j]:End[inew,j]]) %>% 
  ggplot(aes(logc, resp)) + geom_point() + 
  geom_smooth(formula = y~x, method = loess) +
  labs(x=expression(paste("Dose (", log[10], " ", mu, "M)", sep="")), y="Response") +
  theme_minimal() + theme_bw() + theme(aspect.ratio = 0.7) + 
  labs(title=uniq_chnm[i])

gridExtra::grid.arrange(t1, t2, nrow=1)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S2: Histogram of the number of unique doses tested. 
neuro_data %>% group_by(casn, aenm) %>% summarise(n_dose = length(unique(logc))) %>% 
  ggplot() + geom_histogram(aes(n_dose), fill="#666666", col="black", bins = 30) + 
  theme_bw() + 
  geom_vline(aes(xintercept=8), col="red") + 
  scale_x_continuous(breaks=seq(0,22,2)) + 
  theme(aspect.ratio = 0.6, axis.title = element_text(size=20)) + 
  labs(x="Number of unique doses", y="Count of chemical-\nassay endpoint pairs")

# "median number of unique doses is 8"
neuro_dose <- neuro_data %>% group_by(casn, aenm) %>% summarise(n_dose = length(unique(logc))) 
median(neuro_dose$n_dose)

# "30% of them are without replicates" 
neuro_rep <- neuro_data %>% group_by(casn, aenm, logc) %>% summarise(n_rep = n()) 
sum(neuro_rep$n_rep == 1)/nrow(neuro_rep)

rm(neuro_dose, neuro_rep)

#----------------------------------------------------------------------------------------------------------------------------------
# Figure S4: Scatter plot of the responses by assay endpoint. 
data.frame(res=unlist(Y[75:J]), aenm=rep(75:J, apply(End[,75:J], 2, max))) %>% 
  ggplot() + geom_point(aes(as.factor(aenm), res)) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title = element_text(size=20), aspect.ratio = 1/1.8) + 
  labs(y="Response", x="Assay endpoint")
