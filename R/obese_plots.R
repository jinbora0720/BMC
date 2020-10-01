# Plots related to obesity 

# call data 
datalist <- readRDS("~/Documents/GitHub/BMC/datalist.rds")
obese_data <- datalist$neuro_data
hit_vec <- datalist$hit_vec

# arrage data

K_inf = obese_data %>% group_by(aeid, aenm, casn, chnm) %>% 
  summarise(K_ij = n()) %>% as.data.frame()
uniq_aeid = sort(unique(K_inf$aeid))
uniq_aenm = unique(K_inf$aenm)[order(unique(K_inf$aeid))]
J = length(uniq_aeid) # 271
uniq_casn = sort(unique(K_inf$casn))
m = length(uniq_casn) # 30
uniq_chnm = unique(K_inf$chnm)[order(unique(K_inf$casn))]

casn_j = chnm_j = idx_j = list() # list of chemicals for jth assay endpoint 
m_j = rep(0, J) # number of chemicals for jth assay endpoint
for (j in 1:J) {
  chem_j = K_inf %>% filter(aeid==uniq_aeid[j]) %>% dplyr::select(casn, chnm)
  idx = which(uniq_casn %in% chem_j$casn)
  idx_j[[j]] = idx
  casn_j[[j]] = uniq_casn[idx]
  chnm_j[[j]] = uniq_chnm[idx]
  m_j[j] = length(idx)
}

# overall index
K_ij = matrix(0, m, J) # number of measurements
for (j in 1:J) {
  dat = K_inf %>% filter(aeid==uniq_aeid[j], casn %in% casn_j[[j]]) %>% select(K_ij)
  K_ij[idx_j[[j]],j] = dat$K_ij
}

reshape2::melt(K_ij) %>% 
  mutate(countfactor=cut(value, breaks=c(-1,0,1,10,100,500,max(value)), 
                         labels=c("0","1","2-10","11-100","101-500","501-870"))) %>% 
  mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor)))) %>% 
  rename(Chem = Var1, Assay = Var2) %>% 
  ggplot() + geom_tile(aes(Assay, Chem, fill=countfactor), color="grey90", size=0.2) +
  guides(fill=guide_legend(title="Number of \nObservations")) +
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","white")) +
  scale_y_reverse() + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size=20), 
        legend.title = element_text(size=20), legend.text = element_text(size=20)) +
  labs(y="Chemical", x="Assay endpoint") + scale_x_continuous(position = 'top')

ggsave("DataMatWhole2_obese.png")
ggsave("DataMatWhole2_obese.pdf")