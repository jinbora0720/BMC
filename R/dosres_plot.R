# Draw dose-response plots for BMC, ZIPLL and tcpl 

# @param i: chemical index
# @param j: assay endpoint index
# @param data: a list of simdata and truth
# @param result: a list of res_ZIPLL, out, and prc_res
dosres_plot <- function(i, j, data, result) {
  
  # data
  simdata <- data$simdata
  truth <- data$truth
  X <- simdata$X
  orgX <- simdata$orgX
  rxf <- simdata$rxf
  Y <- simdata$orgY
  m_j <- simdata$m_j
  idx_j <- simdata$idx_j
  Start <- simdata$Start
  End <- simdata$End
  
  # ZIPLL
  fit_ZIPLL <- result$res_ZIPLL$fit_ZIPLL
  post.sd_ZIPLL <- fit_ZIPLL$dat[,6]
  post.mn_ZIPLL <- fit_ZIPLL$dat[,5]
  resid_ZIPLL <- fit_ZIPLL$dat[,4] - post.mn_ZIPLL
  
  # BMC
  out <- result$out
  t_ij.postm <- rowMeans(out$t_ij.save, dims=2)
  gamma_ij.save <- out$gamma_ij.save 
  z.save <- lapply(mapply("tcrossprod", out$Lambda.save, out$eta.save, SIMPLIFY = FALSE), function(x) matrix(rnorm(m*J, x, 1), m, J))
  gamma_ij.save[is.na(gamma_ij.save)] <- sapply(z.save, function(x) x[missing_idx]>0)
  gamma_ij.postm <- rowMeans(gamma_ij.save, dim=2, na.rm=TRUE)
  
  # preprocessed BMC results
  prc_res <- result$prc_res
  predY.postL <- prc_res$predY.postL
  predY.postU <- prc_res$predY.postU
  Fit.postm <- prc_res$Fit.postm
  Fit.postL <- prc_res$Fit.postL
  Fit.postU <- prc_res$Fit.postU
    
  # tcpl  
  params <- tcplFit(logc = orgX[[j]][Start[i,j]:End[i,j]],
                    resp = Y[[j]][Start[i,j]:End[i,j]],
                    bmad = (median(Y[[j]]) + sd(Y[[j]]))/3, bidirectional=TRUE)
  aic <- which.min(c(params$cnst_aic, params$hill_aic, params$gnls_aic))
  
  colors <- c("ZIPLL" = "#377EB8", "BMC" = "#E41A1C", "tcpl" = "#4DAF4A") 
  types <- c("ZIPLL" = "dotdash", "BMC" = "longdash", "tcpl" = "dashed")
  
  plotdata <- data.frame(resp = Y[[j]][Start[i,j]:End[i,j]], 
                        respLqt = predY.postL[[j]][Start[i,j]:End[i,j]], 
                        respUqt = predY.postU[[j]][Start[i,j]:End[i,j]], 
                        logc = orgX[[j]][Start[i,j]:End[i,j]], 
                        fit = Fit.postm[[j]][Start[i,j]:End[i,j]],
                        Lqt = Fit.postL[[j]][Start[i,j]:End[i,j]],
                        Uqt = Fit.postU[[j]][Start[i,j]:End[i,j]], 
                        true = rxf[[j]][Start[i,j]:End[i,j]], 
                        fitZ = post.mn_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j], 
                        LqtZ = (post.mn_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j]-2*post.sd_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j]), 
                        UqtZ = (post.mn_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j]+2*post.sd_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j]), 
                        fittcpl = c(rep(0,24), params$hill_modl, params$gnls_modl)[(24*(aic-1)+1):(24*aic)],
                        Lqttcpl = c(rep(0,24)-2*params$cnst_er, params$hill_modl-2*params$hill_er_sd, params$gnls_modl-2*params$gnls_er_sd)[(24*(aic-1)+1):(24*aic)],
                        Uqttcpl = c(rep(0,24)+2*params$cnst_er, params$hill_modl+2*params$hill_er_sd, params$gnls_modl+2*params$gnls_er_sd)[(24*(aic-1)+1):(24*aic)])
  # smoothen true curves and predictive intervals 
  plotdata$respLqt = fitted(loess(respLqt~logc, plotdata))
  plotdata$respUqt = fitted(loess(respUqt~logc, plotdata))
  plotdata$true = fitted(loess(true~logc, plotdata))
  
  p = ggplot(plotdata) + 
    geom_point(aes(logc, resp)) + 
    geom_ribbon(aes(logc, ymin=respLqt, ymax=respUqt), alpha=0.1) +
    geom_line(aes(logc, true), size=1) + 
    geom_ribbon(aes(logc, ymin=LqtZ, ymax=UqtZ, fill="ZIPLL"), alpha=0.3) + 
    geom_ribbon(aes(logc, ymin=Lqttcpl, ymax=Uqttcpl, fill="tcpl"), alpha=0.3) +
    geom_ribbon(aes(logc, ymin=Lqt, ymax=Uqt, fill="BMC"), alpha=0.3) +
    geom_line(aes(logc, fitZ, col= "ZIPLL", linetype="ZIPLL")) +
    geom_line(aes(logc, fittcpl, col= "tcpl", linetype="tcpl")) +
    geom_line(aes(logc, fit, col="BMC", linetype="BMC")) +
    theme_minimal() + theme_bw() + 
    labs(x=expression(paste("Dose (", log[10], " ", mu, "M)", sep="")), y="Response", color="Model", fill="Model", linetype="Model") + 
    scale_color_manual(values = colors) + 
    scale_fill_manual(values = colors) + 
    scale_linetype_manual(values = types) + 
    labs(subtitle=paste("Truth: (Pr(Mean Effect), Pr(Var Effect)) = (", truth$gamma_ij[idx_j[[j]][i],j],",", truth$t_ij[idx_j[[j]][i],j], ")",
                        "\nBMC: (", gamma_ij.postm[idx_j[[j]][i],j],",",t_ij.postm[idx_j[[j]][i],j],")",
                        "\nZIPLL: (", round(fit_ZIPLL$parms[(fit_ZIPLL$parms[,1]==i & fit_ZIPLL$parms[,2]==j),7],3),",", "NA", ")", 
                        "\ntcpl: (", ifelse(aic==1,0,1), ",", "NA", ")")) +
    theme(legend.title=element_text(size=20), legend.text=element_text(size=20))
  
  return(p)
}

# @param i: chemical index
# @param j: assay endpoint index
# @param data: a list of Start and End
# @param result: a list of res_ZIPLL and prc_res
resid_plot <- function(i, j, data, result) {
  
  # data
  Start <- data$Start
  End <- data$End
  
  # ZIPLL
  fit_ZIPLL <- result$res_ZIPLL$fit_ZIPLL
  post.sd_ZIPLL <- fit_ZIPLL$dat[,6]
  post.mn_ZIPLL <- fit_ZIPLL$dat[,5]
  resid_ZIPLL <- fit_ZIPLL$dat[,4] - post.mn_ZIPLL
  
  # preprocessed BMC results
  prc_res <- result$prc_res
  Yhat.postm <- prc_res$Yhat.postm
  res.postm <- prc_res$res.postm
  
  colors <- c("ZIPLL" = "#377EB8", "BMC" = "#E41A1C") 
  shapes <- c("ZIPLL" = 19, "BMC" = 17) 
  
  p = data.frame(fitH = Yhat.postm[[j]][Start[i,j]:End[i,j]], 
             resH = res.postm[[j]][Start[i,j]:End[i,j]], 
             fitZ = post.mn_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j], 
             resZ = resid_ZIPLL[fit_ZIPLL$dat[,1]==i & fit_ZIPLL$dat[,2]==j]) %>% 
    ggplot() + geom_point(aes(fitZ, resZ, color="ZIPLL", shape="ZIPLL"), size=4) + 
    geom_point(aes(fitH, resH, shape="BMC", col="BMC"), size=3) + 
    labs(x="Fitted value", y="Residual", color="Model", shape="Model") +
    scale_color_manual(values = colors) + 
    scale_shape_manual(values = shapes) + 
    theme_bw() + 
    theme(legend.title=element_text(size=20), 
          legend.text=element_text(size=20),
          axis.title = element_text(size=20)) 
  
  return(p)
}

