# Data must have 
## X (a list of spline basis matrix),
## Y (a list of responses), 
## orgX (a list of orginal doses), 
## Start (a matrix of starting indices for chemicals when all measurements are aggregated by each assay endpoint), 
## End (a matrix of end indices), 
## e.g. Start[1:3, 1] = 1, 17, 33, End[1:3, 1] = 16, 32, 40
## The first and the second chemical tested on the first assay endpoint have 16 measurements each, 
## while the third chemical has 8 measurements. 
## idx_j (a list of chemical indices (1 to m) with at least one measurement for each assay endpoint), 
## m_j (a vector of number of chemicals tested for each assay endpoint).

# hyper must have v_d (numeric, normal variance for the amount of heteroscedasticity delta_ij).
# init must have q (numeric, initial or fixed number of factors).
# hetero = TRUE to model heteroscedasticity. 
# adapt = TRUE to adapt the number of factors at each iteration.
# simpler = c("bmc0", "bmci", "bmcj") for simpler structure of gamma_ij.
# verbose = TRUE to show progress bars. 

bmc_new <- function(Data = list(), MCMC = list(thin=1, burnin=0, save=1000), 
                    hyper = list(), init = list(), 
                    hetero = TRUE, adapt = TRUE, simpler = FALSE, verbose = TRUE) {

  ##################
  # to get started #
  ##################
  sp <- any(simpler %in% c("bmc0", "bmci", "bmcj"))
  if (adapt & sp) stop("Adaptive number of factors applies for BMC only.")
  
  ## functions for simpler 
  if (simpler == "bmc0") {
    pi.post = function(gamma_ij, m_j) {
      m = nrow(gamma_ij); J = ncol(gamma_ij)
      RCgam = sum(gamma_ij, na.rm=TRUE)
      return(matrix(rbeta(m*J, 1+RCgam, 1+sum(m_j)-RCgam), m, J))
    }
  } else if (simpler == "bmci") {
    pi.post = function(gamma_ij, idx_j) {
      m = nrow(gamma_ij); J = ncol(gamma_ij)
      J_i = sapply(1:m, function(i) sum(unlist(idx_j)==i))
      RSgam = rowSums(gamma_ij, na.rm=TRUE)
      pi_i = rbeta(m, 1+RSgam, 1+J_i-RSgam)
      return(matrix(rep(pi_i, J), m, J))
    }
  } else {
    pi.post = function(gamma_ij, m_j) {
      m = nrow(gamma_ij); J = ncol(gamma_ij)
      CSgam = colSums(gamma_ij, na.rm=TRUE)
      pi_j = rbeta(J, 1+CSgam, 1+m_j-CSgam)
      return(t(matrix(rep(pi_j, m), J, m)))
    }
  }
  
  ########
  # Data #
  ########
  Y = Data$Y 
  X = Data$X
  orgX = Data$orgX
  p = ncol(X[[1]])
  Start = Data$Start 
  End = Data$End 
  idx_j = Data$idx_j 
  m_j = Data$m_j
  J = length(Y) 
  m = max(unlist(idx_j))
  
  ###################
  # MCMC parameters #
  ###################
  thin = MCMC$thin 
  burnin = MCMC$burnin
  save = MCMC$save
  S = burnin+thin*save
  
  ###################
  # hyperparameters #
  ###################
  if (!sp) {
    ## for factor model: 
    if (is.null(hyper$df)) {
      df = 3                    # gamma hyperparameter (nu) for local shrinkage parameters
    } else {
      df = hyper$df
    }
    if (is.null(hyper$ad1)) {
      ad1 = 2.1                 # gamma hyperparameters for zeta_1
    } else {
      ad1 = hyper$ad1
    }                  
    if (is.null(hyper$bd1)) {
      bd1 = 1
    } else {
      bd1 = hyper$bd1
    } 
    if (is.null(hyper$ad2)) {
      ad2 = 3.1                 # gamma hyperparameters for zeta_2
    } else {
      ad2 = hyper$ad2
    }                  
    if (is.null(hyper$bd2)) {
      bd2 = 1
    } else {
      bd2 = hyper$bd2
    } 
  }

  ## Sigj:
  if (is.null(hyper$a)) {
    a = p+2                     # df parameter for a Wishart distribution
  } else {
    a = hyper$a                   
  }
  if (is.null(hyper$R)) {
    Sigj = sapply(1:J, function(j) {
      covtemp = sapply(1:m_j[j], function(i) {
        coef(lm(Y[[j]][Start[i,j]:End[i,j]]~matrix(X[[j]][Start[i,j]:End[i,j],],ncol=p)-1))
      }) %>% t(.) 
      if (nrow(na.omit(covtemp))>1) { cov(covtemp, use="complete.obs") }
      else {rep(0, p*p)}
    }) 
    Sigj[,colSums(Sigj)==0] = rowSums(Sigj)/(J-sum(colSums(Sigj)==0))
    invSigj = apply(Sigj, 2, function(x) c(ginv(matrix(x, p, p))))
    R = matrix((a-p-1)*rowMeans(Sigj, na.rm=TRUE), p, p)
  } else {
    R = hyper$R                # covariance matrix for a Wishart distribution
  }
                     
  ## normal error:
  if (is.null(hyper$nu0)) {
    nu0 = 1                   # gamma hyperparameters for assay-specific scales 
  } else {
    nu0 = hyper$nu0               
  }
  if (is.null(hyper$sig_sq)) {
    sig0_sq = var(unlist(Y))
  } else {
    sig0_sq = hyper$sig0_sq
  }
  
  if (hetero) {
    ## heteroscedasticity:
    v_d = hyper$v_d             # normal variance for the amount of heteroscedasticity delta_ij 
  }
  
  ##################
  # initial values #
  ##################
  if (!sp) {
    q = init$q                                                  # initial or fixed number of factors
    ## for factor model:
    if (is.null(init$Lambda)) {
      Lambda = matrix(1, nrow = m, ncol = q)                    # factor loadings
    } else {
      Lambda = init$Lambda
    }
    if (is.null(init$eta)) {
      eta = matrix(rnorm(J*q), nrow = J, ncol = q)                   # latent factors
    } else {
      eta = init$eta
    }
    if (is.null(init$delta)) {
      delta = c(rgamma(1,ad1,bd1), rgamma(q-1,ad2,bd2))         # global shrinkage coefficients multipliers
    } else{
      delta = init$delta
    }
    psijh = matrix(rgamma(m*q, df/2, df/2), nrow = m, ncol = q) # local shrinkage parameters 
    tauh = cumprod(delta)                                       # global shrinkage parameters 
    Plam = t(t(psijh) * (tauh))                                 # precision of loadings  
    z_ij = t(mvrnorm(J, mu=rep(0,m), Sigma=tcrossprod(Lambda)+diag(1,m)))
    if (is.null(init$gamma_ij)) {
      gamma_ij = matrix(rbinom(m*J, 1, pnorm(z_ij)), m, J)
      for (j in 1:J) {
        z_ij[-idx_j[[j]],j] = NA
        gamma_ij[-idx_j[[j]],j] = NA
      }
    } else {
      gamma_ij = init$gamma_ij
    }
  } else {
    if (is.null(init$gamma_ij)) {
      gamma_ij = matrix(rbinom(m*J, 1, 0.5), m, J)
      for (j in 1:J) {
        gamma_ij[-idx_j[[j]],j] = NA
      }
    } else {
      gamma_ij = init$gamma_ij
    }
  }
  
  ## Sigj:
  if (is.null(invSigj)) {
    if (is.null(init$invSigj)) {
      invSigj = rWishart(1, df = a, Sigma = solve(R))[,,1]
    } else {
      invSigj = init$invSigj
    }
  }
  if (is.null(Sigj)) {
    if (is.null(init$Sigj)) {
      Sigj = ginv(invSigj)
    } else {
      Sigj = init$Sigj
    }
  }
  
  ## beta_ij:
  if (is.null(init$beta_ij)) {
    ols = lm(unlist(Y)~do.call("rbind",X)-1)
    beta_ij = list()
    for (j in 1:J) {
      beta_ij[[j]] = matrix(coef(ols), nrow=p, ncol=m_j[j])
    }
  } else {
    beta_ij = init$beta_ij
  }

  ## normal error:
  if (is.null(init$sigj_sq)) {
    sigj_sq = 1/rgamma(J, nu0/2, nu0*sig0_sq/2)
  } else {
    sigj_sq = init$sigj_sq
  }
  
  if (hetero) {
    ## heteroscedasticity:
    alpha = c(0,1)
    le = tcrossprod(Lambda, eta) 
    W = alpha[1] + alpha[2]*le
    u_ij = W + matrix(rnorm(m*J), nrow = m, ncol = J)
    t_ij = matrix(rbinom(m*J, 1, pnorm(W)), m, J)
    for (j in 1:J) {
      t_ij[-idx_j[[j]],j] = NA
    }
    if (is.null(init$d_ij)) {
      d_ij = matrix(NA, m, J)
      d_ij[which(t_ij==1, arr.ind=TRUE)] = rnorm(sum(t_ij==1, na.rm=TRUE), 0, sqrt(v_d))
      d_ij[which(t_ij==0, arr.ind=TRUE)] = 0
    } else {
      d_ij = init$d_ij
    }
  }
  
  ################
  # store output #
  ################
  out = list() 
  out$gamma_ij.save = array(NA, dim = c(m,J,save))
  out$beta_ij.save = list()
  for (j in 1:J) {
    out$beta_ij.save[[j]] = array(NA, dim=c(p, m_j[j], save))
  }
  out$Sigj.save = array(NA, dim=c(p*p,J,save))
  out$sigj_sq.save = matrix(NA, J, save)
  if (sp) { 
    out$pi_ij.save = array(NA, dim = c(m,J,save)) 
  } else {
    out$Lambda.save = list()
    out$eta.save = list()
    out$covMean = matrix(0, m, m)
  }
  if (hetero) {
    out$alpha.save = matrix(0, 2, save)
    out$t_ij.save = array(NA, dim = c(m,J,save))
    out$d_ij.save = array(NA, dim = c(m,J,save))
    out$accept = 0
  }
  
  ########
  # MCMC #
  ########
  if (verbose) {
    cat('burnin...\n')
    pb = txtProgressBar(style=3,width=50)
    pb2 = txtProgressBar(style=3,width=50)
  }
  
  for (s in 1:S) {

    if (hetero) {
      # 1. heteroscedasticity
      # update eps
      eps = ProdRes(Y, X, beta_ij, Start, End, m_j)
      
      # update t_ij, d_ij, pi_t, v_d related to heteroscedasticity
      alpha = alpha.post(le, u_ij)
      u_ij = u_ij.post(le, alpha, t_ij)
      td_ij = td_ij.post(eps, orgX, sigj_sq, v_d, t_ij, d_ij, idx_j, m_j, le, alpha, Start, End)
      t_ij = td_ij$t_ij
      d_ij = td_ij$d_ij
      
      # 2. fitting regression parameters
      Y2 = DvidD(Y, orgX, d_ij, idx_j, Start, End, m_j)
      X2 = DvidD_mat(X, orgX, d_ij, idx_j, Start, End, m_j)
    } 
    
    if (sp) {
      if (simpler == "bmci") {
        pi_ij = pi.post(gamma_ij, idx_j)
      } else { pi_ij = pi.post(gamma_ij, m_j) }
    } else {
      # update lambda, eta, z_ij and pi_ij
      fm = linearMGSP(t(z_ij), q, Lambda, delta, tauh, Plam, df, ad1, bd1, ad2, bd2, s, 
                      adapt=adapt, prop=1, epsilon=1e-3)
      Lambda = fm$Lambda
      eta = fm$eta
      le = tcrossprod(Lambda, eta) 
      delta = fm$delta
      tauh = fm$tauh
      Plam = fm$Plam
      Omega = fm$Omega
      if (adapt) { q = fm$k }
      z_ij = z_ij.post(le, gamma_ij)
      for (j in 1:J) {
        z_ij[-idx_j[[j]],j] = NA
      }
      pi_ij = pnorm(le) 
    }
    
    # update Sigj
    Sig = Sigj.post(a, m_j, R, beta_ij)
    invSigj = Sig$invSigj
    Sigj = Sig$Sigj
    
    # update gamma_ij, beta_ij 
    if (hetero) { 
      gbcpp = gbcpp_post(Y2, X2, pi_ij, invSigj, Sigj, sigj_sq, idx_j, m_j, Start, End)
    } else {
      gbcpp = gbcpp_post(Y, X, pi_ij, invSigj, Sigj, sigj_sq, idx_j, m_j, Start, End)
    }
    gamma_ij = gbcpp$gamma_ij
    beta_ij = gbcpp$beta_ij
    
    # 3. normal error
    # update eps
    if (hetero) {
      eps2 = ProdRes(Y2, X2, beta_ij, Start, End, m_j)
    } else {
      eps2 = ProdRes(Y, X, beta_ij, Start, End, m_j)
    }
    
    # update sig_sq 
    sigj_sq = sigj_sq.post(eps2, nu0, sig0_sq)

    # save results 
    if ((s > burnin) & (s-burnin) %% thin==0) {
      out$sigj_sq.save[,(s-burnin)/thin] = sigj_sq
      out$gamma_ij.save[,,(s-burnin)/thin] = gamma_ij
      for (j in 1:J) {
        out$beta_ij.save[[j]][,,(s-burnin)/thin] = beta_ij[[j]]
      }
      out$Sigj.save[,,(s-burnin)/thin] = Sigj
      if (sp) { 
        out$pi_ij.save[,,(s-burnin)/thin] = pi_ij 
      } else {
        out$Lambda.save[[(s-burnin)/thin]] = Lambda
        out$eta.save[[(s-burnin)/thin]] = eta
        out$covMean = out$covMean + Omega/save
      }
      if (hetero) {
        out$alpha.save[,(s-burnin)/thin] = alpha
        out$t_ij.save[,,(s-burnin)/thin] = t_ij
        out$d_ij.save[,,(s-burnin)/thin] = d_ij
        out$accept = out$accept + td_ij$accept/save
      }
    }
    
    if (verbose) {
      setTxtProgressBar(pb, s/burnin) 
      if (s > burnin) {
        if (s == burnin + 1) { 
          close(pb)
          cat('saving...\n')
          setTxtProgressBar(pb2,(s-burnin)/(S-burnin))
        } else setTxtProgressBar(pb2,(s-burnin)/(S-burnin))
      }
    }
  }
  if (verbose) { close(pb2) } 
  return(out)
}
    
