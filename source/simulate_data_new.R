# Functions to create simulated datasets 

simulate_lambda = function(m, q, df=3, ad1=2.1, ad2=3.1) {
  # m: number of chemicals 
  # q: number of factors 
  
  psijh = matrix(rgamma(m*q, df/2, df/2), nrow = m, ncol = q)     
  delta = c(rgamma(1,ad1,1), rgamma(q-1,ad2,1))       
  tauh = cumprod(delta)                                       
  Plam = t(t(psijh) * (tauh))                                    
  return(matrix(rnorm(m*q, mean=0, sd=1/sqrt(Plam)), m, q))
}

generate_data = function(m, J, d=1, seed) {
  # Assumes equal number of measurements (d) for each dose 
  # m: number of chemicals 
  # J: number of assay endpoints 
  # q: number of factors 
  
  set.seed(seed)
  
  # number of chemicals for each assay endpoint 
  m_j = rep(m, J) 
  
  # index for tested chemicals for each assay endpoint 
  idx_j = list()
  for (j in 1:J) {
    idx_j[[j]] = sample(m, m_j[j], replace=FALSE) %>% sort()
  }
  
  # number of measurements for each chemical-assay endpoint combination 
  K_ij = matrix(d*8, m, J)   
  for (j in 1:J) {
    K_ij[-idx_j[[j]], j] = 0
  }
  
  Start = End = matrix(0, m, J)
  for (j in 1:J) {
    Start[1:m_j[j],j] = cumsum(c(1,K_ij[,j]))[idx_j[[j]]]
    End[1:m_j[j],j] = cumsum(K_ij[,j])[idx_j[[j]]]
  }
  
  # Lambda, Eta
  q = 2
  Lambda = simulate_lambda(m, q)
  eta = matrix(rnorm(J*q), J, q)
  
  pi_ij = pnorm(tcrossprod(Lambda, eta))
  gamma_ij = matrix(rbinom(m*J, 1, pi_ij), m, J)
  for (j in 1:J) {
    gamma_ij[-idx_j[[j]], j] = NA
  }
  
  # norm error  
  nu0 = 5; sig0_sq = 0.1
  sigj_sq = 1/rgamma(J, nu0/2, nu0*sig0_sq/2)   
  resid = lapply(1:J, function(j) {
    rnorm(max(End[,j]), 0, sqrt(sigj_sq[j]))
  })
  
  p = 9
  beta_pool = list()
  beta_pool[[1]] = c(0.5, 0.5, 1, 1.3, 1.5, 2, 2.5, 3, 3) # hill
  beta_pool[[2]] = c(0, 0, 0, 0.3, 0.8, 1.3, 2.5, 1.5, 0.8) # gain-loss
  beta_pool[[3]] = c(0, 0, 0, 0, 0, 0, -0.5, -0.7, -2) # decreasing
  
  beta_ij = list()
  beta_ij[1:floor(J/3)] = lapply(1:floor(J/3), function(j) {
    matrix(beta_pool[[1]], p, m_j[j])
  }) 
  beta_ij[(1+floor(J/3)):(2*floor(J/3))] = lapply((1+floor(J/3)):(2*floor(J/3)), function(j) {
    matrix(beta_pool[[2]], p, m_j[j])
  }) 
  beta_ij[(1+2*floor(J/3)):J] = lapply((1+2*floor(J/3)):J, function(j) {
    matrix(beta_pool[[3]], p, m_j[j])
  }) 
  # let some (i,j) pairs to be non-active 
  for (j in 1:J) {
    for (i in 1:m_j[j]) {
      if (gamma_ij[idx_j[[j]][i],j]==0) beta_ij[[j]][,i] = rep(0,p)
    }
  }
  
  # heteroscedasticity 
  t_ij = d_ij = matrix(NA, m, J)
  le = tcrossprod(Lambda, eta)
  alpha = c(-.5, .7)
  W = alpha[1] + alpha[2]*le
  t_ij = matrix(rbinom(m*J, 1, pnorm(W)), m, J)
  for (j in 1:J) {
    t_ij[-idx_j[[j]],j] = NA
  }
  d_ij[which(t_ij==1, arr.ind=TRUE)] = rnorm(sum(t_ij==1, na.rm=TRUE), 2, 0.1)
  d_ij[which(t_ij==0, arr.ind=TRUE)] = 0
  
  # Data
  logc = rep(c(0.301, 0.477, 0.602, 0.845, 1.000, 1.301, 1.602, 2.000), d)
  Xtemp = apply(bs(logc, df=p, intercept=TRUE), 2, function(x) scale(x, center=TRUE, scale=FALSE)) 
  
  orgData = rxf = list()
  for (j in 1:J) {
    orgData[[j]] = matrix(NA, max(End[,j]), 2)
    rxf[[j]] = rep(NA, max(End[,j]))
    for (i in 1:m_j[j]) {
      rxf[[j]][Start[i,j]:End[i,j]] = Xtemp%*%matrix(beta_ij[[j]][,i], ncol=1)
      resp = rxf[[j]][Start[i,j]:End[i,j]]+exp(logc*d_ij[idx_j[[j]][i], j]/2)*resid[[j]][Start[i,j]:End[i,j]]
      orgData[[j]][Start[i,j]:End[i,j],] = cbind(logc, resp)
    }
  }
  orgX = lapply(orgData, function(x) x[,1]) 
  orgY = lapply(orgData, function(x) x[,2])
  
  p1 = 7
  CX = apply(bs(logc, df=p1, intercept=TRUE), 2, function(x) scale(x, center=TRUE, scale=FALSE))
  X = replicate(J, 
                do.call("rbind", replicate(m, CX, simplify = FALSE)), 
                simplify = FALSE)
  
  return(list(simdata = list(orgX=orgX, orgY=orgY, X=X, K_ij=K_ij, m_j=m_j, idx_j=idx_j, Start=Start, End=End),
              truth = list(q=q, Lambda=Lambda, eta=eta, gamma_ij=gamma_ij, rxf=rxf, 
                           t_ij=t_ij, d_ij=d_ij, sigj_sq=sigj_sq)))
}

data_missing = function(simdata, prob_missing=NULL, missing_idx=NULL, seed) {
  set.seed(seed)
  
  m = nrow(simdata$K_ij)
  J = ncol(simdata$K_ij)
  
  # missing pairs  
  if (is.null(missing_idx)) missing_idx = sample(m*J, floor(prob_missing*m*J))
  K_ij = simdata$K_ij
  K_ij[missing_idx] = 0
  
  ## check whether there are chemicals/assay endpoints that are completely missing 
  if (sum(rowSums(K_ij)==0) > 0) {
    is = which(rowSums(K_ij)==0)
    js = rep(0,length(is))
    for (s in 1:length(is)) {
      js[s] = sample(missing_idx[which(missing_idx[,1]==is[s]),2],1)
    }
    K_ij[cbind(is, js)] = simdata$K_ij[cbind(is, js)]
  }
  if (sum(colSums(K_ij)==0) > 0) {
    js = which(colSums(K_ij)==0)
    is = rep(0,length(js))
    for (s in 1:length(js)) {
      is[s] = sample(missing_idx[which(missing_idx[,2]==js[s]),1],1)
    }
    K_ij[cbind(is, js)] = simdata$K_ij[cbind(is, js)]
  }
  
  # index for tested chemicals for each assay endpoint 
  idx_j = apply(K_ij, 2, function(x) which(x>0))
  # number of chemicals for each assay endpoint 
  m_j = apply(K_ij, 2, function(x) length(which(x>0)))
  
  Start = End = matrix(0, m, J)
  for (j in 1:J) {
    Start[1:m_j[j],j] = cumsum(c(1,K_ij[,j]))[idx_j[[j]]]
    End[1:m_j[j],j] = cumsum(K_ij[,j])[idx_j[[j]]]
  }
  
  # training data 
  if (is.null(simdata$rxf)) {
    X = simdata$X
    orgX = simdata$orgX
    orgY = simdata$orgY
    missing_idx2 = which(K_ij==0 & simdata$K_ij>0, arr.ind = TRUE)
    for (s in 1:nrow(missing_idx2)) {
      j = missing_idx2[s,2]
      i = which(simdata$idx_j[[j]] == missing_idx2[s,1])
      X[[j]][simdata$Start[i,j]:simdata$End[i,j],] = NA
      orgX[[j]][simdata$Start[i,j]:simdata$End[i,j]] = NA
      orgY[[j]][simdata$Start[i,j]:simdata$End[i,j]] = NA
    }
    X = lapply(X, na.omit)
    orgX = lapply(orgX, na.omit)
    orgY = lapply(orgY, na.omit)
    
    return(list(missing_idx = missing_idx2, missing_idx_col = missing_idx,
                X=X, orgX=orgX, orgY=orgY, K_ij=K_ij, m_j=m_j, idx_j=idx_j, Start=Start, End=End))
  } else {
    X = simdata$X
    orgX = simdata$orgX
    orgY = simdata$orgY
    rxf = simdata$rxf
    missing_idx2 = which(K_ij==0 & simdata$K_ij>0, arr.ind = TRUE)
    for (s in 1:nrow(missing_idx2)) {
      j = missing_idx2[s,2]
      i = which(simdata$idx_j[[j]] == missing_idx2[s,1])
      X[[j]][simdata$Start[i,j]:simdata$End[i,j],] = NA
      orgX[[j]][simdata$Start[i,j]:simdata$End[i,j]] = NA
      orgY[[j]][simdata$Start[i,j]:simdata$End[i,j]] = NA
      rxf[[j]][simdata$Start[i,j]:simdata$End[i,j]] = NA
    }
    X = lapply(X, na.omit)
    orgX = lapply(orgX, na.omit)
    orgY = lapply(orgY, na.omit)
    rxf = lapply(rxf, na.omit)
    
    return(list(missing_idx = missing_idx2, missing_idx_col = missing_idx,
                X=X, orgX=orgX, orgY=orgY, rxf = rxf, K_ij=K_ij, m_j=m_j, idx_j=idx_j, Start=Start, End=End))
  }
}

# From Evan Poworoznek
# https://github.com/poworoznek/infinitefactor/blob/master/R/jointRot.R

jointRot = function(lambda, eta){
  vari = lapply(lambda, varimax)
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  matches = lapply(loads, msfOUT, piv)
  
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  
  return(list(lambda = lamout, eta = etaout))
}