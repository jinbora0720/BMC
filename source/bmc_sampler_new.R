# ----- Heteroscedasticity ----- #
# alpha
alpha.post = function(le, u_ij, mu_alpha = rep(0, 2), sigsq_alpha = rep(100, 2)) {
  u_ij[is.na(u_ij)] = 0 # purely for computational reason
  W = cbind(1, as.numeric(le))
  V = solve(diag(1/sigsq_alpha) + t(W)%*%W)
  return(V%*%(matrix(mu_alpha/sigsq_alpha, ncol = 1) +
                t(W)%*%matrix(u_ij, ncol = 1)) + t(chol(V))%*%rnorm(2))
}

# u_ij
u_ij.post = function(le, alpha, t_ij) {
  m = nrow(le)
  J = ncol(le)
  lea = alpha[1] + alpha[2]*le
  a = ifelse(t_ij==1 | is.na(t_ij), 0, -Inf)
  b = ifelse(t_ij==1 | is.na(t_ij), Inf, 0)
  u = matrix(runif(m*J, pnorm(a-lea), pnorm(b-lea)), m, J)
  u[is.na(t_ij)] = NA
  return(lea+qnorm(u))
}

# pi_t if heter_simpler = TRUE
pi_t.post = function(t_ij, m_j, a_p, b_p) {
  nt = sum(t_ij==1, na.rm=TRUE)
  return(rbeta(1, a_p+nt, b_p+sum(m_j)-nt))
}

# t_ij and delta_ij
td_ij.post = function(e, orgX, sigj_sq, pi_t, v_d, t_ij, d_ij, idx_j, m_j, Start, End) {
  m = nrow(d_ij)
  J = ncol(d_ij)
  
  # propose t_ij
  t_ij.prop = t_ij
  B = rep(0, J)
  for (j in 1:J) {
    # select random numbers of elements to update
    B[j] = sample(m_j[j], 1)
    # select random indices to update
    rpl = sample(idx_j[[j]], B[j], replace = FALSE)
    # replace 0/1 by 1/0
    t_ij.prop[rpl,j] = ifelse(t_ij[rpl,j]==0,1,0)
  }
  
  # propose d_ij
  d_ij.prop = matrix(NA, m, J)
  d_ij.prop[which(t_ij.prop==1, arr.ind=TRUE)] =
    d_ij[which(t_ij.prop==1, arr.ind=TRUE)]+
    rt(sum(t_ij.prop==1, na.rm=TRUE), df = 4)/sqrt(1/v_d+mean(unlist(orgX)^2)) # sqrt(1/v_d+0.5*sum(unlist(orgX)^2))
  d_ij.prop[which(t_ij.prop==0, arr.ind=TRUE)] = 0
  
  log.r = dbinom(t_ij.prop, 1, pi_t, log=TRUE)-dbinom(t_ij, 1, pi_t, log=TRUE)+
    ifelse(t_ij.prop==0, 0, dnorm(d_ij.prop, 0, sqrt(v_d), log=TRUE))-
    ifelse(t_ij==0, 0, dnorm(d_ij, 0, sqrt(v_d), log=TRUE))+        # p(d_ij=0|t_ij=0)=1, log1 = 0
    loglikecpp(e, d_ij.prop, orgX, sigj_sq, idx_j, Start, End, m_j)-
    loglikecpp(e, d_ij, orgX, sigj_sq, idx_j, Start, End, m_j)
    
  u = runif(1)
  d_ij[which(log(u)<log.r, arr.ind = TRUE)] = d_ij.prop[which(log(u)<log.r, arr.ind = TRUE)]
  t_ij[which(log(u)<log.r, arr.ind = TRUE)] = t_ij.prop[which(log(u)<log.r, arr.ind = TRUE)]
  accept = sum(log(u)<log.r, na.rm=TRUE)/(m*J)
  
  return(list(d_ij = d_ij, t_ij = t_ij, accept = accept))
}

# ----- B-spline Functions ----- #
# Sigj, invSigj
Sigj.post = function(a, m_j, R, beta_ij) {
  J = length(m_j)
  p = nrow(R)
  invSigj = lapply(1:J, function(j) {
    rWishart(n=1, df=a+m_j[j], Sigma=solve(R+lapply(1:m_j[j], function(i) {
      tcrossprod(matrix(beta_ij[[j]][,i], ncol=1))
    }) %>% Reduce('+',.))) %>% matrix(., ncol=1)
  }) %>% do.call(what="cbind")
  Sigj = apply(invSigj, 2, function(x) c(ginv(matrix(x, p, p))))
  return(list(invSigj=invSigj, Sigj=Sigj))
}

# ----- Normal Errors ----- #
# sigj_sq
sigj_sq.post = function(eps, nu0, sig0_sq) {
  J = length(eps)
  return(sapply(1:J, function(j) {
    1/rgamma(1,shape=(nu0+length(eps[[j]]))/2,
             rate=(nu0*sig0_sq+sum(eps[[j]]^2))/2)
  }))
}

# ----- Factor Model ----- #
# linearMGSP modified from Evan's linearMGSP
# https://github.com/poworoznek/infinitefactor/blob/master/R/linearMGSP.R

linearMGSP_simpler = function(Y, k, Lambda, delta, tauh, Plam, 
                      df, ad1, bd1, ad2, bd2, s, adapt=TRUE, prop=1, epsilon=1e-3){
  
  p = ncol(Y)
  n = nrow(Y)
  
  b0 = 1
  b1 = 0.0005
  
  if(is.null(k)) k = floor(log(p)*3)
  
  Y[is.na(Y)] = 0
  
  # update parameters
  ps = rep(1,p)                         # Sigma = diagonal residual covariance, fixed at 1 
  eta = eta_lin_simpler(Lambda, ps, k, n, Y)
  Lambda = lam_lin_simpler(eta, Plam, ps, k, p, Y)
  psijh = psi_mg(Lambda, tauh, ps, k, p, df)
  delta = del_mg(Lambda, psijh, tauh, delta, k, p, ad1, bd1, ad2, bd2)
  tauh = cumprod(delta)
  Plam = plm_mg(psijh, tauh)  
  
  if(adapt){
    # ----- make adaptations ----#
    prob = 1/exp(b0 + b1*s)                    # probability of adapting
    uu = runif(1)
    lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
    vec = lind >= prop
    num = sum(vec)                             # number of redundant columns
    
    if(uu < prob) {
      if((s > 20) & (num == 0) & all(lind < 0.995)) {
        k = k + 1
        Lambda = cbind(Lambda, rep(0,p))
        eta = cbind(eta,rnorm(n))
        psijh = cbind(psijh, rgamma(p,df/2,df/2))
        delta[k] = rgamma(1, ad2,bd2)
        tauh = cumprod(delta)
        Plam = t(t(psijh) * tauh)
      } else {
        if (num > 0) {
          k = max(k - num,1)
          Lambda = Lambda[,!vec, drop = F]
          psijh = psijh[,!vec, drop = F]
          eta = eta[,!vec, drop = F]
          delta = delta[!vec]
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        }
      }
    }
  }
  
  Omega = (tcrossprod(Lambda) + diag(1/c(ps))) 
  
  return(list(Lambda=Lambda, eta=eta, Omega=Omega, 
              delta=delta, tauh=tauh, Plam=Plam, k=k))
}

linearMGSP = function(Y, k, Lambda, delta, tauh, Plam, 
                      df, ad1, bd1, ad2, bd2, s, 
                      alpha, u_ji, adapt=TRUE, prop=1, epsilon=1e-3){
  
  p = ncol(Y)
  n = nrow(Y)
  
  b0 = 1
  b1 = 0.0005
  
  if(is.null(k)) k = floor(log(p)*3)
  
  # purely for computational reason
  Y[is.na(Y)] = 0
  u_ji[is.na(u_ji)] = 0
  
  # update parameters
  ps = rep(1,p)                         # Sigma = diagonal residual covariance, fixed at 1 
  eta = eta_lin(Lambda, ps, k, n, Y, alpha, u_ji)
  Lambda = lam_lin(eta, Plam, ps, k, p, Y, alpha, u_ji)
  psijh = psi_mg(Lambda, tauh, ps, k, p, df)
  delta = del_mg(Lambda, psijh, tauh, delta, k, p, ad1, bd1, ad2, bd2)
  tauh = cumprod(delta)
  Plam = plm_mg(psijh, tauh)  
  
  if(adapt){
    # ----- make adaptations ----#
    prob = 1/exp(b0 + b1*s)                    # probability of adapting
    uu = runif(1)
    lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
    vec = lind >= prop
    num = sum(vec)                             # number of redundant columns
    
    if(uu < prob) {
      if((s > 20) & (num == 0) & all(lind < 0.995)) {
        k = k + 1
        Lambda = cbind(Lambda, rep(0,p))
        eta = cbind(eta,rnorm(n))
        psijh = cbind(psijh, rgamma(p,df/2,df/2))
        delta[k] = rgamma(1, ad2,bd2)
        tauh = cumprod(delta)
        Plam = t(t(psijh) * tauh)
      } else {
        if (num > 0) {
          k = max(k - num,1)
          Lambda = Lambda[,!vec, drop = F]
          psijh = psijh[,!vec, drop = F]
          eta = eta[,!vec, drop = F]
          delta = delta[!vec]
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        }
      }
    }
  }
  
  Omega = (tcrossprod(Lambda) + diag(1/c(ps))) 
  
  return(list(Lambda=Lambda, eta=eta, Omega=Omega, 
              delta=delta, tauh=tauh, Plam=Plam, k=k))
}

# z_ij
z_ij.post = function(le, gamma_ij) {
  m = nrow(le)
  J = ncol(le)
  a = ifelse(gamma_ij==1 | is.na(gamma_ij), 0, -Inf)
  b = ifelse(gamma_ij==1 | is.na(gamma_ij), Inf, 0)
  u = matrix(runif(m*J, pnorm(a-le), pnorm(b-le)), m, J)
  u[is.na(gamma_ij)] = NA
  return(le+qnorm(u))
}

# xi 
xi.post = function(z_ij, m_j, mu_xi = 0, sigsq_xi = 100) {
  var = 1/(1/sigsq_xi + sum(m_j))
  return(var*(mu_xi/sigsq_xi + sum(z_ij, na.rm = TRUE)) + sqrt(var)*rnorm(1))
}
