metadata <- function(data, aeid_coff = NULL) {
  # aeid_coff is a data frame with aeid and coff
  
  ae <- distinct(data[c("aeid","aenm")])
  J <- nrow(ae)
  uniq_aeid <- sort(ae$aeid)
  uniq_aenm <- ae[order(ae$aeid),2]
           
  ch <- distinct(data[c("casn","chnm")])
  m <- nrow(ch)
  uniq_casn <- sort(ch$casn)
  uniq_chnm <- ch[order(ch$casn),2]
  
  if (!is.null(aeid_coff)) {
    coff_ij <- matrix(0, m, J)
  }
  
  K_ij <- reshape2::dcast(data, casn ~ aeid, fun.aggregate = length, value.var = "resp")
  rownames(K_ij) <- K_ij[,1]
  K_ij <- as.matrix(K_ij[uniq_casn, as.character(uniq_aeid)])
  
  casn_j <- unname(apply(K_ij, 2, function(x) rownames(K_ij)[which(x>0)]))
  K_ij <- unname(K_ij)
  idx_j <- unname(apply(K_ij, 2, function(x) as.numeric(which(x>0))))
  m_j <- sapply(idx_j, length)
  
  Start = End <- matrix(0, m, J)
  for (j in 1:J) {
    Start[1:m_j[j],j] <- cumsum(c(1,K_ij[,j]))[idx_j[[j]]]
    End[1:m_j[j],j] <- cumsum(K_ij[,j])[idx_j[[j]]]
  }
  
  p <- 5 # B-spline df 
  orgData <- list()
  X = Y <- list() 
  for (j in 1:J) {
    if (!is.null(aeid_coff)) {
      coff <- aeid_coff[which(aeid_coff$aeid == uniq_aeid[j]),2]
    }
    orgData[[j]] <- matrix(NA, max(End[,j]), 2)
    X[[j]] <- matrix(NA, max(End[,j]), p)
    Y[[j]] <- rep(NA, max(End[,j]))
    for (i in 1:m_j[j]) {
      data_ij <- data %>% 
        filter(aeid == uniq_aeid[j], casn == casn_j[[j]][i]) %>% dplyr::select(logc, resp)
      orgData[[j]][Start[i,j]:End[i,j],] <- cbind(data_ij$logc, data_ij$resp)
      X[[j]][Start[i,j]:End[i,j],] <- apply(bs(data_ij$logc, df=p), 2, function(x) scale(x, center=TRUE, scale=FALSE))
      Y[[j]][Start[i,j]:End[i,j]] <- if (length(data_ij$resp)==1) data_ij$resp else 
        if (var(data_ij$resp)==0) scale(data_ij$resp, center=TRUE, scale=FALSE) else  
          scale(data_ij$resp)
      if (!(is.null(aeid_coff))) {
        coff_ij[idx_j[[j]][i],j] <- if (length(data_ij$resp)==1) coff else 
          if (var(data_ij$resp)==0) abs(coff - mean(data_ij$resp)) else  
            abs((coff - mean(data_ij$resp))/sd(data_ij$resp))
      }
    }
  }
  
  if (is.null(aeid_coff)) {
    return(list(uniq_aeid = uniq_aeid, uniq_aenm = uniq_aenm, 
                uniq_casn = uniq_casn, uniq_chnm = uniq_chnm, 
                casn_j = casn_j, idx_j = idx_j, m_j = m_j, 
                K_ij = K_ij, Start = Start, End = End, 
                orgData = orgData, X = X, Y = Y))
  } else {
    return(list(uniq_aeid = uniq_aeid, uniq_aenm = uniq_aenm, 
                uniq_casn = uniq_casn, uniq_chnm = uniq_chnm, 
                casn_j = casn_j, idx_j = idx_j, m_j = m_j, 
                K_ij = K_ij, Start = Start, End = End, 
                orgData = orgData, X = X, Y = Y, 
                coff_ij = coff_ij))
  }
}
