# ZIPLL 
# modified from https://github.com/AnderWilson/ZIPLL/blob/master/R/simdata.R

simdata <-
  function(nchem, nassay, seed=NULL){
    
    concs <- c(0.301, 0.477, 0.602, 0.845, 1.000, 1.301, 1.602, 2.000)
    nij <- length(concs)
    if(!is.null(seed)) set.seed(seed)
    N <- nchem * nassay
    
    hs <- function(x,theta) scale(theta[1] - (theta[1]-theta[2])/(1+(x/theta[3])^theta[4]), center=TRUE, scale=FALSE)
    
    theta <- cbind(runif(N,0,10), 0, rnorm(N,30,10), runif(N,1,8))
    theta[,1] <- pmax(theta[,1],theta[,2])
    theta[,3] <- pmin(theta[,3],max(concs))
    theta[,3] <- pmax(theta[,3],min(concs))
    
    sigj <- rep(0.1, nassay)
    gamma <- rbinom(N, 1, 0.5)
    
    dat <- NULL
    for(c in 1:nchem){
      for(a in 1:nassay){
        rxf <- gamma[a +nassay*(c-1)]*hs(concs,theta[a +nassay*(c-1),])
        resp <- rxf+rnorm(nij,0,sigj[a])
        dat <- rbind(dat,cbind(c,a,concs,resp))
      }
    }
    
    colnames(dat) <- c("Chemical.Number","Assay.Number","Conc","Resp")
    return(list(dat = dat, truth = list(sigj_sq = sigj^2, gamma_ij = gamma)))
  }
