#######################################################################
### Function to calculate the Hill estimator in exponential margins ###
#######################################################################
# (Also used in our proposed estimation procedure)
eta.hill.exp <- function(data, quant=0.95){
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data.E <- -log(1-data.U)
  
  Q <- apply(data.E,1,min)
  
  u <- quantile(Q,quant)  
  
  Q.ext <- Q[Q>u]
  eta.est <- min(1, mean(Q.ext-u)) 
  
  return(eta.est)
}

####################################
### Function to calculate S_n(k) ###
####################################
# (Used in Peng and Draisma et al. estimators)
Sn <- function(data, k){
  
  n <- nrow(data)
  
  X1.ord <- sort(data[,1],decreasing=T)
  X2.ord <- sort(data[,2],decreasing=T)
  
  v1 <- X1.ord[k]
  v2 <- X2.ord[k]
  
  Sn <- sum(data[,1] >= v1 & data[,2] >= v2)  
  
  return(Sn)  
}

##########################################################
### Function to calculate the estimator of Peng (1999) ###
##########################################################
eta.peng <- function(data, k){
  Sn.k <- Sn(data,k)
  Sn.2k <- Sn(data,(2*k))
  
  eta <- log(2)/(log(Sn.2k/Sn.k))
  
  return(min(eta,1))
}

####################################################################
### Function to calculate the estimator of Draisma et al. (2004) ###
####################################################################
eta.draisma <- function(data, k){
  Sn.j <- sapply(c(1:k), function(i){Sn(data,i)})
  Sn.k <- Sn(data,k)
  
  eta <- sum(Sn.j)/(k*Sn.k - sum(Sn.j))
  
  return(min(eta,1))
}

#########################################################################################
### Function to calculate the Hill estimator for lambda(omega) in exponential margins ###
#########################################################################################
lambda.hill.exp <- function(data, quant=0.95, omegas){
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data.E <- -log(1-data.U)
  
  lambda.est <- NULL
  
  for(iter in 1:length(omegas)){
    omega <- omegas[iter]
    Mw <- apply(cbind(data.E[,1]/omega,data.E[,2]/(1-omega)),1,min)
    
    u <- quantile(Mw,quant)  
    
    Mw.ext <- Mw[Mw>u]
    lambda.est[iter] <- min(1, 1/mean(Mw.ext-u))
  }
      
  return(lambda.est)
}

###############################################################################
### Function to calculate the Hill estimator for tau in exponential margins ###
###############################################################################
tau.hill.exp <- function(data, quant=0.85, deltas){
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data.E <- -log(1-data.U)
  
  tau.est <- NULL
  
  for(iter in 1:length(deltas)){
    delta <- deltas[iter]
    Md.ind <- which(data.E[,2]<=delta*data.E[,1])
    Md <- data.E[Md.ind,1]
    
    u <- quantile(Md,quant)  
    
    Md.ext <- Md[Md>u]
    tau.est[iter] <- min(1, mean(Md.ext-u)) 
  }
  
  return(tau.est)
}


#######################################################################
### Negative log likelihood function for conditional extremes model ###
#######################################################################
nll.HT <- function(par, data){
  
  alph <- par[1]
  beta <- par[2]
  mu <- par[3]
  sig <- par[4]
  
  if(alph < 0 | alph > 1 | beta < 0 | beta > 1 | sig < 0){return(Inf)}
  
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data <- -log(1-data.U)
  
  thresh <- quantile(data[,1],0.95)
  x <- data[data[,1]>thresh,1]
  y <- data[data[,1]>thresh,2]
  
  means <- alph*x + mu*x^beta
  sds <- sig*x^beta
  
  nll <- -sum(dnorm(y, means, sds, log=T))
  
  return(nll)
}

