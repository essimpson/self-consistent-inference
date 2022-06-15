library(evgam)
library(mgcv)

########################################################################################
### Negative log likelihood function for conditional extremes model with fixed alpha ###
########################################################################################
nll.HT.alph <- function(par, data, alph, quantile=0.95){
  beta <- par[1]
  mu <- par[2]
  sig <- par[3]
  
  if(beta<0 | beta>1 | sig<0){return(10^15)}
  
  # Rank transform to exponential margins
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data <- -log(1-data.U)

  # Extract observations with conditioning variable above a high threshold 
  thresh <- quantile(data[,1],quantile)
  x <- data[data[,1]>thresh,1]
  y <- data[data[,1]>thresh,2]
  
  # Calculate the means and sds of the residual process for each observation
  means <- alph*x + mu*x^beta
  sds <- sig*x^beta
  
  # Calculate the negative log likelihood
  nll <- -sum(dnorm(y, means, sds, log=T))
  
  return(nll)
}

###############################################################
### Function to carry out boundary and parameter estimation ###
###############################################################
# 'data' should be a matrix with two columns corresponding to the two variables - can be on the original scale
# 'quant' is the radial quantile to be used for extrapolation
# 'thresh.quant' is the quantile to be used for the GPD threshold - same for the local and smoothed estimators
# 'nei' is the number of observations used to construct each R_w set for the local estimator
# 'len' is the number of w* values to be used
# 'knot.num' is the number of interior knots to be used in the spline functions (equivalent to kappa in the paper)
# 'deltas' is the values for which estimates of tau_i(delta) should be produced
# 'omegas' is the values for which estimates of lambda(omega) should be produced


consistentEstimation <- function(data, quant=0.999, thresh.quant=0.5, nei=100, len=199, knot.num=7, deltas=seq(0.1,1,0.1), omegas=seq(0.1,1,0.1)){
  
  # Rank transform to exponential margins:
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data.E <- -log(1-data.U)
  
  # Calculate pseudo-polar coordinates:
  R.E <- apply(data.E,1,sum) 
  W.E <- data.E[,1]/R.E
  
  # Take log(R) for threshold calculations:
  logR.E <- log(R.E)
  
  # Vector of angles W* at which to estimate the boundary:
  Wstar <- unique(sort(c(quantile(W.E, seq(0,1,length.out=(len-1))),0.5)))
  
  # For each value of W*: calculate the largest angular distance (from W*) in each neighbourhood; 
  # calculate estimates of the radial quantiles from a standard GPD fit with empirical quantile threshold:
  thresh.move <- NULL
  R.W.quant.ind <- NULL
  
  for(sec.num in 1:length(Wstar)){
    sec.W <- Wstar[sec.num]
    eps.W <- sort(abs(W.E - sec.W))[nei]
    nei.W <- which(abs(W.E - sec.W) <= eps.W)
    
    thresh.move[sec.num] <- quantile(R.E[nei.W],thresh.quant)
    fit <- ismev::gpd.fit(R.E[nei.W], threshold=thresh.move[sec.num], npy=length(nei.W))
    R.W.quant.ind[sec.num] <- evd::qgpd((quant-thresh.quant)/(1-thresh.quant), loc=thresh.move[sec.num], scale=fit$mle[1], shape=fit$mle[2])
  }
  
  # Specify placement of knots for splines:
  knots <- c(min(W.E), (min(W.E) + (max(W.E)-min(W.E))*c(1:(knot.num-2))/(knot.num-1)), max(W.E))
  knots[(knot.num+1)/2] <- 0.5
  knots1 <- list(W.E = c(0,knots,1)) # Linear spline knots
  knots2 <- list(W.E = c(0,0,knots,1,1)) # Quadratic spline knots
  knots3 <- list(W.E = c(0,0,0,knots,1,1,1)) # Cubic spline knots
  
  # Use asymmetric Laplace for quantile regression of threshold:
  # Fit on log(R) then back-transform:
  fmla_ald <- paste('logR.E ~ s(W.E, k=',knot.num,', bs="bs", m=1)')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.E,W.E)), family="ald", ald.args=list(tau=0.5), knots=knots1)
  thresh.smooth1 <- exp(predict(m_ald, newdata=list(W.E=W.E))$location) # Estimates at all observed W.E values (for threshold-exceedance calculations)
  thresh.smooth1a <- exp(predict(m_ald, newdata=list(W.E=Wstar))$location) # Estimates at all observed W* values only
  
  fmla_ald <- paste('logR.E ~ s(W.E, k=',knot.num+1,', bs="bs", m=2)')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.E,W.E)), family="ald", ald.args=list(tau=0.5), knots=knots2)
  thresh.smooth2 <- exp(predict(m_ald, newdata=list(W.E=W.E))$location)
  thresh.smooth2a <- exp(predict(m_ald, newdata=list(W.E=Wstar))$location)
  
  fmla_ald <- paste('logR.E ~ s(W.E, k=',knot.num+2,', bs="bs", m=3)')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.E,W.E)), family="ald", ald.args=list(tau=0.5), knots=knots3)
  thresh.smooth3 <- exp(predict(m_ald, newdata=list(W.E=W.E))$location)
  thresh.smooth3a <- exp(predict(m_ald, newdata=list(W.E=Wstar))$location)
  
  # Extract threshold exceedances for subsequent GPD fit:
  exceeds1 <- apply(cbind(R.E,thresh.smooth1),1,function(x){x[1]>x[2]})
  exceeds2 <- apply(cbind(R.E,thresh.smooth2),1,function(x){x[1]>x[2]})
  exceeds3 <- apply(cbind(R.E,thresh.smooth3),1,function(x){x[1]>x[2]})
  
  R.matrix.pred <- as.data.frame(cbind(Wstar, rep(NA,length(Wstar)))) 
  names(R.matrix.pred) <- c("W.E","R.E") 
  
  # Fit the GPD-GAM (with a spline in the scale parameter and constant shape parameter)
  # to the radial components, using the W* values as a covariate in the spline.
  # Calculate the required radial quantile for each W*.
  # Do this for splines of degree 1,2,3 and select the best by comparing to the individual GPD estimates:
  R.matrix <- as.data.frame(cbind(W.E[exceeds1==T], (R.E[exceeds1==T]-thresh.smooth1[exceeds1==T])))
  names(R.matrix) <- c("W.E","R.E") 
  spl=paste('R.E ~ s(W.E, k=',knot.num,', bs="bs", m=1)')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots1)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant1 <- thresh.smooth1a + gpd_pred$`q`
  errors <- sum(abs(R.W.quant1 - R.W.quant.ind))
  
  R.matrix <- as.data.frame(cbind(W.E[exceeds2==T], (R.E[exceeds2==T]-thresh.smooth2[exceeds2==T])))
  names(R.matrix) <- c("W.E","R.E") 
  spl=paste('R.E ~ s(W.E, k=',knot.num+1,', bs="bs", m=2)')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots2)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant2 <- thresh.smooth2a + gpd_pred$`q`
  errors[2] <- sum(abs(R.W.quant2 - R.W.quant.ind))
  
  R.matrix <- as.data.frame(cbind(W.E[exceeds3==T], (R.E[exceeds3==T]-thresh.smooth3[exceeds3==T])))
  names(R.matrix) <- c("W.E","R.E") 
  spl=paste('R.E ~ s(W.E, k=',knot.num+2,', bs="bs", m=3)')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots3)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant3 <- thresh.smooth3a + gpd_pred$`q`
  errors[3] <- sum(abs(R.W.quant3 - R.W.quant.ind))

  # Select R-quantile estimates using spline degree with smallest error  
  degree <- which.min(errors)
  R.W.quant <- R.W.quant1*(degree==1) + R.W.quant2*(degree==2) + R.W.quant3*(degree==3)

  # Transform back to original coordinates to obtain boundary estimates:
  # Only keep estimates **within** range of observed W values:
  boundary <- cbind(R.W.quant*Wstar, R.W.quant*(1-Wstar))
  boundary <- boundary[-c(1,nrow(boundary)),]
  
  # Scale the boundary estimates using Hill estimate of eta then truncating/stretching:
  etaHat <- eta.hill.exp(data)
  
  xstar <- max(apply(boundary,1,min))
  
  bound2 <- etaHat*boundary/xstar
  bound2[bound2[,1]>1,1] <- 1
  bound2[bound2[,2]>1,2] <- 1
  
  bound3 <- bound2
  bound3 <- apply(bound3, 2, function(x){x/max(x)})
  
  bound.scaled <- bound3
  
  # Calculate estimate of eta using result from Nolde (2014):
  eta <- max(apply(bound.scaled,1,min))
  
  # Calculate estimates of alpha using N+W (2020) result:
  alpha1 <- max(bound.scaled[bound.scaled[,1]==1,2])
  alpha2 <- max(bound.scaled[bound.scaled[,2]==1,1])
  # Estimate betas using maximum likelihood estimation with fixed alpha:
  beta1 <- optim(rep(0.5,3), fn=nll.HT.alph, data=data, alph=alpha1)$par[1]
  beta2 <- optim(rep(0.5,3), fn=nll.HT.alph, data=data[,c(2,1)], alph=alpha2)$par[1]
  
  # Calculate estimates of tau1(delta) and tau2(delta) at specified values of delta using N+W (2020) result:
  tau1 <- NULL
  tau2 <- NULL
  
  for(iter in 1:length(deltas)){
    delta <- deltas[iter]
    tau1.inds <- which(bound.scaled[,2]<=delta*bound.scaled[,1])
    tau1.candidates <- bound.scaled[tau1.inds,1]
    if(length(tau1.candidates)==0){
      tau1[iter] <- NA
    }else{
      tau1[iter] <- max(tau1.candidates)
    }
    
    tau2.inds <- which(bound.scaled[,1]<=delta*bound.scaled[,2])
    tau2.candidates <- bound.scaled[tau2.inds,2]
    if(length(tau2.candidates)==0){
      tau2[iter] <- NA
    }else{
      tau2[iter] <- max(tau2.candidates)
    }
  }    
  
  # Calculate estimates of lambda(omega) at specified values of omega using N+W (2020) result:
  lambda <- NULL    
  for(iter in 1:length(omegas)){
    w <- omegas[iter]
    lambda[iter] <- 1/max(apply(bound.scaled,1,function(x){min((x[1]/w),(x[2]/(1-w)))}))
  }   
  
  return(list(boundary=boundary, boundary.scaled=bound.scaled, eta=eta, HT1=c(alpha1,beta1), HT2=c(alpha2,beta2), tau1=tau1, tau2=tau2, lambda=lambda, deg=degree, Wrange=range(W.E)))
  
}
