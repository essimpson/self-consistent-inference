#######################################################################################
#######################################################################################
### This file provides an example for using our self-consistent estimation code ###
### We apply the approach to data simulation from an inverted logistic model ###

### Required packages ###
library(evd) # For data simulation

### Load functions for our method and currently available estimators for comparison ###
source('1 - G-estimation.R')
source('2 - additional-estimators.R')

#######################################################################################
#######################################################################################

# Simulate data from an inverted logistic model 
# (first simulate logistic data, then transform to uniform margins and invert the copula):
  set.seed(1)
  gamma <- 0.5
  data <- rmvevd(10000, dep=gamma, model="log", d=2)
  data.U <- apply(data, 2, function(x){rank(x)/(length(x)+1)})
  data.U <- 1-data.U
  data <- -log(1-data.U)
  
# Run our proposed estimation approach:
# Tuning parameters are selected as proposed in the paper
  estimate <- consistentEstimation(data,  quant=0.999, thresh.quant=0.5, nei=100, len=199, knot.num=7, deltas=seq(0.01,1,0.01), omegas=seq(0.01,1,0.01))
  
  # 'data' should be a matrix with two columns corresponding to the two variables - can be on the original scale or pre-transformed to exponential margins
  # 'quant' is the radial quantile to be used for extrapolation
  # 'thresh.quant' is the quantile to be used for the GPD threshold - same for the local and smoothed estimators
  # 'nei' is the number of observations used to construct each R_w set for the local estimator
  # 'len' is the number of w* values to be used
  # 'knot.num' is the number of interior knots to be used in the spline functions (equivalent to kappa in the paper)
  # 'deltas' is the values for which estimates of tau_i(delta) should be produced
  # 'omegas' is the values for which estimates of lambda(omega) should be produced
  
# The 'consistentEstimation' function returns several results, which we now explain in turn

# First is the estimated boundary set (i.e., with radial values at the specified quantile) before scaling onto [0,1]^2:
  plot(estimate$boundary, type="l", lwd=3, asp=1, xlab=expression(X[1]), ylab=expression(X[2]))

# Then we have our final estimate of the set G:
  plot(estimate$boundary.scaled, type="l", lwd=3, asp=1, xlab=expression(X[1]), ylab=expression(X[2]))
  abline(v=c(0,1), lty=3, col="grey80")
  abline(h=c(0,1), lty=3, col="grey80")
  
# The following give our estimates of the extremal dependence parameters of interest:
  estimate$eta                  # coefficient from Ledford and Tawn (1996)
  estimate$lambda               # coefficients from Wadsworth and Tawn (2013) evaluated at the specified omega values
  estimate$tau1; estimate$tau2  # coefficients from Simpson et al. (2020) evaluated at the specified delta values
  estimate$HT1; estimate$HT2    # coefficients from Heffernan and Tawn (2004) conditioning on X1 and X2

# Finally, the function returns some parameter values selected during the modelling procedure:
  estimate$deg  # the chosen degree of the spline functions
  estimate$Wrange  # the range of angles over which we can provide an estimate of G
  
#######################################################################################
#######################################################################################  
### We also provide code for some existing estimation techniques of individual dependence parameters, for comparison:
  
  eta.hill.exp(data=data, quant=0.95)  # The Hill estimator for eta
  eta.peng(data=data, k=500)           # The Peng estimator for eta
  eta.draisma(data=data, k=500)        # The Draisma et al. estimator for eta
  
  lambda.hill.exp(data=data, quant=0.95, omegas=seq(0.01,1,0.01))  # The Hill estimator for lambda(omega)
    
  tau.hill.exp(data, quant=0.85, deltas=seq(0.01,1,0.01))            # The Hill estimator for tau_1(deltas)
  tau.hill.exp(data[,c(2,1)], quant=0.85, deltas=seq(0.01,1,0.01))   # The Hill estimator for tau_2(deltas)

  optim(rep(0.5,4), fn=nll.HT, data=data)$par[1:2]            # MLEs for (alpha1, beta1)
  optim(rep(0.5,4), fn=nll.HT, data=data[,c(2,1)])$par[1:2]   # MLEs for (alpha2, beta2)    
  
