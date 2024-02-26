#File: EMM_Assignment1.R
#Project: Assignment 1 for Econometrics Methods (Extremum Estimators)
#Author: Bingxue Li
#Date: 2024-02-03

# Setup -------------------------------------------------------------------
### CLEAN UP MY ENVIRONMENT
cat("\014")
install.packages("scatterplot3d")
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library("scatterplot3d") # Load required package for 3D plotting
rm(list=ls()) 		# Clear workspace
set.seed(1) 		# Set seed for random number generator
# My working directory
setwd('/Users/renee/Desktop/Assignment1') 

### DATA GENERATING PROCESS
n = 300 			# Set the sample size
alpha1 = 0    # Set the parameter value
alpha2 = 1    # Set the parameter value
beta = 1      # Set the parameter value
theta = c(alpha1,alpha2,beta)	# True parameter value vector
k = length(theta) # Dimension of parameter space
# Data generating process
x = rnorm(n,0.5,1) 	# regressors
u = rnorm(n,0,1)		# error term
y_star = x * beta + u	 # latent "utility"
y = rep(0,n)				# observed outcome
is_one = as.logical((y_star>alpha1) * (y_star<=alpha2))
y[is_one] = rep(1,sum(is_one))
is_two = (y_star>alpha2)
y[is_two] = rep(2,sum(is_two))

# Load the log-likelihood
source("OP_LL.R")
# Compute the Maximum Log-likelihood estimator
# optimization without user-specified gradient (Broyden–Fletcher–Goldfarb–Shanno algorithm)
result_b <- optim(par = theta, OP_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_b$par

# Load the objective function for NLS
source("OP_NLS.R")
# optimization without user-specified gradient
result_d <- optim(par = theta, OP_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_d$par

# Load the objective function for GMM given in (e)
source("OP_GMMe.R")
# optimization without user-specified gradient
result_f <- optim(par = theta, OP_GMMe, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_f$par

# Load the objective function for GMM given in (h)
source("OP_GMMh.R")
# optimization without user-specified gradient
result_i <- optim(par = theta, OP_GMMh, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
result_i$par
