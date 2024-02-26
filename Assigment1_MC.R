#File: EMM_Assignment1.R
#Project: Assignment 1 for Econometrics Methods (Extremum Estimators)
#Author: Bingxue Li
#Date: 2024-02-03

# Setup -------------------------------------------------------------------
### CLEAN UP MY ENVIRONMENT
cat("\014")
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
rm(list=ls()) 		# Clear workspace
set.seed(123) 		# Set seed for random number generator
# My working directory
setwd('/Users/renee/Desktop/Assignment1') 
library(tictoc)
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian

tic()

rm(list=ls()) 		# Clear workspace

set.seed(1) 		# Set seed for random number generator

n = 300 			# Set the sample size

alpha1 = 0
alpha2 = 1
beta = 1

theta = c(alpha1,alpha2,beta) 	# True parameter value

k = length(theta) 

# Load the log-likelihood, its derivative, and the hessian
source("OP_LL.R")
source("OP_NLS.R")
source("OP_GMMe.R")
source("OP_GMMh.R")
source("OP_Var_GMMh.R")

num = 1000			# Number of Monte Carlo iterations

theta_hat_ML_vec = matrix(0,num,k)
theta_hat_NLS_vec = matrix(0,num,k)
theta_hat_GMMe_vec = matrix(0,num,k)
theta_hat_GMMh_vec = matrix(0,num,k)

inside_N_ML = rep(0,num)
inside_N_NLS = rep(0,num)
inside_N_GMMe = rep(0,num)
inside_N_GMMh = rep(0,num)

J_2_sum = matrix(rep(0,9),3)
Var_hat_NLS_sum = matrix(rep(0,9),3)
Var_hat_GMM_sum = matrix(rep(0,9),3)

epsilon = 0.1
J_2_sum = matrix(rep(0,9),3)
for (it in 1:num) {
  # Data generating process
  x = rnorm(n,0.5,1) 	# regressor
  
  u = rnorm(n,0,1)							 	# error term
  
  y_star = x * beta + u						# latent "utility"
  
  y = rep(0,n)				# observed outcome
  
  is_one = as.logical((y_star>alpha1) * (y_star<=alpha2)) # element-wise logical operations, indicating whether or not in the range alpha1 and 2
  
  y[is_one] = rep(1,sum(is_one))
  
  is_two = (y_star>alpha2)
  
  y[is_two] = rep(2,sum(is_two))
  
  # ML
  # call the OP_LL function to obtrain ml estimators 
  result <- optim(par = theta, OP_LL, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_ML = result$par
  # store ml estimator for simulation it
  theta_hat_ML_vec[it,1:k] = theta_hat_ML
  # count if the estimator is close enough to true para
  if (sqrt(sum((theta_hat_ML - theta)^2)) < epsilon) {
    inside_N_ML[it] = 1
  }
  # k by k matrix for variance
  J_2_sum = J_2_sum + solve(result$hessian)/n

  # NLS
  result <- optim(par = theta, OP_NLS, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_NLS = result$par
  
  theta_hat_NLS_vec[it,1:k] = theta_hat_NLS
  
  if (sqrt(sum((theta_hat_NLS - theta)^2)) < epsilon) {
    inside_N_NLS[it] = 1
  }

  # GMMe
  result <- optim(par = theta, OP_GMMe, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMMe = result$par
  
  theta_hat_GMMe_vec[it,1:k] = theta_hat_GMMe
  
  if (sqrt(sum((theta_hat_GMMe - theta)^2)) < epsilon) {
    inside_N_GMMe[it] = 1
  }
  
  # GMMh
  result <- optim(par = theta, OP_GMMh, y = y, x = x, method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)
  theta_hat_GMMh = result$par
  
  theta_hat_GMMh_vec[it,1:k] = theta_hat_GMMh
  
  if (sqrt(sum((theta_hat_GMMh - theta)^2)) < epsilon) {
    inside_N_GMMh[it] = 1
  }
  Var_hat_GMM_sum = Var_hat_GMM_sum + OP_Var_GMMh(y,x,theta_hat_GMMh)*(1/n)
  
}

# Averages
colMeans(theta_hat_ML_vec)
colMeans(theta_hat_NLS_vec)
colMeans(theta_hat_GMMe_vec)
colMeans(theta_hat_GMMh_vec)

# Variances
var(theta_hat_ML_vec)
var(theta_hat_NLS_vec)
var(theta_hat_GMMe_vec)
var(theta_hat_GMMh_vec)

# (Estimated) Probability that theta_hat lies inside a neigborhood around theta_0
mean(inside_N_ML)
mean(inside_N_NLS)
mean(inside_N_GMMe)
mean(inside_N_GMMh)

# # Plot the histogram of the distribution of the estimator over Monte Carlo iterations
par(mfrow=c(2,2))
hist(theta_hat_ML_vec[1:num,3])
hist(theta_hat_NLS_vec[1:num,3])
hist(theta_hat_GMMe_vec[1:num,3])
hist(theta_hat_GMMh_vec[1:num,3])

J_2_sum/num
Var_hat_GMM_sum/num

toc()


