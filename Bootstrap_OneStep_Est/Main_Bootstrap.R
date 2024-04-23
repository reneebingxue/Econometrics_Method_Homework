#File: EMM_Final_Project.R
#Project: Final_Project: Boundary Issue with ARCH Model
#Author: Bingxue Li, Cong Minh Nguyen
#Date: 2024-03-23
# Setup -------------------------------------------------------------------
### CLEAN UP MY ENVIRONMENT
cat("\014")
library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
library(Rsolnp) # Optimisation with constraint
library(corpcor)
library(tictoc)
rm(list=ls()) 		# Clear workspace
# My working directory
setwd('...')
source("DGP.R")
source("ARCH_LL.R")
source("ARCH_LL_constrained.R")
source("score_gradient.R")
source("DGP.R")
source("J_1.R")
source("onestep_est.R")
source("J_hessian.R")
source("datafunction_boot.R")
source("stat_boot.R")
load("cv_value.RData")

##########################################
# Part.I: Data Generating Process for RCP#
##########################################
n <- 152  # The first two values are auto zero to initiate the series
num <- 1000 # iterations
B= 199
alpha1 <- 0  # ARCH(1) parameter
alpha2 <- 0  # ARCH(2) parameter
omega <- 1  # constant term parameter
theta <- c(alpha1, alpha2, omega)
k=length(theta)                      # Dimension of parameter space
alpha2_null <- 0 # The null hypothesis
reject_t = matrix(0,num,1)
reject_LR_boot_restricted = matrix(0,num,1)

##########################################################################################
# Part.II: Bootstrapping                                                                 #
##########################################################################################

## a) Cavaliere et al. (2017) : RESTRICTED BOOTSTRAP ##

set.seed(5) 		# Set seed for random number generator
for (it in 1:num) {
  # Data generating process
  X <- simulate_arch2(n, theta)$X
  result <- solnp(pars = c(0.5,0.5,0.5), 
                  fun = log_likelihood_arch2, 
                  LB = c(0,0,0), 
                  UB = c(1,1,Inf),
                  control = c(delta= 1e-7, tol= 1e-8, trace=0))
  
  result_constrained <-  solnp(pars =  c(0.5, 0.5), 
                               fun = log_likelihood_arch2_constrained,
                               LB = c(0, 0), 
                               UB = c(1, Inf),
                               control = c(delta= 1e-7, tol= 1e-8, trace=0))
  
  theta_hat_MLE = result$par
  theta_hat_MLE_constrained = c(result_constrained$par[1],alpha2_null,result_constrained$par[2])  
  LR = -2*(result$values[length(result$values)] -  result_constrained$values[length(result_constrained$values)])
  sigma_squared_b <- rep(0, n-2)
  z_b <- rep(0, n-2)
  for (t in 3:n) {
    sigma_squared_b[t-2] <- theta_hat_MLE_constrained[3] + theta_hat_MLE_constrained[1] * X[t - 1]^2 + theta_hat_MLE_constrained[2] * X[t - 2]^2
    z_b[t-2] <- X[t] / sqrt(sigma_squared_b[t-2])
  }
  z_bs = (z_b - mean(z_b))/sqrt(mean((z_b-mean(z_b))^2))   # recentered and standarsdized for E(z_b) = 0, V(z_b) = 1
  LR_boot_restricted = rep(0,B)
  for (b in 1:B) {
    z_boot = sample(z_bs, size = n, replace = TRUE)
    X <- simulate_arch2_boot(n, z_boot, theta_hat_MLE_constrained)
    result_b <- solnp(pars = c(0.5,0.5,0.5), 
                      fun = log_likelihood_arch2, 
                      LB = c(0,0,0), 
                      UB = c(1,1,Inf),
                      control = c(delta= 1e-7, tol= 1e-8, trace=0))
    
    result_constrained_b <-  solnp(pars =  c(0.5, 0.5), 
                                   fun = log_likelihood_arch2_constrained,
                                   LB = c(0, 0), 
                                   UB = c(1, Inf),
                                   control = c(delta= 1e-7, tol= 1e-8, trace=0))
    LR_boot = -2*(result_b$values[length(result_b$values)] - result_constrained_b$values[length(result_constrained_b$values)])
    LR_boot_restricted[b] = LR_boot
  }
  LR_boot_restricted <- sort(LR_boot_restricted)
  cv_b_restricted <- LR_boot_restricted[floor(B*0.975)]
  
  if (is.nan(LR) == 0) {
    if (LR > cv_b_restricted) {
      reject_LR_boot_restricted[it] = 1
    }
  }
  print(it)
}
mean(reject_LR_boot_restricted)
  
## b) Wald-test and t-test (block length = { T^{1/3}, T^{1/4}, T^{1/5}}) bootstrapping based on the one-step estimator ##
set.seed(5) 		# Set seed for random number generator
for (it in 1:num) {
  # Data generating process
  X <- simulate_arch2(n, theta)$X
  result <- solnp(pars = c(0.5,0.5,0.5), 
                  fun = log_likelihood_arch2, 
                  LB = c(0,0,0), 
                  UB = c(1,1,Inf),
                  control = c(delta= 1e-9, tol= 1e-7, trace=0))
  one_step_MLE <- onestep_est(result$par)
  temp = one_step_MLE[2] - alpha2_null
  t = sqrt(n-2)*(temp/sqrt(pseudoinverse(J_1(one_step_MLE))[2,2]))
  X_estimation <- X
  t_boot = rep(0,B)
  for (b in 1:B) {
    X <- stat.boot(X_estimation, (n-2)^(1/3))
    # Restricted model in each boot
    result_b <- solnp(pars = c(0.5,0.5,0.5), 
                                fun = log_likelihood_arch2, 
                                LB = c(0,0,0), 
                                UB = c(1,1,Inf),
                                control = c(delta= 1e-9, tol= 1e-7, trace=0))
    one_step_MLE_b <- onestep_est(result_b$pars) 
    temp_b = one_step_MLE_b[2] - one_step_MLE[2]
    t_boot[b] = sqrt(n-2)*(temp_b/sqrt(pseudoinverse(J_1(one_step_MLE_b))[2,2]))
  }
  t_boot <- sort(t_boot)
  cv_b_t <- t_boot[floor(B*0.95)]
  if (is.nan(t) == 0) {
    if (t > cv_b_t) {
      reject_t[it] = 1
    }
  }
  print(it)
}
mean(reject_t)

