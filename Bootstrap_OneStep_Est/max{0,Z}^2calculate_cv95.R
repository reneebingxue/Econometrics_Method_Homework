library(numDeriv) 	# loads the functions grad and hessian which numerically evaluate the gradient and hessian
rm(list=ls()) 		# Clear workspace
#Calculate critical value at 95% for max{0,Z}^2
num <- 1000 # iterations
cv_list <- matrix(0,num,1)
# Derive the CORRECT asymptotic CV
# Simulate n standard normal random variables
for (it in 1:num) {
  Z <- rnorm(n=1e7)
  # Apply the transformation to get Y = (max{0,Z})^2
  Y <- (pmax(0, Z))^2
  # Find the 95th percentile of the simulated distribution
  cv_list[it] <- quantile(Y, .95)
}
cv <- mean(cv_list)
save(cv, file = "cv_value.RData")
