onestep_est_ARCH3 <- function(theta_hat) {
  
  T <- length(X)
  
  sumscore <- score_gradient_ARCH3(theta_hat)
  
  theta_tilde <- theta_hat - solve(J_hessian_ARCH3(theta_hat))%*%t(sumscore)
  
  return(theta_tilde)  
}
