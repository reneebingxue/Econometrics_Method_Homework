onestep_est <- function(theta_hat) {
  
  T <- length(X)
  
  sumscore <- score_gradient(theta_hat)
  
  theta_tilde <- theta_hat - pseudoinverse(J_hessian(theta_hat))%*%t(sumscore)
  
  return(theta_tilde)  
}
