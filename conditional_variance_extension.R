#------------------------------------------------------------------------------
#
#       condtional variance
#
#		This function calculates the conditional variance of a factor process with
#   diagonal drift matrix. (see Corollary 5.1)
#
#------------------------------------------------------------------------------


conditional_variance<- function(theta_diag, mu_diag, d, Delta, a, alpha, x) {
  
  # initializing values 
  v1 <- matrix(rep(0, d * d), nrow = d) 
  v2 <- matrix(rep(0, d * d), nrow = d) 
  v3 <- matrix(rep(0, d * d), nrow= d)
  f <- rep(0, d)
  v <- matrix(rep(0, d * d), nrow= d)
  num1 <- rep(0, d)
  num2 <- rep(0, d)
  denom <- rep(0, d)
  
  epsilon <- 0.00001
  # calculating the conditional variance according to Corollary 1 
  # note that a and b have to be determined first 
  for (i in 1:d) {
    for (j in 1:d) {
      if ( abs(diag(theta_diag)[i] + diag(theta_diag)[j]) < epsilon) {
        v1[i, j] <- Delta 
        v2[i, j] <- a[i, j] + t(alpha[, i, j]) %*% mu_diag
      } else {
        v1[i, j] <- ((-1 + exp((diag(theta_diag)[i] + diag(theta_diag)[j]) * Delta)) 
                     / (diag(theta_diag)[i] + diag(theta_diag)[j]))
        v2[i, j] <- a[i, j] + t(alpha[, i, j]) %*% mu_diag
      }
      for (k in 1:d) {
        if( abs(diag(theta_diag)[i] + diag(theta_diag)[j] - diag(theta_diag)[k]) < epsilon) {
          f[k] <- (Delta * exp((diag(theta_diag)[i] + diag(theta_diag)[j]) * Delta) * 
                     (alpha[k, i, j] * (x[k] - mu_diag[k]))) 
        } else {
          num1[k] <- (exp(diag(theta_diag)[k] * Delta) - 
                        exp((diag(theta_diag)[i] + diag(theta_diag)[j]) * Delta))
          denom[k] <- - diag(theta_diag)[i] - diag(theta_diag)[j] + diag(theta_diag)[k]
          num2[k] <- alpha[k, i, j] * (x[k] - mu_diag[k])
          f[k] <- (num1[k] / denom[k]) * num2[k]
        }
      }
      v3[i,j] <- sum(f)
    }
  }
  cond_variance <- v1 * v2 + v3
  return(cond_variance)
}




