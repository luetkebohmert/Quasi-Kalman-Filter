#------------------------------------------------------------------------------
#
#       diagonalize drift
#
#		This file contains a function which diagonalizes the drift matrix. 
#
#------------------------------------------------------------------------------

# In the theoretical part we made the assumption of the existence of a non-singular
# matrix L, that diagonalizes the matrix theta_P. A sufficient condition for that 
# assumption is that the eigenvalues of theta_P are distinct. Furthermore the
# positivity of the eigenvalues is required for staionarity and we restrict ourself
# to real eigenvalues. Those conditions of well-posedness are checked in this
# function and if the restrictions are fulfilled the matrix L is computed and a list
# of the matrix L, the diagonalized drift and a logical value is returned.


diagonalize_drift <- function(theta_P) {
  
  # check if it is possible to diagonalize (the condition is too strong)
  if (anyDuplicated(eigen(theta_P)$values) == 0 &&
     Im(eigen(theta_P)$values) == rep(0,d) && any(Re(eigen(theta_P)$values) < 0) == FALSE) {
    L_inverse <- eigen(theta_P)$vectors
    L_inverse <- cbind(L_inverse[, c(dim(theta_P)[1]:1)])
    L <- solve(L_inverse)
    theta_P_diagonalized <- round(L %*% theta_P %*% L_inverse, digits = 3)
    return(list(L, theta_P_diagonalized, TRUE))
  } else {
    list(FALSE, FALSE, FALSE)
  }
}



