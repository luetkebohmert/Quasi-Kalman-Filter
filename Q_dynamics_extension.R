#------------------------------------------------------------------------------
#
#       Q_dynamics
#
#		This file contains the function which converts the parameters of the 
#   the factor process solving the stochastic differential equation under P
#   in to the parameters according to the riskneutral measure.
#
#------------------------------------------------------------------------------

# In the theoretical part we define the market price of risk proportional to 
# diffusion (see Assumption 1). Now we are able to convert the parameters 
# (theta_P, sigma, G, mu_P, g0) to the corresponding parameters 
# (theta_Q, sigma, G, mu_Q, g0). In addition the function requires the dimension of
# the factor process and the parameter lambda determining the market price of
# risk. The output of the functin is a list of theta_Q and mu_Q. 



Q_dynamics <- function(theta_P, sigma, lambda, G, mu_P, g0, d) {
  # initializing matrices
  H <- matrix(rep(0, (d * d)), d)
  # check for right dimension of the parameters 
  if(d > 1 && dim(theta_P) == dim(sigma) && dim(sigma) == dim(G) && dim(sigma) == c(d, d) 
     && length(lambda) == dim(G)[1] && length(g0) == d) {
    for(i in 1:d) {
      H[i, ] <- lambda[i] * G[, i]
    }
    theta_Q <- theta_P - sigma %*% H
    J <- lambda * g0
    thetamu_Q <- theta_P %*% mu_P + sigma %*% J
    return(list(theta_Q, thetamu_Q))
  } else if (d == 1 && length(theta_P) == length(sigma) && length(sigma) == length(G) 
             && length(sigma) == d && length(lambda) == d && length(g0) == d) {
    theta_Q <- theta_P - sigma * lambda * G 
    J <- lambda * g0
    thetamu_Q <- theta_P * mu_P + sigma * J
    return(list(theta_Q, thetamu_Q))
  } else {
    print("wrong dimension")
  }
}
