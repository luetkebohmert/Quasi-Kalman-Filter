#------------------------------------------------------------------------------
#
#       objective function
#
#		This file provides the input function for the Nelder-Mead optimization. 
#
#------------------------------------------------------------------------------

# We want to use the Nelder-Mead algorithm to optimize the Log-Likelihood.
# The Nelder-Mead function is prefered because one is able to pre-define intervals for 
# the parameter search. Some of the bounds arise from the classification of Dai  and 
# Singleton. Other restrictions are chosen by our own considerations in order to receive
# meaningful parameter values. 
# The Nelder-Mead function minimizes a function according to a single numeric vector argument
# returning a numeric scalar. In our case this function is the objective function. 
# Therefore the input parameters are a vector x, the dimensions d and m, the data yt, 
# the size of the discretization given by intermediate.steps, tau and the starting values 
# of the Kalman filter x0 and P0. The objective function will be minimized according to 
# the vector x, the other input variables are fixed. Inside the function the vector x defines
# the required parameters for the quasi_kalman_filter function.  
# The output of the objective function is the negative value of the Log-Likelihood. 


objective_extension <- function(x, d, m, yt, intermediate.steps, tau_I, tau_A, x0, P0) {
  theta_P <- define_theta_P_extension(x, d, m)
  mu_P <- define_mu_P_extension(x, d, m)
  sigma <- define_sigma_extension(x, d, m)
  g0 <- define_g0_extension(x, d, m)
  G <- define_G_extension(x, d, m)
  delta_I <- define_delta_I_extension(x, d, m)
  c_I <- define_c_I_extension(x, d, m)
  lambda <- define_lambda_extension(x, d, m)
  rrt <- define_rrt_extension(x, d, m)
  delta_A <- define_delta_A_extension(x, d, m)
  c_A <- define_c_A_extension(x, d, m)
  e <- define_e_extension(x, d, m)
  f <- define_f_extension(x, d, m)


  LL <- c(quasi_kalman_filter_extension(
    theta_P, sigma, lambda, rrt, G, mu_P, 
    g0, delta_I, c_I, delta_A, c_A, e, f, 
    intermediate.steps, x0, P0, d, m, yt,
    tau_I, tau_A)[[1]])
  return(-LL)
}

# This auxilliary function serves to tackle a small programming issue. The pre-implemented
# Nelder-Mead function in R confuses the "initial guess" with the starting values of the
# Kalman iteration process. We therefore propose to use that extra function.
objective.1_extension <- function(x, d, m, yt, intermediate.steps, tau_I, tau_A, xx0, P0) {
  objective_extension(x = x, d = d, m = m, yt = yt, 
                      intermediate.steps = intermediate.steps, tau_I = tau_I,
                      tau_A = tau_A, x0 = xx0, P0 = P0)
}
