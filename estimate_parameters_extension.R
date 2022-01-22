#------------------------------------------------------------------------------
#
#       parameter estimation 
#
#		This file combines previous functions for the parameter estimation of the 
#   multi currency affine model based on the Quasi-Kalman-filter. 
#
#
#------------------------------------------------------------------------------


# Based on the Quasi-Maximum-Likelihood approach, the following function estimates 
# the parameters of an affine model according to the definiton in Chapter 8.
# The factor process belongs to the canonical (d,m) class defined by Dai and Singleton. 
# The domestic short rate as well as the foreign short rate process is an affine function
# of the factor process. Moreover the exchange rate is exponentially affine in the 
# factor process. 







estimate_parameters_extension <-
  function(yt, tau_I, tau_A, d, m, intermediate.steps, num_combi, 
           discretization, sample_theta_diag, sample_theta_nondiag, 
           sample_theta_tri, sample_mu, sample_G, sample_c, sample_delta,
           sample_lambda, sample_rrt, sample_e, sample_f, x0, P0){ 

    
    
    Terminal <- length(yt[1,]) / 261
    
    X_best <-
      select_parameters_extension(d = d, m = m, yt = yt, 
                                  intermediate.steps = 261,
                                  tau_I = tau_I, tau_A= tau_A, x0 = 0, P0= 1,
                                  num_combi = num_combi, discretization = discretization, 
                                  sample_theta_diag = sample_theta_diag,
                                  sample_theta_nondiag = sample_theta_nondiag, 
                                  sample_theta_tri = sample_theta_tri,
                                  sample_mu = sample_mu, sample_G = sample_G, 
                                  sample_c= sample_c, sample_delta = sample_delta,
                                  sample_lambda = sample_lambda, sample_rrt = sample_rrt,
                                  sample_e, sample_f)
    
    
    
    lower = define_lower_bounds_extension(d = d, m = m , sample_theta_diag = sample_theta_diag,
                                          sample_theta_nondiag = sample_theta_nondiag, 
                                          sample_theta_tri = sample_theta_tri, sample_mu = sample_mu,
                                          sample_G = sample_G, sample_c = sample_c, sample_delta = sample_delta, 
                                          sample_lambda = sample_lambda, sample_rrt = sample_rrt,
                                          sample_e = sample_e, sample_f = sample_f)
    
    
    upper = define_upper_bounds_extension(d = d, m = m , sample_theta_diag = sample_theta_diag,
                                          sample_theta_nondiag = sample_theta_nondiag, 
                                          sample_theta_tri = sample_theta_tri, sample_mu = sample_mu,
                                          sample_G = sample_G, sample_c = sample_c, sample_delta = sample_delta, 
                                          sample_lambda = sample_lambda, sample_rrt = sample_rrt,
                                          sample_e = sample_e, sample_f = sample_f)
    
    x <- X_best[1, ]
    parameter_fitted_1 <- neldermead(x = x, fn = objective.1_extension, lower = lower,
                                     upper = upper, d = d, m = m, yt = yt,
                                     intermediate.steps = intermediate.steps,
                                     tau_I = tau_I, tau_A = tau_A, xx0 = NA, P0= NA)
    
    #x <- X_best[2,]
    #parameter_fitted_2 <- neldermead(x, objective, d = d, lower = lower,
    #                                  upper = upper, m = m, yt = yt,
    #                                  intermediate.steps =intermediate.steps,
    #                                  tau = tau)$par
    
    #x <- X_best[3,]
    #parameter_fitted_3 <- neldermead(x, objective, d = d, lower = lower,
    #                                 upper = upper, m = m, yt = yt,
    #                                 intermediate.steps =intermediate.steps,
    #                                 tau = tau)$par
    
    
    return(list(parameter_fitted_1, 2, 3))
    
    
} 










