#------------------------------------------------------------------------------
#
#       parameter restriction 
#
#		This file contains two function which fix the lower and upper boundaries
#   for the parameter search. 
#
#
#------------------------------------------------------------------------------

# Dai and Singleton formulate some parameter restrictions. The following functions 
# constrain the parameter search of the Nelder-Mead-Simplex to these restricitons. 
# We postulate the same restrictions on the foreign short rate model as on the 
# domestic one.
# The input variables for the functions are meaningful boundaries for the parameters 
# and returns a vector. 
# Of course it depends on the context which boundaries are meaningful.
# One needs to be careful with the function define_lower_bounds since the function 
# automatically multiples the input variable by -1.


define_lower_bounds_extension <- function(d, m, sample_theta_diag, sample_theta_nondiag, 
                                          sample_theta_tri, sample_mu, sample_G, sample_c,
                                          sample_delta, sample_lambda, sample_rrt,
                                          sample_e, sample_f) {
  n <- d - m
  if(m != 0) {
    lower <- c(-sample_theta_diag, rep(-sample_theta_nondiag, d - 1))
    if (d > 1){
      for (i in 2:d) {
        lower1 <- c(rep(-sample_theta_nondiag, (i - 1)), -sample_theta_diag, rep(-sample_theta_nondiag, d - i)) 
        lower <- c(lower, lower1)
      }
      without <- (d * m + 1):(d * m + m) 
      for (i in 1:(n - 1)) {
        without1 <- (d * (m + i) + 1) : (d * (m + i) + m) 
        without <- c(without, without1)
      }
      lower <- lower[-without]
    } else {
      lower <- lower
    }
    
    lower <- c(lower, rep(0, m), rep(0, m * n), -sample_c, rep(-sample_delta, m), rep(0, n), 
               rep(-sample_lambda, d), 0.0000000001, -sample_c, rep(-sample_delta, m), rep(0, n), 
               -sample_e, rep(-sample_f, d))
  } else {
    lower <- c(-sample_theta_diag, rep(-sample_theta_diag, d - 1)) 
    if (d > 1) {
      for (i in (d - 1):1) {
        lower1 <- c(-sample_theta_diag, rep(-sample_theta_diag, i - 1)) 
        lower <- c(lower, lower1)
      }
    }
    lower <- c(lower, rep(0, m), rep(0, m * n), -sample_c, rep(-sample_delta, m), rep(0, n), 
               rep(-sample_lambda, d), 0.0000000001, -sample_c, rep(-sample_delta, m), rep(0, n),
               -sample_e, rep(-sample_f, d))
  }
  return(lower)
}



define_upper_bounds_extension <- function(d, m, sample_theta_diag, sample_theta_nondiag, 
                                          sample_theta_tri, sample_mu, sample_G, sample_c,
                                          sample_delta, sample_lambda, sample_rrt,
                                          sample_e, sample_f) {
  n <- d - m
  if (m != 0) {
    upper <- c(sample_theta_diag, rep(0, d - 1)) 
    if (d > 1) {
      for (i in 2:d) {
        if (i <= m) {
          upper1 <- c(rep(0, (i - 1)), sample_theta_diag, rep(0, d - i)) 
          upper <- c(upper, upper1)
        } else {
          upper1 <- c(rep(sample_theta_nondiag, (i - 1)), sample_theta_diag,
                      rep(sample_theta_nondiag, d - i)) 
          upper <- c(upper, upper1) 
        }
      }
      without <- (d * m + 1) : (d * m + m) 
      for (i in 1:(n - 1)) {
        without1 <- (d * (m + i) + 1):(d * (m + i) + m) 
        without <- c(without, without1)
      }
      upper <- upper[-without]
    } else {
      upper <- upper
    }
    
    #sample_c restricted for negativ interest
    upper <- c(upper, rep(sample_mu, m), rep(sample_G, m * n), 0, rep(sample_delta, d),
               rep(sample_lambda, d), sample_rrt, 0, rep(sample_delta, d),  0,
               rep(sample_f, d))
  } else {
    upper <- c(sample_theta_diag, rep(sample_theta_diag, d - 1)) 
    if (d > 1) {
      for (i in (d - 1):1) {
        upper1 <- c(sample_theta_diag, rep(sample_theta_diag, i - 1)) 
        upper <- c(upper, upper1)
      }
    }
    upper <- c(upper, rep(sample_mu, m), rep(sample_G, m * n), 0, rep(sample_delta, d),
               rep(sample_lambda, d), sample_rrt, 0, rep(sample_delta, d),  0,
               rep(sample_f, d))
  }
  return(upper)
}


restrict_parameters_extension<- function(d, m, sample_theta_diag, sample_theta_nondiag, 
                                         sample_theta_tri, sample_mu, sample_G, sample_c,
                                         sample_delta, sample_lambda, sample_rrt,
                                         sample_e, sample_f) {
  lower <- define_lower_bounds_extension(d, m, sample_theta_diag, sample_theta_nondiag, 
                                         sample_theta_tri, sample_mu, sample_G, sample_c,
                                         sample_delta, sample_lambda, sample_rrt,
                                         sample_e, sample_f)
  upper <- define_upper_bounds_extension(d, m, sample_theta_diag, sample_theta_nondiag, 
                                         sample_theta_tri, sample_mu, sample_G, sample_c,
                                         sample_delta, sample_lambda, sample_rrt,
                                         sample_e, sample_f)
  return(list(lower, upper))
}



