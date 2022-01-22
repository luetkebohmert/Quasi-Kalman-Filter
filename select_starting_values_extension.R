#------------------------------------------------------------------------------
#
#       starting values of the parameter combinations
#
#		This file contains serveral functions for selecting the starting points
#   of parameter combinations.
#
#------------------------------------------------------------------------------


# We randomly sample points of parameter combinations and evaluate the 
# Likelihood-function at each point. The points with the best values 
# are selected as the starting points for the parameter search using 
# the Nelder-Mead method. 





# This function calculates the number of parameters, which have to be 
# estimated. The input are the parameters d and m concerning the 
# classificaiton of the factor process of Dai and Singleton. 
number_of_parameters_extension <- function(d, m) {
  if (m != 0) {
    n <- d - m
    num <- d * d - m * n + m + m * n + 1 + d + d + 1 + d + 1 + d + 1 
  } else {
    num <- sum(1:d) + 1 + d + d + 1 + d + 1 + d + 1 
  }
  return(num)
}


# This function randomly generates parameter combinations in specified 
# intervals. Some of the bounds are fixed by the classification. 
# The function returns a matrix which contains a parameter combination in 
# each row. The input is the dimension d and m und the bounds for the 
# parameter intervals. We choose the same parameter restrictions for 
# domestic and the foreign short rate. 
generate_parameters_extension <- function(d, m, num_combi, discretization, 
                                          sample_theta_diag, sample_theta_nondiag,
                                          sample_theta_tri, sample_mu, sample_G, 
                                          sample_c, sample_delta, sample_lambda,
                                          sample_rrt, sample_e, sample_f) {
  
  num_para <- number_of_parameters_extension(d, m)
  num_para1 <- number_of_parameters_extension(d, m) - (2 * d) - 2 
  # initializing
  # X is the matrix which contains a parameter combination in each row 
  X <- matrix(rep(0, num_para * num_combi), ncol = num_para)
  
  # auxilliary matrix
  x <- matrix(rep(0, num_para * discretization), ncol = num_para)
  
  n <- d - m
  
  if (m != 0) {
    # random generator for theta 
    # it would be good to make a difference between the off diagonal elements 
    # and the elements on the the diagonal 
    # conditioned on non negativ values on the diagonal
    if (d == 1) {
      x[, 1] <- seq(0, sample_theta_diag, length.out = discretization)
    } else {
      i <- 1
      for (j in 1:(m * d)) {
        if(i != (d + 2) && i != 1) {
          x[, j] <- seq(-sample_theta_nondiag, 0, length.out = discretization)
          i <- i + 1
        } else {
          x[, j] <- seq(0, sample_theta_diag, length.out = discretization)
          i <- 2
        }
      }
    }
    
    i <- 1 
    if((m * d + 1) < (d * d - m * n)) {
      for (j in (m * d + 1):(d * d - m * n)) {
        if(i != (n + 2) && i != 1) {
          x[, j] <- seq(-sample_theta_nondiag, sample_theta_nondiag, length.out = discretization)
          i <- i + 1
        } else {
          x[, j] <- seq(0, sample_theta_diag, length.out = discretization)
          i <- 2
        } 
      }
    }
    
    # random generator for mu
    for (j in (d * d - m * n + 1):(d * d - m * n + m)) {
      x[, j] <- seq(0, sample_mu, length.out = discretization)
    }
    
    # random generator for G
    if ((d * d - m * n + m + 1) < (d * d + m )) {
      for (j in (d * d - m * n + m + 1):(d * d + m )) {
        x[, j] <- seq(0, sample_G, length.out = discretization)
      }
    }
    
    # random generator for c_I
    x[, d * d + m + 1] <- seq(-sample_c, 0, length.out = discretization)
    
    # random generator for delta_I
    if (d == 1) {
      x[, 4] <- seq(0, sample_delta, length.out = discretization)
    } else {
      for (j in (d * d + m + 2):(d * d + m + 2 + (m - 1))) {
        x[, j] <- seq(-sample_delta, sample_delta, length.out = discretization)
      }
      if ((d * d + m + m + 2) < (d * d + m + 2 + d)) {
        for (j in (d * d + m + m + 2):(d * d + m + 1 + d)) {
          x[, j] <- seq(0, sample_delta, length.out = discretization)
        } 
      }
    }
    
    # random generator for lambda
    for (j in (num_para1 - d):(num_para1 - 1)) {
      x[, j] <- seq(-sample_lambda, sample_lambda, length.out = discretization)
    }
    
    # random generator for rrt 
    x[, num_para1] <- seq(0.00000001, sample_rrt, length.out = discretization)
    
    # random generator for c_A
    x[, num_para1 + 1] <- seq(-sample_c, 0, length.out = discretization)
    
    # random generator for delta_A
    if (d == 1) {
      x[, num_para1 + 2 ] <- seq(0, sample_delta, length.out = discretization)
    } else {
      for (j in (num_para1 + 2):(num_para1 + 2 + (m - 1))) {
        x[, j] <- seq(-sample_delta, sample_delta, length.out = discretization)
      }
      if ((num_para1 + m + 2) < (num_para1 + 2 + d)) {
        for (j in (num_para1 + m + 2):(num_para1 + 1 + d)) {
          x[, j] <- seq(0, sample_delta, length.out = discretization)
        } 
      }
    }
    
    # random generator for e
    x[, num_para - d] <- seq(-sample_e, 0, length.out = discretization)
    
    # random generator for f
    for (j in (num_para + 1 - d):num_para) {
      x[, j] <- seq(-sample_f, sample_f, length.out = discretization)
    }
    
  } else {
    # random generator for theta 
    for (j in 1:sum(1:d)) {
      x[, j] <-  seq(0, sample_theta_diag, length.out = discretization)
    }
    
    # random generator for c_I
    x[, (sum(1:d) + 1)] <-  seq(-sample_c, 0, length.out = discretization)
    
    # random generator for delta_I
    for (j in (sum(1:d) + 2):(sum(1:d) + 1 + d)) {
      x[, j] <-  seq(0, sample_delta, length.out = discretization)
    }
    
    # random generator for lambda
    for (j in (num_para1 - d):(num_para1 - 1)){
      x[, j] <- seq(-sample_lambda, sample_lambda, length.out = discretization)
    }
    
    # random generator for rrt 
    x[, num_para1] <- seq(0.00000001, sample_rrt, length.out = discretization)
  }
  
  # random generator for c_I
  x[, num_para1 + 1] <-  seq(-sample_c, 0, length.out = discretization)
  
  # random generator for delta_I
  for (j in (num_para1 + 2):(num_para1 + 1 + d)) {
    x[, j] <-  seq(0, sample_delta, length.out = discretization)
  }
  
  # random generator for e
  x[, num_para1 + d + 2] <-  seq(-sample_e, 0, length.out = discretization)
  
  # random generator for f
  for (j in (num_para + 1 - d):num_para) {
    x[, j] <-  seq(-sample_f, sample_f, length.out = discretization)
  }
  
  for (i in 1:num_para) {
    X[, i] <- sample(x[, i], num_combi, replace = TRUE)
  }
  
  # check the remaining parameter restrictions
  j <- 0
  k <- 0
  
  for (i in 1:num_combi) {
    x <- X[i, ]
    theta_P <- define_theta_P_extension(x, d, m)
    mu_P <- define_mu_P_extension(x, d, m)
    sigma <- define_sigma_extension(x, d, m)
    g0 <- define_g0_extension(x, d, m)
    G <- define_G_extension(x, d, m)
    delta <- define_delta_I_extension(x, d, m)
    c <- define_c_I_extension(x, d, m)
    lambda <- define_lambda_extension(x, d, m)
    if (diagonalize_drift(theta_P)[[3]] == TRUE && Im(eigen(theta_P)$values) == rep(0,d) &&
        any(Re(eigen(Q_dynamics(theta_P, sigma, lambda, G, mu_P, g0, d)[[1]])$values) < 0) == FALSE
        && Im(eigen(Q_dynamics(theta_P, sigma, lambda, G, mu_P, g0, d)[[1]])$values) == rep(0,d)){
      if (if(m > 1) {
        any(rowSums(t(t(theta_P) * mu_P)[1:m, 1:m]) > 0)
      } else if (m == 1) {
        t(t(theta_P) * mu_P)[1:m, 1:m] > 0
      } else {
        help1 = 1
      } ){
        if (is.complex(diagonalize_drift(theta_P)[[1]]) == FALSE){
          k <- c(k, i)
        } else{
          j <- c(j, i)
        }
      } else {
        j <- c(j, i)
      }
    } else 
    {
      j <- c(j, i)
    }
  }
  
  j <- j[-1]
  
  if(length(j) == num_combi){
    print("The random generator produced invalid combindations. Try again.")
  } else if (length(j) == 0) {
    X <- X
  } else {
    X <- X[-j, ]
  }
  return(X)
}


# This function return the 5 best parameter combinations concerning the value 
# of the Log-Likelihood-function. The input of the function are the dimensions 
# d and m, as well as the observation data and tau, x0 and P0 and all the bounds 
# for the parameter intervals. Because inside the function the parameter-generator
# provides the parameter combinations for which the Log-Likelihood-function will
# be evaluated. 
select_parameters_extension <- function(d, m, yt, intermediate.steps, tau_I, tau_A, x0, P0,
                                        num_combi, discretization,  sample_theta_diag,
                                        sample_theta_nondiag, sample_theta_tri, sample_mu,
                                        sample_G, sample_c, sample_delta, sample_lambda,
                                        sample_rrt, sample_e, sample_f) {
  
  num_para <- number_of_parameters_extension(d, m)
  
  X <- generate_parameters_extension(d, m, num_combi, discretization, 
                                     sample_theta_diag, sample_theta_nondiag,
                                     sample_theta_tri, sample_mu, sample_G, 
                                     sample_c, sample_delta, sample_lambda,
                                     sample_rrt, sample_e, sample_f)
  # initializing the vector for saving the results of the Log-Likelihood-function
  L <- rep(0, dim(X)[1])
  
  # combuting the Likelihood-Function for each parameter combination and evaluating 
  # which one is the best 
  pb <- txtProgressBar(min = 0, max = dim(X)[1], style = 3)
  for(i in 1:dim(X)[1]){
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
    x <- X[i, ]
    L[i] <- objective_extension(x, d, m, yt, intermediate.steps, tau_I, tau_A, x0, P0)
  }
  L[is.na(L)] <- 0
  LL <- sort(L, decreasing = FALSE)
  # matrix containing the  five "best" parameter combination 

  X_best <- X[c(which(L == LL[1]), which(L == LL[2]), which(L == LL[3]), 
                which(L == LL[4]), which(L == LL[5]),  which(L == LL[6]),  
                which(L == LL[7]),  which(L == LL[8]),  which(L == LL[9]),
                which(L == LL[10])),]
  return(X_best)
}




