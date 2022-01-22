#------------------------------------------------------------------------------
#
#       define parameters 
#
#		This file contains several small auxiliary functions which convert a single 
#   numeric vector x to the parameters of the extended affine model according 
#   to Chapter 8. 
#
#
#------------------------------------------------------------------------------


# For optimization purposes it is important that the objective function depends on 
# a single numeric vector. Therefore we define several functions converting a single 
# numeric vector to the required parameters. The vector x contains the values in the 
# following order:
#     c(theta_P(column-wise), mu_P,  G (column-wise), c_I, delta_I, lambda, rrt,
#         c_A, delta_A, e, f)

define_theta_P_extension <- function(x, d, m) {
  if (d > m && m != 0) {
    n <- d - m
    m1 <- matrix(x[1:(d * m)], nrow = d)
    m2 <- matrix(x[(d * m + 1):((d * m) + n ^ 2)], nrow = n)
    m3 <- matrix(rep(0, m * n), nrow = m)
    m4 <- rbind(m3, m2)
    theta_P <- cbind(m1, m4)
    return(theta_P)
  } else if (d > m && m == 0) {
    theta_P <- matrix(rep(0, d * d), nrow = d)
    theta_P[, 1] <- x[1:d]
    if(d >= 3) {
      theta_P[, 2] <- c(0, x[(d + 1):(2 * d - 1)])
      for (i in 3:d) {
        theta_P[, i] <- c(rep(0, i - 1),
                          x[((i - 1) * d + sum((- 1):(- i + 2)) + 1):
                              (i * d + sum(- 1:(1 - i)))])
      }
    } else if(d == 2) {
      theta_P[, 2] <- c(0, x[(d + 1):(2 * d - 1)])
    } else {
      theta_P <- theta_P[, 1]
    }
    return(as.matrix(theta_P, nrow = d))
  } else {
    theta_P <- matrix(x[1:(d * d)], nrow = d)
    return(theta_P)
  }
}

define_mu_P_extension <- function(x, d, m) {
  n <- d - m
  mu_P<- rep(0, d)
  if(m != 0 ) {
    mu_P[1:m] <- x[(d * d - m * n + 1):(d * d - m * n + m)]
  } else {
    mu_P <- mu_P 
  }
  return(mu_P)
}

define_sigma_extension <- function(x, d, m) {
  n <- d - m
  sigma <- diag(d)
  return(sigma)
}

define_g0_extension <- function(x, d, m) {
  g0 <- rep(1, d)
  if(m != 0 ) {
    g0[1:m] <- 0
  } else {
    g0 <- g0
  }
  return(g0)
}

define_G_extension <- function(x, d, m) {
  if(m != 0 && d >= 2) {
    n <- d - m
    m1 <- diag(m)
    m2 <- matrix(x[((d * d) - (m * n) + m + 1):((d * d) + m)], nrow = m)
    m3 <- cbind(m1, m2)
    m4 <- matrix(rep(0, d * n), nrow = n)
    G <- rbind(m3, m4)
    return(G)
  } else if (m != 0 && d == 1) {
    G <- diag(1)
    return(G)
  } else {
    G <- matrix(rep(0, (d * d)), nrow = d)
    return(G)
  }
}

define_c_I_extension <- function(x, d, m) {
  if(m != 0) {
    c_I <- x[(d * d) + m + 1]
    return(c_I)
  } else {
    c_I <- x[(sum(1:d) + 1)]
    return(c_I)
  }
}

define_delta_I_extension <- function(x, d, m) {
  if( m != 0 ) {
    delta_I <- x[((d * d)  + m + 2):((d * d) + m + d + 1)]
    return(delta_I)  
  } else {
    delta_I <- x[(sum(1:d) + 2):(sum(1:d) + d + 1)]
    return(delta_I)  
  } 
}

define_lambda_extension <- function(x, d, m) {
  if ( m != 0 ) {
    lambda <- x[((d * d) + m + d + 2):((d * d) + m + d + d + 1)]
    return(lambda)
  } else {
    lambda <- x[(sum(1:d) + d + 2):(sum(1:d) + 1 + d + d)]
    return(lambda)
  }
}

define_rrt_extension <- function(x, d, m) {
  if ( m != 0 ) {
    rrt <- x[(d * d) + m + d + d + 2]
    return(rrt)
  } else {
    rrt <- x[(sum(1:d) + 1 + d + d + 1)]
    return(rrt)
  }
}

define_c_A_extension <- function(x, d, m) {
  if(m != 0) {
    c_A <- x[(d * d) + m + d + d + 3]
    return(c_A)
  } else {
    c_A <- x[(sum(1:d) + 1 + d + d + 2)]
    return(c_A)
  }
}

define_delta_A_extension <- function(x, d, m) {
  if( m != 0 ) {
    delta_A <- x[((d * d) + m + d + d + 4):((d * d) + m + d + d + 3 + d)]
    return(delta_A)  
  } else {
    delta_A <- x[(sum(1:d) + 1 + d + d + 3):(sum(1:d) + 1 + d + d + 2 + d)]
    return(delta_A)  
  } 
}

define_e_extension <- function(x, d, m) {
  if(m != 0) {
    e <- x[((d * d) + m + d + d + 3 + d) + 1]
    return(e)
  } else {
    e <- x[(sum(1:d) + 1 + d + d + 2 + d + 1)]
    return(e)
  }
}

define_f_extension <- function(x, d, m) {
  if( m != 0 ) {
    f <- x[((d * d) + m + d + d + 3 + d + 2):((d * d) + m + d + d + 3 + d + 1 + d)]
    return(f)  
  } else {
    f <- x[(sum(1:d) + 1 + d + d + 3 + d + 1):(sum(1:d) + 1 + d + d + 3 + d +  d)]
    return(f)  
  } 
}


