#------------------------------------------------------------------------------
#
#       Quasi-Kalman-Filter 
#
#		This file is an implementation of the Quasi-Kalman-Filter for 
#   the affine model of two currency areas following Chapter 8.
#
#------------------------------------------------------------------------------


# Similar to the Quasi-Kalman-Filter of the affine short rate model of one country, 
# this function implements the Quasi-Kalman-Filter according to the equations
# (KF1) to (KF8) of the theoretical part. Keep in mind that the aim is to estimate
# the parameters of the extended model. Therefore this function fully 
# depends on the parameters of the factor process under the measure P 
# (theta_P, sigma, lambda, rrt, G, mu_P, g0), the parameters of the short rate 
# delta_d, delta_f, c_d and c_f,  and the parameters of the exchange rates e and f, as
# well as the proposed starting values of the recursion (x0, P0) (note that in this
# implementation the starting values of the filter are defined within the function.
# Hence it does not matter to choose correctly a pripori), the dimension of the
# factor process according to the classification of Dai and Singleton (d, m),
# the observation data (yt) including the domestic yields to 
# the times to maturity tau_d and the foreign yields to the times to maturity
# tau_f and the exchange rate. The function returns a list of components: 
# list(LL, s, p, xtt, xt, vt). The value of the Log-Likelihood, which could 
# be maximized for parameter estimation purposes. The two numbers: p and s. 
# Due to rounding errors Ft could be non-singular and non-symmetric. Theoretically 
# this can not be the case. But if this phenomenon arises, we igonore the iteration 
# step and raise the value of p respectively s by one. In addition the filtered and 
# predicted states as well as the forecasting error are returned.

# Besides the straightforward camputation of the filter we concentrate on the 
# classifictaion of Dai and Singelton. Of course we restrict the intervals of the 
# parameters according to the classification. But some assumptions like strictly 
# positiveness and real eigenvalues of the drift matrix and the assumption that 
# the sum of the i-th row of theta_P multiplied with the i-th entry of mu_P is 
# greater zero for 0 < i < (m+1) are difficult to implement as a priori restrictions.
# Hence we check those restrictions in this function and if at least one is not fulfilled, 
# we assign a high negative value to the Likelihood-function, to prevent us from 
# choosing that parameter constellation.

quasi_kalman_filter_extension <- function(theta_P, sigma, lambda, rrt, G, mu_P, 
                                          g0, delta_I, c_I, delta_A, c_A, e, f, 
                                          intermediate.steps, x0, P0, d, m, yt,
                                          tau_I, tau_A) {
  
  # counts how often Ft is not symmetric
  s <- 0
  # counts how often Ft is not singular 
  p <- 0
  # theoretically both counts are equal to 0
  
  # auxilliary variable serves to tackle a small programming issue
  sfactor <- 10^10
  
  
  # dimension of the factor process
  d <- d 
  # depth of the observation series
  q <- (length(tau_I) + length(tau_A) + 1)
  # length of the observation series
  l <- length(yt[1, ])
  
  # initializing a and b for conditional variance 
  a <- matrix(rep(0, d * d), nrow = d)
  alpha <- array(rep(0, d * d * d), dim = c(d, d, d))
  
  
  
  # initializing matrices and arrays for the kalman filter: 
  # matrix containing the intercept of the transition equation 
  ff <- matrix(rep(0, d), nrow = d)
  # matrix containing the factor of the transition equation
  FF <- matrix(rep(0, d * d), nrow = d)
  # matrix containing the factor of the measurement equation
  HH <- matrix(rep(0, q * d), nrow = q)
  # matrix containing the intercept of the measurement equation
  hh <- matrix(rep(0, q, nrow = q))
  # matrix containing the variance of the disturbances of the measurement equation
  RR <- matrix(rep(0, q * q), nrow = q)
  # array containing the variance of the innovations of the transition equation 
  QQt<- array(rep(0, d * d * l), dim= c(d, d, l))
  
  # initialzinig matrices and arrays for the results of the kalman filter 
  # d x (l + 1) - matrix containing the predicted states 
  xt <- matrix(rep(0, d * (l + 1)), nrow = d)  
  # d x l - matrix containing the filtered state variables, 
  xtt <- matrix(rep(0, d * l), nrow = d)
  # d x d x (l + 1) - array containing the mean squared error matrizes of xt
  Pt <- array(rep(0) , dim=c(d, d, l + 1))
  # d x d x l  - array containing mean squared error matrizes of xtt
  Ptt <- array(rep(0) , dim=c(d, d, l))
  # q x l - matrix containing the prediction errors
  vt <- matrix(rep(0, q * l), nrow = q)
  # q x q x l  - array containing the variances of vt
  Ft <- array(rep(0), dim=c(q, q, l))
  # d x q x l - array containing the “Kalman gain”
  Kt <- array(rep(0), dim=c(d, q, l))
  # vecotor containing the log-likelihood values
  LL <- rep(0, l) 
  
  # size of discetization 
  Delta <- 1 / intermediate.steps
  
  
  # check if the restrictions of the canonical model of 
  # Dai an Singleton are fulfilled
  if (if (m > 1) {
    any(rowSums(t(t(theta_P) * mu_P)[1:m, 1:m]) <= 0)
  } else if(m == 1) {
    t(t(theta_P) * mu_P)[1:m, 1:m] <= 0
  } else {
    1 == 0
  }) {
    LL <- -100000
  } else {
    if (diagonalize_drift(theta_P)[[3]] == TRUE && Im(eigen(theta_P)$values) == rep(0,d) &&
        any(Re(eigen(Q_dynamics(theta_P, sigma, lambda, G, mu_P, g0, d)[[1]])$values) < 0) == FALSE
        && Im(eigen(Q_dynamics(theta_P, sigma, lambda, G, mu_P, g0, d)[[1]])$values) == rep(0,d)) {
      
      # parameter after diagonalizing theta
      # using the function diagonalize_drift
      L <- diagonalize_drift(theta_P)[[1]]
      if (is.complex(L)) {
        LL <- -100000
      } else {
        theta_diag <- diagonalize_drift(theta_P)[[2]]
        sigma_diag <- L
        G_diag <- matrix(rep(0, d * d), d)
        for (i in 1:d) {
          G_diag[, i] <- solve(t(L)) %*% G[, i]
        }
        mu_diag <- L %*% mu_P
        g0_diag <- g0
        delta_diag <- solve(t(L)) %*% delta_I 
        c_diag <- c_I
        
        
        a <- L %*% diag(g0_diag,d) %*% t(L)
        
        for (i in 1:d) {
          for (j in 1:d) {
            for (k in 1:d) {
              alpha[k, i, j] <- sum(G_diag[k, ] * L[i, ] * L[j, ]) 
            }
          }
        }
        
        
        # all matrices according to the transition equation
        # using the results after diagonalization
        ff <- matrix(- solve(L) %*% diag(exp( - diag(theta_diag) * Delta), d) %*% mu_diag + mu_P, d)
        FF <- solve(L) %*% diag(exp( - diag(theta_diag) * Delta), d) %*% L
        
        
        
        # all matrices according to the observation equation
        # using the Riccati_solver and the Riccati_solver_modification function 
        
        HH[1:length(tau_I), ] <- - (as.matrix(Riccati_solver(
          theta_P, sigma, lambda, G, mu_P, g0, delta_I, c_I, d, tau_I)
          [1 : length(tau_I), 3 : (d + 2)], ncol = d) 
          / matrix(rep(tau_I, each = d), ncol=d, byrow = T))
        
        HH[(length(tau_I) + 1):(length(tau_I) +length(tau_A)), ] <-
          ((matrix(rep(f, length(tau_A)), ncol = d, byrow = TRUE) - 
              as.matrix(Riccati_solver_modification(theta_P, sigma, lambda, G, mu_P, g0,
                                                    delta_A, c_A, d, tau_A, c(0, f))
                        [1 : length(tau_A), 3 : (d + 2)], ncol = d)) 
           / matrix(rep(tau_A, each = d), ncol=d, byrow = T))
                                                    
        HH[(length(tau_I) + length(tau_A) + 1), ] <- f
        
        
        
        hh[1:length(tau_I), ] <- - as.matrix(Riccati_solver(
          theta_P, sigma, lambda, G, mu_P, g0, delta_I, c_I, d, tau_I)
          [1 : length(tau_I), 2], ncol = 1) / tau_I
        
        hh[(length(tau_I) + 1):(length(tau_I) +length(tau_A)), ] <- 
          - as.matrix(Riccati_solver_modification(theta_P, sigma, lambda, G, mu_P,
                                                  g0, delta_A, c_A, d, tau_A, c(0, f))
                      [1 : length(tau_A), 2], ncol = 1) / tau_A
        
        hh[(length(tau_I) + length(tau_A) + 1), ] <- as.numeric(e) 
        
        
        RR <- diag(rrt, q)
        
        # defining the starting values of the filter 
        xt[, 1] <- mu_P
        
        epsilon <- 0.00000000001
        for (i in 1:d) {
          for (j in 1:d) {
            if ( abs(diag(theta_diag)[i] + diag(theta_diag)[j]) < epsilon) {
              Pt[i, j, 1] <- 1
            } else {
              Pt[i, j, 1] <- ((a[i, j] + t(alpha[, i, j]) %*% mu_diag)/ 
                                (diag(theta_diag)[i] + diag(theta_diag)[j])) 
            }
          }
        }
        
        Pt[, , 1] <- solve(L) %*% Pt[, , 1] %*% solve(t(L))
        
        # kalman iteration:
        for (i in 1:l) {
          vt[, i] = yt[, i] - hh - HH %*% xt[, i]
          Ft[, , i] = HH %*% Pt[, , i] %*% t(HH) + RR
          # if (is.singular.matrix(Ft[,,i])) {
          # print(list(theta_P, sigma, lambda, rrt, G, mu_P, g0))}
          if (is.non.singular.matrix (sfactor^2 * Ft[, , i], tol = 1e-30)) {
            Kt[, , i] = Pt[, , i] %*% t(HH) %*% solve(sfactor * Ft[, , i]) * sfactor
            xtt[, i] = xt[, i] + Kt[, , i] %*% vt[, i]
            Ptt[, , i] = Pt[, , i] -  Kt[, , i] %*% HH %*% Pt[, , i]
            # HHt depends on the value of the filtered state variables att
            QQt[, , i] <- array( solve(L) %*% 
                                   conditional_variance(theta_diag, mu_diag, d, Delta, 
                                                        a, alpha, L %*% xtt[, i])
                                 %*% t(solve(L)), dim = c(d, d, 1))
            xt[, i + 1] = ff + FF %*% xtt[, i]
            Pt[, , i + 1] = FF %*% Ptt[, , i] %*% t(FF) + QQt[, , i]
            mu <- as.vector(hh + HH %*% xt[, i])
            sigmaLL <- Ft[, , i]
            # due to rounding errors there could arise a problem with
            # the symmetry of sigma, which is theoretically symmetric
            if (isSymmetric(sigmaLL, tol = sqrt(.Machine$double.eps))) {
              LL[i] <- log(dmvnorm(yt[, i], mu, sigmaLL))
            } else {
              LL[i] <- -Inf
              s <- s + 1 
            }
          } else {
            p <- p + 1 
            LL[i] <- -Inf 
          }
        }
        LL <- sum(LL[LL != -Inf])
      }
    } else {
      LL <- -100000 
    }
  }
  if (LL != -Inf && is.nan(LL) == FALSE && is.na(LL) == FALSE) {
    return(list(LL, s, p, xtt, xt, vt))
  } else {
    # this should not be the case
    LL <- -100000 
    return(list(LL, s, p, xtt, xt, vt))
  }
}



