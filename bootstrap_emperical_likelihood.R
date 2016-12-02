##########################
### BEL implementation ###
##########################
# Implementation by Kristian Davidsen

#' BEL confidence interval on CV
#'
#' @export
#' @description Calculate the bootstrapped emperical likelihood based confidence interval for the coefficient of variance.
#' @keywords BEL, EL, bootstrap
#' @param X List of samples
#' @param boot_samples Number of bootstraps to generate the correction factor
#' @param prob Probability level for the confidence interval
#' @return List with two elements, lower and upper interval bounds
#' @examples
#' BEL_CVconfint(X = rnorm(n = 100, mean = 1, sd = 1), boot_samples = 1000, prob = 0.95)
#' BEL_CVconfint(X = rlnorm(n = 100, meanlog = 1, sdlog = 1), boot_samples = 1000, prob = 0.95)
BEL_CVconfint <- function(X = X, boot_samples = 1000, prob = 0.95) {
  k_hat <- CV(X)  # Estimate of the CV
  # Array for bootstrap EL ratio evaluated in k_hat
  ell <- array(data = NA, dim = boot_samples)
  # Run the bootstrap:
  for (i in 1:boot_samples) {
    Xi <- sample(x = X, size = length(X), replace = T)
    ell[i] <- EL_ratio(Xi, k = k_hat)
  }
  # Calculate the bootstrap estimated scale constant:
  c_star <- 1 / mean(ell, na.rm = T)
  # Now find the CI:
  BEL <- EL_CVconfint(X = X, k_hat = k_hat, c = c_star, prob = prob)
  return(BEL)
}


#' EL confidence interval on CV
#'
#' @description Internat function. Calculates the emperical likelihood based confidence interval
#' for the coefficient of variance. Calls the EL_ratio function to get the likelihood ratio.
#' This is primarily a function used by the BEL_CVconfint function.
#' @param X List of samples
#' @param k_hat Estimate of the CV
#' @param c Scaling constant, typically calculated by bootstrapping
#' @param prob Probability for the confidence interval
#' @return List with two elements, lower and upper interval bounds
EL_CVconfint <- function(X, k_hat, c, prob) {
  chisq_cut <- qchisq(p = prob, df = 1) / c
  CI_obj <- function(ki) {EL_ratio(X, ki) - chisq_cut}

  find_BEL <- function (r) {
    # There must be at least two roots:
    if (length(r) < 2) return(NULL)
    # Look fro the true roots:
    for (i in 2:length(r)) {
      if (r[i-1] < k_hat && r[i] > k_hat) {
        BEL <- c(r[i-1], r[i])
        return(BEL)
        break
      }
    }
    return(NULL)
  }

  lb <- k_hat * 0.1
  ub <- k_hat * 5

  BEL <- NULL
  ns <- c(100, 1000, 10000)  # intervals to try
  i <- 1
  # Try more intervals if root finding fails:
  while (is.null(BEL) && i <= length(ns)) {
    r <- uniroot.all(CI_obj, interval = c(lb, ub), n = ns[i])
    BEL <- find_BEL(r)
    i <- i + 1
  }
  return(BEL)
}


#' Emperical likelihood ratio
#'
#' @description Internat function. Finds CV EL ratio e.g. in a bootstrap approach.
#' @keywords emperical likelihood ratio, CV
#' @param Xi Sample to test (these could be bootstrapped)
#' @param k Input coefficient of variance
EL_ratio <- function(Xi, k) {
  # Calculate for the EL ratio for a single k:
  EL_ratio_single <- function(Xi, k) {
    W_hat <- sd(Xi) - k * Xi
    n <- length(Xi)
    # Objective function to find lambda roots:
    lambda_obj <- function(lambda) {
      # lapply over the input lambda to allow a list of lambda's
      # as input to the objective function:
      fl <- lapply(X = lambda, FUN = function(l) {sum(W_hat / (1 + l * W_hat))})
      return(as.numeric(fl))
    }
    # Boundaries suggested in Owen 1988:
    ub <- 1 / (2*n * min(W_hat)) - 1 / min(W_hat)
    lb <- 1 / (2*n * max(W_hat)) - 1 / max(W_hat)

    lambda <- NA
    ns <- c(10, 100, 1000, 10000)  # intervals to try
    i <- 1
    # Try more intervals if root finding fails:
    while ((is.na(lambda) || min(lW) < -1) && i < length(ns)) {
      lambda <- uniroot.all(lambda_obj, interval = c(lb, ub), n = ns[i], nroots = 1)[1]
      lW <- lambda * W_hat
      # Max 100 intervals for searching for a root to make a lW bigger than -1:
      if (min(lW) < -1 && ns[i] >= 100) break
      i <- i + 1
    }
    # Return 'NA' if a lambda root cannot be found:
    if (is.na(lambda)) {
      # print('Cant find any lambda')
      return(NA)
    }
    # Correction if the root was bad:
    if (min(lW) < -1) {
      return(9999)
    } else {
      ell_i <- 2 * sum(log(1 + lW))
      return(ell_i)
    }
  }  # EL_ratio

  # To allow list input:
  ELR <- lapply(k, FUN = function(ki) {EL_ratio_single(Xi, ki)})
  return(as.numeric(ELR))
}


#' Uniroot extented
#'
#' @export
#' @description This function is an extension, to the uniroot in base, that allows finding multiple roots.
#' The function have been adapted from the rootSolve package.
#' @keywords root, root finding, all roots, never using uniroot again
#' @param f Function to find roots for
#' @param interval Interval to search in
#' @param lower Lower bound in interval
#' @param upper Upper bound in interval
#' @param tol Tolerance of the return values
#' @param n Number of intervals to do root finding
#' @param nroots Number of root to find before stopping (-1 means all)
#' @return All roots in the required interval or a limited number specified
#' and starting from the lower interval.
#' @examples
#' Find all roots of the cosine function in the interval between -10 and 10:
#' uniroot.all(f = function(x) {cos(x)}, interval = c(-10, 10), n = 1000)
#'
#' Find the first 10 roots of the positive cosine function:
#' uniroot.all(f = function(x) {cos(x)}, interval = c(0, 1000), n = 1000, nroots = 10)
uniroot.all <- function (f, interval, lower= min(interval),
                         upper= max(interval), tol = .Machine$double.eps^0.2,
                         n = 100, nroots = -1, ... ) {

  # Error checking as in uniroot:
  if (!missing(interval) && length(interval) != 2) {
    stop("'interval' must be a vector of length 2")
  } else if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper) {
    stop("lower < upper  is not fulfilled")
  }

  # Subdivide interval in n subintervals and estimate the function values:
  xseq <- seq(lower, upper, len = n+1)
  mod  <- f(xseq,...)

  # Some function values may already be 0:
  Equi <- xseq[which(mod==0)]

  ss   <- mod[1:n]*mod[2:(n+1)]  # interval where functionvalues change sign
  ii   <- which(ss < 0)

  for (i in ii) {
    Equi <- c(Equi, uniroot(f, lower = xseq[i], upper = xseq[i+1], tol = tol ,...)$root)
    # Return when the required number of roots are found:
    if (length(Equi) == nroots) {
      return(Equi)
    }
  }
  return(Equi)
}

