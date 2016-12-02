#########################
### Utility functions ###
#########################

#' Calculate the coefficient of variance
#'
#' @export
#' @param x list of samples
#' @return Coefficient of variance
CV <- function(x) {sd(x) / mean(x)}


#' CV to DF in a chisq dist.
#'
#' @export
#' @description Convert coefficient of variance to degrees of freedom in the chisq distribution.
#' @param coeffv Input coefficient of variance
#' @return chisq Degrees of freedom
chiqs_CV2df <- function(coeffv) {2 / coeffv^2}


#' Sample gamma by mean and std
#'
#' @export
#' @description Sample from the gamma distribution using mean and standard deviation as input.
#' @param ns Number of samples to draw
#' @param mu Mean of the sample
#' @param std Standard deviation of the sample
sample_gamma <- function(ns = 100, mu = 1, std = 1) {
  x <- rgamma(ns, shape = mu^2 / std^2, scale = std^2 / mu)
  return(x)
}


#' Sample log norm by mean and std
#'
#' @export
#' @description Sample from the log normal distribution using mean and standard deviation as input:
#' @param ns Number of samples to draw
#' @param mu Mean of the sample
#' @param std Standard deviation of the sample
sample_lognorm <- function(ns = 100, mu = 1, std = 1) {
  x <- rlnorm(ns, meanlog = log(mu) - 0.5 * log(std^2 / mu^2 + 1), sdlog = sqrt(log(std^2 / mu^2 + 1)))
  return(x)
}


#' Sample chisq by the mean
#'
#' @export
#' @description Sample from the chi squared distribution using mean as input:
#' @param ns Number of samples to draw
#' @param mu Mean of the sample
sample_chisq_fixed_mu <- function(ns = 100, mu = 1) {
  x <- rchisq(n = ns, df = mu)
  return(x)
}


#' Sample chisq by the std
#'
#' @export
#' @description Sample from the chi squared distribution using standard deviation as input:
#' @param ns Number of samples to draw
#' @param std Standard deviation of the sample
sample_chisq_fixed_sd <- function(ns = 100, std = 1) {
  x <- rchisq(n = ns, df = std^2 / 2)
  return(x)
}


# Internal function to read a user input,
# and yield the statistics function requested:
yield_mystatistic <- function(x) {
  # Read in the statistics to use:
  if (x == "mean")   mystatistic <- mean
  if (x == "median") mystatistic <- median
  if (x == "sd")     mystatistic <- sd
  if (x == "CV")     mystatistic <- CV
  return(mystatistic)
}



################################
### End of utility functions ###
################################





###############################
### Sampling for simulation ###
###############################

#' Draw one sample for CV simulation
#'
#' @export
#' @description Draw one sample data for CV simulation purpose from different distributions
#' using mean and standard deviation as input:
#' @param n_data1 Number of draws to make for the sample
#' @param mu1 Mean of the sample
#' @param std1 Standard deviation of the sample
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @return A single sample with size and distribution as specified.
#' @examples
#' draw_one_sampleCV_data(n_data1 = 100, mu1 = 1, std1 = 1, simdist = 'norm')
#' draw_one_sampleCV_data(n_data1 = 100, mu1 = 1, std1 = 1, simdist = 'lognorm')
draw_one_sampleCV_data <- function(n_data1, mu1, std1, simdist = 'norm') {
  coeffv1 <- std1 / mu1
  
  onesim <- list()
  if (simdist == 'norm') {
    onesim$X1 <- rnorm(n = n_data1, mean = mu1, sd = std1)
  } else if (simdist == 'lognorm') {
    onesim$X1 <- sample_lognorm(ns = n_data1, mu = mu1, std = std1)
  } else if (simdist == 'gamma') {
    onesim$X1 <- sample_gamma(ns = n_data1, mu = mu1, std = std1)
  } else if (simdist == 'chisq') {
    onesim$X1 <- rchisq(n = n_data1, df = chiqs_CV2df(coeffv = coeffv1))
  } else {
    stop('Could not understand simdist parameter.')
  }
  
  return(onesim)
}


#' Draw two samples for CV simulation
#'
#' @export
#' @description Draw two sample data for CV simulation purpose from different distributions
#' using mean and standard deviation as input:
#' @param n_data1 Number of draws to make for the first sample
#' @param n_data2 Number of draws to make for the second sample
#' @param mu1 Mean for the first sample
#' @param mu2 Mean for the second sample
#' @param std1 Standard deviation for the first sample
#' @param std2 Standard deviation for the second sample
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @return Two samples with size and distribution as specified.
#' @examples
#' draw_two_sampleCV_data(n_data1 = 100, n_data2 = 50, mu1 = 1, mu2 = 1, std1 = 1, std2 = 1, simdist = 'lognorm')
#' draw_two_sampleCV_data(n_data1 = 100, n_data2 = 50, mu1 = 1, mu2 = 1, std1 = 1, std2 = 1, simdist = 'chisq')
draw_two_sampleCV_data <- function(n_data1, n_data2, mu1, mu2, std1, std2, simdist = 'norm') {
  coeffv1 <- std1 / mu1
  coeffv2 <- std2 / mu2
  
  twosim <- list()
  if (simdist == 'norm') {
    twosim$X1 <- rnorm(n = n_data1, mean = mu1, sd = std1)
    twosim$X2 <- rnorm(n = n_data2, mean = mu2, sd = std2)
  } else if (simdist == 'lognorm') {
    twosim$X1 <- sample_lognorm(ns = n_data1, mu = mu1, std = std1)
    twosim$X2 <- sample_lognorm(ns = n_data2, mu = mu2, std = std2)
  } else if (simdist == 'gamma') {
    twosim$X1 <- sample_gamma(ns = n_data1, mu = mu1, std = std1)
    twosim$X2 <- sample_gamma(ns = n_data2, mu = mu2, std = std2)
  } else if (simdist == 'chisq') {
    twosim$X1 <- rchisq(n = n_data1, df = chiqs_CV2df(coeffv = coeffv1))
    twosim$X2 <- rchisq(n = n_data2, df = chiqs_CV2df(coeffv = coeffv2))
  } else {
    stop('Could not understand simdist parameter.')
  }
  
  return(twosim)
}

######################################
### End of sampling for simulation ###
######################################






###################################
### Comparing bootstrap methods ###
###################################

#' Bootstrap on single sample with different methods
#'
#' @description Single sample bootstrap confidence interval on different statistics with different bootstrapping methods.
#' @keywords comparing bootstraps, one sample
#' @import simpleboot
#' @param x input data
#' @param statistic the statistic to run [mean, median, sd, CV]
#' @param nsim Number of bootstrap replicates
#' @param student should the bootstrap be studentized?
#' @param nstudent if bootstrap is studentized how many replicates should this be done with?
#' @param BEL should Bootstrap emperical likelihood confidence interval be calculated? (only available for CV)
#' @param prob probability level for the confidence interval
#' @return list of confidence intervals for the different bootstrapping methods
#' @export
#' @examples
#' bootCI_single(x = rlnorm(n = 100, meanlog = 1, sdlog = 1), statistic = 'CV')
#' bootCI_single(x = rnorm(n = 100, mean = 1, sd = 1), statistic = 'sd')
bootCI_single <- function(x, statistic = "mean", nsim = 1000, student = TRUE,
                          nstudent = 100, BEL = TRUE, prob = 0.95) {
  # Read in the statistics to use:
  mystatistic <- yield_mystatistic(statistic)
  
  # Run the bootstrapping:
  bootres <- one.boot(x, FUN = mystatistic, R = nsim, student = student, M = nstudent)
  # If studentized act on this:
  bootci_list <- list()
  if (student) {
    bootcis <- boot.ci(bootres, conf = prob)
    bootci_list <- list("norm" = bootcis$normal[2:3], "basic" = bootcis$basic[4:5], "student" = bootcis$student[4:5], "perc" = bootcis$percent[4:5], "bca" = bootcis$bca[4:5])
  } else {
    btypes <-  c("norm","basic", "perc", "bca")
    bootcis <- boot.ci(bootres, type = btypes, conf = prob)
    bootci_list <- list("norm" = bootcis$normal[2:3], "basic" = bootcis$basic[4:5], "perc" = bootcis$percent[4:5], "bca" = bootcis$bca[4:5])
  }
  
  if (statistic == "CV" && BEL) {
    bootci_list$BEL <- BEL_CVconfint(X = x, boot_samples = nsim, prob = prob)
  }
  bootci_list['estimate'] <- mystatistic(x)
  return(bootci_list)
}


#' Bootstrap on two samples with different methods
#'
#' @description Compare two samples by different NON-par bootstrapping methods
#' by finding confidence interval on different statistics with different bootstrapping methods.
#' @keywords comparing bootstraps, two samples
#' @import simpleboot
#' @param x input data from sample 1
#' @param y input data from sample 2
#' @param statistic the statistic to run [mean, median, sd, CV]
#' @param nsim Number of bootstrap replicates
#' @param student should the bootstrap be studentized?
#' @param nstudent if bootstrap is studentized how many replicates should this be done with?
#' @param prob probability level for the confidence interval
#' @return List of confidence intervals for the different bootstrapping methods
#' @export
#' @examples
#' bootcompare(x = rlnorm(n = 100, meanlog = 1, sdlog = 1), y = rlnorm(n = 100, meanlog = 1, sdlog = 1), statistic = 'CV')
#' bootcompare(x = rnorm(n = 50, mean = 1, sd = 1), y = rnorm(n = 100, mean = 1, sd = 1), statistic = 'CV')
bootcompare <- function(x, y, statistic = "mean", nsim = 1000, student = TRUE, nstudent = 100, prob = 0.95) {
  # Read in the statistics to use:
  mystatistic <- yield_mystatistic(statistic)
  
  # Run the bootstrapping:
  bootres <- two.boot(x, y, FUN=mystatistic, R=nsim, student = student, M = nstudent)
  # If studentized act on this:
  bootci_list <- list()
  if (student) {
    bootcis <- boot.ci(bootres, conf = prob)
    bootci_list <- list("norm" = bootcis$normal[2:3], "basic" = bootcis$basic[4:5], "student" = bootcis$student[4:5], "perc" = bootcis$percent[4:5], "bca" = bootcis$bca[4:5])
  } else {
    bootcis <- boot.ci(bootres, type = c("norm","basic", "perc", "bca"), conf = prob)
    bootci_list <- list("norm" = bootcis$normal[2:3], "basic" = bootcis$basic[4:5], "perc" = bootcis$percent[4:5], "bca" = bootcis$bca[4:5])
  }
  
  bootci_list['estimate'] <- mystatistic(x) - mystatistic(y)
  return(bootci_list)
}

##########################################
### End of comparing bootstrap methods ###
##########################################

