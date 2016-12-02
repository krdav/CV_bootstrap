#################################################################
### Functions to find coverage probabilities of CI's for CV's ###
#################################################################
# Implementation by Kristian Davidsen

#' Bootstrap coverage probability on one sample with different methods
#'
#' @description Coverage probabilities for one sample CI of CV calculated with different NON-par bootstrapping methods.
#' @keywords comparing bootstraps, one sample, coverage probabilities
#' @import foreach
#' @import doMC
#' @param n_data1 Number of data points in the sample to test
#' @param mu1 Mean of the simulated distribution
#' @param std1 Standard deviation of the simulated distribution
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @param prob Probability level for the confidence interval
#' @param np Number of processors to use. Default set to 4 because this is a fairly normal number these days.
#' @param boot_samples Number of bootstraps to generate
#' @param reps Replicates to make the coverage probability
#' @param flavors Choice of bootstrapping flavors [norm, basic, perc, bca, student, BEL]. Default all.
#' @return List of coverage probabilities for the different bootstrapping methods
#' @export
#' @examples
#' coverage_prob_single_sample(n_data1 = 100, mu1 = 1, std1 = 1,
#'                             simdist = 'norm', boot_samples = 100, reps = 100)
#'
#' This will take some time (approx. 5 min):
#' coverage_prob_single_sample(n_data1 = 100, mu1 = 1, std1 = 1, simdist = 'norm',
#'                             boot_samples = 100, reps = 1000, flavors = c('BEL'))
coverage_prob_single_sample <- function(n_data1, mu1, std1, simdist = 'norm', np = 4,
                                        boot_samples = 200, reps = 1000, prob = 0.95,
                                        flavors = c("norm", "basic", "perc", "bca", "student", "BEL")) {
  registerDoMC(np)
  CIarray_of_lists <- foreach (i = 1:reps) %dopar% {
    coeffv <- std1 / mu1

    onesim <- draw_one_sampleCV_data(n_data1 = n_data1, mu1 = mu1, std1 = std1, simdist = simdist)
    X <- onesim[['X1']]

    allCI <- bootCI_single(x = X, statistic = 'CV', nsim = boot_samples, prob = prob, student = T, nstudent = 100, BEL = T)

    res <- list()
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      CI <- allCI[[flavor]]
      if (length(CI) != 2) {
        res[flavor] <- 'NA'
      } else if (coeffv < CI[1]) {
        res[flavor] <- 'L'
      } else if (coeffv > CI[2]) {
        res[flavor] <- 'H'
      } else {
        res[flavor] <- 'C'
      }
    }
    return(res)
  }

  CI_list_of_arrays <- list()
  for (i in 1:length(flavors)) {
    flavor <- flavors[i]
    CI_list_of_arrays[[flavor]] <- array(data = NA, dim = length(CIarray_of_lists))
  }

  for (j in 1:length(CIarray_of_lists)) {
    CI_list <- CIarray_of_lists[[j]]
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      CIres <- CI_list[[flavor]]
      CI_list_of_arrays[[flavor]][j] <- CIres
    }
  }

  allres <- list()
  for (i in 1:length(flavors)) {
    flavor <- flavors[i]
    CIarray <- CI_list_of_arrays[[flavor]]
    CP <- sum(CIarray == 'C') / sum(CIarray != 'NA')
    if (sum(CIarray == 'L') > 0) {
      skew <- sum(CIarray == 'H') / sum(CIarray == 'L')
    } else {
      skew <- NA
    }
    res <- list()
    res[['NAs']] <- sum(CIarray == 'NA')
    res[['CP']] <- CP
    res[['skew']] <- skew
    res[['CIarray']] <- CIarray
    allres[[flavor]] <- res
  }
  return(allres)
}


#' Bootstrap coverage probability on one sample with different methods and different distributions
#'
#' @description Coverage probabilities for one sample CI of CV calculated with different NON-par bootstrapping methods,
#' and on multiple distributions. Also a list of different CV's can be provided. Results are printed in a matrix format.
#' Resulsts include the coverage probabilities but also the skewness of this and how many NA's returned in a run.
#' @keywords comparing bootstraps, one sample, coverage probabilities
#' @import foreach
#' @import doMC
#' @param n_data1 Number of data points in the sample to test
#' @param mu1 Mean of the simulated distribution
#' @param std1 Standard deviation of the simulated distribution
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @param prob Probability level for the confidence interval
#' @param np Number of processors to use. Default set to 4 because this is a fairly normal number these days.
#' @param boot_samples Number of bootstraps to generate
#' @param reps Replicates to make the coverage probability
#' @param flavors Choice of bootstrapping flavors [norm, basic, perc, bca, student, BEL]. Default all.
#' @return None. All results are printed to comma separated text files.
#' @export
#' @examples
#' This one takes a bit of time to run through (approx. 15 min),
#' so either grab a coffee, be patient or run with fewer reps:
#' multiple_coverage_prob_single_sample(mu1 = 1, std1 = 1,
#'                                      boot_samples = 200,
#'                                      reps = 300, prob = 0.95,
#'                                      n_data_list = list('30' = 30),
#'                                      sim_distributions = c('lognorm', 'gamma'))
#'
#' Alternatively to above this will run faster (approx. 2 min):
#' multiple_coverage_prob_single_sample(mu1 = 1, std1 = 1,
#'                                      boot_samples = 100,
#'                                      reps = 200, prob = 0.95,
#'                                      n_data_list = list('30' = 30),
#'                                      sim_distributions = c('gamma'),
#'                                      flavors = c("perc"))
multiple_coverage_prob_single_sample <- function(mu1, std1, boot_samples = 200, reps = 1000,
                                                 prob = 0.95, n_data_list = list('16' = 16, '30' = 30),
                                                 flavors = c("norm", "basic", "perc", "bca", "student", 'BEL'),
                                                 sim_distributions = c('norm', 'lognorm', 'gamma', 'chisq'),
                                                 k_list = list('0.2' = 0.2, '0.5' = 0.5), write_dir = '.') {

  # Set the write directory:
  setwd(write_dir)

  # Make the list of matrices to store the results:
  all_CP_mat <- list()
  all_skew_mat <- list()
  all_NAs_mat <- list()
  for (d in 1:length(sim_distributions)) {
    simdist <- sim_distributions[[d]]
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      all_CP_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
      all_skew_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
      all_NAs_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
    }
  }


  # Make the run:
  for (d in 1:length(sim_distributions)) { # all the distributions to test
    simdist <- sim_distributions[[d]]
    for (i in 1:length(n_data_list)) {  # rows, n_data_list
      for (j in 1:length(k_list)) {     # cols, k_list
        res <- coverage_prob_single_sample(n_data1 = n_data_list[[i]][1], mu1 = mu1, std1 = k_list[[j]][1],
                                           simdist = sim_distributions[[d]], boot_samples = boot_samples,
                                           reps = reps, prob = prob, flavors = flavors)
        # Asign the results for different bootstrapping flavors:
        for (f in 1:length(flavors)) {
          flavor <- flavors[f]
          all_CP_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['CP']]
          all_skew_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['skew']]
          all_NAs_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['NAs']]
        }
      }
    }
  }

  # Print the results:
  for (d in 1:length(sim_distributions)) {
    simdist <- sim_distributions[[d]]
    # Coverage probabilities:
    fnam_CP <- paste('CP_mat_', simdist, '.csv', sep = '')
    write(x = 'Coverage probabilities:\n', file = fnam_CP)
    # Skew:
    fnam_skew <- paste('skew_mat_', simdist, '.csv', sep = '')
    write(x = 'Skew matrix:\n', file = fnam_skew)
    # NA's (usefull diagnosis for BEL):
    fnam_NA <- paste('NAs_mat_', simdist, '.csv', sep = '')
    write(x = 'NAs matrix:\n', file = fnam_NA)
    for (f in 1:length(flavors)) {
      flavor <- flavors[f]

      # Write the CP matrix:
      CP_mat <- all_CP_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_CP, append = T)
      suppressWarnings(
        write.table(x = CP_mat, file = fnam_CP, quote = F, sep = ',', row.names = rownames(CP_mat),
                    col.names = colnames(CP_mat), append = T)
      )
      write(x = '', file = fnam_CP, append = T)
      # Write the skew matrix:
      skew_mat <- all_skew_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_skew, append = T)
      suppressWarnings(
        write.table(x = skew_mat, file = fnam_skew, quote = F, sep = ',', row.names = rownames(skew_mat),
                    col.names = colnames(skew_mat), append = T)
      )
      write(x = '', file = fnam_skew, append = T)
      # Write the NA matrix:
      NAs_mat <- all_NAs_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_NA, append = T)
      suppressWarnings(
        write.table(x = NAs_mat, file = fnam_NA, quote = F, sep = ',', row.names = rownames(NAs_mat),
                    col.names = colnames(NAs_mat), append = T)
      )
      write(x = '', file = fnam_NA, append = T)
    }
  }
}


#' Bootstrap coverage probability on two samples with different methods
#'
#' @description Coverage probabilities for two samples CI of CV calculated with different NON-par bootstrapping methods.
#' @keywords comparing bootstraps, two sample, coverage probabilities
#' @import foreach
#' @import doMC
#' @param n_data1 Number of data points in the first sample to test
#' @param n_data2 Number of data points in the second sample to test
#' @param mu1 Mean of the simulated distribution for the first sample
#' @param mu2 Mean of the simulated distribution for the second sample
#' @param std1 Standard deviation of the simulated distribution for the first sample
#' @param std2 Standard deviation of the simulated distribution for the second sample
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @param prob Probability level for the confidence interval
#' @param np Number of processors to use. Default set to 4 because this is a fairly normal number these days.
#' @param boot_samples Number of bootstraps to generate
#' @param reps Replicates to make the coverage probability
#' @param flavors Choice of bootstrapping flavors [norm, basic, perc, bca, student]. Default all.
#' @return List of coverage probabilities for the different bootstrapping methods
#' @export
#' @examples
#' coverage_prob_two_sample(n_data1 = 100, n_data2 = 100, mu1 = 1, mu2 = 1,
#'                             std1 = 1, std2 = 1, simdist = 'norm', boot_samples = 500, reps = 4)
#'
#' This will take some time (approx. 5 min):
#' coverage_prob_two_sample(n_data1 = 100, n_data2 = 100, mu1 = 1, mu2 = 1, std1 = 1, std2 = 1,
#'                             simdist = 'norm', boot_samples = 250, reps = 200, flavors = c('norm'))

coverage_prob_two_sample <- function(n_data1, n_data2, mu1, mu2, std1, std2, simdist = 'norm',
                                     boot_samples = 200, reps = 1000, prob = 0.95,
                                     flavors = c("norm","basic", "perc", "bca", "student")) {
  registerDoMC(np)
  CIarray_of_lists <- foreach (i = 1:reps) %dopar% {
    coeffv1 <- std1 / mu1
    coeffv2 <- std2 / mu2
    CVdiff <- coeffv1 - coeffv2

    twosim <- draw_two_sampleCV_data(n_data1 = n_data1, n_data2 = n_data2, mu1 = mu1, mu2 = mu2,
                                     std1 = std1, std2 = std2, simdist = simdist)
    X1 <- twosim[['X1']]
    X2 <- twosim[['X2']]

    allCI <- bootcompare(x = X1, y = X2, statistic = 'CV', nsim = boot_samples, prob = prob, student = T, nstudent = 100)

    res <- list()
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      CI <- allCI[[flavor]]
      if (length(CI) != 2) {
        res[flavor] <- 'NA'
      } else if (CVdiff < CI[1]) {
        res[flavor] <- 'L'
      } else if (CVdiff > CI[2]) {
        res[flavor] <- 'H'
      } else {
        res[flavor] <- 'C'
      }
    }
    return(res)
  }

  CI_list_of_arrays <- list()
  for (i in 1:length(flavors)) {
    flavor <- flavors[i]
    CI_list_of_arrays[[flavor]] <- array(data = NA, dim = length(CIarray_of_lists))
  }

  for (j in 1:length(CIarray_of_lists)) {
    CI_list <- CIarray_of_lists[[j]]
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      CIres <- CI_list[[flavor]]
      CI_list_of_arrays[[flavor]][j] <- CIres
    }
  }

  allres <- list()
  for (i in 1:length(flavors)) {
    flavor <- flavors[i]
    CIarray <- CI_list_of_arrays[[flavor]]
    CP <- sum(CIarray == 'C') / sum(CIarray != 'NA')
    if (sum(CIarray == 'L') > 0) {
      skew <- sum(CIarray == 'H') / sum(CIarray == 'L')
    } else {
      skew <- NA
    }
    res <- list()
    res[['NAs']] <- sum(CIarray == 'NA')
    res[['CP']] <- CP
    res[['skew']] <- skew
    res[['CIarray']] <- CIarray
    allres[[flavor]] <- res
  }
  return(allres)
}


#' Bootstrap coverage probability on two samples with different methods and different distributions
#'
#' @description Coverage probabilities for two sample CI of CV calculated with different NON-par bootstrapping methods,
#' and on multiple distributions. Also a list of different CV's can be provided. Results are printed in a matrix format.
#' Resulsts include the coverage probabilities but also the skewness of this and how many NA's returned in a run.
#' @keywords comparing bootstraps, two samples, coverage probabilities
#' @import foreach
#' @import doMC
#' @param n_data_list List of data point to test
#' @param mu1 Mean of the simulated distribution for sample 1
#' @param mu2 Mean of the simulated distribution for sample 2
#' @param k_list List of standard deviations to test
#' @param simdist The distribution to sample from [norm, lognorm, gamma, chisq]
#' @param prob Probability level for the confidence interval
#' @param np Number of processors to use. Default set to 4 because this is a fairly normal number these days.
#' @param boot_samples Number of bootstraps to generate
#' @param reps Replicates to make the coverage probability
#' @param flavors Choice of bootstrapping flavors [norm, basic, perc, bca, student]. Default all.
#' @return None. All results are printed to comma separated text files.
#' @export
#' @examples
#' This example will take around 5 min to run and is just a prove of concept
#' that multiple CV's and sample sizes can easily be run. Doing this with
#' enough replicates to make the coverage probability robust will take
#' very long time and in such cases running it on a cluster is necessary.
#' # Combinations of CV's to try:
#' k_list <- list()
#' k_list[['0.2_0.2']] <- c(0.2, 0.2)
#' k_list[['0.5_1']] <- c(0.5, 1)
#' k_list[['0.5_0.5']] <- c(0.5, 0.5)
#' k_list[['1_1']] <- c(1, 1)
#' k_list[['0.2_0.5']] <- c(0.2, 0.5)
#' k_list[['0.5_0.2']] <- c(0.5, 0.2)
#' k_list[['0.2_1']] <- c(0.2, 1)
#' k_list[['1_0.2']] <- c(1, 0.2)
#'
#' # Combinations of sample sizes to try:
#' n_data_list <- list()
#' n_data_list[['30_30']] <- c(30, 30)
#' n_data_list[['30_30']] <- c(30, 30)
#' n_data_list[['60_60']] <- c(60, 60)
#' n_data_list[['90_90']] <- c(90, 90)
#'
#' multiple_coverage_prob_two_sample(mu1 = 1, mu2 = 1, boot_samples = 250, reps = 4,
#'                                   prob = 0.95, n_data_list = n_data_list,
#'                                   sim_distributions = c('gamma'),
#'                                   k_list = k_list)
multiple_coverage_prob_two_sample <- function(mu1, mu2, boot_samples = 200, reps = 1000,
                                              prob = 0.95, n_data_list,
                                              flavors = c("norm", "basic", "perc", "bca", "student"),
                                              sim_distributions = c('norm', 'lognorm', 'gamma', 'chisq'),
                                              k_list, write_dir = '.') {
  # Make the list of matrices to store the results:
  all_CP_mat <- list()
  all_skew_mat <- list()
  all_NAs_mat <- list()
  for (d in 1:length(sim_distributions)) {
    simdist <- sim_distributions[[d]]
    for (i in 1:length(flavors)) {
      flavor <- flavors[i]
      # Following three lines are such a lousy hack,
      # but for unknown reasons it must be used to initialize
      # the correct list structure...
      all_CP_mat[[simdist]][[flavor]] <- data.frame(NA)
      all_skew_mat[[simdist]][[flavor]] <- data.frame(NA)
      all_NAs_mat[[simdist]][[flavor]] <- data.frame(NA)
      all_CP_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
      all_skew_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
      all_NAs_mat[[simdist]][[flavor]] <- matrix(data = NA, nrow = length(n_data_list), ncol = length(k_list), dimnames = list(names(n_data_list), names(k_list)))
    }
  }

  # Make the run:
  for (d in 1:length(sim_distributions)) { # all the distributions to test
    simdist <- sim_distributions[[d]]
    for (i in 1:length(n_data_list)) {  # rows, n_data_list
      for (j in 1:length(k_list)) {     # cols, k_list
        res <- coverage_prob_two_sample(n_data1 = n_data_list[[i]][1], n_data2 = n_data_list[[i]][2],
                                        mu1 = mu1, mu2 = mu2, std1 = k_list[[j]][1], std2 = k_list[[j]][2],
                                        simdist = sim_distributions[[d]], boot_samples = boot_samples,
                                        reps = reps, prob = prob, flavors = flavors)
        # Asign the results for different bootstrapping flavors:
        for (f in 1:length(flavors)) {
          flavor <- flavors[f]
          all_CP_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['CP']]
          all_skew_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['skew']]
          all_NAs_mat[[simdist]][[flavor]][i, j] <- res[[flavor]][['NAs']]
        }
      }
    }
  }

  # Print the results:
  for (d in 1:length(sim_distributions)) {
    simdist <- sim_distributions[[d]]
    # Coverage probabilities:
    fnam_CP <- paste('CP_mat_', simdist, '.csv', sep = '')
    write(x = 'Coverage probabilities:\n', file = fnam_CP)
    # Skew:
    fnam_skew <- paste('skew_mat_', simdist, '.csv', sep = '')
    write(x = 'Skew matrix:\n', file = fnam_skew)
    # NA's (usefull diagnosis for BEL):
    fnam_NA <- paste('NAs_mat_', simdist, '.csv', sep = '')
    write(x = 'NAs matrix:\n', file = fnam_NA)
    for (f in 1:length(flavors)) {
      flavor <- flavors[f]

      # Write the CP matrix:
      CP_mat <- all_CP_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_CP, append = T)
      suppressWarnings(
        write.table(x = CP_mat, file = fnam_CP, quote = F, sep = ',', row.names = rownames(CP_mat),
                    col.names = colnames(CP_mat), append = T)
      )
      write(x = '', file = fnam_CP, append = T)
      # Write the skew matrix:
      skew_mat <- all_skew_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_skew, append = T)
      suppressWarnings(
        write.table(x = skew_mat, file = fnam_skew, quote = F, sep = ',', row.names = rownames(skew_mat),
                    col.names = colnames(skew_mat), append = T)
      )
      write(x = '', file = fnam_skew, append = T)
      # Write the NA matrix:
      NAs_mat <- all_NAs_mat[[simdist]][[flavor]]
      write(x = flavor, file = fnam_NA, append = T)
      suppressWarnings(
        write.table(x = NAs_mat, file = fnam_NA, quote = F, sep = ',', row.names = rownames(NAs_mat),
                    col.names = colnames(NAs_mat), append = T)
      )
      write(x = '', file = fnam_NA, append = T)
    }
  }
}

