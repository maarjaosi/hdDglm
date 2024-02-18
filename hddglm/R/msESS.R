#' Monotone sequence estimator for ESS
#'
#'Monotone sequence estimator of the effective sample size (ESS) for an MCMC chain.
#' @param chain Vector containing values produced by an MCMC algorithm.
#'
#' @return The montone sequence estimator of the chain.
#' @export
msESS <- function(chain) {

  autocorrelations <- coda::autocorr(coda::mcmc(chain), lags = 1:length(chain))
  nr_autocorr <- length(autocorrelations)

  nr_of_seq_values <- nr_autocorr/2 - 0.5^(nr_autocorr %% 2)

  seq_values <- lapply(1:nr_of_seq_values,
                       function(s) {
                         autocorrelations[2*s] + autocorrelations[2*s + 1]
                         }) %>%
                  unlist()

  #index up to which the values of the sequence are positice
  index_positive <- which(seq_values < 0)[1] - 1
  index_positive <- ifelse(is.na(index_positive), nr_of_seq_values, index_positive)

  #index up to which the values of the sequence are monotone
  index_monotone <- 1
  while(seq_values[index_monotone] >= seq_values[index_monotone + 1]) {
    index_monotone <- index_monotone + 1
  }

  index_ESS <- 2*min(index_positive, index_monotone) +1

  ESS <- 1/(1 + 2*sum(autocorrelations[1:index_ESS]))

  return(ESS)
}

