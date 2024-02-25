#' DGLM Moment estimates.
#'
#' Estimates the mean and covariance of the regression parameter \eqn{\boldsymbol{beta}_t}
#' from values produced by an MCMC algorithm.
#'
#' @param chain  The MCMC chain of the regression parameters of a DGLM.
#' @param p Number of parameters.
#' @param t Time step for which the moments are calculated
#' @param K Number of runs of the MCMC algorithm .
#' @param burn_in The burn-in period for the MCMC chain. An integer.
#'
#' @return A named list. The first element \code{mean} is the mean and the second element \code{cov} is the covariance of the DGLM regression parameter
#'at time t
#'@export
estimateMoments <- function(chain,
                            p,
                            t,
                            K,
                            burn_in = 500) {

  #Data frame where the i-th column contains the values of the simulation for the
  #beta_{ti}
  beta_df <- purrr::map(.x  = 1:p,
                        .f = function(j) {
                          index_covariate <- seq(j, (K-1)*p + j , p)
                          covariate_name <- paste("beta", j, sep = "_")
                          df <- data.frame(chain[t, index_covariate][burn_in:K])
                          return(df)
                        }

  ) %>% purrr::list_cbind()

  #estimate the posterior moments
  cov <- cov(beta_df)
  mean  <- apply(beta_df, 2, mean)

  return(list(mean = mean, cov = cov))
}






