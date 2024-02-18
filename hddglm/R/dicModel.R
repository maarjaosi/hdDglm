#' Variable selection using DIC.
#'
#' Select variables according to the DIC for a dynamic binary logistic model based
#' on a solution path obtained by applying the penalized credible regions method.
#'
#'
#' @param beta_chain A chain of values for the regression parameters of a DGLM produced by an MCMC algorithm.
#' @param solution_path Sequence of models produced by penalized credible regions method.
#' A \eqn{p \times p} matrix whose (i,j) element is TRUE
#' when the j-th parameter belongs into the i-th model else FALSE.
#' @param mean Mean of the regression parameter at time t. A vector of length p.
#' @param x Covariance vector at time t. A vector of length t.
#' @param y Value of the response variable at time t.
#' @param t Time step for which the models are evaluated.
#' @param burn_in The burn in period for the MCMC chain.
#'
#' @return The selected variables according to DIC. A vector containing the selected variables.
#'
#' @export
dicModel <- function(beta_chain,
                     solution_path,
                     mean,
                     x,
                     y,
                     t,
                     burn_in = 500) {

  beta_df <- purrr::map(.x  = 1:p,
                        .f = function(j) {
                          index_covariate <- seq(j, (K-1)*p + j , p)
                          covariate_name <- paste("beta", j, sep = "_")
                          df <- data.frame(beta_chain[t+1, index_covariate][burn_in:K])
                          return(df)
                        }

  ) %>% purrr::list_cbind()

  dic_scores <- rep(0,(nrow(solution_path)))

  for (i in 1:nrow(solution_path)) {
    cur_coef <- solution_path[i,]
    cur_df <- beta_df
    cur_df[, !cur_coef] <- 0
    mean_deviance <- mean(apply(cur_df,  1, logisticDeviance, y,  x))
    mean_cur <- mean
    mean_cur[!cur_coef] <- 0
    dic_scores[i] = 2*mean_deviance - logisticDeviance(mean_cur, y, x)
  }

  selected_model_index <- which(dic_scores == min(dic_scores, na.rm = TRUE))
  selected_variables <- (1:p)[solution_path[selected_model_index,]]

  return(selected_variables)
}


#' Calculate deviance of the binary logistic model.
#'
#' @param beta Regression parameters. Vector of length p.
#' @param y Value of the response at time t.
#' @param x covariate vector; vector of length p
#'
#' @return Value of deviance.
logisticDeviance <- function(beta, y, x) {
  eta <- x%*%(beta)
  q <- exp(eta)/(1 + exp(eta))
  value <- - 2*log((q^y)*((1-q)^(1- y)))
  return(value)
}


