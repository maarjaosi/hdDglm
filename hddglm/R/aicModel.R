#' Variable selection using AIC.
#'
#' Select variables according to the AIC for a dynamic binary logistic model based
#' on a solution path obtained by applying the penalized credible regions method.
#'
#' @param solution_path Sequence of models produced by penalized credible regions method.
#' A \eqn{p \times p} matrix whose (i,j) element is TRUE
#' when the j-th parameter belongs into the i-th model else FALSE.
#' @param coef_values Values of the coefficients obtained by creating the model sequence.
#' @param mean Mean of the regression parameter at time t. A vector of length p.
#' @param x Covariance vector at time t. A vector of length t.
#' @param y Value of the response variable at time t.
#' @param t Time step for which the models are evaluated.
#'
#' @return The selected variables according to AIC. A vector containing the selected variables.
#' @export
aicModel <- function(solution_path,
                     coef_values,
                     mean,
                     x,
                     y,
                     t) {

  AIC <- rep(0,(nrow(solution_path)))

  for (i in 1:nrow(solution_path)) {
    eta_i <- x%*%(coef_values[i, ]*mean^2)
    q_i <- exp(eta_i)/(1 + exp(eta_i))
    AIC[i] = 2*(i)  - 2*log((q_i^y)*((1-q_i)^(1- y)))
  }

  selected_model_index <- which(AIC == min(AIC, na.rm = TRUE))

  selected_variables <- (1:p)[solution_path[selected_model_index,]]
  return(selected_variables)

}
