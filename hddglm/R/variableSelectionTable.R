#' Variable selection.table
#'
#' Creates a table to evaluate variable selection by calculating TPR and FPR values.
#'
#' @param selected_model Model selected by some method. A vector containing the variables included into the model.
#' @param true_model The true model. A vector containing the variables included into the model.
#' @param p  Number of parameters.
#' @param t Time step at which the model was selected
#'
#' @return Data frame containing the columns \code{p}, \code{T}, \code{TPR} and \code{FPR} with the respective values
#' for the selected model.
#' @export
variableSelectionTable <- function(selected_model, true_model, p , t) {

  TP <- sum( true_model %in% selected_model)
  FP <- sum(!(selected_model %in% true_model))
  FN <- length(true_model) - TP
  TN <- (p - length(true_model)) - FP
  TPR <- TP/length(true_model)
  FPR <- FP/(p - length(true_model))


  table <- data.frame(p = p,
                      T = t,
                      TPR = TPR,
                      FPR = FPR)

  return(table)

}
