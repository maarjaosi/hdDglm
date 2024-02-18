#' Variable ranking table.
#'
#' Creates a variable ranking table based on a solution path. The table shows
#' for each model in the sequence which variable was newly included into this model.
#'
#' @param solution_path Sequence of models produced by penalized credible regions method.
#' A \eqn{p \times p} matrix whose (i,j) element is TRUE
#' when the j-th parameter belongs into the i-th model else FALSE.
#' @return A data frame with rows \code{model_index} and \code{variable}.
#' The row \code{model_index} contains the index of the model in the sequence.
#' The row \code{variable} contains the variable (index).
#' @export
variableRanking <- function(model_sequence) {
  p <- ncol(model_sequence)
  variable_ranking <- rep(0, p)

  for (i in 1:p) {
    i_appears <- which(model_sequence[,i])[1]
    variable_ranking[i_appears] <- i
  }

  variable_ranking_table <- rbind(model_index = 1:p,
                                  variable = variable_ranking)

  return(variable_ranking_table)
}
