#' Calculate AUC of a ROC curve
#'
#' @param data data frame containing columns FPR and TPR that the ROC curve is based on
#'
#' @return AUC value
calculateAUC <- function(data) {
  data_aggregated <- data %>% dplyr::group_by(FPR) %>%
    dplyr::summarize(TPR = max(TPR)) %>%
    dplyr::group_by(TPR) %>%
    dplyr::summarize(min_FPR = min(FPR),
                     max_FPR = max(FPR)) %>%
    dplyr::mutate(auc = (max_FPR - min_FPR)*TPR)
    AUC <- sum(data_aggregated$auc)
  return(AUC)

}



#' ROC plots and AUC of the solution paths.
#'
#' Creates a ROC plot for a solution path and calculates its AUC.
#'
#' @param solution_path Sequence of models produced by penalized credible regions method.
#' A \eqn{p \times p} matrix whose (i,j) element is TRUE
#' when the j-th parameter belongs into the i-th model else FALSE.wise false.
#' @param true_model Vector of length p indicating which variables belong into the true model.
#' The i-th element is \code{TRUE} if the i-th element belongs into the true model else
#' it is \code{FALSE}.
#' @param method Name of the method used. Will be included into the title of the plot.
#'
#' @return A named list. First element \code{plot} is the ROC-curve plot and the second element
#' \code{AUC} is the AUC value for the ROC-curve.
#' @export
rocPlot <- function(solution_path, true_model, method) {
  n <- nrow(solution_path)
  p <- ncol(solution_path)

  false_positive_rate <- rep(0,n)
  true_positive_rate <- rep(0,n)
  for (i in 1:n) {
    false_positive_rate[i] <- sum(solution_path[i, ] & !(true_model))/sum(!true_model)
    true_positive_rate[i] <- sum(solution_path[i, ] & true_model)/sum(true_model)
  }



  data <- data.frame(FPR = false_positive_rate,
                     TPR = true_positive_rate)

  AUC <- round(calculateAUC(data), digits = 3)

  p <- ggplot2::ggplot(data, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(colour = "red") +
    ggplot2::geom_line(ggplot2::aes(y = FPR),linetype = "dashed") +
    ggplot2::labs(title = paste(method, ", p = ", p,
                                ", T = ", T, "\n", "ROC-curve of the solution path \n",
                                "AUC = ", AUC, sep = "")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) #title aligned in the middle


  return(list(plot = p, AUC = AUC))
}
