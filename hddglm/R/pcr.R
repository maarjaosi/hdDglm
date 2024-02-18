#' Credible Regions Method
#'
#' Calcualtes a solution path by applying the credible regions method.
#'
#' @param mean Mean of the regression parameter. Vector of length p.
#' @param cov Covariance of the regression parameter. A \eqn{p \times p} matrix.
#'
#' @return A named list containing elements \code{solution_path} and
#' \code{coef_values}.
#'
#' The element \code{solution_path} and is a \eqn{p \times p} matrix whose (i,j) element is TRUE
#' when the j-th parameter belongs into the i-th model else FALSE. The element \code{coef_values}
#' is a \eqn{p \times p} matrix containing the corresponding parameter values obtained
#' by solving the related optimization problem.
#'
#' @export
pcr <- function(mean, cov) {

  p <- length(mean)

  #artifical data for a linear model
  D <- diag(mean^2)
  X <- solve(expm::sqrtm(cov))%*%D
  Y <- solve(expm::sqrtm(cov))%*%mean

  #apply lasso to solve the optimization problems
  coef.lars <- lars::lars(X,
                          Y,
                          intercept = FALSE,
                          use.Gram=F,
                          normalize = FALSE)

  coef_values <- coef(coef.lars)
  included_coef <- coef(coef.lars) != 0

  #the first model is the null model which we remove
  return(list(solution_path = included_coef[-1,],
              coef_values = coef_values[-1,]))
}
