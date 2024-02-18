#' Linear Bayes Method for the dynamic logistic binary model
#'
#' Calculates estimates of the posterior means and covariance matrices
#' of regression parameters of a dynamic logistic binary model using the Linear Bayes method.
#'
#' @param y Realizations of the response variables. A vector of length T.
#' @param x Covariates. A \eqn{T \times p} matrix.
#' @param m0 Prior mean of \eqn{\boldsymbol{\beta}_0}. A vector of length p.
#' @param C0 Prior covariance of \eqn{\boldsymbol{\beta}_0}. A \eqn{p \times p} matrix.
#' @param G Evolution matrices. A function of time.
#' @param W Error covariance.  A \eqn{p \times p} matrix.
#'
#' @return A named list. First element \code{m.list} is a list containing the posterior means
#' of the regression parameters. Second element \code{C.list} is a list containing the posterior
#' covariance matrices of the regression parameters.
#' @export
logisticLinearBayes <- function(y,
                                x,
                                m0,
                                C0,
                                G,
                                W) {
  T <- length(y)
  m.list <- vector("list",length=T)
  C.list <- vector("list",length=T)

  m <- m0
  C <- C0

  for (t in 1:T) {
    xt <- x[t,]
    a <- G(t)%*%m
    R <- G(t)%*%C%*%t(G(t)) + W

    f_prior <- as.numeric(t(xt)%*%a)
    q_prior <- as.numeric(t(xt)%*%R%*%xt)
    u <- as.numeric(R%*%xt)

    s_prior <- q_prior^(-1)*(1 + exp(f_prior))
    r_prior <- q_prior^(-1)*(1 + exp(-f_prior))


    s_post <- s_prior + y[t]
    r_post <- r_prior + 1 - y[t]

    f_post <- log(s_post) - log(r_post)
    q_post <- (1/s_post) + (1/r_post)


    m <- a + u*(f_post -f_prior)/q_prior
    C <- R - u%*%t(u)*(1 - q_post/q_prior)/q_prior


    m.list[[t]] = m
    C.list[[t]] = C
  }

  return(list(m.list = m.list, C.list = C.list))

}
