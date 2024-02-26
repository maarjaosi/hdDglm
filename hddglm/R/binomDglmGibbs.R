#'Updates the values of the rate parameter of the error covariance component.
#'
#' @param gamma Values of the gamma parameter. Vector of length p.
#' @param beta Values of the regression parameters; \eqn{T \times p} matrix.
#' @param G Evolution matrices. An \code{R} function \code{G(t)} that returns the \eqn{p \times p} evolution matrix \eqn{G_t} dependent on the time argument \eqn{t}.
#'
#' @return Updated value of the gamma parameters. A vector of length p.
updateError <- function(gamma,
                         beta,
                         G) {
  T <- nrow(beta) - 1
  p <- ncol(beta)
  beta_evolution <- rlist::list.rbind(purrr::map(.x = 1:T,
                        .f = function(t) {
                              return(t(G(t)%*%beta[t, ]))
                        }
  ))
  gamma_updated <- gamma + 0.5*t(colSums((beta[2:(T+1),] - beta_evolution)^2))
  return(gamma_updated)
}

#'Sample latent polya-gamma variables
#'
#' @param b parameter of the polya-gamma distribution; vector of length T
#' @param x covariates; a Txp matrix
#' @param beta values of the regression parameters; (T)xp matrix
#'
#' @return T samples from a Polya-Gamma distribution, where the t-th sample is
#' from a PG(b[t], x_t^Tbeta_t + offset) distribution
sampleLatentPg <- function(b,
                           x,
                           beta) {
  T <- length(b)
  omega <- rep(0,T)
  for (t in 1:T) {
    eta <- t(x[t,])%*%beta[t+1,]
    omega[t] <- BayesLogit::rpg(1, b[t], eta)
  }

  return(omega)
}

#' Gibbs sampler for dynamic binomial models.
#'
#' Samples from the posterior distribution of the regression parameters of a dynamic binomial model using a Gibbs sampler based on data augmentation.
#'
#' @param model A string indicating which observation distribution to use. Value can be either "binom" for the binomial or "nbinom" for the negative-binomial distributio.
#' the default value is "binom".
#' @param n Number of trials. A parameter of the (negative)-binomial distribution.
#' @param K Number of iterations of the Gibbs sampler. An integer.
#' @param y Realizations of the response variables. A vector of length T.
#' @param x Covariates. A \eqn{T \times p} matrix.
#' @param m0 Prior mean of \eqn{\boldsymbol{\beta}_0}. A vector of length p
#' @param C0 Prior covariance of \eqn{\boldsymbol{\beta}_0}. A \eqn{p \times p} matrix
#' @param G Evolution matrices. An \code{R} function \code{G(t)} that returns the \eqn{p \times p} evolution matrix \eqn{G_t} dependent on the time argument \eqn{t}.
#' @param omega Initial value for the latent variables. A vector of length T.
#' @param alpha Initial values of the shape parameter for the distribution of the elements of the error covariance. Vector of length p.
#' @param gamma Initial values of the rate parameter for the distribution of the elements of the error covariance. Vector of length p.
#'
#' @return K samples of the regression parameters and the covariance matrix from the posterior \eqn{p(\boldsymbol{\beta}_{1:T}, W \: \vert \: \boldsymbol{y}_{1:T})}
#' where \eqn{\boldsymbol{\beta}_t} are the regression parameters of the dynamic binomial model \eqn{(G_t, W)_{t=1}^T}.
#' A named list containing two matrices. The first matrix \code{W_chain} contains
#' the diagonal values of sampled values for the error covariance matrix. The k-th column of
#' \code{W_chain} contains the diagonal for the k-th sample.
#' The second matrix \code{beta_chain} contains the sampled values for the regression parameters \eqn{\boldsymbol{\beta}_t}
#' The element in the in \eqn{t}-th row and \eqn{p\cdot(k-1) + i}-th column of \code{beta_chain} contains the K-th sample of the regression parameter \eqn{\beta}_{ti}.
#' @export
binomDglmGibbs <- function(model = "binom",
                           K,
                           y,
                           x,
                           G,
                           m0,
                           C0,
                           n,
                           omega,
                           alpha,
                           gamma) {

  T <- length(y) #number of time steps
  p <- ncol(x) #number of parameters

  #initialize the error covariance matrix
  tau_cur <- rgamma(p, shape = alpha, rate = gamma)
  W_cur <- diag(1/tau_cur)

  #initialize parameters
  W_all <- matrix(0, nrow = p, ncol = K)
  alpha_cur <- alpha
  gamma_cur <- gamma
  omega_cur <- omega
  beta_all <- matrix(0, nrow = T + 1, ncol =  p*K)
  beta_cur <- matrix(0, nrow = T + 1, ncol = p)


  for (k in 1:K) {
    #create the pseudo data
    if (model == "binom") {
      z <- (y - n/2)/omega_cur
    }
    else if (model == "nbinom") {
      z <- (y - (y + n)/2)/omega_cur - log(n)
    }


    #sample the DGLM parameters
    beta_cur <- FFBS(y = z,
                     x = x,
                     m0 = m0,
                     C0 = C0,
                     G = G,
                     W = W_cur,
                     V = 1/omega_cur)

    #update error covariance
    alpha_cur <- alpha + T/2
    gamma_cur <-  updateError(gamma,
                               beta_cur,
                               G)

    tau_cur <-  rgamma(p, shape = alpha_cur, rate = gamma_cur)
    W_cur <- diag(1/tau_cur, ncol = p, nrow = p)

    #sample the latent variables
    if (model == "binom") {
      omega_cur <- sampleLatentPg(rep(n, T), x, beta_cur)
    }
    else if (model == "nbinom") {
      omega_cur <- sampleLatentPg(y + rep(n, T), x, beta_cur)
    }


    #gather the current parameters
    W_all[, k] <- matrix(diag(W_cur), nrow = p)
    beta_all[, (p*(k-1) + 1):(p*(k-1) + p)] <- beta_cur
  }
  return(list(W_chain = W_all, beta_chain = beta_all))
}

