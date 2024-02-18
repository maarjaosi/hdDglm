#' Forward filtering backward sampling
#'
#' Applies the forward filtering backward sampling algorithm sample the from the posterior distribution of regression parameters of a DLM.
#' T denotes the number of time steps and p the number of covariates / parameters of the DLM.
#'
#' @param y Realizations of the response variables. A vector of length \eqn{T}.
#' @param x Covariate vectors. A \eqn{T \times p} matrix.
#' @param m0 Prior mean of \eqn{\boldsymbol{\beta}_0}. A vector of length \eqn{p}.
#' @param C0 Prior covariance of \eqn{\boldsymbol{\beta}_0}. A \eqn{p \times p} matrix.
#' @param G Evolution matrices. A function of time.
#' @param W Error covariance. A \eqn{p \times p} matrix.
#' @param V Variances of the observation variables. A vector of length \eqn{T}.
#'
#' @return A sample from the posterior distribution \eqn{p(\boldsymbol{\beta}_{1:T} \: \vert \: \boldsymbol{y}_{1:T})} for the
#' DLM \eqn{(G_t, W, 1/V_t)_{t=1}^T}. A \eqn{(T+1) \times p} matrix whose \eqn{(t,i)} element
#' is the sampled value for \eqn{\beta_{t,i}}.
#' @export

FFBS <- function(y, x, m0, C0, G, W,V) {


  #Initialize
  T = length(y)
  p = ncol(x)
  m.list = vector("list",length=T + 1)
  C.list = vector("list",length=T + 1)
  R.list = vector("list",length=T + 1)
  a.list = vector("list",length=T + 1)

  m <- m0
  C <- C0


  m.list[[1]] <- m0
  C.list[[1]] <- C0

  #Forward Filter
  for(t in 1:T) {

    xt = as.vector(x[t,])

    a = G(t)%*%m
    R = G(t)%*%C%*%t(G(t)) + W

    f = as.numeric(xt%*%a)
    q = as.numeric(t(xt)%*%R%*%xt) + V[t]

    m = a + R%*%xt*(y[t]-f)/q
    C = R - R%*%xt%*%t(xt)%*%R/q

    #Collect m,C,R,a
    m.list[[t + 1]] = m
    C.list[[t + 1]] = C
    R.list[[t + 1]] = R
    a.list[[t + 1]] = a

  }

  #Backward Sampling
  h = m.list[[T + 1]]
  H = C.list[[T + 1]]
  theta.t1 = as.vector(t(chol(H))%*%rnorm(p) + h)
  a.t1 = as.vector(a.list[[T + 1]])
  back.samp = rbind(theta.t1,NULL)

  for(t in (T-1):1){

    m <- m.list[[t + 1]]
    C <- C.list[[t + 1]]
    R <- chol2inv(chol(R.list[[t+2]]))
    a.t1 <- as.vector(a.list[[t+2]])

    h <- m + C%*%G(t+1)%*%R%*%(theta.t1-a.t1)
    H <- C - C%*%t(G(t+1))%*%t(R)%*%G(t+1)%*%t(C)
    theta.t1 <- as.vector(t(chol(H))%*%rnorm(p) + h)
    back.samp <- rbind(theta.t1,back.samp)

  }
  R <- chol2inv(chol(R.list[[2]]))
  a.t1 <- as.vector(a.list[[2]])

  h <- m0 + C%*%G(1)%*%R%*%(theta.t1-a.t1)
  H <- C0 - C%*%t(G(1))%*%t(R)%*%G(1)%*%t(C0)

  theta.t1 <- as.vector(t(chol(C0))%*%rnorm(p) + m0)
  back.samp<- rbind(theta.t1, back.samp)

  return(back.samp)

}
