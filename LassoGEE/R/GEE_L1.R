#' @title Function to fit penalized GEE by I-CGD algorithm.
#' @description  This function fits a $L_{1}$ penalized GEE model to longitudinal
#'  data by I-CGD algorithm.
#' @param X A design matrix of dimension  \code{(nm) * p}.
#' @param y A response vector of length  \code{m * n}.
#' @param id A vector for identifying subjects/clusters.
#' @param family A family object: a list of functions and expressions
#' for defining link and variance functions. Families supported here is same as
#'   in \pkg{PGEE} which are binomial, gaussian, gamma and poisson.
#' @param lambda A numerical value for the penalization parameter.
#' @param corstr A character string, which specifies the type of correlation
#' structure. Structures supported in PGEE are "AR-1","exchangeable",
#' "independence", and "unstructured". The default corstr type is "independence".
#' @param beta.ini User specified initial values for regression parameters. The default value is NULL.
#' @param R User specified correlation matrix. The default value is \code{NULL}.
#' @param scale.fix A logical variable; if true, the scale parameter is
#' fixed at the value of scale.value. The default value is TRUE.
#' @param scale.value If  \code{scale.fix = TRUE}, this assignes a
#' numeric value to which the scale parameter should be fixed.
#' The default value is 1.
#' @return A list containing the following components:
#'   \item{betaest}{return final estimation}
#'   \item{beta_all_step}{return estimate in each iteration}
#'   \item{inner.count}{iterative count in each stage}
#'   \item{outer.iter}{iterate number of outer loop}
#' @references Li, Y., Gao, X., and Xu, W. (2020). Statistical consistency for
#' generalized estimating equation with $L_{1}$ regularization.
#' @export
#' @importFrom PGEE mycor_gee1
#' @importFrom stats binomial
#' @import MASS
#' @import RcppArmadillo
#' @import Rcpp evalCpp
#' @useDynLib LassoGEE, .registration = TRUE
#' @import mvtnorm
#' @examples
#' \dontrun{
#' set.seed(123)
#' p <- 256
#' s <- ceiling(p^{1/3})
#' n <- ceiling(10 * s * log(p))
#' m <- 4
#' # covariance matrix of p number of continuous covariates
#' X.sigma <- matrix(0, p, p)
#' {
#'   for (i in 1:p)
#'     X.sigma[i,] <- 0.5^(abs((1:p)-i))
#' }
#'
#' # generate matrix of covariates
#' X <- as.matrix(rmvnorm(n*m, mean = rep(0,p), X.sigma))
#'
#' # true regression parameter associated with the covariate
#' bt <- runif(s, 0.05, 0.5) # = rep(1/s,s)
#' beta.true <- c(bt,rep(0,p-s))
#' # intercept
#' beta_intercepts <- 0
#' # unstructure
#' tt <- runif(m*m,-1,1)
#' Rtmp <- t(matrix(tt, m,m))%*%matrix(tt, m,m)+diag(1,4)
#' R_tr <- diag(diag(Rtmp)^{-1/2})%*%Rtmp%*%diag(diag(Rtmp)^{-1/2})
#' diag(R_tr) = round(diag(R_tr))
#'
#' # library(SimCorMultRes)
#' # simulation of clustered binary responses
#' simulated_binary_dataset <- rbin(clsize = m, intercepts = beta_intercepts,
#'                                  betas = beta.true, xformula = ~X, cor.matrix = R_tr,
#'                                  link = "probit")
#' lambda <- 0.2* s *sqrt(log(p)/n)
#' data = simulated_binary_dataset$simdata
#' y = data$y
#' X = data$X
#' id = data$id
#'
#' ptm <- proc.time()
#' nCGDfit = LassoGEE(X = X, y = y, id = id, family = binomial("probit"),
#'                  lambda = lambda, corstr = "unstructured")
#' proc.time() - ptm
#' betaest <- nCGDfit$betaest
#' }
LassoGEE <- function(X, y, id, family = binomial("probit"), lambda,
                   corstr = "independence",  beta.ini = NULL,
                   R = NULL, scale.fix = TRUE, scale.value = 1)  {


  N<-length(unique(id))
  nx=ncol(X)

  avec <- as.integer(unlist(lapply(split(id, id), "length")))
  maxclsz <-max(avec)
  maxcl <- maxclsz
  nt<-avec
  nobs<-sum(nt)
  Mv <- NULL
  Mv <- as.integer(Mv)


  if (!is.null(beta.ini)) {
    betaest <- matrix(beta.ini, ncol = 1)
    if(nrow(betaest) != nx) {stop("Dimension of beta != ncol(X)!")}
    #message("user\'s initial regression estimate")

  } else {
    betaest = c(rep(0,nx))
  }

  aindex=cumsum(nt)
  index=c(0,aindex[-length(aindex)])


  diff<-1
  iter<-0
  maxiter <- 30
  count <- c()
  beta_all_step <- list()

  while(iter < maxiter) {

    R.fi.hat <- PGEE::mycor_gee1( N, nt, y, X, family, beta_new = betaest,
                                  corstr = corstr, Mv = Mv, maxclsz = maxclsz,
                                  R = R, scale.fix = scale.fix,
                                  scale.value = scale.value)
    Rhat <- R.fi.hat$Ehat
    fihat <- R.fi.hat$fi

    eta <- drop(X%*%betaest)
    mu=family$linkinv(eta)
    mu_eta = family$mu.eta(eta)
    vari = family$variance(mu)
    # S.H.E.M.val = .Call("SHM", X, y, mu, mu_eta, vari, nt, index, Rhat, N,
    #                     fihat, PACKAGE = "LassoGEE")
    S.H.E.M.val = SHM(X = X, y = y, mu = mu, mu_eta = mu_eta, vari = vari,
                      nt = nt, index = index, Rhat = Rhat,
                      N = N, fihat = fihat)
    S<-S.H.E.M.val$S
    v<-S.H.E.M.val$H
    u<- v%*%betaest + S

    # betaest1<-WLreglass(v,u,rep(lambda*N,nx))
    inner.fit <- ac_prox_grad(u = u, v = v, lambda = rep(lambda*N,nx))
    betaest1 <- inner.fit$beta_k
    diff<-sum(abs(betaest-betaest1))
    betaest<-betaest1

    iter<-iter+1
    count[iter] <- inner.fit$k
    beta_all_step[[iter]] <- inner.fit$beta_inner_step
    # if (silent==0) cat("iter",iter,"beta_new",beta_new,"diff",diff,"\n")
    if (diff <= 0.0001) {
      cat("iter: ",iter, "diff: ",diff,"\n")
      break
    }

  }

  return(list(betaest = c(betaest), beta_all_step = beta_all_step,
              inner.count = count, outer.iter = iter))
}





# proximal of L1 norm
prox_L1 = function(x, lambda){
  return(sign(x) * pmax(0, abs(x) - lambda))
}

##
ac_prox_grad <- function(u, v, lambda) {

  # initilization
  L <- abs((eigen(v + 0.0001*diag(ncol(v)))$values)[1])
  ## step size
  # nu <- 0.5             ## line search for t paramter

  beta_last <-solve(v)%*%u          ## initilization of x
  z_last <- beta_last
  t_last <- 1

  k <- 0
  maxiter <- 100           ## number of iterations
  tol <- 1e-4
  beta_inner_step <- beta_last

  while(k < maxiter) {

    need_project <- z_last-(v%*%z_last)/L + u/L
    beta_new <- prox_L1(need_project, lambda/L)
    distance <- beta_new-beta_last
    k <- k+1
    cat(k, fill = TRUE)
    beta_inner_step <- cbind(beta_inner_step, beta_new)

    if (sqrt(sum(distance^2)) <= tol * sqrt(sum(beta_last^2))) {
      break
    }
    t_new <- (1 + sqrt(1 + 4 * t_last * t_last)) / 2
    z_new <- beta_new + (t_last - 1) * distance / t_new
    beta_last <- beta_new
    t_last <- t_new
    z_last <- z_new
  }

  return(list(beta_k = beta_new, beta_inner_step = beta_inner_step, k = k))

}



