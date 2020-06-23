#' The Multi-measure INdividual Deconvolution (MIND) algorithm
#'
#' It calculates the empirical Bayes estimates of subject- and cell-type-specific gene expression, via a computationally efficient EM algorithm.
#'
#' @param X bulk gene expression (gene x subject x measure).
#' @param W subject-specific cell type fraction (subject x measure x cell type).
#' @param maxIter maximum number of iterations for the EM algorithm.
#' @param tol tolerance level of absolute relative change of the log-likelihood to stop the EM algorithm.
#' @param verbose logical, to print the detailed information for each iteration: iter (the iteration number), logLike_change, sigma2_e, mean(diag(Sigma_c))).
#' @param ncore number of cores to run in parallel
#' 
#' @return A list containing the output of the EM deconvolution algorithm
#' \item{A}{the deconvolved cell-type-specific gene expression (gene x cell type x subject).}
#' \item{mu}{the estimated profile matrix (gene x cell type).}
#' \item{iter}{the number of iterations used in the EM algorithm.}
#' \item{Sigma_c}{the covariance matrix for the deconvolved cell-type-specific expression (cell type x cell type).}
#' \item{sigma2_e}{the error variance.}
#' \item{loglike}{the log-likelihood for each EM iteration.}
#' \item{var_A}{the posterior covariance matrix for A (vectorized covariance matrix by subject).}
#'
#' @references Wang, Jiebiao, Bernie Devlin, and Kathryn Roeder. "Using multiple measurements of tissue to estimate subject-and cell-type-specific gene 
#' expression." Bioinformatics 36.3 (2020): 782-788. https://doi.org/10.1093/bioinformatics/btz619
#'
#' @examples
#'
#' data(example)
#'
#' deconv = mind(X = example$X, W = example$W, ncore = 2)
#'
#' @export mind


mind = function(X, W, maxIter = 100, tol = 0.001, verbose = F, ncore = 4) {
  
  package.check('doParallel')
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  
  deconv = deconv_xx(X, W, maxIter = maxIter, tol = tol, verbose = verbose)
  A = get_A(X, W, deconv)
  
  var_A = sapply(deconv$tmp, function(x) x$var_A_x)
  
  closeAllConnections()
  
  return(list(A = A, iter = deconv$iter, Sigma_c = deconv$Sigma_c, sigma2_e = deconv$sigma2_e, loglike = deconv$loglike,
              mu = deconv$mu, var_A = var_A))
}


# three intermediate steps: get_xx, deconv_xx, get_A

get_xx = function(X, W, mu) {
  # mu: gene x cell type
  p = nrow(X)
  n = ncol(X)
  t = dim(X)[3]
  
  XX = foreach(i = 1:n) %dopar% {

    # demeaned
    X[,i,] = X[,i,] - mu %*% t(matrix(W[i,,], nrow = t)) # gene x measure
    
    miss = which(is.na(X[1,i,])) # missing/uncollected measures
    ti = t - length(miss)
    
    if(t != ti) {
      Xi = matrix(X[,i,-miss], ncol = ti)
    } else {
      Xi = matrix(X[,i,], ncol = t)
    }
    
    XX = matrix(0, ti, ti)
    for(j in 1:p) XX = XX + Xi[j,] %*% t(Xi[j,])
    return(XX)
  }
  
  return(XX)
}


deconv_xx = function(X, W, maxIter = 100, tol = 0.001, verbose = FALSE) {
  
  p = nrow(X)
  n = ncol(X)
  t = dim(X)[3]
  k = dim(W)[3]
  
  ## starting values
  Sigma_c = diag(k)
  sigma2_e = 1
  mu = matrix(apply(X, 1, mean, na.rm = T), p, k)
  
  iter = 0
  likelihood = NULL
  cond = TRUE
  
  while (cond) {
    
    iter = iter + 1
    
    Sigma_c_inv = solve(Sigma_c)
    
    ls = list()
    XX_all = get_xx(X, W, mu) # this slows down the algorithm
    
    WW = matrix(0, k, k)
    WX = matrix(0, k, p)
    
    for(i in 1:n) {
      
      XX = XX_all[[i]]
      
      miss = which(is.na(X[1,i,])) # missing/uncollected measures
      if(length(miss) > 0) {
        Wi = W[i,,][-miss,,drop = F]
        Xi = matrix(X[,i,-miss], ncol = t - length(miss))
      } else {
        Wi = matrix(W[i,,], nrow = t)
        Xi = matrix(X[,i,], ncol = t)
      }
      
      I0 = diag((t - length(miss)))
      Sigma = Wi %*% Sigma_c %*% t(Wi) + sigma2_e * I0
      Sigma_inv = solve(Sigma)
      
      obs_likelihood = - 0.5 * (p * determinant(x = Sigma, logarithm = TRUE)$modulus + sum(diag(Sigma_inv %*% XX)))
      
      ## E step
      
      var_A_x = solve(t(Wi)%*%Wi / sigma2_e + Sigma_c_inv)
      
      ## M step
      
      # for Sigma_c
      # A * A
      SS_sigma_c = var_A_x * p + var_A_x %*% t(Wi) %*% XX %*% Wi %*% var_A_x / sigma2_e^2
      
      var_e_x = sigma2_e * I0 - sigma2_e^2 * Sigma_inv
      
      SS_sigma2_e = sum(diag(var_e_x)) * p + sigma2_e^2 * sum(diag(Sigma_inv %*% XX %*% Sigma_inv)) #Sigma_inv %*% XX:repeated in obs_likelihood
      
      ls[[i]] = (list(obs_likelihood = obs_likelihood,
                      SS_sigma_c = SS_sigma_c,
                      SS_sigma2_e = SS_sigma2_e,
                      var_A_x = var_A_x
      ))
      
      WW = WW + t(Wi) %*% Sigma_inv %*% Wi
      
      WX = WX + t(Wi) %*% Sigma_inv %*% t(Xi)
      
    }
    
    mu = t(solve(WW) %*% WX)
    
    obs_likelihood = sum(sapply(ls, function(x) x$obs_likelihood))
    SS_sigma_c = Reduce('+', lapply(ls, function(x) x$SS_sigma_c))
    SS_sigma2_e = sum(sapply(ls, function(x) x$SS_sigma2_e))
    
    sigma2_e = as.numeric(SS_sigma2_e / (p * sum(!is.na(X[1,,]))))
    
    Sigma_c = SS_sigma_c / (n*p)
    
    # log-likelihood
    likelihood = c(likelihood, obs_likelihood)
    
    # relative change of log-likelihood
    likelihood_change = ifelse(iter >= 2, (likelihood[iter] - likelihood[iter - 1])/abs(likelihood[iter - 1]), NA)
    cond = (iter < maxIter & ifelse(iter >= 2, abs(likelihood_change) > tol, TRUE))
    
    if (verbose) {
      print(round(c(iter = iter, logLike_change = likelihood_change, sigma2_e = sigma2_e,
                    Sigma_c = mean(diag(Sigma_c))), log10(1/tol)))
    }
  }
  
  rownames(mu) = rownames(X)
  
  return(list(iter = iter, Sigma_c = Sigma_c, sigma2_e = sigma2_e, loglike = likelihood, tmp = ls, mu = mu))
}


get_A = function(X, W, deconv) {
  
  p = nrow(X)
  n = ncol(X)
  t = dim(X)[3]
  k = dim(W)[3]
  
  sigma2_e = deconv$sigma2_e
  
  ls = array(NA, dim = c(p, k, n))
  for(i in 1:n) {
    
    miss = which(is.na(X[1,i,])) # missing/uncollected measures
    ti = t - length(miss)
    
    if(t != ti) {
      Wi = W[i,,][-miss,,drop = F]
      Xi = matrix(X[,i,-miss], ncol = ti)
    } else {
      Wi = matrix(W[i,,], nrow = t)
      Xi = matrix(X[,i,], ncol = t)
    }
    
    A = matrix(0, p, k)
    
    var_A_x = deconv$tmp[[i]]$var_A_x
    for(j in 1:p) A[j,] = var_A_x %*% t(Wi) %*% (Xi[j,] - Wi %*% deconv$mu[j,]) / sigma2_e + deconv$mu[j,]
    
    ls[,,i] = A
  }
  
  rownames(ls) = rownames(X)
  colnames(ls) = dimnames(W)[[3]]
  dimnames(ls)[[3]] = colnames(X)
  
  return(ls)
}



#' Estimating cell type fractions with a signature matrix using non-negative least squares (NNLS)
#'
#' It calls the nnls package to estimate cell type fractions of bulk data using a pre-estimated signature matrix. It is recommended to 
#' keep the row and column names of the input data.
#'
#' @param sig signature matrix (marker gene x cell type).
#' @param bulk bulk data that need to be deconvolved (gene x tissue sample).
#'
#' @return A matrix containing the estimated cell type fractions (tissue sample x cell type). Row sums have been normalized to be 1 per sample.
#'
#' @export est_frac

est_frac = function (sig, bulk) {
  sig = as.matrix(sig)
  bulk = as.matrix(bulk[rownames(sig),])
  package.check("nnls")
  nls <- apply(bulk, 2, function(b) nnls(sig, b)$x)
  nls = t(apply(nls, 2, function(x) x/sum(x)))
  colnames(nls) = colnames(sig)
  rownames(nls) = colnames(bulk)
  print(round(colMeans(nls), 2))
  return(nls)
}



# use this function to check if each package is on the local machine
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package.check <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}


#' A data example
#'
#' A data list for demonstration.
#'
#'
#' @name example
#' @docType data
#' @return A list containing
#' \item{X}{bulk gene expression (gene x subject x measure).}
#' \item{W}{subject-specific cell type fraction (subject x measure x cell type).}
#' @examples
#'
#' data(example)
#'
NULL
