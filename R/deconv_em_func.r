#' The Multi-measure INdividual Deconvolution (MIND) algorithm
#'
#' It calculates the empirical Bayes estimates of cell-type-specific gene expression, per individual, via a computationally efficient EM algorithm.
#'
#' @param X bulk gene expression (gene x individual x measure).
#' @param W subject-specific cell type fraction (individual x measure x cell type).
#' @param maxIter maximum number of iterations for the EM algorithm.
#' @param tol tolerance level of absolute relative change of the log-likelihood to stop the EM algorithm.
#' @param verbose logical, to print the detailed information for each iteration: iter (the iteration number), logLike_change, sigma2_e, Sigma_c.
#'
#' @return A list containing the output of the EM deconvolution algorithm
#' \item{alpha}{the deconvolved cell-type-specific gene expression (gene x cell type x individual).}
#' \item{iter}{the number of iterations used in the EM algorithm.}
#' \item{Sigma_c}{the covariance matrix for the deconvolved cell-type-specific expression (cell type x cell type).}
#' \item{sigma2_e}{the error variance.}
#' \item{loglikelihood}{the log-likelihood for each EM iterations.}
#'
#' @references Wang, Jiebiao, Bernie Devlin, Kathryn Roeder.
#' Using multiple measurements of tissue to estimate individual- and cell-type-specific gene expression via deconvolution. Submitted.
#'
#' @examples
#'
#' data(example)
#'
#' deconv = mind(X = example$X, W = example$W)
#'
#' get_network(alpha = deconv$alpha, W = example$W, cell_type = 3, cor_cutoff = 0.7)
#'
#'
#' @export mind

mind = function(X, W, maxIter = 100, tol = 0.001, verbose = F) {

  # all individuals deconvolved together
  XX1 = get_xx(X)
  deconv = deconv_xx2(X, XX1, W, maxIter = maxIter, tol = tol, verbose = verbose)
  alpha = get_alpha(X, W, deconv)

  return(list(alpha = alpha, iter = deconv$iter, Sigma_c = deconv$Sigma_c, sigma2_e = deconv$sigma2_e, loglikelihood = deconv$loglikelihood))
}



# three steps: get_xx, deconv_xx, get_alpha

get_xx = function(X) {

  p = nrow(X)
  n = ncol(X)
  t = dim(X)[3]

  XX = list()

  for(i in 1:(n)) {

    miss = which(is.na(X[1,i,])) # missing/uncollected measures
    ti = t - length(miss)

    if(t != ti) {
      Xi = matrix(X[,i,-miss], ncol = ti)
    } else {
      Xi = matrix(X[,i,], ncol = t)
    }

    XX[[i]] = matrix(0, ti, ti)

    for(j in 1:p) XX[[i]] = XX[[i]] + Xi[j,] %*% t(Xi[j,])
  }

  return(XX)
}

# allow subject-specific W (individual x measure x cell type)

deconv_xx2 = function(X, XX_all, W, maxIter = 100, tol = 0.001, verbose = FALSE) {

  ## for balanced data first

  p = nrow(X)
  n = ncol(X)
  t = dim(X)[3]
  k = dim(W)[3]

  ## starting values

  Sigma_c = diag(k)

  sigma2_e = 1

  iter = 0
  likelihood = NULL
  cond = TRUE

  while (cond) {

    iter = iter + 1

    Sigma_c_inv = solve(Sigma_c)

    ls = list()

    for(i in 1:(n)) {

      XX = XX_all[[i]]

      miss = which(is.na(X[1,i,])) # missing/uncollected measures
      if(length(miss) > 0) {
        Wi = W[i,,][-miss,,drop = F]
      } else {
        Wi = matrix(W[i,,], nrow = t)
      }

      I0 = diag((t - length(miss)))
      Sigma = Wi %*% Sigma_c %*% t(Wi) + sigma2_e * I0
      Sigma_inv = solve(Sigma)

      obs_likelihood = - 0.5 * (p * determinant(x = Sigma, logarithm = TRUE)$modulus + sum(diag(Sigma_inv %*% XX)))

      ## E step

      ### make a list to save
      var_alpha_x = solve(t(Wi)%*%Wi / sigma2_e + Sigma_c_inv)

      ## M step

      # for Sigma_c

      SS_sigma_c = var_alpha_x * p + var_alpha_x %*% t(Wi) %*% XX %*% Wi %*% var_alpha_x / sigma2_e^2

      var_e_x = sigma2_e * I0 - sigma2_e^2 * Sigma_inv

      SS_sigma2_e = sum(diag(var_e_x)) * p + sigma2_e^2 * sum(diag(Sigma_inv %*% XX %*% Sigma_inv))

      ls[[i]] = (list(obs_likelihood = obs_likelihood,
                      SS_sigma_c = SS_sigma_c,
                      SS_sigma2_e = SS_sigma2_e,
                      var_alpha_x = var_alpha_x
      ))

    }

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
                    Sigma_c = mean(diag(Sigma_c))), 3))
    }
  }

  return(list(iter = iter, Sigma_c = Sigma_c, sigma2_e = sigma2_e, loglikelihood = likelihood, tmp = ls))
}

get_alpha = function(X, W, deconv) {

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

    alpha = matrix(0, p, k)

    var_alpha_x = deconv$tmp[[i]]$var_alpha_x
    for(j in 1:p) alpha[j,] = var_alpha_x %*% t(Wi) %*% Xi[j,] / sigma2_e

    ls[,,i] = alpha
  }

  rownames(ls) = rownames(X)
  colnames(ls) = dimnames(W)[[3]]

  return(ls)
}


#' Making a network plot for co-expressed genes
#'
#' It makes a circle network plot for co-expressed genes based on the weighted correlation of the deconvolved cell-type-specific expression. The correlation is
#' weighted by the average cell type fraction per person.
#'
#' @param alpha the deconvolved cell-type-specific gene expression (gene x cell type x individual).
#' @param W subject-specific cell type fraction (individual x measure x cell type).
#' @param cell_type a numeric index indicating which cell type to be plotted.
#' @param cor_cutoff the threshold of weighted correlation for linkages to be plotted.
#' @param gene_color the color to depict genes in the network. When a single value is provided, all genes are plotted in the same color.
#' If a vector of length of genes is provided, each gene will be plotted with the provided color.
#'
#' @return An interactive circle network plot made by edgebundleR to visualize the weighted co-expression of genes based on deconvolved expression.
#'
#' @export get_network

get_network = function(alpha, W, cell_type, cor_cutoff = 0.9, gene_color = "#CD2836") {

  wcor_exp = array(NA, dim = c(nrow(alpha), nrow(alpha), ncol(alpha)))
  for(k in 1:ncol(alpha)) {
    frac_mean = rowMeans(W[,,k], na.rm = T)
    wcor_exp[,,k] = cov.wt(t(alpha[,k,]), frac_mean, cor = T)$cor
  }

  rownames(wcor_exp) = colnames(wcor_exp) = rownames(alpha)
  g = graph_from_adjacency_matrix(wcor_exp[,,cell_type] > cor_cutoff, mode = 'undirected')
  V(g)$color = gene_color
  eb = edgebundle(g)
  print(eb)
}


#' A data example
#'
#' A data list for demonstration.
#'
#'
#' @name example
#' @docType data
#' @return A list containing
#' \item{X}{bulk gene expression (gene x individual x measure).}
#' \item{W}{subject-specific cell type fraction (individual x measure x cell type).}
#' @examples
#'
#' data(example)
#'
NULL
