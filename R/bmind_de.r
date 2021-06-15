#' The bMIND algorithm that considers Bayesian testing and covariates in the deconvolution model
#'
#' It calculates the Bayesian estimates of sample- and cell-type-specific (CTS) gene expression, via MCMC.
#'
#' @param bulk bulk gene expression (gene x sample).
#' @param frac sample-specific cell type fraction (sample x cell type). If not specified (NULL), it will be estimated by non-negative least squares (NNLS) by 
#' providing signature matrix or Bisque by providing single-cell reference.
#' @param sample_id sample/subject ID vector. The default is that sample ID will be automatically provided for sample-level bMIND analysis, otherwise 
#' subject ID should be provided for subject-level bMIND analysis. Note that the subject ID will be sorted in the output and different sample_id would 
#' produce slightly different results in MCMCglmm.
#' @param ncore number of cores to run in parallel for providing sample/subject-level CTS estimates. The default is all available cores.
#' @param profile prior profile matrix (gene by cell type). Gene names should be in the same order of bulk, and cell type names should be in the same order
#' as frac.
#' @param profile_co prior profile matrix (gene by cell type) for controls. 
#' @param profile_ca prior profile matrix (gene by cell type) for cases. 
#' @param covariance prior covariance array (gene by cell type by cell type). Gene names should be in the same order of bulk, and cell type names should be 
#' in the same order as frac. The default is 0.5 * Identity matrix for covariance of fixed effects.
#' @param covariance_co prior covariance array (gene by cell type by cell type) for controls. 
#' @param covariance_ca prior covariance array (gene by cell type by cell type) for cases. 
#' @param nu hyper-parameter for the prior covariance matrix. The larger the nu, the higher the certainty about the information in covariance, and the more 
#' informative is the distribution. The default is 50.
#' @param nitt number of MCMC iterations.
#' @param thin thinning interval for MCMC.
#' @param burnin burn-in iterations for MCMC.
#' @param y binary (0-1) outcome/phenotype vector for CTS DE analysis (0 for controls, 1 for cases). Should be the same 
#' length and order as sample_id or sort(unique(sample_id)) and row names of covariate.
#' @param covariate matrix for covariates to be adjusted in deconvolution model.
#' @param covariate_bulk colnames of covariate denoting variables that affect bulk expression
#' @param covariate_cts colnames of covariate denoting variables that affect CTS expression
#' @param np option to use non-informative prior
#' @param noRE option to not calculate sample-level CTS estimates
#' @param max_samp max number of posterior samples to generate in testing. An adaptive procedure is used to increase nitt for those genes with p-values = 
#' 1/number of posterior samples.
#' @param frac_method method to be used for estimating cell type fractions, either 'NNLS' or 'Bisque'. 
#' **All arguments starting from this one will be used to estimate cell-type fractions only, if those fractions are not pre-estimated.**
#' @param sc_count sc/snRNA-seq raw count as reference for Bisque to estimate cell type fractions.
#' @param sc_meta meta data frame for sc/snRNA-seq reference. A binary (0-1) column of 'case' is expected to indicate case/control status.
#' @param signature signature matrix for NNLS to estimate cell type fractions. Log2 transformation is recommended.
#' @param signature_case signature matrix from case samples for NNLS to estimate cell type fractions. Log2 transformation is recommended. If this is 
#' provided, signature will be treated as signature matrix for unaffected controls.
#' @param case_bulk case/control status vector for bulk data when using case/control reference to estimate the cell type fractions for case/control subjects
#' separately.
#' 
#' @return A list containing the output of the bMIND algorithm (some genes with error message in MCMCglmm will not be outputted, 
#' e.g., with constant expression)
#' \item{A}{the deconvolved cell-type-specific gene expression (gene x cell type x sample).}
#' \item{SE}{the standard error of cell-type-specific gene expression (gene x cell type x sample).}
#' \item{coef}{the estimated coefficients matrix (gene x variables).}
#' \item{frac}{the estimated cell type fractions (sample x cell type).}
#' \item{pval}{the p-values of CTS-DE testing (gene x cell type).}
#' \item{qval}{the q-values of CTS-DE testing by BH FDR adjustment (gene x cell type).}
#'
#' @references Wang, Jiebiao, Kathryn Roeder, and Bernie Devlin. "Bayesian estimation of cell-type-specific gene expression per bulk sample with prior 
#' derived from single-cell data." bioRxiv (2020).
#' @export bMIND2
bMIND2 = function(bulk, frac = NULL, sample_id = NULL, ncore = NULL, profile = NULL, covariance = NULL, 
                  profile_co = NULL, covariance_co = NULL, profile_ca = NULL, covariance_ca = NULL, y = NULL, covariate = NULL, 
                  covariate_bulk = NULL, covariate_cts= NULL, noRE = T, np = F, 
                  nu = 50, nitt = 1300, burnin = 300, thin = 1, max_samp = 1e6,
                  frac_method = NULL, sc_count = NULL, sc_meta = NULL, signature = NULL, signature_case = NULL, case_bulk = NULL) {
  
  # check if bulk has genes with constant expression, exclude them, together with those constant genes in profile and covariance
  
  # estimate cell type fractions
  if(is.null(frac)) est_frac = TRUE else est_frac = FALSE
  if(est_frac) frac = est_frac_sc(bulk, sc_count, signature, signature_case, frac_method, case_bulk, sc_meta)
  
  if(is.null(ncore)) ncore = detectCores()
  if(is.null(y)) cts_est = bmind_all(X = bulk, W = frac, sample_id = sample_id, ncore = ncore, mu = profile, var_fe = covariance, 
                                     covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts, np = np, noRE = noRE,
                                     nu = nu, nitt = nitt, burnin = burnin, thin = thin) else
                                       cts_est = bmind_all(X = bulk, W = frac, sample_id = sample_id, ncore = ncore, mu = cbind(profile_co, profile_ca), 
                                                           var_fe_co = covariance_co, var_fe_ca = covariance_ca, y = y, max_samp = max_samp,
                                                           covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts, 
                                                           np = np, noRE = noRE,
                                                           nu = nu, nitt = nitt, burnin = burnin, thin = thin)
                                     if(est_frac) cts_est$frac = frac
                                     
                                     return(cts_est)
}

# bmind with input of 2-dimensional data: X (gene x sample), W (sample x cell type); sample should be in the same order
# with covariate*y interaction

bmind_all = function(X, W, y = NULL, mu = NULL, var_fe = NULL, var_fe_co = NULL, var_fe_ca = NULL, max_samp = 1e6, 
                     covariate = NULL, covariate_bulk = NULL, covariate_cts = NULL,
                     V_re = NULL, sample_id = NULL, nu = 50, nitt = 1300, burnin = 300, thin = 1, noRE = T, np = F, ncore) {
  
  K = ncol(W)
  colnames(W)[1] = gsub("[^0-9A-Za-z///' ]", "", colnames(W)[1], ignore.case = TRUE)
  colnames(W)[1] = gsub(' ', '_', colnames(W)[1])
  
  if(is.null(rownames(X))) rownames(X) = 1:nrow(X)
  # if(is.null(rownames(mu))) rownames(mu) = rownames(X)
  if(is.null(var_fe)) {
    var_fe = array(NA, dim = c(nrow(X), K, K))
    for(i in 1:nrow(X)) var_fe[i,,] = diag(.5, K)
  }
  if(is.null(var_fe_co)) {
    var_fe_co = array(NA, dim = c(nrow(X), K, K))
    for(i in 1:nrow(X)) var_fe_co[i,,] = diag(.5, K)
  }
  if(is.null(var_fe_ca)) {
    var_fe_ca = array(NA, dim = c(nrow(X), K, K))
    for(i in 1:nrow(X)) var_fe_ca[i,,] = diag(.5, K)
  }
  
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  getDoParWorkers()
  ng = nrow(X)
  
  if(is.null(y)) {
    res = foreach(j = 1:ng, .errorhandling = 'pass') %dopar% {
      
      if(np) bmind1(X[j,], W, np = np, nu = nu, nitt = nitt, burnin = burnin, thin = thin, noRE = noRE,
                    covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts) else
                      bmind1(X[j,], W, mu[j,], var_fe[j,,], np = np, nu = nu, nitt = nitt, burnin = burnin, thin = thin, noRE = noRE,
                             covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts)
    }
    names(res) = rownames(X)
    nerr_id = which(sapply(res, length) != 2)
    err_id = which(sapply(res, length) == 2)
    if(length(err_id) > 0) {
      print(paste(length(err_id), 'errors'))
      print(str(unique(res[err_id])))
    }
    res = res[nerr_id]
    stopCluster(cl)
    
    coef = do.call(rbind, lapply(res, function(x) x$coef))
    rownames(coef) = names(res)
    
    if(noRE) return(list(coef = coef)) else {
      res1 = get_A_se(res, X)
      return(list(coef = coef, A = res1$A, SE = res1$SE))
    }
    
  } else {
    res = foreach(j = 1:ng, .errorhandling = 'pass') %dopar% {
      
      if(np) bmind1_y(X[j,], W, y, max_samp = max_samp, np = np, nu = nu, noRE = noRE, 
                      covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts, nitt = nitt, burnin = burnin, thin = thin) else
                        bmind1_y(X[j,], W, y, mu[j,], var_fe_co[j,,], var_fe_ca[j,,], max_samp = max_samp, np = np, nu = nu, noRE = noRE, 
                                 covariate = covariate, covariate_bulk = covariate_bulk, covariate_cts = covariate_cts, 
                                 nitt = nitt, burnin = burnin, thin = thin)
    }
    stopCluster(cl)
    names(res) = rownames(X)
    nerr_id = which(sapply(res, length) != 2)
    err_id = which(sapply(res, length) == 2)
    if(length(err_id) > 0) {
      print(paste(length(err_id), 'errors'))
      print(str(unique(res[err_id])))
    }
    res = res[nerr_id]
    
    coef = do.call(rbind, lapply(res, function(x) x$coef))
    pval = do.call(rbind, lapply(res, function(x) x$pval))
    qval = apply(pval, 2, function(x) p.adjust(x, method = 'fdr'))
    rownames(coef) = rownames(pval) = rownames(qval) = names(res)
    
    if(noRE) return(list(qval = qval, coef = coef, pval = pval)) else {
      res1 = get_A_se(res, X)
      return(list(qval = qval, coef = coef, pval = pval, A = res1$A, SE = res1$SE))
    }
  }
}


# convert a list of results for each gene to 3D array
get_A_se = function(mind1, X) {
  P = length(mind1)
  N = dim(mind1[[1]]$A)[2]
  
  deconv1_A = array(NA, dim = c(P, nrow(mind1[[1]]$A), N))
  rownames(deconv1_A) = names(mind1)
  colnames(deconv1_A) = rownames(mind1[[1]]$A)
  dimnames(deconv1_A)[[3]] = colnames(mind1[[1]]$A)
  SE = deconv1_A
  for(i in names(mind1)) {
    deconv1_A[i,,] = mind1[[i]]$A
    SE[i,,] = mind1[[i]]$se
  }
  
  deconv1_A[deconv1_A < min(X)] = min(X)
  deconv1_A[deconv1_A > max(X)] = max(X)
  
  return(list(A = deconv1_A, SE = SE))
}



# with y and case-control priors
bmind1_y = function(x, W, y = NULL, mu = NULL, var_fe_co, var_fe_ca, max_samp = 1e6, covariate = NULL, covariate_bulk = NULL, covariate_cts = NULL,
                    V_re = NULL, sample_id = NULL, nu = 50, nitt = 1300, burnin = 300, thin = 1, noRE = T, np = F) {
  
  if(is.null(sample_id)) sample_id = rownames(W)
  if(!is.null(y)) {
    y = as.numeric(y)
    co = as.numeric(y == 0)
    ca = as.numeric(y == 1)
  }
  K = ncol(W)
  
  cell = colnames(W)
  random = as.formula(paste('~us(', paste(cell, collapse = '+'), '):sample_id'))
  if(is.null(V_re) & np == F) V_re = (var_fe_co + var_fe_co) / 2
  
  miss = which(is.na(x))
  if(length(miss) > 0) {
    x = x[-miss]
    W = W[-miss,]
    sample_id = sample_id[-miss]
  }
  
  set.seed(1)
  
  if(!is.null(y) & is.null(covariate)) {
    df = data.frame(W, sample_id, x, co, ca)
    fe_formula = as.formula(paste('x ~ -1 +', paste(c(paste0(cell, ':co'), paste0(cell, ':ca')), collapse = '+')))
    
    n_fe = 2*K
    if(!np) V_fe = bdiag(var_fe_co, var_fe_ca) # assumes block diagonal
  }
  if(!is.null(y) & !is.null(covariate)) {
    df = data.frame(W, sample_id, x, co, ca, covariate)
    variables = c(paste0(cell, ':co'), paste0(cell, ':ca'))
    if(!is.null(covariate_cts)) variables = c(as.vector(sapply(covariate_cts, function(x) paste0(cell, ':', x))), variables)
    if(!is.null(covariate_bulk)) variables = c(covariate_bulk, variables)
    fe_formula = as.formula(paste('x ~ -1 +', paste(variables, collapse = '+')))
    
    n_fe = (length(covariate_cts) + 2) * K + length(covariate_bulk)
    if(!np) {
      V_fe = bdiag(diag(n_fe - 2*K) * 1e10, var_fe_co, var_fe_ca)
      mu = c(rep(0, n_fe - 2*K), mu)
    }
  }
  
  # number of posterior samples/Bayesian iterations
  nsamp = (nitt - burnin)/thin
  repeat {
    
    if(noRE) {
      
      if(np) lme = MCMCglmm(fe_formula, data = df, verbose = F, nitt = nitt, burnin = burnin, thin = thin) else 
        lme = MCMCglmm(fe_formula, data = df, verbose = F, prior = list(B = list(mu = mu, V = V_fe)), nitt = nitt, burnin = burnin, thin = thin)
      
    } else {
      
      if(np) lme = MCMCglmm(fe_formula, random, data = df, verbose = F, pr = T, nitt = nitt, burnin = burnin, thin = thin) else 
        lme = MCMCglmm(fe_formula, random, data = df, verbose = F, pr = T, prior = list(B = list(mu = mu, V = V_fe), G = list(G1 = list(nu = nu, V = V_re))), 
                       nitt = nitt, burnin = burnin, thin = thin)
    }
    
    for(k in 1:K) {
      if(any(colnames(lme$Sol) == paste0('co', ":", cell[k]))) colnames(lme$Sol)[colnames(lme$Sol) == paste0('co', ":", cell[k])] = paste0(cell[k], ':co')
      if(any(colnames(lme$Sol) == paste0('ca', ":", cell[k]))) colnames(lme$Sol)[colnames(lme$Sol) == paste0('ca', ":", cell[k])] = paste0(cell[k], ':ca')
    }
    pval = apply(lme$Sol[, paste0(cell, ':ca')] - lme$Sol[, paste0(cell, ':co')], 2, function(x) 2*min(mean(x > 0), mean(x < 0)))
    pval[pval == 0] = 1/nsamp
    if(!any(pval == 1/nsamp) | nsamp == max_samp) break
    nsamp = nsamp*10
    nitt = nsamp*thin + burnin
  }
  names(pval) = cell
  
  if(noRE) return(list(pval = pval, coef = (summary(lme)$solutions)[, 1], res_mcmcglmm = lme$Sol)) else {
    N = length(unique(sample_id))
    
    ## lme$Sol's column in the order of cell1.sample1:N ... cellK.sample1:N
    # RE: subject x cell
    # re2 = matrix(colMeans(lme$Sol)[-(1:n_fe)], ncol = K)
    # rownames(re2) = sapply(matrix(colnames(lme$Sol)[-(1:n_fe)], ncol = K)[,1], function(x) unlist(strsplit(x, '[.]'))[3]) # sample name
    # colnames(re2) = cell
    
    # # Sigma_c
    # Sigma_c = matrix(colMeans(lme$VCV)[-ncol(lme$VCV)], K, K)
    # sigma2_e = colMeans(lme$VCV)[ncol(lme$VCV)]
    # rownames(Sigma_c) = colnames(Sigma_c) = cell
    
    # 3d array for CTS estimates: sample x cell x Bayesian iterations (note that sample ID will be sorted by characters)
    cts_est1 = array(NA, dim = c(N, K, nsamp))
    for(k in 1:K) cts_est1[,k,] = t(lme$Sol[,n_fe+N*(k-1)+(1:N)])
    for(i in 1:N) if(y[i] == 0) cts_est1[i,,] = cts_est1[i,,] + t(lme$Sol[, paste0(cell, ':co')]) else 
      cts_est1[i,,] = cts_est1[i,,] + t(lme$Sol[, paste0(cell, ':ca')])
    
    se = apply(cts_est1, 2:1, sd) # cell x sample, as A or t(re2)
    re = apply(cts_est1, 2:1, mean)
    
    sample_name = sapply(colnames(lme$Sol)[n_fe+(1:N)], function(x) unlist(strsplit(x, '[.]'))[3])
    colnames(se) = colnames(re) = sample_name
    rownames(se) = rownames(re) = cell
    
    return(list(pval = pval, coef = (summary(lme)$solutions)[, 1], res_mcmcglmm = lme$Sol, A = re[,sample_id], se = se[,sample_id]))
  }
}


# without y
bmind1 = function(x, W, mu = NULL, var_fe, covariate = NULL, covariate_bulk = NULL, covariate_cts = NULL,
                  V_re = NULL, sample_id = NULL, nu = 50, nitt = 1300, burnin = 300, thin = 1, noRE = T, np = F) {
  
  if(is.null(sample_id)) sample_id = rownames(W)
  K = ncol(W)
  
  cell = colnames(W)
  random = as.formula(paste('~us(', paste(cell, collapse = '+'), '):sample_id'))
  if(is.null(V_re) & np == F) V_re = var_fe
  
  miss = which(is.na(x))
  if(length(miss) > 0) {
    x = x[-miss]
    W = W[-miss,]
    sample_id = sample_id[-miss]
  }
  
  set.seed(1)
  
  if(is.null(covariate)) {
    df = data.frame(W, sample_id, x)
    fe_formula = as.formula(paste('x ~ -1 +', paste(cell, collapse = '+')))
    
    n_fe = K
    if(!np) V_fe = var_fe # assumes block diagonal
  } else {
    df = data.frame(W, sample_id, x, covariate)
    variables = cell
    if(!is.null(covariate_bulk)) variables = c(variables, covariate_bulk)
    if(!is.null(covariate_cts)) variables = c(variables, as.vector(sapply(covariate_cts, function(x) paste0(cell, ':', x))))
    fe_formula = as.formula(paste('x ~ -1 +', paste(variables, collapse = '+')))
    
    n_fe = (length(covariate_cts) + 1) * K + length(covariate_bulk)
    if(!np) {
      V_fe = bdiag(var_fe, diag(n_fe - K) * 1e10) # check parameters' order
      mu = c(mu, rep(0, n_fe - K))
    }
  }
  
  # number of posterior samples/Bayesian iterations
  nsamp = (nitt - burnin)/thin
  if(noRE) {
    
    if(np) lme = MCMCglmm(fe_formula, data = df, verbose = F, nitt = nitt, burnin = burnin, thin = thin) else 
      lme = MCMCglmm(fe_formula, data = df, verbose = F, prior = list(B = list(mu = mu, V = V_fe)), nitt = nitt, burnin = burnin, thin = thin)
    
  } else {
    
    if(np) lme = MCMCglmm(fe_formula, random, data = df, verbose = F, pr = T, nitt = nitt, burnin = burnin, thin = thin) else 
      lme = MCMCglmm(fe_formula, random, data = df, verbose = F, pr = T, prior = list(B = list(mu = mu, V = V_fe), G = list(G1 = list(nu = nu, V = V_re))), 
                     nitt = nitt, burnin = burnin, thin = thin)
  }
  
  if(noRE) return(list(coef = (summary(lme)$solutions)[, 1], res_mcmcglmm = lme$Sol, noRE = noRE)) else {
    N = length(unique(sample_id))
    
    ## lme$Sol's column in the order of cell1.sample1:N ... cellK.sample1:N
    # RE: subject x cell
    # re2 = matrix(colMeans(lme$Sol)[-(1:n_fe)], ncol = K)
    # rownames(re2) = sapply(matrix(colnames(lme$Sol)[-(1:n_fe)], ncol = K)[,1], function(x) unlist(strsplit(x, '[.]'))[3]) # sample name
    # colnames(re2) = cell
    
    # # Sigma_c
    # Sigma_c = matrix(colMeans(lme$VCV)[-ncol(lme$VCV)], K, K)
    # sigma2_e = colMeans(lme$VCV)[ncol(lme$VCV)]
    # rownames(Sigma_c) = colnames(Sigma_c) = cell
    
    # 3d array for CTS estimates: sample x cell x Bayesian iterations (note that sample ID will be sorted by characters)
    cts_est1 = array(NA, dim = c(N, K, nsamp))
    for(k in 1:K) cts_est1[,k,] = t(lme$Sol[,n_fe+N*(k-1)+(1:N)])
    for(i in 1:N) cts_est1[i,,] = cts_est1[i,,] + t(lme$Sol[, cell])
    
    se = apply(cts_est1, 2:1, sd) # cell x sample, as A or t(re2)
    re = apply(cts_est1, 2:1, mean)
    
    sample_name = sapply(colnames(lme$Sol)[n_fe+(1:N)], function(x) unlist(strsplit(x, '[.]'))[3])
    colnames(se) = colnames(re) = sample_name
    rownames(se) = rownames(re) = cell
    
    return(list(coef = (summary(lme)$solutions)[, 1], res_mcmcglmm = lme$Sol, A = re[,sample_id], se = se[,sample_id]))
  }
}
