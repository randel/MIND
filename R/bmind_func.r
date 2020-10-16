#' The bMIND algorithm to estimate sample-level cell-type-specific expression and conduct CTS differential expression (DE) analysis
#'
#' It calculates the Bayesian estimates of sample- and cell-type-specific (CTS) gene expression, via MCMC. For all input, dim names are recommended if applicable.
#'
#' @param bulk bulk gene expression (gene x sample). We recommend log2-transformed data for better performance, except when using Bisque to estimate 
#' cell type fractions, raw count is expected for Bisque. If the max(bulk) > 50, bulk will be transformed to log2(count per million + 1) 
#' before running bMIND.
#' @param frac sample-specific cell type fraction (sample x cell type). If not specified (NULL), it will be estimated by non-negative least squares (NNLS) by 
#' providing signature matrix or Bisque by providing single-cell reference.
#' @param frac_method method to be used for estimating cell type fractions, either 'NNLS' or 'Bisque'. 
#' @param sc_count sc/snRNA-seq raw count as reference for Bisque to estimate cell type fractions.
#' @param sc_meta meta data frame for sc/snRNA-seq reference. A binary (0-1) column of 'case' is expected to indicate case/control status.
#' @param signature signature matrix for NNLS to estimate cell type fractions. Log2 transformation is recommended.
#' @param signature_case signature matrix from case samples for NNLS to estimate cell type fractions. Log2 transformation is recommended. If this is 
#' provided, signature will be treated as signature matrix for unaffected controls.
#' @param case_bulk case/control status vector for bulk data when using case/control reference to estimate the cell type fractions for case/control subjects
#' separately.
#' @param sample_id sample/subject ID vector. The default is that sample ID will be automatically provided for sample-level bMIND analysis, otherwise 
#' subject ID should be provided for subject-level bMIND analysis. Note that the subject ID will be sorted in the output and different sample_id would 
#' produce slightly different results in MCMCglmm.
#' @param ncore number of cores to run in parallel for providing sample/subject-level CTS estimates. The default is all available cores.
#' @param profile prior profile matrix (gene by cell type). Gene names should be in the same order of bulk, and cell type names should be in the same order
#' as frac. If not specified (NULL), the bulk mean will be supplied.
#' @param covariance prior covariance array (gene by cell type by cell type). Gene names should be in the same order of bulk, and cell type names should be 
#' in the same order as frac. If not specified (NULL), bulk variance / sum(colMeans(frac)^2) will be supplied.
#' @param nu hyper-parameter for the prior covariance matrix. The larger the nu, the higher the certainty about the information in covariance, and the more 
#' informative is the distribution. The default is 50.
#' @param V_fe hyper-parameter for the covariance matrix of fixed-effects. The default is 0.5 * Identity matrix.
#' @param nitt number of MCMC iterations.
#' @param thin thinning interval for MCMC.
#' @param burnin burn-in iterations for MCMC.
#' @param y binary (0-1) outcome/phenotype vector for CTS DE analysis. It can also be a factor with two levels. Should be the same 
#' length and order as sample_id or sort(unique(sample_id)) and row names of covariate.
#' @param covariate matrix for covariates to be adjusted in CTS differential testing.
#' 
#' @return A list containing the output of the bMIND algorithm (some genes with error message in MCMCglmm will not be outputted, 
#' e.g., with constant expression)
#' \item{A}{the deconvolved cell-type-specific gene expression (gene x cell type x sample).}
#' \item{Sigma_c}{the covariance matrix for the deconvolved cell-type-specific expression (gene x cell type x cell type).}
#' \item{mu}{the estimated profile matrix (gene x cell type).}
#' \item{frac}{the estimated cell type fractions (sample x cell type).}
#' \item{pval}{the p-values of CTS-DE testing (cell type x gene).}
#' \item{qval}{the q-values of CTS-DE testing by MANOVA and BH FDR adjustment (cell type x gene).}
#'
#' @references Wang, Jiebiao, Kathryn Roeder, and Bernie Devlin. "Bayesian estimation of cell-type-specific gene expression per bulk sample with prior derived from single-cell data." 
#' bioRxiv (2020).
#'
#' @examples
#'
#' data(example)
#' bulk = t(na.omit(apply(example$X, 1, as.vector)))
#' frac = na.omit(apply(example$W, 3, as.vector))
#' colnames(bulk) = rownames(frac) = 1:nrow(frac)
#'
#'# with provided cell type fractions
#' deconv1 = bMIND(bulk, frac = frac, y = rbinom(n = nrow(frac), size = 1, prob = .5), ncore = 2)
#' 
#' set.seed(1)
#' data(signature)
#' bulk = matrix(rnorm(300*ncol(bulk), 10), ncol = ncol(bulk))
#' rownames(bulk) = rownames(signature)[1:nrow(bulk)]
#' colnames(bulk) = 1:ncol(bulk)
#' 
#' # without provided cell type fractions
#' deconv2 = bMIND(bulk, signature = signature[,-6], y = rbinom(n = nrow(frac), size = 1, prob = .5), 
#' ncore = 2)
#'
#' @export bMIND
#' 
bMIND = function(bulk, frac = NULL, frac_method = NULL, sc_count = NULL, sc_meta = NULL, signature = NULL, signature_case = NULL, case_bulk = NULL,
                 sample_id = NULL, ncore = NULL, profile = NULL, covariance = NULL, nu = 50, V_fe = NULL, nitt = 1300, burnin = 300, thin = 1, 
                 y = NULL, covariate = NULL) {
  
  # check if bulk has genes with constant expression, exclude them, together with those genes in profile and covariance
  
  # estimate cell type fractions
  if(is.null(frac)) est_frac = TRUE else est_frac = FALSE
  if(est_frac) frac = est_frac_sc(bulk, sc_count, signature, signature_case, frac_method, case_bulk, sc_meta)
  
  if(is.null(ncore)) ncore = detectCores()
  if(is.null(sample_id)) sample_id = rownames(frac)
  if(max(bulk) > 50) bulk = log2(apply(bulk, 2, function(x) x/sum(x)*1e6) + 1)
  cts_est = bmind(bulk, frac, sample_id, ncore, profile, covariance, nu, nitt, burnin, thin, V_fe)
  if(est_frac) cts_est$frac = frac
  
  if(is.null(y)) return(cts_est) else {
    bmind_test = test(cts_est$A, y, covariate)
    return((c(cts_est, bmind_test)))
  }
  
}

# bmind with input of 2-dimensional data: X (gene x sample), W (sample x cell type); sample should be in the same order

bmind = function(X, W, sample_id, ncore = 30, profile = NULL, covariance = NULL, nu = 50, nitt = 1300, burnin = 300, thin = 1, V_fe = NULL) {
  
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  getDoParWorkers()
  
  K = ncol(W)
  
  if(is.null(profile)) profile = matrix(1, nrow(X), K) * apply(X, 1, mean)
  if(is.null(V_fe)) V_fe = diag(.5, ncol(profile))
  colnames(W) = gsub(' ', '.', colnames(W))
  
  if(is.null(covariance)) {
    covariance = array(NA, dim = c(nrow(X), K, K))
    for(i in 1:nrow(X)) covariance[i,,] = diag(K) * var(X[i,]) / sum(colMeans(W)^2)
    # summary(covI[1:nrow(X),1,1])
  }
  
  if(is.null(rownames(X))) rownames(X) = 1:nrow(X)
  if(is.null(rownames(profile))) rownames(profile) = rownames(X)
  if(is.null(rownames(covariance))) rownames(covariance) = rownames(X)
  if(is.null(colnames(profile))) colnames(profile) = colnames(W)
  if(is.null(colnames(covariance))) colnames(covariance) = colnames(W)
  if(any(rownames(X) != rownames(profile)) | any(rownames(profile) != rownames(covariance))) 
    print(('Warning: check gene names of bulk data and prior'))
  if(any(colnames(W) != colnames(profile)) | any(colnames(profile) != colnames(covariance))) 
    print(('Warning: check cell type names of fraction and prior'))
  
  mind_mc_ls = foreach(i = rownames(X), .errorhandling = 'pass') %dopar% {
    return(lme_mc2(x = X[i,], W = W, sample_id, mu = profile[i,], V_fe = V_fe, V_re = covariance[i,,], nu = nu, nitt = nitt, burnin = burnin, 
                   thin = thin))
  }
  names(mind_mc_ls) = rownames(X)
  nerr_id = which(sapply(mind_mc_ls, length) != 2)
  err_id = which(sapply(mind_mc_ls, length) == 2)
  if(length(err_id) > 0) {
    print(paste(length(err_id), 'errors'))
    print(str(unique(mind_mc_ls[err_id])))
  }
  res = allto1(mind_mc_ls[nerr_id])
  # rownames(res$A) = rownames(X)[nerr_id]
  # dimnames(res$A)[[3]] = unique(sample_id)
  # colnames(res$A) = colnames(res$mu) = colnames(W)
  # res$lme = mind_mc_ls[[1]]$lme
  res$A[res$A < min(X)] = min(X)
  res$A[res$A > max(X)] = max(X)
  
  if(nrow(W) == length(sample_id)) res$A = res$A[,,unique(sample_id)]
  
  stopCluster(cl)
  return(res)
}


# for one gene for 2D data
lme_mc2 = function(x, W, sample_id, mu, V_fe, V_re, nu = 50, nitt = 1300, burnin = 300, thin = 1) {
  
  K = ncol(W)
  
  # if(is.null(colnames(W))) cell = paste0('cell', 1:K) else 
  cell = colnames(W)
  random = as.formula(paste('~us(', paste(cell, collapse = '+'), '):sample_id'))
  
  miss = which(is.na(x))
  if(length(miss) > 0) {
    x = x[-miss]
    W = W[-miss,]
    sample_id = sample_id[-miss]
  }
  
  set.seed(1)
  fe_formula = as.formula(paste('x ~ -1 + ', paste(cell, collapse = '+')))
  lme2 = MCMCglmm(fe_formula, random, data = data.frame(W, sample_id, x), 
                  verbose = F, pr = T, prior = list(B = list(mu = mu, V = V_fe), 
                                                    G = list(G1 = list(nu = nu, V = V_re))), 
                  nitt = nitt, burnin = burnin, thin = thin) # check the order of subjects' output
  
  N = length(unique(sample_id))
  
  # RE: subject x cell
  re2 = matrix(colMeans(lme2$Sol)[-(1:K)], ncol = K)
  rownames(re2) = sapply(matrix(colnames(lme2$Sol)[-(1:K)], ncol = K)[,1], function(x) unlist(strsplit(x, '[.]'))[3])
  colnames(re2) = cell
  
  mu = colMeans(lme2$Sol)[1:K]
  names(mu) = cell
  
  # Sigma_c
  D2 = matrix(colMeans(lme2$VCV)[-ncol(lme2$VCV)], K, K)
  
  sigma2_e = colMeans(lme2$VCV)[ncol(lme2$VCV)]
  rownames(D2) = colnames(D2) = cell
  
  # 3d array for CTS estimates: sample x cell x Bayesian iterations
  cts_est1 = array(NA, dim = c(N, K, nitt - burnin))
  for(k in 1:K) cts_est1[,k,] = t(lme2$Sol[,k] + lme2$Sol[,K+N*(k-1)+(1:N)])
  se = apply(cts_est1, 2:1, sd) # cell x sample, as A/re2
  
  return(list(A = t(re2) + mu, sigma2_e = sigma2_e, Sigma_c = D2, mu = mu, se = se)) # , lme = lme2
}


# convert a list of results for each gene to 3D array
allto1 = function(mind1) {
  P = length(mind1)
  K = ncol(mind1[[1]]$Sigma_c)
  N = dim(mind1[[1]]$A)[2]
  
  mu = t(sapply(mind1, function(x) x$mu))
  
  deconv1_A = array(NA, dim = c(P, nrow(mind1[[1]]$A), N))
  deconv1_cov = array(NA, dim = c(P, K, K))
  rownames(deconv1_A) = rownames(deconv1_cov) = rownames(mu) = names(mind1)
  colnames(deconv1_A) = colnames(deconv1_cov) = dimnames(deconv1_cov)[[3]] = rownames(mind1[[1]]$A)
  dimnames(deconv1_A)[[3]] = colnames(mind1[[1]]$A)
  SE = deconv1_A
  for(i in names(mind1)) {
    deconv1_A[i,,] = mind1[[i]]$A
    SE[i,,] = mind1[[i]]$se
    deconv1_cov[i,,] = mind1[[i]]$Sigma_c
  }
  return(list(A = deconv1_A, Sigma_c = deconv1_cov, mu = mu, SE = SE))
}


est_frac_sc = function(bulk, sc_count = NULL, signature = NULL, signature_case = NULL, frac_method, case_bulk = NULL, sc_meta = NULL) {
  
  if(!is.null(signature)) {
    frac_method = 'NNLS'
    signature = signature[rownames(signature) %in% rownames(bulk),]
  }
  
  # NNLS
  if(frac_method == 'NNLS') {
    
    if(max(bulk) > 50) bulk = log2(apply(bulk, 2, function(x) x/sum(x)*1e6) + 1)
    if(max(signature) > 50) signature = log2(signature + 1)
    
    if(is.null(signature_case)) frac = est_frac(signature, bulk) else {
      if(max(signature_case) > 50) signature_case = log2(signature_case + 1)
      
      frac0 = est_frac(signature, bulk[, case_bulk == 0])
      frac1 = est_frac(signature_case, bulk[, case_bulk == 1])
      frac = rbind(frac0, frac1)[colnames(bulk),]
    }
  }
  
  # Bisque
  if(frac_method == 'Bisque') {
    
    package.check("BisqueRNA")
    package.check("Biobase")
    bulk_eset = ExpressionSet(assayData = bulk)
    # (Expects read counts for both datasets, as they will be converted to counts per million (CPM))
    sc_eset = ExpressionSet(assayData = as.matrix(sc_count), 
                            phenoData = methods::new("AnnotatedDataFrame", data = sc_meta, 
                                                     varMetadata = data.frame(labelDescription = colnames(sc_meta), 
                                                                              row.names = colnames(sc_meta))))
    
    if(is.null(sc_meta$case)) frac = t(ReferenceBasedDecomposition(bulk_eset, sc_eset, cell.types = 'cell_type',
                                                                   subject.names = 'subject', use.overlap = F)$bulk.props) 
    if(!is.null(sc_meta$case)) {
      frac0 = t(ReferenceBasedDecomposition(bulk_eset[, case_bulk == 0], sc_eset[, sc_meta$case == 0], 
                                            cell.types = 'cell_type', subject.names = 'subject', use.overlap = F)$bulk.props)
      frac1 = t(ReferenceBasedDecomposition(bulk_eset[, case_bulk == 1], sc_eset[, sc_meta$case == 1], 
                                            cell.types = 'cell_type', subject.names = 'subject', use.overlap = F)$bulk.props)
      frac = rbind(frac0, frac1)[colnames(bulk),]
    }
  }
  
  return(frac)
}


############################### association testing

meta_sample2sub = function(meta_sample, sub_id) {
  N = length(unique(meta_sample[,sub_id]))
  col_keep = sub_id
  for(i in (1:ncol(meta_sample))[-sub_id]) {
    if(length(unique(paste(meta_sample[,i], meta_sample[,sub_id]))) <= N) col_keep = c(col_keep, i)
  }
  meta_sub = unique(meta_sample[,col_keep])
  rownames(meta_sub) = meta_sub[,1]
  meta_sub = meta_sub[,-1]
  meta_sub = meta_sub[sort(rownames(meta_sub)),]
  meta_sub = meta_sub[, apply(meta_sub, 2, function(x) length(unique(x))) != 1]
  meta_sub = meta_sub[, colMeans(is.na(meta_sub)) != 1]
  return(meta_sub)
}


# pval for one gene in case some cell types have missing pval
get_pval = function(pval, cell_type, K) {
  pval0 = rep(NA, K)
  names(pval0) = cell_type
  names = intersect(names(pval), cell_type)
  pval0[names] = pval[names]
  return(pval0)
}


test = function(A, y, covariate = NULL) {
  
  if(dim(A)[3] != length(y)) print('CTS estimates and y have different length')
  if(!is.null(covariate)) if(dim(A)[3] != nrow(covariate)) print('CTS estimates and covariate have different number of samples/subjects') else {
    if(!is.null(rownames(covariate)) & any(rownames(covariate) != dimnames(A)[[3]])) covariate = covariate[dimnames(A)[[3]],]
  }
  
  K = ncol(A)
  cell_type = colnames(A)
  if(is.null(covariate)) pval = apply(A, 1, function(x) {
    pval = coef(summary(glm(y ~ ., data = data.frame(t(x)), family = 'binomial')))[,4]
    return(get_pval(pval, cell_type, K))
  }) else
    pval = apply(A, 1, function(x) {
      pval = coef(summary(glm(y ~ ., data = data.frame(t(x), covariate), family = 'binomial')))[,4]
      return(get_pval(pval, cell_type, K))
    })
  
  qval = pval2qval(pval, A, y, covariate)
  # rownames(qval) = rownames(pval) = substring(rownames(pval), 5)
  return(list(qval = qval, pval = pval))
}


# MANOVA; pval: K x ngene
pval2qval = function(pval, A, y, covariate = NULL) {
  ng = nrow(A)
  # pval for each gene
  if(is.null(covariate)) pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y))$stats[1, "Pr(>F)"], silent = T)) else 
    pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y + covariate))$stats[1, "Pr(>F)"], silent = T))
  pval = pval[,!is.na(as.numeric(pval1))]
  pval1 = na.omit(as.numeric(pval1))
  qval1 = p.adjust(pval1, 'fdr')
  # hist(pval1)
  # print(min(qval1))
  qval = pval
  K = ncol(A)
  for(i in 1:ncol(pval)) {
    qval[,i] = 1
    if(min(pval[,i], na.rm = T) < .05/K) qval[,i][which.min(pval[,i])] = qval1[i]
  }
  return(qval)
}

                  

########################################################################## get prior CTS profile and covariance matrix from single-cell data
## input
# sc: single-cell count matrix, gene x cell
# meta_sc: data.frame for meta of cells (cell x features, including sample (ID), cell_type)

## output
# profile: CTS profile matrix (gene x cell type), in log2(CPM + 1) scale
# cov: CTS covariance matrix (gene x cell type x cell type)
                   
get_prior = function(sc, meta_sc) {
  
  meta_sc$sample = as.character(meta_sc$sample)
  sample = unique(meta_sc[, c('sample')])
  
  library(data.table)
  library(edgeR)
  
  cell_type = sort(unique(meta_sc$cell_type))
  K = length(cell_type)
  
  cts = array(NA, dim = c(nrow(sc), length(sample), K))
  rownames(cts) = rownames(sc)
  colnames(cts) = sample
  dimnames(cts)[[3]] = cell_type
  for(j in sample) {
    for(k in dimnames(cts)[[3]]) cts[,j,k] = rowMeans(sc[, meta_sc$sample == j & meta_sc$cell_type == k])
    id = which(colMeans(is.na(cts[,j,])) == 0)
    cts[,j,id] = log2(cpm(cts[,j,id]) + 1) # make it log2 CPM + 1
  }
  
  library(matrixcalc)
  cov = array(NA, dim = c(nrow(cts), K, K))
  rownames(cov) = rownames(sc)
  colnames(cov) = dimnames(cov)[[3]] = cell_type
  for(i in 1:nrow(cov)) {
    cov[i,,] = cov(cts[i,,], use = 'pairwise')
  }
  
  cov <- cov[apply(cov, 1, is.positive.definite),,]
  
  profile = matrix(NA, nrow(sc), K)
  rownames(profile) = rownames(sc)
  colnames(profile) = cell_type
  for(i in cell_type) {
    profile[,i] = log2(cpm(rowMeans(sc[, meta_sc$cell_type == i])) + 1)
  }
  
  return(list(cov = cov, ctsExp = cts, profile = profile))
}
