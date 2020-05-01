
# bmind for 2-dimensional data: X (gene x sample), W (sample x cell)
# id: sample subject ID (note that the subject ID will be sorted in the output); 
# different subject ID order would produce slightly diff results

bmind = function(X, W, id, ncore = 30, profile = NULL, covariance = NULL, nu = 50, nitt = 1300, burnin = 300, thin = 1, V_fe = NULL, 
                 pmedian = F) {
  
  library(doParallel)
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  getDoParWorkers()
  
  K = ncol(W)
  if(is.null(V_fe)) V_fe = diag(.5, ncol(profile))
  colnames(W) = gsub(' ', '.', colnames(W))
  
  if(is.null(profile)) profile = matrix(1, nrow(X), K) * apply(X, 1, mean)
  
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
    print(('Warning: check cell names of fraction and prior'))
  
  library(foreach)
  mind_mc_ls <- foreach(i = rownames(X), .errorhandling = 'pass') %dopar% {
    
    source('https://raw.githubusercontent.com/randel/MIND/master/R/bmind_func.r')
    
    return(lme_mc2(x = X[i,], W = W, id, mu = profile[i,], V_fe = V_fe, V_re = covariance[i,,], nu = nu, nitt = nitt, burnin = burnin, 
                   thin = thin,
                   pmedian = pmedian))
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
  # dimnames(res$A)[[3]] = unique(id)
  # colnames(res$A) = colnames(res$mu) = colnames(W)
  # res$lme = mind_mc_ls[[1]]$lme
  res$A[res$A < min(X)] = min(X)
  res$A[res$A > max(X)] = max(X)
  
  if(nrow(W) == length(id)) res$A = res$A[,,rownames(W)]
  
  return(res)
}


# for one gene for 2D data
lme_mc2 = function(x, W, id, mu, V_fe, V_re, nu = 50, nitt = 1300, burnin = 300, thin = 1, pmedian) {
  
  K = ncol(W)
  
  # if(is.null(colnames(W))) cell = paste0('cell', 1:K) else 
  cell = colnames(W)
  random = as.formula(paste('~us(', paste(cell, collapse = '+'), '):id'))
  
  miss = which(is.na(x))
  if(length(miss) > 0) {
    x = x[-miss]
    W = W[-miss,]
    id = id[-miss]
  }
  
  library(MCMCglmm)
  set.seed(1)
  fe_formula = as.formula(paste('x ~ -1 + ', paste(cell, collapse = '+')))
  lme2 <- MCMCglmm(fe_formula, random, data = data.frame(W, id, x), 
                   verbose = F, pr = T, prior = list(B = list(mu = mu, V = V_fe), 
                                                     G = list(G1 = list(nu = nu, V = V_re))), 
                   nitt = nitt, burnin = burnin, thin = thin) # check the order of subjects' output
  
  N = length(unique(id))
  
  # RE: subject x cell
  if(pmedian) re2 = matrix(apply(lme2$Sol, 2, median)[-(1:K)], ncol = K) else re2 = matrix(colMeans(lme2$Sol)[-(1:K)], ncol = K)
  rownames(re2) = sapply(matrix(colnames(lme2$Sol)[-(1:K)], ncol = K)[,1], function(x) unlist(strsplit(x, '[.]'))[3])
  colnames(re2) = cell
  
  if(pmedian) mu = apply(lme2$Sol, 2, median)[1:K] else mu = colMeans(lme2$Sol)[1:K]
  names(mu) = cell
  
  # Sigma_c
  if(pmedian) D2 = matrix(apply(lme2$VCV, 2, median)[-ncol(lme2$VCV)], K, K) else D2 = matrix(colMeans(lme2$VCV)[-ncol(lme2$VCV)], K, K)
  
  if(pmedian) sigma2_e = apply(lme2$VCV, 2, median)[ncol(lme2$VCV)] else sigma2_e = colMeans(lme2$VCV)[ncol(lme2$VCV)]
  rownames(D2) = colnames(D2) = cell
  
  return(list(A = t(re2) + mu, sigma2_e = sigma2_e, Sigma_c = D2, mu = mu)) # , lme = lme2
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
  for(i in names(mind1)) {
    deconv1_A[i,,] = mind1[[i]]$A
    deconv1_cov[i,,] = mind1[[i]]$Sigma_c
  }
  return(list(A = deconv1_A, Sigma_c = deconv1_cov, mu = mu))
}


est_frac_sc = function(bulk, sc = NULL, sig = NULL, sig_case = NULL, method, case_bulk = NULL, sc_meta = NULL, marker = NULL) {
  
  # NN LS
  if(method == 'NNLS') {
    if(is.null(sig_case)) frac = est_frac(sig, bulk) else {
      frac0 = est_frac(sig, bulk[, case_bulk == 0])
      frac1 = est_frac(sig_case, bulk[, case_bulk == 1])
      frac = rbind(frac0, frac1)[colnames(bulk),]
    }
  }
  
  # Bisque
  if(method == 'Bisque') {
    
    package.check("BisqueRNA")
    package.check("Biobase")
    bulk_eset <- ExpressionSet(assayData = bulk)
    # using only markers? no, all genes (Expects read counts for both datasets, as they will be converted to counts per million (CPM))
    sc_eset <- ExpressionSet(assayData = as.matrix(sc), 
                             phenoData = methods::new("AnnotatedDataFrame", data = sc_meta, 
                                                      varMetadata = data.frame(labelDescription = colnames(sc_meta), 
                                                                               row.names = colnames(sc_meta))))
    
    if(is.null(sc_meta$case)) frac = t(ReferenceBasedDecomposition(bulk_eset, sc_eset, markers = marker, cell.types = 'cell_type',
                                                                   subject.names = 'subject', use.overlap = F)$bulk.props) 
    if(!is.null(sc_meta$case)) {
      frac0 = t(ReferenceBasedDecomposition(bulk_eset[, case_bulk == 0], sc_eset[, sc_meta$case == 0], markers = marker, 
                                            cell.types = 'cell_type', subject.names = 'subject', use.overlap = F)$bulk.props)
      frac1 = t(ReferenceBasedDecomposition(bulk_eset[, case_bulk == 1], sc_eset[, sc_meta$case == 1], markers = marker, 
                                            cell.types = 'cell_type', subject.names = 'subject', use.overlap = F)$bulk.props)
      frac = rbind(frac0, frac1)[colnames(bulk),]
    }
  }
  
  return(frac)
}


# NNLS
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

get_pval = function(pval, cell_type, K) {
  pval0 = rep(NA, K)
  names(pval0) = cell_type
  names = intersect(names(pval), cell_type)
  pval0[names] = pval[names]
  return(pval0)
}

test = function(A, y, covariate = NULL) {
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
