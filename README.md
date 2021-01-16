bMIND: Bayesian estimation of cell-type-specific (CTS) gene expression for each tissue sample and CTS differential expression analysis
===============================================================

![](man/bMIND.png)

`bMIND` is a Bayesian deconvolution method to integrate bulk and scRNA-seq data. With a prior derived from scRNA-seq data, we estimate sample-level cell-type-specific (CTS) expression from bulk tissue expression via MCMC.

## Installation

Installation requires the `devtools` package.

``` r
devtools::install_github('randel/MIND')
```
## Example

<!-- end list -->

``` r
library(MIND)
data(example)
bulk = t(na.omit(apply(example$X, 1, as.vector)))
frac = na.omit(apply(example$W, 3, as.vector))
colnames(bulk) = rownames(frac) = 1:nrow(frac)
y = rbinom(n = nrow(frac), size = 1, prob = 0.5)
colnames(frac) = gsub(' ', '.', colnames(frac))
covariate = data.frame(c1 = rnorm(length(y)), c2 = rnorm(length(y)))

deconv = bMIND2(bulk, frac = frac, y = y, covariate = covariate, covariate_bulk = 'c1', covariate_cts = 'c2', ncore = 12, np = T, noRE = F)
```

For detailed tutorial, please see the [PDF
manual](https://github.com/randel/MIND/blob/master/MIND-manual.pdf).

The cell type fraction can be pre-estimated using 1) non-negative least squares (NNLS), which requires a
signature matrix derived from reference samples of single-cell RNA-seq data; 2) Bisque, which requires raw single-cell data.

## Reference

**bMIND**: Wang, Jiebiao, Kathryn Roeder, and Bernie Devlin. "Bayesian estimation of cell-type-specific gene expression per bulk sample with prior derived from single-cell data." [bioRxiv](https://doi.org/10.1101/2020.08.05.238949) (2020).

MIND (frequentist method): see https://github.com/randel/MIND/blob/master/MIND.md