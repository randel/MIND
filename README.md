Estimating sample/subject-level cell-type-specific gene expression
===============================================================

## bMIND: Bayesian estimation of cell-type-specific gene expression for each tissue sample
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
data(signature_brain_darmanis)

set.seed(1)
bulk = matrix(rnorm(300*100, 10), ncol = 100)
rownames(bulk) = rownames(signature)[1:nrow(bulk)]
colnames(bulk) = 1:ncol(bulk)
y = rbinom(n = ncol(bulk), size = 1, prob = .5)

deconv = bMIND(bulk, signature = signature[,-6], y = y)
```

For detailed tutorial, please see the [PDF
manual](https://github.com/randel/MIND/blob/master/MIND-manual.pdf).

The cell type fraction can be pre-estimated using 1) non-negative least squares (NNLS), which requires a
signature matrix derived from reference samples of single-cell RNA-seq data; 2) Bisque, which requires raw single-cell data.

## Reference

**bMIND**: Wang, Jiebiao, Kathryn Roeder, and Bernie Devlin. "Bayesian estimation of sample- and cell-type-specific gene expression, illustrated by differential expression analysis." *Submitted* (2020).



# MIND: using multiple measurements of tissue to estimate subject-and cell-type-specific gene expression
![](man/MIND.png)


## Example

Estimate subject- and cell-type-specific gene expression (saved as
`deconv$A` below) using bulk gene expression data and pre-estimated cell
type fractions:

  - `X`: bulk gene expression (gene x subject x measure)
  - `W`: cell type fraction (subject x measure x cell type)

<!-- end list -->

``` r
library(MIND)

data(example)

deconv = mind(X = example$X, W = example$W)
```

For details, please see the [PDF
manual](https://github.com/randel/MIND/blob/master/MIND-manual.pdf).

The cell type fraction can be pre-estimated using `est_frac()` (based on
non-negative least squares, see an example [here](http://rpubs.com/randel/est_frac)) or another standard deconvolution method. It requires a
signature matrix derived from reference samples of single-cell RNA-seq
data.

## Reference

**MIND**: Wang, Jiebiao, Bernie Devlin, and Kathryn Roeder. "Using multiple measurements of tissue to estimate subject-and cell-type-specific gene expression." *Bioinformatics* 36.3 (2020): 782-788. https://doi.org/10.1093/bioinformatics/btz619
