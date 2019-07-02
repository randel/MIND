MIND (Multi-measure INdividual Deconvolution)
================

## Using Multiple Measurements of Tissue to Estimate Subject- And Cell-Type-Specific Gene Expression

`MIND` is a method to glean insights from bulk gene expression. It
borrows information across multiple measurements of the same tissue per
subject, such as multiple regions of the brain, using an empirical Bayes
approach to estimate subject- and cell-type-specific gene expression via
deconvolution.

## Installation

Installation requires the `devtools` package.

``` r
devtools::install_github('randel/MIND')
```

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

Jiebiao Wang, Bernie Devlin, Kathryn Roeder. Using multiple measurements
of tissue to estimate subject- and cell-type-specific gene expression.
bioRxiv 379099; doi: <https://doi.org/10.1101/379099>
