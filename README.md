MIND (Multi-measure INdividual Deconvolution)
=====

### Using Multiple Measurements of Tissue to Estimate Subject- And Cell-Type-Specific Gene Expression via Deconvolution

A method to glean more insights from bulk gene expression. It borrows information across multiple measurements of the same tissue per subject, such as multiple regions of the brain, using an empirical Bayes approach to estimate subject- and cell-type-specific gene expression via deconvolution.

### Installation

```r
devtools::install_github('randel/MIND') # requiring the `devtools` package
```

### Example

```r
library(MIND)

data(example)

deconv = mind(X = example$X, W = example$W)
# deconv$A # deconvolved subject- and cell-type-specific gene expression
```

For details, please see the [PDF manual](https://github.com/randel/MIND/blob/master/MIND-manual.pdf).

The cell type fraction can be pre-estimated using standard deconvolution method (e.g., [CIBERSORT](https://cibersort.stanford.edu)) and reference samples of purified cells (e.g., [NeuroExpresso](https://pavlab.msl.ubc.ca/data-and-supplementary-information/supplement-to-mancarci-et-al-neuroexpresso/) for brain tissue) or single-cell RNA-seq data.


### Reference
Jiebiao Wang, Bernie Devlin, Kathryn Roeder. Using multiple measurements of tissue to estimate subject- and cell-type-specific gene expression. bioRxiv 379099; doi: https://doi.org/10.1101/379099
