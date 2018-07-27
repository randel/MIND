MIND (Multi-measure INdividual Deconvolution)
=====

### Using Multiple Measurements of Tissue to Estimate Individual- And Cell-Type-Specific Gene Expression via Deconvolution

A method to glean more insights from bulk gene expression. It borrows information across multiple measurements of the same tissue per individual, such as multiple regions of the brain, using an empirical Bayes approach to estimate individual- and cell-type-specific gene expression via deconvolution.

### Installation

```r
devtools::install_github('randel/MIND') # requiring the `devtools` package
```

### Example

```r
library(MIND)

data(example)

deconv = mind(X = example$X, W = example$W)
# deconv$alpha # deconvolved individual- and cell-type-specific gene expression

# downstream analysis: cell-type-specific network
get_network(alpha = deconv$alpha, W = example$W, cell_type = 3, cor_cutoff = 0.7)
```

For details, please see the [PDF manual](https://github.com/randel/MIND/blob/master/MIND-manual.pdf).

The cell type fraction can be pre-estimated using standard deconvolution method (e.g., [CIBERSORT](https://cibersort.stanford.edu)) and reference samples of purified cells (e.g., [NeuroExpresso](https://pavlab.msl.ubc.ca/data-and-supplementary-information/supplement-to-mancarci-et-al-neuroexpresso/) for brain tissue) or single-cell RNA-seq data.


### Reference
Jiebiao Wang, Bernie Devlin, Kathryn Roeder. Using multiple measurements of tissue to estimate individual- and cell-type-specific gene expression via deconvolution. bioRxiv 379099; doi: https://doi.org/10.1101/379099
