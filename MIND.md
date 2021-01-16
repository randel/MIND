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
