# Proteus

Proteus is an R package to analyse proteomics data from MaxQuant. It starts with the evidence file, can create peptide and protein tables and perform differential expression.

At this stage only unlabelled data are accepted.

## Installation

Proteus can be installed directly from GitHub

```r
devtools::install_github("bartongroup/Proteus", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
```

In order to run examples or vignette code, additional packages with example data need to be installed

```r
devtools::install_github("bartongroup/proteusUnlabelled", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
devtools::install_github("bartongroup/proteusTMT", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
```


## Tutorial

Proteus contains tutorial vignettes for unlabelled and TMT analysis

```r
vignette("unlabelled", package="proteus")
vignette("TMT", package="proteus")
```
