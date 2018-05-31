# Proteus

*Proteus* is an R package for downstream analysis of *MaxQuant* output. The input for *Proteus* is the evidence file. Evidence data are aggregated into peptides and then into proteins. *Proteus* offers many visualisation and data analysis tools both at peptide and protein level. In particular it allows simple differential expression using *limma*.

## Installation

Proteus can be installed directly from GitHub. First, you need to install BioConductor and limma:

```r
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
```

You also need devtools:

```r
install.packages("devtools")
```

In order to run examples or vignette code, additional packages with example data need to be installed:

```r
devtools::install_github("bartongroup/proteusUnlabelled", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
devtools::install_github("bartongroup/proteusTMT", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
devtools::install_github("bartongroup/proteusSILAC", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e")
```

Finally, you can install proteus:

```r
devtools::install_github("bartongroup/Proteus", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e", build_vignettes=TRUE)
```

## Tutorial

Proteus contains tutorial vignettes. We suggest starting with the comprehensive tutorial for unlabelled proteomics:

```r
vignette("unlabelled", package="proteus")
```

There are additional, shorter vignettes, showing the specifics of using *Proteus* with TMT and SILAC data:

```r
vignette("TMT", package="proteus")
vignette("SILAC", package="proteus")
```
