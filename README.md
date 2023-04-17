# Proteus

*Proteus* is an R package for downstream analysis of *MaxQuant* output. The input for *Proteus* is the evidence file. Evidence data are aggregated into peptides and then into proteins. *Proteus* offers many visualisation and data analysis tools both at peptide and protein level. In particular it allows simple differential expression using *limma*.

## Warning

*Proteus* is no longer under active development, and we believe that some of its features have become outdated. Specifically, we advise against using peptide and protein aggregation from the event file. Instead, we recommend importing the `proteinGroups` file directly into R and utilizing the *maxLFQ* normalization method. We would like to highlight that alternative protein quantification tools, such as *Proteome Discoverer* and *Spectronaut*, are available and offer normalized protein intensity data.

## Installation

Proteus can be installed directly from GitHub. First, you need to install BioConductor and limma:

```r
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("limma")
```

You also need devtools:

```r
install.packages("devtools")
```

In order to run examples or vignette code, additional packages with example data need to be installed:

```r
devtools::install_github("bartongroup/proteusLabelFree")
devtools::install_github("bartongroup/proteusTMT")
devtools::install_github("bartongroup/proteusSILAC")
```

Finally, you can install proteus:

```r
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
```

Note: use `build_vignettes = FALSE` if you run into problems with vignettes installation.

## Tutorial

Proteus contains tutorial vignettes. We suggest starting with the comprehensive tutorial for label-free proteomics:

```r
vignette("proteus", package="proteus")
```

There are additional, shorter vignettes, showing the specifics of using *Proteus* with TMT and SILAC data:

```r
vignette("TMT", package="proteus")
vignette("SILAC", package="proteus")
```

If, for some reason, you cannot build vignettes, copies are available online:

- [Proteus](http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/proteus.html)
- [SILAC](http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/SILAC.html)
- [TMT](http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/TMT.html)

## Application note

The article about this package can be found on [BioRxiv](https://www.biorxiv.org/content/early/2018/09/20/416511).
