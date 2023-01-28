# RiboBayes: Detection and analysis of ribosome pause sites from Ribo-seq data

![](https://github.com/amberluo1/RiboBayes/blob/logo_testing/RiboBayes%2520Logo.png)

## Capabilities

The primary functionality of RiboBayes is to detect the location and expression level of all ribosome pause sites in Ribo-seq (ribosome profiling) datasets. In addition to this function, RiboBayes provides utility for further data exploration and visualization as detailed below. RiboBayes provides functions to:

-   Assign each ribosome pause site a significance score based on its expression change across user-specified conditions

    -   Classify ribosome pause sites as constant, downregulated, and upregulated

-   Visualize the distribution of constant, upregulated, and downregulated pause sites

-   Visualize the localization of constant, upregulated, and downregulated pause sites along mRNA transcripts (e.g. reveal that upregulated pause sites are more common closer to the start codon)

-   Visualize the conservation of pause sites across experiments (e.g. a certain pause site is called as a pause site in 5/6 replicates)

Each of these functions is explored in the vignette below.

## Installation

RiboBayes R package can be installed using devtools as shown below:

```{r}
devtools::install_github("amberluo/RiboBayes")
```

Dependencies include RiboR and edgeR, which must be installed from Bioconductor.
