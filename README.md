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
devtools::install_github("amberluo1/RiboBayes")
```

## Tutorial

RiboBayes uses the RiboR environment to process ribosome profiling data. First, we download a sample .ribo file:

```{bash}
! wget https://github.com/ribosomeprofiling/ribo_manuscript_supplemental/raw/master/sidrauski_et_al/ribo/without_coverage/all.ribo

```

This data is from HEK293 cells (GEO accession number [GSE65778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65778)) and is publicly available.

```{r}
#generates the 'ribo' class object 
ribo <- Ribo(PASTE PATH TO FILE HERE, rename = rename_default )
```

To familiarize yourself with the RiboR environment or to explore the .ribo file, see [this vignette](https://ribosomeprofiling.github.io/ribor/ribor.html).

For this vignette, we can select the top 500 highly expressed mRNA transcripts to detect pause sites on. Most mRNA transcripts with low expression are too noisy to reliably detect pause sites.

```{r}

# Get all viable read lengths for analysis based on UTR expression (higher = less viable)
lengths = get_read_lengths(ribo)
lengths = as.numeric(lengths)

rc_CDS <- get_region_counts(ribo.object    = ribo,
                            range.lower = lengths[[1]],
                            range.upper = lengths[[length(lengths)]],
                            tidy       = TRUE,
                            alias      = TRUE,
                            transcript = FALSE,
                            normalize=TRUE,
                            region     = "CDS",
                            compact    = FALSE)

region_lengths <- get_internal_region_lengths(ribo.object = ribo, alias = TRUE)
cds=region_lengths%>%select(transcript, CDS)
rc_CDS = rc_CDS %>% left_join(cds)%>%mutate(count=count/CDS)

rc_CDS_w = dcast(rc_CDS[,-5], transcript ~ experiment)
high_exp = rowSums( cpm(rc_CDS_w[,-1]) > 353) > 1
sum(high_exp)

high_transcripts=rc_CDS_w$transcript[high_exp]
```

Next, we find all the pause sites in this experiment with the `get_pause_sites()` function:

```{r}
sites = get_pause_sites(ribo, transcripts = high_transcripts)
```

This function may take some time to run. Optionally, you can specify the number of cores to run `get_pause_sites()` on:

```{r}
n_cores = detectCores() - 1
sites = get_pause_sites(ribo, transcripts = high_transcripts, cores = n_cores)
```
