# RiboBayes: Detection and analysis of ribosome pause sites from Ribo-seq data

::: {style="text-align: center;"}
![](RiboBayes%20Logo.png){style="text-align: center" width="387"}
:::

Credit for figuring out how to center logo: Jack Liu (MIT Computer Science and Math '26, IChO '22, USNCO camper '21, USNCO High Honors '20, USAMO '22 '21, Carmel Scibowl Captain '22, Boyfriend of Amber Luo)

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
! wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158374/suppl/GSE158374_manuscript.HEK293.ribo.hdf5
```

This data is from a study on how expression of the SARS-CoV-2 NSP1 and NSP2 proteins affect translation in HEK293 cells (GEO accession number [GSE158374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158374)).

```{r}
# generates the 'ribo' class object 
ribo <- Ribo("./GSE158374_manuscript.HEK293.ribo.hdf5", rename = rename_default)
ribo
```

To familiarize yourself with the RiboR environment, see [this vignette](https://ribosomeprofiling.github.io/ribor/ribor.html).

The Ribo object has 3 experiments (NSP1 expression, NSP2 expression, wild-type), with 3 replicates per experiment. In this vignette, we evaluate the effect of NSP1 expression on ribosome pausing by comparing it to wild-type.

```{r}
# selecting the relevant experiments for analysis 
all_experiments = get_experiments(ribo)
experiments = all_experiments[c(1:3, 7:9)]
```

For this vignette, we can select the top 500 most highly expressed mRNA transcripts to detect pause sites on. Most mRNA transcripts with low expression are too noisy to reliably detect pause sites.

```{r}
# Get all viable read lengths for analysis based on UTR expression (higher = less viable)
lengths <- get_read_lengths(ribo)
lengths <- as.numeric(lengths)

# Find the transcripts with the highest expression, normalized for coding sequence length 
rc_CDS <- get_region_counts(
  ribo.object = ribo,
  range.lower = lengths[[1]],
  range.upper = lengths[[length(lengths)]],
  tidy = TRUE,
  alias = TRUE,
  transcript = FALSE,
  normalize = TRUE,
  region = "CDS",
  compact = FALSE
)
region_lengths <- get_internal_region_lengths(ribo.object = ribo, alias = TRUE)
cds <- region_lengths %>% dplyr::select(transcript, CDS)
rc_CDS <- rc_CDS %>%
  left_join(cds) %>%
  mutate(count = count / CDS)

rc_CDS_w <- dcast(rc_CDS[, -5], transcript ~ experiment)
rc_CDS_w <- rc_CDS_w %>% mutate(summed = rowSums(cpm(rc_CDS_w[, -1])))
rc_CDS_w <- rc_CDS_w %>% arrange(-summed)
high_transcripts <- rc_CDS_w[1:500, ]$transcript
```

### Pause Site Detection on a Single Transcript

We can first run RiboBayes on a single transcript to visualize its functions for pause site detection.

### Pause Site Detection on 500 Highly Expressed Transcripts

Next, we find all the pause sites in this experiment with the `get_pause_sites()` function:

```{r}
sites = get_pause_sites(ribo, transcripts = high_transcripts)
```

This function may take some time to run. Optionally, you can specify the number of cores to run `get_pause_sites()` on:

```{r}
n_cores = detectCores() - 1
sites = get_pause_sites(ribo, transcripts = high_transcripts, cores = n_cores)
```
