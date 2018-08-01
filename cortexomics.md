Cortexomics
================
Dermot Harnett

``` r
Sys.time()
```

    ## [1] "2018-08-01 12:09:18 CEST"

## Intro

  - Importance of understanding gene regulation

  - Importance of understanding gene regulation at multiple levels

  - Specific importance to the brain
    
      - complexity of translational regulation in the brain
      - importance of complexes

  - Distinction between intra and inter gene variation, discussion of
    it’s importance

  - Lack of studies in living tissue that examine Gene expression at
    multiple levels

  - 
## Refs to include

## results

### Figure 1 - Data Quality

Rnaseq, Riboseq, Mass spec correlations, outlier robust correlations. -
Number, size, mapped reads, correlation between replicates - Enzymes
used, for RIboseq, rationale - Likely effects this has on the data -
e.g. cyclohexamine (see sup figure for sequence
preference)

``` r
include_graphics("plots/tmp.pdf")
```

<embed src="plots/tmp.pdf" title="caption" alt="caption" width="900px" height="900px" type="application/pdf" />

### Figure 2 Riboseq

``` r
include_graphics("plots/tmp.pdf")
```

<embed src="plots/tmp.pdf" title="caption" alt="caption" width="900px" height="900px" type="application/pdf" />

Basic discussion of Riboseq - Periodicity in our Riboseq - Number of
ORFs recovered - New ORFS recoovered - e.g. - Translational Regulation
Tests -
Ribodiff

### Figure 3 Mass Spec Integration

``` r
include_graphics("plots/tmp.pdf")
```

<embed src="plots/tmp.pdf" title="caption" alt="caption" width="900px" height="900px" type="application/pdf" />

  - Observation that many genes don’t have perceptible reaction to
    changing riboseq levels

  - calculation of degredation
rates?

  - 
### Figure 4 - clustering/sequence analysis

``` r
include_graphics("plots/tmp.pdf")
```

<embed src="plots/tmp.pdf" title="caption" alt="caption" width="900px" height="900px" type="application/pdf" />

  - possible clustering methods
      - timeless
      - self organizing maps
  - Sequence analysis

## Discussion

Resource for study of cortical development For comparison to single cell
data

## Supplementary Figures

## Supp Figure 1

## Supp Figure 2

  - Correlation between Riboseq Replicates

  - Sequence biases in the Riboseq replicates

  - Attemps to compensate? ‘We leveraged replicate measurements to
    assess data quality and its dependence on several parameters,
    including alignment strategy, mRNA enrichment method, PCR artifacts,
    gene length normalization, and batch effects (Supplemental Fig.
    S2A–E; Supplemental Methods).’’

  - 
## Supp Figure 3

## References

Javonavic et al 2015
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4506746/>

Cenik et al 2015 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617958/>
- Lots of influence on translation by variants in start codons, ORFs. -
Inter individual variation

Yuzwa et al 2017 cell rep Developmental Emergence of Adult Neural Stem
Cells as Revealed by Single-Cell Transcriptional Profiling.
<https://www.ncbi.nlm.nih.gov/pubmed/29281841>

Cho et al 2015 memory and translational reg Methods should be useful,
note they measure enrichment of reads in coding sequences.

Zahr et al 2018 <https://www.ncbi.nlm.nih.gov/pubmed/29395907>

## Methods
