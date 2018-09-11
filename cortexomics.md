Cortexomics
================
Dermot Harnett

``` r
Sys.time()
```

    ## [1] "2018-09-11 16:29:46 CEST"

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
In this work, we have carried out a comprehensive analysis of gene
regulation in the developing murine neocortex, assaying gene regulation
at the level of steady state RNA levels, Ribosome bound RNA, and steady
state protein.

First, we demonstrate how by accounting for the unique technical biases
of each assay, a more accurate account of gene regulation throughout the
central dogma (intro to dogma quant here) can be achieved.

Secondly, our analysis provides a unique window into the relative
magnitude and prevalence of transcriptional, translational and post
translational gene regulation mechanisms over the course of development
in a well studied tissue, allowing us to compare variance between and
within gene’s expression levels over time.

Thirdly, our analysis builds on existing literature (\_ demonstrating
the importance of translational regulatory mechanisms to in cortical
development.

### Refs to include

[Fujii et al 2007](https://www.ncbi.nlm.nih.gov/pubmed/28195124)
FACS+Riboseq translational regulation of lots of cell signalling
molecules - Shh wnt hippo pi3k mapk uORFs prevalent. confirm importance
of uORF for Ptch1 signaling. lots of mRNAs have translational but no
transcriptional reg. - let’s look at their processing, did they do p
sites? I worry about uorfs and offsets - exclude the last 5 and first 15
codons… this is going to kill power, but I see where they’re coming from
with CHX. - Base a lot of theri analysis on absolute TE levels…. there
are a LOT of confounders in that not sure if they accounted for all -
mapping, splicing, ramping, elongation speed, stalling…

[](https://www.cell.com/cell/abstract/S0092-8674\(17\)30582-2) \#\#
results other attached packages: \[1\] gdtools0.1.7 compiler==3.4.3
Rcpp==0.12.18 svglite==1.2.1 \>

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

### Figure 3 Mass Spec Integration (jpeg file)

``` r
include_graphics("plots/tmp.jpeg")
```

<img src="plots/tmp.jpeg" title="caption" alt="caption" width="900px" height="900px" />
- Observation that many genes don’t have perceptible reaction to
changing riboseq levels - calculation of degredation rates?
-

### Figure 4 - clustering/sequence analysis (jpeg directly made)

``` r
plot(3)
```

<img src="cortexomics_files/figure-gfm/f4-1.jpeg" title="caption" alt="caption" width="900px" height="900px" />

``` r
# include_graphics("plots/tmp.jpeg")
```

  - possible clustering methods
      - timeless
      - self organizing maps
  - Sequence analysis

## Discussion

Resource for study of cortical development For comparison to single cell
data

## Supplementary Figures

## Supp Figure 1

Data quality Duplication rates

![The multiqc report](pipeline/multiqc_report.html) shows various
summary statistics for the libraries.

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

We processed data with a pipeline constructed using [Snakemake](_link),
and a set of scripts written in R and Python. Code is available on
github at [github link](_link). Riboseq reads were first trimmed of
illumina adaptors and (2xNNNN) UMI barcodes, before being filtered for
tRNA and rRNA. tRNA and rRNA was removed by aligning to tRNA and
(spliced) rRNA sequences from gencode, version M12 using bowtie2. Any
Riboseq read aligning to these sequences, allowing for up to 1 mismatch,
was eliminated from subsequent analysis. we then used [STAR](_link) with
the following parameters (insert later) randomly assigning multimapping
reads to one of their best mapping positions and eliminating reads which
mapped to 20 positions or more. [samtools](_link), [bamtools](_link),
[Picard](_link) and [multiqc](_link) were then used to produce an
aggregated QC report for all RNAseq and RIboseq reads. \[we’ll want some
stuff about the final results of our QC here\] We then applied
\[\_\_\_\_either ORFik or Riboqc\_\_\_, assuming we go this route\] to
obtain Psite counts from our Riboseq data, and used these, along with
Rsubread feature counts, to obtain gene-level RNAseq and P-site counts
for our data. [Ribodiff](_link) was then used to assess differential
translational efficiency. Joint modeling of riboseq and rnaseq data was
carried out by transforming all data with a variance stabilizing
transformation (which resulted in approximately heteroskedastic data,
see figure \_ ) and then fitting a linear model with [Limma](_link).

(Notes on the above - right now I’m assumign we’ll sort out he
periodicity issue and then get p-sites, if not we might just end up
counting reads and not using ORFik or riboqc. I’m leaving out SaTAnn and
the location of ORF.)
