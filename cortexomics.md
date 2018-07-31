Cortexomics
================
Dermot Harnett
July 30, 2018

<!-- ## R Markdown -->

<!-- This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>. -->

\#\#Figure One - generate in
doc

``` r
plot(1)
```

<img src="cortexomics_files/figure-gfm/fig1-1.svg" title="caption" alt="caption" width="85%" />
\#\#Figure Two - include\_graphics svg

``` r
include_graphics("plots/tmp.svg")
```

<img src="plots/tmp.svg" title="caption" alt="caption" width="85%" />
\#\# include\_graphics
pdf

``` r
include_graphics("plots/tmp.pdf")
```

<embed src="plots/tmp.pdf" title="caption" alt="caption" width="85%" type="application/pdf" />

\#\#Figure Three - svg display direct in markdown
![imagename](plots/tmp.svg)

\#\#Figure Three - pdf display direct in markdown
![imagename](plots/tmp.pdf) 

<iframe src="https://github.com/zslastman/cortexomics/blob/master/plots/tmp.pdf" style="width:718px; height:700px;" frameborder="0"></iframe>



<!-- ``` -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- ## Including Plots -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo=FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
