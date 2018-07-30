Cortexomics
================
Dermot Harnett
July 30, 2018

<!-- ## R Markdown -->
<!-- This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>. -->
Figure One - generate in doc
----------------------------

``` r
plot(1)
```

<img src="cortexomics_files/figure-markdown_github/fig1-1.svg" alt="caption" width="85%" />
<p class="caption">
caption
</p>

Figure Two - include\_graphics svg
----------------------------------

``` r
include_graphics("plots/tmp.svg")
```

<img src="plots/tmp.svg" alt="caption" width="85%" />
<p class="caption">
caption
</p>

include\_graphics pdf
---------------------

``` r
include_graphics(file.path(root,"tmp.pdf")%T>%{normalizePath(.,must=T)%>%cat})
```

    ## /home/zslastman/bih_cluster/fast/groups/ag_ohler/dharnet_m/cortexomics/plots/tmp.pdf

<img src="plots//tmp.pdf" alt="caption" width="85%" />
<p class="caption">
caption
</p>

Figure Three - svg display direct in markdown
---------------------------------------------

![imagename](/home/zslastman/tmp.svg)

Figure Three - pdf display direct in markdown
---------------------------------------------

![imagename](/home/zslastman/tmp.pdf) <!-- ##Figure Three --> <!-- ```{r fig3, out.width = "85%", fig.cap = "caption"} --> <!-- include_graphics(file.path(root,"tmp.pdf")) --> <!-- ``` --> <!-- ##Figure Four --> <!-- ```{r fig4, out.width = "85%", fig.cap = "caption"} --> <!-- include_graphics(file.path(root,"tmp.pdf"), auto_pdf = T) --> <!-- ``` -->

<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- ## Including Plots -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo=FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
