library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

args = R.utils::commandArgs(trailingOnly=TRUE,asValues=TRUE,defaults =list(
  bedfile='test.bed',
  phastconsfile='ext_data/mm10.60way.phastCons.bw',
  output='phastcons.out.tsv'
))
if(interactive()) args=args[-1]
attach(args)
# stopifnot(str_detect(bedfile,'fasta$'))

# gr = GRanges('chr10:1000000-1000650')
# names(gr)<-'foo'
# export(gr,bedfile)

gr = rtracklayer::import(bedfile)

if(all(is.na(gr$name))) gr$name <- seq_along(gr)
# gr2 = GRanges(c('chr10:1000000-1000600','chr10:1000700-1000800'))

# gr2%>%coverage%>%rtracklayer::export('test.bw')

pcdata = rtracklayer::import(phastconsfile,which=gr,as='Rle')

pcofbeds = pcdata[gr]

pcofbeds = mean(pcofbeds)

names(pcofbeds) <- gr$name

enframe(pcofbeds,'bedname','mean_phastcons')%>%write_tsv(output)