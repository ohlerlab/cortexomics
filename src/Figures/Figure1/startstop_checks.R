################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
## Use user specified psite mapping offset 
library(RiboDiPA)
offsets = read_tsv('ext_data/offsets_manual.tsv')%>%{list(qwidth = .$length,psite = .$offset)}
# offsets <- cbind(qwidth = c(28, 29, 30, 31, 32), 
#     psite = c(18, 18, 18, 19, 19))
bam_file_list = Sys.glob('pipeline/star/data/*/*.bam')%>%str_subset('ribo_\\d+.bam')%>%str_subset('E13|E175')
gtf_file = gtf
data.psite2 <- psiteMapping(bam_file_list, gtf_file, psite.mapping = offsets)


## Merge the P-site data into bins with a fixed or an adaptive width
data.binned <- dataBinning(data = data.psite2$coverage, bin.width = 1, 
    zero.omit = FALSE, bin.from.5UTR = TRUE, cores = 4)

## Perform differential pattern analysis
## Perform differential pattern analysis
result.pst <- diffPatternTest(data = data.binned, 
    classlabel = classlabel, method=c('gtxr', 'qvalue'))


ebp1_startstop = read_tsv("../Ebp1_ribo/tables/start_change_tbl.tsv")

