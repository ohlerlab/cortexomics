################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
	load('data/1_integrate_countdata.R')
}

besttrs = iso_tx_countdata$abundance%>%
	as.data.frame%>%
	rownames_to_column('transcript_id')%>%
	left_join(tx2genemap%>%set_colnames(c('transcript_id','gene_id')))%>%
	pivot_longer(-one_of('gene_id','transcript_id'))%>%
	separate(name,c('time','assay','rep'))%>%
	group_by(gene_id,transcript_id)%>%summarise(value=median(value))%>%
	group_by(gene_id)%>%
	filter(gene_id %in% highcountgenes)%>%
	slice(which.max(value))%>%.$transcript_id


## Use user specified psite mapping offset 
library(RiboDiPA)
offsets = read_tsv('ext_data/offsets_manual.tsv')%>%{list(qwidth = .$length,psite = .$offset)}
# offsets <- cbind(qwidth = c(28, 29, 30, 31, 32), 
#     psite = c(18, 18, 18, 19, 19))
bam_file_list = Sys.glob('pipeline/star/data/*/*.bam')%>%str_subset('ribo_\\d+.bam')%>%str_subset('E13|E175')
classlabel = data.frame(condition = bam_file_list%>%str_extract('E13|E175'))%>%
		mutate(comparison = ifelse(condition %in% c('E175'),2,1))

gtf_file = gtf
data.psite2 <- psiteMapping(bam_file_list, gtf_file, psite.mapping = offsets)

## Merge the P-site data into bins with a fixed or an adaptive width
data.binned <- dataBinning(data = data.psite2$coverage, bin.width = 1, 
    zero.omit = FALSE, bin.from.5UTR = TRUE, cores = 2)

stopifnot(all(highcountgenes%in%names(data.binned)))
{
mod_databin = data.binned




startseclen = 100

mod_databin = mod_databin[highcountgenes]

## Perform differential pattern analysis
## Perform differential pattern analysis

stbin = data.binned%>%lapply(.%>%.[,1:min(ncol(.),startseclen)])

nullmat = matrix(0,nrow=length(bam_file_list))

truncbin = data.binned%>%lapply(possibly(o=nullmat,.%>%.[,(startseclen+1):ncol(.)]%>%rowSums%>%matrix%>%set_colnames(startseclen+1)))

mod_databin = map2(stbin,truncbin,.f=cbind)

# for(gn in highcountgenes){
# 	mod_databin[[gn]] %<>% .[,1:min(ncol(.),50)]
# }
highc = mod_databin%>%map_dbl(sum)%>%which.max
mod_databin[[highc]] = cbind(mod_databin[[highc]],nullmat)
mod_databin[[highc]]%<>%{colnames(.)[ncol(.)]=startseclen+2;.}

}


mod_databin[[highc]]%>%dim
mod_databin[[highc-1]]%>%dim

result.pst <- diffPatternTest(data = data.binned[(highc-300):(highc+300)], classlabel = classlabel, method=c('gtxr', 'qvalue'))
result.pst <- diffPatternTest(data = mod_databin[(highc-300):(highc+300)], classlabel = classlabel, method=c('gtxr', 'qvalue'))

ebp1_startstop = read_tsv("../Ebp1_ribo/tables/start_change_tbl.tsv")

#now plot
plotfile<- here(paste0('plots/','ribodipa','.pdf'))
pdf(plotfile)
plotTest(result = result.pst, genes.list = NULL, threshold = 0.05)%>%lapply(function(p)
	p+xlim(c(1,100))
)
message(normalizePath(plotfile))
dev.off()
 

