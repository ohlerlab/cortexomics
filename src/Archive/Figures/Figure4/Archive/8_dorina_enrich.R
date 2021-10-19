################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}
if(!exists("long_exons")) {
	base::source("src/Figures/Figure4/7_motseqfiles.R")
}

{
library(liftOver)
chainpath = here('ext_data/mm9ToMm10.over.chain')
ch = import.chain(chainpath)
seqlevelsStyle(cur) = "UCSC"  # necessary
cur19 = liftOver(cur, ch)

relevantprots = 'src/dorina_relevant_regulators.txt'%>%readLines

relevantprotreg = relevantprots%>%unique%>%paste0(collapse='|')%>%regex

rnaexpfiles = c('ext_data/CLIP-Seq_mm9_all','ext_data/HITS-CLIP_mm9_all','ext_data/PAR-CLIP_mm9_all')
rnaexpfile = rnaexpfiles[[1]]

allrnaexpdata = lapply(rnaexpfiles,function(rnaexpfile){
	rnaexpdata =  import(rnaexpfile,format='bed')
	rnaexpdata%<>%subset(name%>%str_detect(relevantprotreg))
	rnaexpdata
})
allrnaexpdata = allrnaexpdata%>%do.call(what=c,.)


allrnaexpdata = liftOver(allrnaexpdata,ch)%>%unlist
allrnaexpdatalist = allrnaexpdata%>%split(.,.$name)

reggrlist = list(
	tr = long_exons,
	cds = long_cds,
	fputr = long_fputr,
	tputr = long_tputr
)
}
allrnaexpdatalist%>%names%>%str_subset('Zhang')


dorina_enrich_df = map_df(.id='region',names(regseqlist)%>%setNames(.,.),function(regnm){
	map_df(.id='set',names(setlist)%>%setNames(.,.),function(setnm){
		#
		set = setlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
		grs = reggrlist[[regnm]][set]
		controlset = ctlsetlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
		controlgrs = reggrlist[[regnm]][controlset]
		#
		if(length(controlset)!=0){
			rnaexpdataxptests = lapply(allrnaexpdatalist,function(rnaexpdata){
				setovnums = grs%>%overlapsAny(rnaexpdata)%>%{c(sum(!.),sum(.))}
				ctlovnums = controlgrs%>%overlapsAny(rnaexpdata)%>%{c(sum(!.),sum(.))}
				tidy(fisher.test(matrix(c(ctlovnums,setovnums),ncol=2)))
			})
		}else{
			rnaexpdataxptests = lapply(allrnaexpdatalist,function(rnaexpdata){
				setov = grs%>%overlapsAny(rnaexpdata)
				if(!any(setov)) return(NULL)
				split(order(names(grs)),setov) %>% {wilcox.test(.[[1]],.[[2]])}%>%tidy
			})
		}
		rnaexpdataxptests%>%bind_rows(.id='dataset')
	})
})

dorina_enrich_df%>%
		filter(region=='tputr')%>%
		filter(set %in% 'teupgeneslmatch')%>%
		filter(p.value<0.05)%>%
		select(region,set,estimate,conf.low,conf.high,dataset)%>%
		arrange(estimate)%>%as.data.frame

