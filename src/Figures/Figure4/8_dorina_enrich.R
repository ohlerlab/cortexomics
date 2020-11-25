################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
if(!exists("long_exons")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure4/7_motseqfiles.R")
}
#mm9ToMm10.over.chain.gz
#wget http://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/mm9ToMm10.over.chain.gz
#mv  mm9ToMm10.over.chain.gz ext_data/mm9ToMm10.over.chain.gz
#gunzip ext_data/mm9ToMm10.over.chain.gz


# wget https://dorina.mdc-berlin.de/api/v1.0/download/regulator/mm9/CLIP-Seq_mm9_all
# wget https://dorina.mdc-berlin.de/api/v1.0/download/regulator/mm9/HITS-CLIP_mm9_all
# wget https://dorina.mdc-berlin.de/api/v1.0/download/regulator/mm9/PAR-CLIP_mm9_all
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

stop()

regnm='tputr'
setnm='tedowngeneslmatch'

dorina_enrich_df = map_df(.id='region',names(regseqlist)%>%setNames(.,.),function(regnm){
	map_df(.id='set',names(setlist)%>%setNames(.,.),function(setnm){
		# regnm <<- regnm
		# setnm <<- setnm

		set = setlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
		grs = reggrlist[[regnm]][set]
		controlset = ctlsetlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
		controlgrs = reggrlist[[regnm]][controlset]

		# rnaexpdata = allrnaexpdatalist[[1]]
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
		# satb2rnaexptests = lapply(allrnaexpdatalist,function(rnaexpdata){
			# setovnums = grs['Satb2']%>%overlapsAny(rnaexpdata)
		# })
		rnaexpdataxptests%>%bind_rows(.id='dataset')
	})
})

dorina_enrich_df%>%
		filter(region=='tputr')%>%
		filter(set %in% 'teupgeneslmatch')%>%
		filter(p.value<0.05)%>%
		select(region,set,estimate,conf.low,conf.high,dataset)%>%
		arrange(estimate)%>%as.data.frame


# t.test(grs%>%width%>%sum%>%log2%>%.[between(.,6,14)],controlgrs%>%width%>%sum%>%log2%>%.[between(.,6,14)])
# >%sum%>%log2%>%.[between(.,6,14)])



grs%>%width%>%sum%>%log2%>%.[between(.,6,14)]%>%txtdensity
controlgrs%>%width%>%sum%>%log2%>%.[between(.,6,14)]%>%txtdensity