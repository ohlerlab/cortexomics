


################################################################################
########Let's also check relative occupancy within protein secondary structure
################################################################################
#FT domain?
#  [1] "Entry"                      "Entry name"
#  [3] "Status"                     "Protein names"
#  [5] "Gene names"                 "Organism"
#  [7] "Coiled coil"                "Compositional bias"
#  [9] "Domain [CC]"                "Domain [FT]"
# [11] "Repeat"                     "Beta strand"
# [13] "Helix"                      "Turn"
# [15] "Cross-reference (InterPro)" "Signal peptide"
# [17] "Transit peptide"            "Transmembrane"
# [19] "Length"
#read in the swissprot info on secondary strcuture
{
sec_struct_df = fread('ext_data/uniprot-yourlist_M202101185C475328CEF75220C360D524E9D456CE0D1394W.tab.gz')

structcols = c(
 # "Domain [CC]"                ,
 "Domain [FT]",
 "Repeat"                     ,"Beta strand",
 "Helix"                      ,"Turn",
 # "Cross-reference (InterPro)" ,
 "Signal peptide",
 "Transit peptide",
 "Transmembrane"
)
#
sec_domain_df <- map_df(structcols,(function(structcol){
	sec_struct_df%>%
		# select(Entry,rcol = !!sym(structcol))%>%
		select(swissprot_id=Entry,rcol = !!structcol)%>%
		mutate(rcol = str_extract_all(rcol,'\\d+\\.\\.\\d+'))%>%
		unnest(rcol)%>%
		separate(rcol,c('start','end'),sep='\\.\\.')%>%
		mutate(domain=structcol)
}))
#also the swissprot tr_id to swissprot mapping
tr_id2sw_id = 'ext_data/gencode.vM12.metadata.SwissProt.gz'%>%fread(header=F)%>%select(tr_id=V1,swissprot_id=V3)
tr_id2sw_id$tr_id%<>%str_replace('\\.\\d+$','')
#get sw_id to tr_id match, only in cases where the length is correct
sw2tr_id_lmatch <- sec_struct_df%>%select(aalen=Length,swissprot_id=Entry)%>%
	inner_join(tr_id2sw_id)%>%
	inner_join(cdsgrl%>%width%>%sum%>%enframe('tr_id','length'))%>%
	mutate(lenmatch=(length-3)==(aalen*3))%>%
	filter(lenmatch)%>%
	select(swissprot_id,tr_id)
#
sec_domain_df <- sec_domain_df%>%inner_join(sw2tr_id_lmatch)
#now convert from AA space, to codon, to transcript
sec_domain_df$start %<>% as.numeric %>% multiply_by(3) %>% subtract(-2)
sec_domain_df$end %<>% as.numeric %>% multiply_by(3)
sec_domain_df%<>%filter(tr_id %in% ribocovtrs)
sec_domain_df$start %<>% add(cdsstarts[sec_domain_df$tr_id]-1)
sec_domain_df$end %<>% add(cdsstarts[sec_domain_df$tr_id]-1)
}

{

#convert to a granges object
sec_domain_gr = sec_domain_df%>%
	rename('seqnames' := tr_id)%>%
	data.frame%>%GRanges
#harmonize the seqnames
domribocovtrs = ribocovtrs%>%intersect(sec_domain_gr@seqnames)
sec_domain_gr=keepSeqlevels(sec_domain_gr,domribocovtrs,pruning='coarse')
# sec_domain_gr%<>%resize(width(.)+30,'center')
seqinfo(sec_domain_gr) = trseqinfo[domribocovtrs]
#mark out the start and end as 2ndary domains
names(trcds) = seqnames(trcds)
sec_domain_gr%<>%c(trcds%>%keepSeqlevels(domribocovtrs,'coarse')%>%resize(30,'start')%>%{.$domain='start';.})
sec_domain_gr%<>%c(trcds%>%keepSeqlevels(domribocovtrs,'coarse')%>%resize(30,'end')%>%{.$domain='end';.})
sec_domain_gr = c(sec_domain_gr,sec_domain_gr%>%gaps%>%{.$domain='no_domain';.})
domaintypes=sec_domain_gr%>%.$domain%>%unique
#now compute the mean relative density inside each of the domains
mainsamps%<>%setNames(.,.)
sample=mainsamps[1]
domname='start'

trseqinfo['ENSMUST00000212864']
seqinfo(sec_domain_gr)['ENSMUST00000212864']
'ENSMUST00000212864'%in%domribocovtrs
seqinfo(sec_domain_gr)@seqnames

}

if(!file.exists(here('data/struct_dens_tbls.rds'))){

	domnms=sec_domain_gr$domain%>%unique%>%setNames(.,.)
	# domnms=c('Transmembrane','no_domain','start','end')%>%setNames(.,.)
	struct_dens_tbls <-  lapply(domnms,function(domname){
		mclapply(mc.cores=4,mainsamps,function(sample){
			domranges = sec_domain_gr%>%subset(domain==domname)
			domaindens = psitecovnorm[[sample]][domranges]%>%mean
			trmean = psitecovnorm[[sample]]%>%mean
			(domaindens)/(as(trmean,"List")[domranges@seqnames])
		})
	})
	saveRDS(struct_dens_tbls,here('data/struct_dens_tbls.rds'))

}else{
	struct_dens_tbls<-readRDS(here('data/struct_dens_tbls.rds'))
}




#format table for plotting, taking mean 
#over all domain ranges
struct_dens_tbl = struct_dens_tbls%>%
	map_depth(2,.%>%mean%>%mean(na.rm=T))%>%
	map(enframe,'sample','mean_rel_dens')%>%
	bind_rows(.id='domain')%>%
	unnest(mean_rel_dens)
#
#now plot
plotfile<- here(paste0('plots/','sec_struct_dens','.pdf'))
pdf(plotfile)
struct_dens_tbl%>%
	filter(!domain%in%c('start','end'))%>%
	mutate(stage=sample%>%str_replace('_\\w+_\\d+$',''))%>%
	group_by(domain)%>%mutate(mean_rel_dens = mean_rel_dens/mean(mean_rel_dens[stage=='E13']))%>%
	ggplot(.,aes(y=mean_rel_dens,color=stage,fill=stage,x=stage))+
	# stat_identity(geom='bar',position='dodge')+
	geom_point(geom='bar',position='dodge')+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_color_manual(name='stage',values=stagecols)+
	scale_y_continuous(paste0('Mean Relative RPF Density'))+
	facet_wrap(domain ~ .,scale='free_y')+
	ggtitle(paste0('Density vs Secondary Structure'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

################################################################################
########metageneplots for signal pep genes
################################################################################
	
hassigpep = sec_domain_df%>%{setNames(.$domain%>%str_detect('Signal'),.$tr_id)}

sigpeptrs = longcdstrs%>%intersect(names(hassigpep)[hassigpep])
nosigpeptrs = longcdstrs%>%intersect(names(hassigpep)[!hassigpep])


metasignaldf_stgrp <-  bind_rows(.id='fraction',list(
		sig_pep=get_metasignaldf(bindatamats[mainsamps],sigpeptrs),
		no_sig_pep=get_metasignaldf(bindatamats[mainsamps],nosigpeptrs)
	))%>% 
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
{
library(rlang)
plotfile<-'plots/sigpepmetaplot.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.012))
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

################################################################################
########Also for low/high TE
################################################################################

tevals = iso_tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(dataset,val,-tr_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	group_by(time,assay,tr_id)%>%summarise(val=mean(val))%>%
	spread(assay,val)%>%
	mutate(te = log2(ribo) - log2(total))%>%
	select(tr_id,time,ribo,total,te)
tevals%<>%filter(is.finite(te))
te_topqtrs = tevals%>%filter(tr_id%in%longcdstrs,time=="E13")%>%filter(te>=quantile(te,0.75,na.rm=T))%>%.$tr_id
te_lowqstrs = tevals%>%filter(tr_id%in%longcdstrs,time=="E13")%>%filter(te <=quantile(te,0.25,na.rm=T))%>%.$tr_id
te_midqstrs = setdiff(longcdstrs,c(te_topqtrs,te_lowqstrs))
message(length(te_topqtrs))
message(length(te_midqstrs))
message(length(te_lowqstrs))
metasignaldf_stgrp <-  bind_rows(.id='fraction',list(
		te_topq=get_metasignaldf(bindatamats[mainsamps],te_topqtrs),
		te_midqs =get_metasignaldf(bindatamats[mainsamps],te_midqstrs),
		te_lowqs=get_metasignaldf(bindatamats[mainsamps],te_lowqstrs)
	))%>% 
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
metasignaldf_stgrp$fraction%<>%factor(levels=c('te_topq','te_midqs','te_lowqs'))
{
library(rlang)
plotfile<-'plots/te_quart_metaplot.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot()+facet_grid(scale='free')
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}


################################################################################
########And for high/low monosomal enrichment
################################################################################

fr_alllimmares = 'data/fr_alllimmares.rds'%>%readRDS
monovals = 'data/fr_alllimmares.rds'%>%readRDS
gnm2trid = setNames(longcdstrs,trid2gnm[[longcdstrs]])
monovals%<>%mutate(tr_id = gnm2trid[gene_name])
monovals%<>%filter(!is.na(tr_id))
#
monotopqtrs = monovals%>%filter(tr_id%in%longcdstrs,contrast=="fraction80S")%>%filter(logFC>=quantile(logFC,0.75,na.rm=T))%>%.$tr_id
monolowqstrs = monovals%>%filter(tr_id%in%longcdstrs,contrast=="fraction80S")%>%filter(logFC <=quantile(logFC,0.25,na.rm=T))%>%.$tr_id
monomidqstrs = setdiff(longcdstrs,c(monotopqtrs,monolowqstrs))
#
metasignaldf_stgrp <-  bind_rows(.id='fraction',list(
		mono_topq=get_metasignaldf(bindatamats[mainsamps],monotopqtrs),
		mono_midqs =get_metasignaldf(bindatamats[mainsamps],monomidqstrs),
		mono_lowqs=get_metasignaldf(bindatamats[mainsamps],monolowqstrs)
	))%>% 
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
metasignaldf_stgrp$fraction%<>%factor(levels=c('mono_topq','mono_midqs','mono_lowqs'))
{
library(rlang)
plotfile<-'plots/mono_quart_metaplot.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.0075))+facet_grid()
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}
tevals%<>%mutate(gene_name = trid2gnm[[tr_id]])
fr_alllimmares%>%head
tevals%>%head
#
fr_alllimmares%>%filter(contrast=='fraction80S')%>%left_join(tevals%>%filter(time=='E13'))%>%
	filter(te>-12)%>%
	filter(is.finite(te),is.finite(logFC))%>%
	{quicktest(.$logFC,.$te)}


################################################################################
########And high/low transcriptional change
################################################################################

dtxnvals = iso_tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(dataset,val,-tr_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	group_by(time,assay,tr_id)%>%summarise(val=mean(val))
dtxnvals = dtxnvals%>%filter(assay=='total')%>%group_by(tr_id)%>%summarise(dtxn = val[time=='P0']/val[time=='E13'])
dtxnvals%<>%filter(is.finite(dtxn))
dtxn_topqtrs = dtxnvals%>%filter(tr_id%in%longcdstrs)%>%filter(dtxn>=quantile(dtxn,0.75,na.rm=T))%>%.$tr_id
dtxn_lowqstrs = dtxnvals%>%filter(tr_id%in%longcdstrs)%>%filter(dtxn <=quantile(dtxn,0.25,na.rm=T))%>%.$tr_id
dtxn_midqstrs = setdiff(longcdstrs,c(dtxn_topqtrs,dtxn_lowqstrs))
message(length(dtxn_topqtrs))
message(length(dtxn_midqstrs))
message(length(dtxn_lowqstrs))
metasignaldf_stgrp <-  bind_rows(.id='fraction',list(
		dtxn_topq=get_metasignaldf(bindatamats[mainsamps],dtxn_topqtrs),
		dtxn_midqs =get_metasignaldf(bindatamats[mainsamps],dtxn_midqstrs),
		dtxn_lowqs=get_metasignaldf(bindatamats[mainsamps],dtxn_lowqstrs)
	))%>% 
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
metasignaldf_stgrp$fraction%<>%factor(levels=c('dtxn_topq','dtxn_midqs','dtxn_lowqs'))
{
library(rlang)
plotfile<-'plots/dtxn_quart_metaplot.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot()+facet_grid(scale='free')
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

################################################################################
########
################################################################################
	

nonmainsamps = names(bindatamats)%>%setdiff(mainsamps)

metasignaldf_stgrp <- get_metasignaldf(bindatamats[nonmainsamps],longcdstrs) %>% 
	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
	mutate(stage = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
	group_by(stage,section,start,sample,fraction)%>%
	summarise(signal=mean(signal))
{
library(rlang)
plotfile<-'plots/figures/figure1/fig1c_myribowaltz_frac_metaplots.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.006))
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

################################################################################
########Check for ramp effect on the scale of weinberg.
################################################################################
	##Let's also check weinbergs start ramp effect
wbergstartwinds = trcds[ribocovtrs]%>%subset(width>1200)%>%resize(width(.),'end')%>%resize(1200)
seqinfo(wbergstartwinds) = trseqinfo
wbergstartwinds = wbergstartwinds[!is_out_of_bounds(wbergstartwinds)]
#
rampprofiles <- mclapply(mc.cores=4,mainsamps,function(sample){
	psitecovnorm[[sample]][wbergstartwinds]%>%as.matrix%>%colMeans
})
#
crampprofiles <- lapply(mainsamps,function(sample){
	rampprofiles[[sample]][-400]%>%matrix(.,nrow=3)%>%colSums%>%
		enframe('position','val')
})
#
#now plot
plotfile<- here(paste0('plots/','wberg_ramp_plot','.pdf'))
pdf(plotfile)
crampprofiles%>%bind_rows(.id='sample')%>%
		mutate(stage=sample%>%str_replace('_\\w+_\\d+$',''))%>%
	ggplot(.,aes(x=position,y=val,color=stage,group=sample))+
	geom_line()+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_color_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('codon'))+
	scale_y_continuous(paste0('Mean RPF density'))+
	ggtitle(paste0('Ramp Effect on RPF density '))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


	

################################################################################
########Mono enrichment vs Ramp
################################################################################
#statically at one timepoint, does monosomeal poly, enrich for the start effect?	


'data/fr_alllimmares.rds'%>%readRDS






