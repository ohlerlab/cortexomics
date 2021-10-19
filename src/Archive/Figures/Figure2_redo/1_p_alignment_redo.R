if(TRUE){
# if(!exists('fpcovlist')) base::source(here::here("src/Figures/Figure2/1_load_pos_data.R"))
if(!exists('fafile')) base::source(here::here("src/Figures/load_annotation.R"))
if(is.environment(fafile)){ fafileob = fafile} else{fafileob = Rsamtools::FaFile(fafile)}

offsets <- here('ext_data/offsets_manual.tsv')%>%read_tsv
# samp1=names(fpcovlist)[[1]]
samp2='E175_ribo_2'
samp1=samp2
# rl = names(fpcovlist[[samp1]])[[1]]
rl='26'
rln=as.numeric(rl)
rloffset=offsets%>%filter(length==rln)%>%.$offset
#
#we probably want to train our model on all samples together
ribocovtrs<-readRDS('data/ribocovtrs.rds')
#filter for only long enough ones, if we're doing windows
trspacecds = pmapToTranscripts(cdsgrl[ribocovtrs],exonsgrl[ribocovtrs])
trspacecds%<>%unlist
#
cdsstarts = trspacecds%>%start%>%setNames(names(trspacecds))
cdsprestops = trspacecds%>%end%>%`-`(2+3)%>%setNames(names(trspacecds))
#


addphase <- function(gr,cdsstarts){
	gr$phase = unlist((start(gr) - cdsstarts[as.vector(seqnames(gr))])%%3)
	gr
}
mainrls=as.character(25:31)%>%setNames(.,.)
string2onehot<-function(scol,acgt=c('A','C','T','G')){
	egval = rep(1,length(scol))
	acgt%>%
		vapply(function(lev){as.numeric(scol==lev)},egval)%>%
		matrix(ncol=4,dimnames=list(NULL,acgt))
}
#
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	nbp = ncol(mat)/2
	n=1
	posmats = lapply(1:ncol(mat),function(n){
		string2onehot(mat[,n,drop=FALSE])%>%
		set_colnames(paste0(pre,n-nbp-1,'_',colnames(.)))
	})%>%
	purrr::reduce(cbind)
	posmats
}
prop <- function(x,rnd=3) round(x/sum(x),rnd)
enddist <- function(gr){
	end(gr) - seqlengths(gr)[as.vector(seqnames(gr))]
}

width1grs <- function(gr){
	stopifnot(Negate(is.unsorted)(gr))
	isw1 <- width(gr)==1
	broad <- gr[!isw1]
	#vector of integers - 1,2,3 for range 1-3
	narrowstarts <- unlist(as(broad@ranges,'IntegerList'))
	narrow <- {GRanges(
			rep(seqnames(broad),width(broad)),
			IRanges(narrowstarts,w=1)
		)}
	mcols(narrow) <- mcols(broad)[rep(seq_along(broad),width(broad)),,drop=F]
	sort(c(gr[isw1],narrow))
}

add_flank_bases <- function(reads,exonsgrl,fafileob,nbp=2){
	#
	#this gets the flanking bases for reads,
	stopifnot(all(unique(seqnames(reads))%in%names(exonsgrl)))
	stopifnot(nbp>0)
	stopifnot((nbp%%1)==0)
	stopifnot(width(head(reads))!=1)
	#
	alltrs = names(exonsgrl)
	nbp=2
	#calculate how much the ends need to be extended
	maxtpext = reads%>%enddist%>%
		enframe('tr_id','tpext')%>%
		group_by(tr_id)%>%
		summarise_at(vars(tpext),max)
	maxtpext$tpext = maxtpext$tpext+nbp
	maxtpext$tpext %<>% pmax(0)
	maxtpext %<>% {setNames(.[[2]],.[[1]])}
	def = rep(0,length(alltrs))%>%setNames(alltrs)
	missing = setdiff(names(def),names(maxtpext))
	maxtpext = c(maxtpext)%>%{c(.,def[missing])}
	#calculate how much the starts need to be extended
	maxfpext = reads%>%{setNames(start(.),seqnames(.))}%>%
		enframe('tr_id','tpext')%>%
		group_by(tr_id)%>%
		summarise_at(vars(tpext),min)
	maxfpext$tpext = -(maxfpext$tpext-1-2)
	maxfpext$tpext %<>% pmax(0)
	maxfpext %<>% {setNames(.[[2]],.[[1]])}
	def = rep(0,length(alltrs))%>%setNames(alltrs)
	missing = setdiff(names(def),names(maxfpext))
	maxfpext = c(maxfpext)%>%{c(.,def[missing])}
	#expanded exons grl with room for the flanks of our reads, and any
	#3' ends that align off the edge.
	read_trs <- seqnames(reads)%>%unique%>%as.vector
	exp_exonsgrl <- exonsgrl[read_trs]%>%sort_grl_st%>%
		{resize_grl(.,sum(width(.))+maxfpext[read_trs],fix='end')}%>%
		{resize_grl(.,sum(width(.))+maxtpext[read_trs],fix='start')}
	#and it's sequence
	exp_exonseq <- extractTranscriptSeqs(exp_exonsgrl,x=fafileob)
	#expanded seqinfo object
	expseqinfo = Seqinfo(seqnames=names(exp_exonsgrl),seqlengths=sum(width(exp_exonsgrl)))
	#define 5' ends in expanded exons spaces
	fpflank = reads%>%resize(nbp,'start')%>%
		shift(maxfpext[as.vector(seqnames(.))])%>%
		keepSeqlevels(.,read_trs)%>%
		{seqinfo(.)<-expseqinfo;.}%>%
		resize(nbp+nbp,'end')
	exp_exonseq[fpflank]->mat
	fpseqmat = exp_exonseq[fpflank]%>%dnaseq2onehot('fp_')
	#reads%>%head(1)%>%resize(4)%>%shift(-2)%>%getSeq(.,x=exonseq)
	#define 3' ends in expanded exon space
	tpflank = reads%>%
		keepSeqlevels(read_trs)%>%
		{seqinfo(.)<-expseqinfo;.}%>%
		resize(nbp,'end')%>%
		shift(maxfpext[as.vector(seqnames(.))])%>%
		resize(nbp+nbp,'start')
	#and for 3' end
	tpseqmat = exp_exonseq[tpflank]%>%dnaseq2onehot('tp_')
	#combine the sequence 1 hot matrices
	seqmat <- cbind(fpseqmat,tpseqmat)
	mcols(reads)%<>%.[,colnames(.)%>%setdiff(colnames(seqmat))]
	mcols(reads)%<>%cbind(seqmat)
	reads
}
setstrand<-function(gr,str='+'){strand(gr)<-str;gr}
startwinds = trspacecds%>%setstrand%>%resize(1)%>%resize(3+6)%>%resize(width(.)+6,'end')%>%
	.[!is_out_of_bounds(.)]
stopwinds = trspacecds%>%setstrand%>%resize(1,'end')%>%resize(3+6,'end')%>%resize(width(.)+6)%>%
	.[!is_out_of_bounds(.)]


refends=c('start','stop')
ref_end=refends[[2]]
if(!file.exists(here('data/inbams.rds'))){
	library(GenomicAlignments)
	trbams = Sys.glob('pipeline/star_transcript/data/*/*.sort.bam')%T>%{stopifnot(has_length(.))}
	trbams%T>%{stopifnot(file.exists(.))}
	trbams%<>%setNames(.,basename(dirname(.)))
	inbams <- mclapply(mc.cores=10,trbams,ribocovtrs,F=function(trbam,ribocovtrs){
		inbam = readGAlignments(trbam,
			param=ScanBamParam(what=c('mapq')))
		mcols(inbam)$length<-qwidth(inbam)
		inbam%<>%renameSeqlevels(seqlevels(inbam)%>%str_extract('\\w+'))
		stopifnot(mean(ribocovtrs%in%seqnames(inbam))>.9)
		inbam%>%keepSeqlevels(ribocovtrs,pruning='coarse')
	})
	saveRDS(inbams,here('data/inbams.rds'))
}else{
	inbams<-readRDS(here('data/inbams.rds'))
}
# fpcov = fpcovlist[[samp1]][[rl]]

#read objects
readlens=(25:31)%>%setNames(.,.)

rln=29
cdsstarts_=cdsstarts

if(!file.exists(here('data/samp_rl_reads.rds'))){
		
	samp_rl_reads = mclapply(mc.cores=10,inbams,function(ireads,cdsstarts_=cdsstarts){
		lapply(readlens,function(rln){
			#
			if(ref_end=='start'){tr_target_pos=cdsstarts_}else{tr_target_pos=cdsprestops}
			ireads = ireads%>%subset(length==rln)%>%
				subset(seqnames%in%ribocovtrs)
			ireads = ireads[,NULL]%>%
				as("GRanges")%>%
				{.$targetpos=tr_target_pos[as.vector(seqnames(.))];.}%>%
				addphase(cdsstarts_)%>%
				{.$cor_offset=(.$targetpos+.$phase)-start(.);.}
			# return(length(ireads))
			#
			allpos = expand.grid(phase=0:2,cor.offset=0:rln)%>%
				filter((cor.offset%%3)==0)%>%as.data.frame
			frame_stat_df=ireads%>%
					subset(cor_offset<(rln-6))%>%
					subset(cor_offset>5)%>%
					{names(.)<-NULL;.}%>%
					as.data.frame%>%
					# mutate_at(vars(score),replace_na,0)%>%
					left_join(allpos,.)%>%
					group_by(phase,cor_offset)%>%
					# tally(wt=score)%>%
					tally()%>%
					# {(.$n*.$score)%>%sum}
					# group_by(phase)%>%
					ungroup%>%
					arrange(desc(cor_offset))%>%
					mutate(twind=n+lag(n)+lead(n))
					# group_by(cor_offset)%>%
					# filter(n()==3)%>%
					# mutate(twind=sum(n))
			frame_stat_df%>%as.data.frame
			phasevect = 
				frame_stat_df%>%
				slice(((which.max(twind)-1):(which.max(twind)+1)))%>%
				{setNames(.$cor_offset,.$phase)}
			message(capture.output(phasevect)%>%paste0(collapse='\n'))
			#
			ireads%<>%{.$offset <- phasevect[as.character(.$phase)];.}
			ireads%<>%{.$error=(.$targetpos+.$phase)-(start(.)+.$offset);.}
			ireads$cor_offset=NULL
			ireads$targetpos=NULL
			list(ireads,frame_stat_df)
	})
	})
	saveRDS(samp_rl_reads,here('data/samp_rl_reads.rds'))
}else{
	samp_rl_reads<-readRDS(here('data/samp_rl_reads.rds'))
}
rm(inbams)
gc()


if(!file.exists(here('data/offsetstatdf.rds'))){
	offsetstatdf = samp_rl_reads%>%map_depth(2,2)
	offsetstatdf%<>% map_df(.id='sample',~bind_rows(.id='length',.))
	saveRDS(offsetstatdf,here('data/offsetstatdf.rds'))
}else{
	offsetstatdf<-readRDS(here('data/offsetstatdf.rds'))

}

lphaseoffsetdf <- offsetstatdf%>%
	# filter(sample=='E13_ribo_1',length==29)%>%
	group_by(sample,length)%>%
	dplyr::slice(((which.max(twind)-1):(which.max(twind)+1)))


samp_rl_reads <- samp_rl_reads%>%map_depth(2,1)

#a ten gigabyte object.... sub optimal.
if(!file.exists(here('data/seqshiftmodels.rds'))){
	seqshiftmodels <- samp_rl_reads%>%mclapply(mc.cores=10,function(samplreadlist){
		samplreadlist%>%lapply(function(ireads){
			trs = seqnames(ireads)%>%unique
			traintrs = trs%>%sample(.,floor(length(.)*0.99))
			testtrs = trs%>%setdiff(traintrs)
			# ireads%>%subset(seqnames%in%traintrs)%>%subset(phase==iphase)%>%subset(between(error,-6,6)) -> ireads
			rfdf <- ireads%>%
				subset(seqnames%in%traintrs)%>%
				subset(between(error,-6,6))%>%
				add_flank_bases(exonsgrl,fafileob)%>%
				mcols%>%as.data.frame%>%
				# select(-phase,-mapq,-length,-offset)
				select(-offset)
			rfdf$error%<>%as_factor
			rfdf%<>%select(error,phase,everything())
			nvars = length(colnames(rfdf))-2
			selweights = c(0,rep(1/nvars,nvars))
			message('training probabalistic random forrest')
			seqshiftmodel <- ranger::ranger(
				formula= error ~ . ,
				data=rfdf,
				# importance='permutation',
				probability=TRUE,
				split.select.weights=selweights
			)
			seqshiftmodel
		})
	})
	saveRDS(seqshiftmodels,here('data/seqshiftmodels.rds'))
}else{
	seqshiftmodels<-readRDS(here('data/seqshiftmodels.rds'))
}

seqshiftmodels%>%object.size()%>%divide_by(1e9)

ireads<- samp_rl_reads[[1]][['28']]
seqshiftmodel<-seqshiftmodels[[1]][['28']]
prob=0.5
library(ranger)

get_probforrest_preds <- function(ireads,seqshiftmodel,prob=0.5){
	require(ranger)
	ureads <- unique(ireads)
	testrfdf <- ureads[,c('phase')]%>%
		add_flank_bases(exonsgrl,fafileob)%>%
		mcols%>%
		as.data.frame
	utestrfdf<-testrfdf%>%distinct
	colnames(utestrfdf)[1]<-'phase'
	predmat = predict(seqshiftmodel,utestrfdf)$predictions
	#re-arrange the prob matrix so 0 is first, and col.max can bias us towards inaction
	predmat <- cbind(predmat[,'0',drop=FALSE],predmat[,setdiff(colnames(predmat),'0')])
	ismat = predmat > prob
	bestcols = max.col(ismat,ties.method='first')
	predictions <- as.numeric(colnames(ismat)[bestcols])
	#go from unique sequence matrix rows, to the offset for each unique read	
	predictions <- predictions[vctrs::vec_group_id(testrfdf)]
	#and now for each actual read
	predictions <- predictions[match(ireads,ureads)]
	predictions
}

}

# toptrs
toptrs <- ribocovtrs
# top_samp_rl_reads <- samp_rl_reads%>%map_depth(2,~subset(.,seqnames %in% toptrs))
top_samp_rl_reads <- samp_rl_reads
if(!file.exists(here('data/top_sample_rl_reads.rds'))){
	seqshifts <- mclapply(mc.cores=20,names(top_samp_rl_reads)[1]%>%setNames(.,.),function(isample){
		lapply(names(samp_rl_reads[[isample]]['29'])%>%setNames(.,.),function(rl){
		# for(rl in '28'){
			seqshift<-get_probforrest_preds(
				top_samp_rl_reads[[isample]][[rl]],
				seqshiftmodels[[isample]][[rl]]
			)
			cat('.')
			seqshift
		})
	})

	for(isample in names(seqshifts)){
		for(rl in names(seqshifts[[isample]])){
			cat('.')
			top_samp_rl_reads[[isample]][[rl]]$seqshift <- seqshifts[[isample]][[rl]]
		}
	}

	top_samp_rl_reads%>%saveRDS(here('data/top_samp_rl_reads.rds'))
 }else{
 	top_sample_rl_reads<-readRDS(here('data/top_sample_rl_reads.rds'))
 }

top_samp_rl_reads <- top_samp_rl_reads[1]%>%map(~.['29'])


# top_samp_rl_reads<-readRDS(here('data/top_samp_rl_reads.rds'))

phaseshift<-function(ireads){
	ireads$phase%>%
		table%>%.[as.character(c(0:2))]%>%
		multiply_by(-1)%>%rank%>%
		subtract(1)%>%
		add(0:-2)
}

# tr_target_pos=cdsstarts
# pprops<-ireads%>%
# 	as("GRanges")%>%
# 	shift(.,.$offset+.$seqshift+phaseshift(.)[as.character(.$phase)])%>%
# 	{.$targetpos=tr_target_pos[as.vector(seqnames(.))];.}%>%
# 	addphase(cdsstarts)%>%
# 	.$phase%>%table%>%prop
# pprops


get_utr_exts <- function(trspacecds,n_fp_ext,n_tp_ext){
	trlens <- seqlengths(trspacecds)
	seqnms = names(trlens)
	fputrlens = start(trspacecds[seqnms])-1
	fputrext = pmax(0,n_fp_ext - fputrlens)%>%setNames(seqnms)
	#
	tputrlens = trlens - end(trspacecds[seqnms])
	tputrext = pmax(0,n_tp_ext - tputrlens)%>%setNames(seqnms)
	#
	list(fputrext,tputrext)
}
ext_cov <- function(cov,fputrext,tputrext){
	gr = as(cov,"GRanges")
	seqnms = names(fputrext)
	stopifnot(identical(seqnms,names(tputrext)))
	stopifnot(identical(seqnms,names(cov)))
	seqlengths(gr)[seqnms]%<>%add(fputrext+tputrext)
	gr %<>% shift(fputrext[as.character(seqnames(gr))])
	coverage(gr,weight='score')
}

get_phase_seqshift_psites <-	function(sampreads,reduce_rls=TRUE){
		psites <- sampreads%>%lapply(
			.%>%
			subset(seqnames%in%names(ltrspacecds))%>%
			shift(.,.$offset+.$seqshift+phaseshift(.)[as.character(.$phase)])%>%
			# addphase(cdsstarts)%>%subset(phase!=0)%>%
			resize(1)%>%
			coverage
		)
		if(reduce_rls) psites = purrr::reduce(psites,`+`)
		psites
}

get_cds_bin_counts <- function(
    sampreads,
    trspacecds,
    window_dims=list(
    	n_start_inner = 6,
		n_stop_inner = 6,
		n_fp_ext = 3,
		n_tp_ext = 3,
		psitefunc = get_phase_seqshift_psites
    )
  ){
  	#
	#
	wd=window_dims
	# attach(window_dims)
	STARTWINDSIZE = window_dims$n_start_inner+window_dims$n_fp_ext
	STOPWINDSIZE = window_dims$n_stop_inner+window_dims$n_tp_ext
	MINCDSSIZE = window_dims$n_start_inner+window_dims$n_stop_inner+1
	TOTBINS = (STARTWINDSIZE)+(STOPWINDSIZE)+1
	stopifnot(seqlevels(trspacecds)%in%seqlevels(sampreads[[1]]))
	#
	ltrspacecds = trspacecds
	longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
	ltrspacecds = ltrspacecds[longcdstrs]
	#see how many we eliminate for length reasons
	n_inp_trs <- length(trspacecds)
	n_lon_trs <- length(ltrspacecds)
	message(str_interp('started with ${n_inp_trs} of which ${n_lon_trs} long enough'))
	#creating shifted psite coverage
	# browser()
	# sampreads[[1]]%>%subset(seqnames%in%names(ltrspacecds))%>%
		# shift(.,phaseshift(.)[as.character(.$phase)])%>%
		# addphase(cdsstarts)%>%
		# .$phase%>%table%>%prop
	psitecov = sampreads%>% psitefunc
	#
    message('extending coverage rles for the plot windows')
	ltrspacecds%<>%keepSeqlevels(unique(seqnames(ltrspacecds)))
	c(fputrext,tputrext) %<-% get_utr_exts(ltrspacecds,wd$n_fp_ext,wd$n_tp_ext)
  # 	epsitecovlist  <- psitecovlist%>% 
		# lapply(.%>%.[longcdstrs])%>%
		# lapply(ext_cov,fputrext,tputrext)
  	epsitecovlist  <- psitecov%>%.[longcdstrs]%>%ext_cov(fputrext,tputrext)
  	#
    ltrspacecds = trspacecds
    longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
    ltrspacecds = ltrspacecds[longcdstrs]
    #
    eltrspacecds <- ltrspacecds[longcdstrs]
    seqlengths(eltrspacecds)[longcdstrs] %<>% add(fputrext+tputrext)
    eltrspacecds %<>% shift(fputrext)
    eltrspacecds%<>%resize(width(.)+wd$n_fp_ext,'end')%>%resize(width(.)+wd$n_tp_ext,'start')
    #
    message('gathering window counts')
    midwind = eltrspacecds%>%
      resize(width(.)-(STOPWINDSIZE),'end')%>%
      resize(width(.)-(STARTWINDSIZE),'start')
    midmat = epsitecovlist[midwind]%>%
      sum%>%#compress these variable length Rles down to 1 value per gene
      matrix#make a 1 column matrix
    stwind = eltrspacecds%>%resize(STARTWINDSIZE,ignore.strand=TRUE)
    startmat = epsitecovlist[stwind]%>%as.matrix
    endwind = eltrspacecds%>%resize(STOPWINDSIZE,fix='end',ignore.strand=TRUE)
    endmat = epsitecovlist[endwind]%>%as.matrix
    out = cbind(startmat,midmat,endmat)
    rownames(out) <- longcdstrs
    out
  }
#
#
library(abind)
library(tidyverse)
displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(displaystagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols <- displaystagecols%>%setNames(names(stageconv))
roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}
rownorm <-function(x) x %>%sweep(.,MARGIN=1,F='/',STAT=rowSums(.)%>%{pmax(.,min(.[.>0])/10)})
get_metasignaldf<-function(cds_bin_counts,genelist=TRUE,window_dims){
	wd = window_dims
	n_st_wind = wd$n_fp_ext+wd$n_start_inner
	n_stp_wind = wd$n_tp_ext+wd$n_stop_inner
	metasignaldf = cds_bin_counts%>%
		map_df(.id='sample',.%>%.[genelist,]%>%
			rownorm%>%
			# {./sum(.)}%>%
			colMeans(na.rm=T)%>%
			enframe('start','signal'))%>%
		group_by(sample)%>%
		mutate(section=c(rep('AUG',n_st_wind),'middle',rep('stop',n_stp_wind)))
	#assign stage
	metasignaldf$stage <- metasignaldf$sample%>%str_extract('[^_]+')
	# metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',STARTWINDSIZE+1,TOTBINS- STOPWINDSTART ))
	metasignaldf$start = c(
		-wd$n_fp_ext:(wd$n_start_inner-1),
		NA,
		(-(wd$n_stop_inner-3)):(wd$n_tp_ext+3-1)
	)%>%rep(n_distinct(metasignaldf$sample))
	# metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',1+wd$n_fp_ext,1+wd$n_fp_ext+wd$n_start_inner+wd$n_stop_inner-2 ))
	metasignaldf
}
get_metaplot <- function(metasignaldf,window_dims,ylims=NULL){
	if(!'fraction'%in%colnames(metasignaldf)) metasignaldf$fraction='total'
	if(!'sample'%in%colnames(metasignaldf)) metasignaldf$sample=metasignaldf$stage
	#
	wd=window_dims
	startlims = c(-(wd$n_fp_ext+1),((wd$n_start_inner+1)-1))
	stoplims = c((-((wd$n_stop_inner+1)-3)),(wd$n_tp_ext+3-1+1))
	metasignaldf%>%
	filter(section!='middle')%>%
	split(.,.$section)%>%map( .%>%
		# slice_by(sample,8)%>%
		# filter(position%%3 == 0)%>%
		{
			isfirst = .$section[1]=='AUG'
			qplot(data=.,group=sample,x=start,y=signal,geom='blank')+
			geom_line(aes(color='stage'))+
			stat_identity(geom='bar',aes(y=signal,color=as_factor(start%%3),fill=as_factor(start%%3)),width=I(0.5))+
			scale_x_continuous(name='position (bp)',
				limits=if(.$section[1]=='AUG') startlims else stoplims ,
				minor_breaks=number_ticks,breaks=partial(number_ticks,n=12)
			)+
			# facet_grid( fraction ~ section,scale='free_x')+
			facet_grid( sample ~ section,scale='free_x')+
			scale_color_manual(values=stagecols)+
			scale_y_continuous(name='Mean Psite Count / CDS Total',limits=ylims)+
			theme_bw()+
			# theme_minimal()+
			theme(panel.grid = element_line(color=I('grey')))
		}
	)%>%
	ggarrange(plotlist=.,ncol=2,common.legend=T)
}

# sampreads = top_samp_rl_reads[1]%>%map(~.['28']) 
sampreads = top_samp_rl_reads
# toptrs = samp_rl_reads[[1]]%>%lapply(coverage)%>%lapply(sum)%>%.[[1]]%>%.[names(trspacecds)]%>%sort%>%names%>%tail(500)
window_dims=list(
	n_start_inner = 60,
	n_stop_inner = 60,
	n_fp_ext = 36,
	n_tp_ext = 36
)


if(!file.exists(here('data/psitecov.rds'))){
	psitecov <- top_samp_rl_reads%>%mclapply(mc.cores=10,get_phase_seqshift_psites,reduce_rls=FALSE)
	saveRDS(psitecov,here('data/psitecov.rds'))
}else{
	psitecov<-readRDS(here('data/psitecov.rds'))
}

if(!file.exists(here('data/cds_bin_counts.rds'))){
	cds_bin_counts <- sampreads%>%mclapply(mc.cores=10,get_cds_bin_counts,trspacecds[toptrs],window_dims)
	saveRDS(cds_bin_counts,here('data/cds_bin_counts.rds'))
}else{
	cds_bin_counts<-readRDS(here('data/cds_bin_counts.rds'))
}

metasignaldf_stgrp <- get_metasignaldf(cds_bin_counts,TRUE,window_dims) %>% 
	filter(str_detect(sample,'ribo'))%>%
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal))
{
'plots/psite_redo/'%>%dir.create(showWarn=F,rec=T)
library(rlang)
plotfile<-'plots/psite_redo/framemetaplots.pdf'%T>%pdf(h=20,w=12)
	rwplot <- metasignaldf_stgrp%>%get_metaplot(window_dims)
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

nonmainsamps<-names(cds_bin_counts)%>%str_subset(neg=T,'ribo')
fracmetasignaldf_stgrp <- get_metasignaldf(cds_bin_counts[nonmainsamps],TRUE,window_dims) %>% 
	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
	mutate(stage = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))

{
library(rlang)
plotfile<-'plots/psite_redo/frac_framemetaplots.pdf'%T>%pdf(h=12,w=12)
rwplot <- metasignaldf_stgrp%>%filter(fraction=='80S')%>%get_metaplot(window_dims)
rwplotpoly <- metasignaldf_stgrp%>%filter(fraction=='Poly')%>%get_metaplot(window_dims)
print(ggarrange(nrow=2,plotlist=list(rwplot,rwplotpoly)))
dev.off()
normalizePath(plotfile)%>%message
}


{
library(rlang)
plotfile<-'plots/psite_redo/poly_test.pdf'%T>%pdf(h=6,w=12)
print(rwplot)
dev.off(

)
normalizePath(plotfile)%>%message
}

{
library(rlang)
plotfile<-'plots/psite_redo/frac_test.pdf'%T>%pdf(h=12,w=12)
rwplot <- metasignaldf_stgrp%>%filter(fraction=='80S')%>%get_metaplot_stagecol(window_dims)
rwplotpoly <- metasignaldf_stgrp%>%filter(fraction=='Poly')%>%get_metaplot_stagecol(window_dims)
print(ggarrange(nrow=2,plotlist=list(rwplot,rwplotpoly)))
dev.off()
normalizePath(plotfile)%>%message
}

{
library(rlang)
plotfile<-'plots/psite_redo/poly_test.pdf'%T>%pdf(h=6,w=12)
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

sampfpcov=psitecov[[1]]
if(!file.exists(here('data/ptrsums.rds'))){
	ptrsums <- 	mcmapply(mc.cores=8,SIMPLIFY=F,psitecov,names(psitecov),FUN=function(sampfpcov,sampname){
				sampfpcov = sampfpcov
				sampfpcov%>%map(regScoreSums,innercds)%>%purrr::reduce(.,`+`)
			})
	saveRDS(ptrsums,here('data/ptrsums.rds'))
}else{
	ptrsums<-readRDS(here('data/ptrsums.rds'))
	stopifnot(names(ptrsums)==names(psitecov))
	stopifnot(names(ptrsums[[1]])==names(innercds))
}
psitecodonmats <- 	mcmapply(mc.cores=4,SIMPLIFY=F,psitecov[1:10],names(psitecov)[1:10],FUN=function(samp_pcov,sampname){
	ptrsum = ptrsums[[sampname]]
	samp_pcov%>%lapply(safely(function(rlfpcov){
		rlfpcov = rlfpcov[names(ptrsum)]
		rlfpcov = rlfpcov/(ptrsum)
		('.')
		out = rlfpcov[allcodlistnz]
		stopifnot(out%>%max(na.rm=T)%>%max(na.rm=T)%>%is.finite)
		out = out %>%split(cods)%>%lapply(as.matrix)
		out = lapply(out,Matrix::Matrix,sparse=TRUE)
		out
	}))
})
psitecodonmats%<>%map_depth(2,'result')
saveRDS(psitecodonmats,'data/psitecodonmats.rds')
psitecodonmats<-readRDS('data/psitecodonmats.rds')

row_gnames <- seqnames(allcodlistnz)%>%split(cods)
matlist=psitecodonmats[[1]][[1]]
.x='AAC'
get_glist_profiles<-function(glist,row_gnames,psitecodonmats){
	matlist = psitecodonmats[[1]][[1]]
	.x='AAA'
	ucods <- unique(cods)%>%setNames(.,.)
	psitecodonmats%>%map_depth(2,function(matlist)map(ucods,function(cod){
		matlist[[cod]][as.vector(row_gnames[[cod]])%in%unique(glist),]%>%colMeans(na.rm=T)
	}))
}
glist=dtselgenelist[['nochange']]
get_codon_prof_df <- function(glist,row_gnames,psitecodonmats){
	codonprofiledat <- get_glist_profiles(glist,row_gnames,psitecodonmats)	
	codonprofiledat <- codonprofiledat%>%
		map_depth(3,.%>%
			enframe('position','count'))%>%
			map_df(.id='sample',.%>%
				map_df(.id='readlen',.%>%
					bind_rows(.id='codon')))
	# codonprofiledat <- codonprofiles%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
	codonprofiledat %<>% mutate(position = position - 1 - (FLANKCODS*3))
	# codonprofiledat %<>% group_by(subgroup,sample,readlen,codon)%>%mutate(count= count / median(count))
	codonprofiledat %<>% group_by(sample,readlen,codon)%>%
		# mutate(count= count / median(count))
		identity
	codonprofiledat %<>% filter(!codon %in% c('TAG','TAA','TGA'))
	codonprofiledat$readlen%<>%str_replace('^(\\d)','rl\\1')
	codonprofiledat
}
sub_p_profilelist <- dtselgenelist['allhigh']%>%lapply(get_codon_prof_df,row_gnames,psitecodonmats)
# subfpprofilelist[['nochange']]




# subfpprofilelist%>%names
geneset='allhigh'
codonvarprofiles <-  
	sub_p_profilelist[[geneset]]%>%
# rustprofiledat%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	# mutate(count = count / mean(count))%>%
	# filter(count==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	# filter(sample%>%str_detect('ribo_1'))%>%
	filter(!is.nan(count))%>%
	# summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(rep,time,readlen,position)%>%
	# mutate(count = count / mean(count))%>%
	summarise(sdsig=sd(count,na.rm=T))%>%
	group_by(readlen,time,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	filter(-18<position,position < 12)%>%
	# filter(position> -numreadlen+6,position< -6)%>%
	filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)
#
codonvaroffsets = codonvarprofiles%>%
	filter(-9<position,position < 9)%>%
	group_by(readlen,time)%>%slice(which.max(sdsig))%>%
	mutate(offset=position+3)%>%
	as.data.frame
#
codonvaroffsets = sub_p_profilelist[[geneset]]%>%
	ungroup%>%
	distinct(sample)%>%
	separate(sample,c('time','assay','rep'),remove=F)%>%
	left_join(codonvaroffsets)%>%
	select(sample,readlen,numreadlen,offset)%>%
	separate(sample,c('time','assay','rep'),remove=F)


{
# codonvaroffsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_p_pos_vs_codon_variance.pdf'
pdf(plotfile,w=12,h=3*n_distinct(codonvarprofiles$readlen))
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	filter(time%>%is_in(names(stagecols)))%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(codonvaroffsets,readlen%in%.$readlen),aes(xintercept= offset),color=I('blue'),linetype=2)+
		geom_vline(data=filter(codonvaroffsets,readlen%in%.$readlen),aes(xintercept= offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}

{
offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_p_pos_vs_codon_variance_zoom.pdf'
pdf(plotfile,w=6,h=3*1)
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	filter(time%>%is_in(names(stagecols)))%>%
	filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29'))%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		# geom_vline(data=codonvaroffsets2plot%>%filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29')),aes(xintercept= -offset),color=I('blue'),linetype=2)+
		# geom_vline(data=codonvaroffsets2plot%>%filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29')),aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}



################################################################################
########## Now check against tRNA levels
################################################################################
if(!exists('allcodsig_isomerge')) base::source(here('src/Figures/Figure3/3_tRNA_array_analysis.R'))
allcodsig_isomerge%<>%addcodon
trna_ab_df_samp = allcodsig_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]
#get data on the trnas
trna_dat <- trna_ab_df_samp%>%
	select(-sample)%>%
    select(fraction,time,rep,codon,abundance,availability)%>%
    group_by(fraction,time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean),na.rm=T)
#
clean_fr_sampnames<-function(x) x%>%str_replace('.*_(Poly|80S)(.*)_()','\\2_\\1ribo_\\3')
rloffsets<-offsets%>%select(readlen,length,offset)%>%mutate(readlen=str_replace(length,'^(\\d)','rl\\1'))
lexp = 3+3 #include positions for positions corresponding to bigger offsets
rexp = 3#include positions for positions corresponding to smaller offsets
# codonprofiles%>%
codondata <- 
	sub_p_profilelist[['allhigh']]%>%
	# codonprofiledat%>%
	# filter(subgroup=='nochangehighe')%>%
	# filter(subgroup=='all')%>%
# codondata <- rustprofiledat%>%
	# filter(sample%>%str_detect('ribo'))%>%
	# mutate_at(vars(sample),clean_fr_sampnames)%>%
    group_by(sample)%>%
    # mutate(readlen=paste0('rl',length))%>%
    # safe_left_join(rloffsets)%>%
    safe_left_join(codonvaroffsets%>%mutate_at('sample',clean_fr_sampnames))%>%
    separate(sample,c('time','assay','rep'))%>%
    # filter(readlen=='rl29')%>%
    # filter(readlen%in%c('rl29','rl30','rl28','rl27'))%>%
    # filter(readlen%in%c('rl29','rl31','rl30','rl27','rl26'))%>%
    # mutate(offset=0)%>%#for the shifted data
    # filter(position <= -(offset-rexp))%>%
    # filter(position >= -(offset+lexp))%>%
    # filter(position == -offset)%>%
    filter(between(position, -offset-5,-offset))%>%
    # filter(between(position, -offset-5,-offset-3))%>%
    group_by(assay,time,codon,rep)%>%
    summarise(
    	inf_p_site_occ= sum(count[position%>%between(-2,0)]),
    	p_site_occ = sum(count[position%>%between(-offset-2,-offset)]),
    	a_site_occ = sum(count[position%>%between(-offset-3,-offset-3)])
    )%>%
    arrange(assay!='ribo')
#seperate poly and total tRNA info
total_trnadat <- trna_dat%>%filter(fraction=='Total')%>%select(time,codon,abundance,availability)
poly_trnadat <- trna_dat%>%filter(fraction=='Poly')%>%
	select(time,codon,poly_abundance=abundance,poly_availability=availability,inf_p_site_occ)
#add tRNA info, 
#AA and aa corrected dwell time
codondata%<>%
    left_join(total_trnadat,by=c('time','codon'))%>%
    left_join(poly_trnadat,by=c('time','codon'))%>%
    mutate(AA=GENETIC_CODE[codon])%>%
    group_by(assay,time,rep,AA)%>%
    arrange(assay!='ribo')%>%
    mutate(aacor_p_site_occ=ifelse(n()==1,NA,p_site_occ-mean(p_site_occ)))
#summarise the replicates
repsumcodondata<-codondata%>%
	group_by(assay,time,codon,AA)%>%
	summarise_at(vars(abundance,availability,poly_abundance,poly_availability,a_site_occ,p_site_occ,aacor_p_site_occ,inf_p_site_occ),mean)%>%
    arrange(assay!='ribo')
#check all codons there
stopifnot(codondata%>%group_by(assay,time,rep)%>%tally%>%.$n%>%`==`(61)%>%all)
stopifnot(repsumcodondata%>%group_by(assay,time)%>%tally%>%.$n%>%`==`(61)%>%all)

# codondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$availability)}
repsumcodondata%>%filter(assay=='ribo',time=='E13')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$a_site_occ,.$availability)}
repsumcodondata%>%filter(assay=='ribo',time=='E13')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$availability)}
# codondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='Polyribo')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$poly_abundance)}


codondata%>%lm(data=.,a_site_occ~AA)%>%aov%>%summary%>%as.list%>%.[[1]]%>%.[[2]]%>%{.[1]/sum(.)}
codondata%>%lm(data=.,p_site_occ~AA)%>%aov%>%summary%>%as.list%>%.[[1]]%>%.[[2]]%>%{.[1]/sum(.)}
codondata%>%lm(data=.,inf_p_site_occ~AA)%>%aov%>%summary%>%as.list%>%.[[1]]%>%.[[2]]%>%{.[1]/sum(.)}

# best_samp_rl_reads
# 	group_by(sample,length)%
# 	slice(((which.max(twind)-1):(which.max(twind)+1)))%>%
# 	group_by(length,phase)%>%
# 	nest()%>%
# 	summarise(
# 		offset=map(data,~table(.$cor_offset)%>%enframe('val','freq'))
# 	)%>%
# 	unnest(offset)%>%
# 	as.data.frame%>%
# 	arrange(length,phase)

		# psitecov <- ireads%>%resize(1)%>%{shift(.,.$offset)}
		# psitecov%<>%{.$error=(.$targetpos+.$phase)-start(.);.}
		# phs_psitecov%>%coverage(weight='score')

		# # iphase=1
		# trs = seqnames(ireads)%>%unique
		# traintrs = trs%>%sample(.,floor(length(.)*0.8))
		# testtrs = trs%>%setdiff(traintrs)
		# # ireads%>%subset(seqnames%in%traintrs)%>%subset(phase==iphase)%>%subset(between(error,-6,6)) -> ireads
		# rfdf <- ireads%>%subset(seqnames%in%traintrs)%>%
		# 	# subset(phase==iphase)%>%
		# 	subset(between(error,-6,6))%>%
		# 	add_flank_bases(exonsgrl,fafileob)%>%
		# 	mcols%>%as.data.frame%>%
		# 	# select(-phase,-mapq,-length,-offset)
		# 	select(-offset)
		# rfdf$error%<>%as_factor
		# rfdf%<>%select(error,phase,everything())
		# nvars = length(colnames(rfdf))-2
		# selweights = c(0,rep(1/nvars,nvars))
		# message('training probabalistic random')
		# seqshiftmodel <- ranger::ranger(
		# 	formula= error ~ . ,
		# 	data=rfdf,
		# 	# importance='permutation',
		# 	probability=TRUE,
		# 	split.select.weights=selweights
		# )
		# testrds<- ireads%>%subset(seqnames%in%testtrs)%>%
		# 	# subset(phase==iphase)%>%
		# 	identity
		# testrfdf <- testrds%>%
		# 	add_flank_bases(exonsgrl,fafileob)%>%
		# 	mcols%>%as.data.frame
		# # testrds$seqshift <- predict(seqshiftmodel,testrfdf)$predictions
		# predmat = predict(seqshiftmodel,testrfdf)$predictions
		# predmat <- cbind(predmat[,'0',drop=FALSE],predmat[,setdiff(colnames(predmat),'0')])
		# predmat%>%colnames
		# ismat = predmat > 0.5
		# predictions <- as.numeric(colnames(ismat)[max.col(ismat,ties.method='first')])
		# testrds$seqshift <- predictions
		# #convert from factor back to integer
		# testrds$seqshift%<>%as.character%<>%as.numeric
		# #
		# testrds%>%subset(between(error,-6,6))%>%mcols%>%as.data.frame%>%
		# 	transmute(offset,error,cor_offset=offset+error,seqoffset=offset+seqshift)%>%
		# 		mutate(aligned=cor_offset==offset)%>%
		# 		mutate(seqaligned=cor_offset==seqoffset)%>%
		# 		group_by(aligned,seqaligned)%>%tally

		# preseqshifttots=testrds%>%resize(1)%>%shift(.$offset)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# preseqshifttots%>%prop%>%txtplot(c((-6 : -1),0:8),.,ylim=0:1)
		# postseqshifttots=testrds%>%resize(1)%>%shift(.$offset+.$seqshift)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# postseqshifttots%>%prop%>%txtplot(c((-6 : -1),0:8),.,ylim=0:1)

		# ireads


		# testrds%>%subset(between(error,-6,6))%>%mcols%>%as.data.frame%>%
		# 	transmute(cor_offset,offset,error,seqoffset=offset+seqshift)%>%
		# 	mutate(aligned=cor_offset==offset)%>%
		# 	mutate(seqaligned=cor_offset==seqoffset)%>%
		# 	group_by(aligned,seqaligned)%>%tally


		# psitecov$error%>%table%>%{.[between(as.numeric(names(.)),-rln,rln)]}%>%{./sum(.)}%>%round(3)
		# psitecov%>%subset(phase==0)%>%.$error%>%table%>%{.[between(as.numeric(names(.)),-rln/2,rln/2)]}%>%{./sum(.)}%>%round(3)
		# psitecov%>%subset(phase==1)%>%.$error%>%table%>%{.[between(as.numeric(names(.)),-rln/2,rln/2)]}%>%{./sum(.)}%>%round(3)
		# psitecov%>%subset(phase==2)%>%.$error%>%table%>%{.[between(as.numeric(names(.)),-rln/2,rln/2)]}%>%{./sum(.)}%>%round(3)


		# trspacecds%<>%setstrand
		# startwinds%<>%setstrand
		# psitecov%<>%subset(strand=='+')
		# psitecov%>%subset(error==0)%>%overlapsAny(trspacecds%>%setstrand())%>%table

		# testrange = phs_psitecov%>%subset(error==0)%>%subsetByOverlaps(invert=TRUE,trspacecds)%>%head(1)

		# #so why this weird profile around my start windows.
		# possums = psitecov%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# possums %>% txtplot(c((-6 : -1),0:8),.)

		# errzeropossums = psitecov%>%subset(error==0)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# errzeropossums %>% txtplot(c((-6 : -1),0:8),.)

		# errthreepossums = psitecov%>%subset(error==3)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# errthreepossums %>% txtplot(c((-6 : -1),0:8),.)


		# #okay so lots of other things contribute to the coverage in these start windows - I'm confused about something.
		# noterrzeropossums = psitecov%>%subset(error!=0)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# noterrzeropossums %>% txtplot(c((-6 : -1),0:8),.)

		# #okay so this overlaps a start window even though it souldn't....
		# weirdread = psitecov%>%subsetByOverlaps(startwinds)%>%tail(5)%>%head(1)
		# startwinds%>%subsetByOverlaps(weirdread)


		# possums = phs_psitecov%>%resize(1)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
		# # possums = phs_psitecov%>%subset(cor_offset<phasevect[as.character(.$phase)])%>%resize(1)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums

		# startreads = psitecov%>%resize(1)%>%subsetByOverlaps(startwinds)
		# startreads

		# testread = startreads%>%subset(error>30)%>%head(1)

		# psitecov%>%coverage()



# tputr_trs=(tputrs%>%width%>%sum > 30 )%>%names(.)[.]
# fputr_trs=(fputrs%>%width%>%sum > 30 )%>%names(.)[.]
# train_trs = intersect(tputr_trs,fputr_trs)%>%intersect(ribocovtrs)

# trainreads = fpcov[train_trs]%>%
# 		GRanges%>%subset(score!=0)%>%
# 		resize(rln)
# 		# .[!is_out_of_bounds(.)]
# trainreads%>%resize()



# #checks out

# # trainreads%>%subset(seqnames=='ENSMUST00000020909')%>%start%>%min
# get_exp_exonseq




# # trainreads%>%head(1)%>%resize(nbp*2,'end')%>%shift(nbp)%>%getSeq(.,x=exonseq)
# # tpseqmat%>%head(1)

# phs_psitecov%>%add_flank_bases(exonsgrl,fafileob,nbp=2)

# stopifnot('fp_0' %in% mcols(phs_psitecov))

# phs_psitecov%<>%setNames(refends)
# #
# phs_psitelist %<>% lapply(.%>%{Reduce(`+`,.)})
# phs_psitelist %<>% lapply(.%>%{ . / (sum(.)+1e-9)})


# startwinds = trspacecds%>%resize(1)%>%resize(3+6)%>%resize(width(.)+6,'end')%>%
# 	.[!is_out_of_bounds(.)]
# stopwinds = trspacecds%>%resize(1,'end')%>%resize(3+6,'end')%>%resize(width(.)+6)%>%
# 	.[!is_out_of_bounds(.)]
# #align over prev MORE than start!
# phs_psitecov[['start']][['26']]%>%.[startwinds]%>%as.matrix%>%colSums%>%{txtplot(c(-6:8),.)}
# #align over prev MORE than start!
# phs_psitecov[['start']][['27']]%>%.[startwinds]%>%as.matrix%>%colSums%>%{txtplot(c(-6:8),.)}
# #align over both start and prev, though start better
# phs_psitecov[['start']][['28']]%>%.[startwinds]%>%as.matrix%>%colSums%>%{txtplot(c(-6:8),.)}
# #align over start only
# phs_psitecov[['start']][['29']]%>%.[startwinds]%>%as.matrix%>%colSums%>%{txtplot(c(-6:8),.)}




# phs_psitecov = fpcov[ribocovtrs]%>%
# 	GRanges%>%subset(score!=0)%>%
# 	{.$cdsstart=cdsstarts[as.vector(seqnames(.))];.}%>%
# 	{.$cor_offset=start(.)-.$cdsstart;.}%>%
# 	addphase(cdsstarts)
# # phs_psitecov%>%coverage(weight='score')

# # phasevect[phs_psitecov$phase]%>%length
# # phs_psitecov$phase%>%table
# # phasevect[as.character(phs_psitecov$phase)]
# # 	shift(rln-1)%>%
# # 	trim()%>%
# # 	coverage(weight='score')

# tpcov = fpcov[ribocovtrs]%>%GRanges%>%subset(score!=0)%>%
# 	{.$cdsstart=cdsstarts[seqnames(gr)];.}%>%
# 	{.$cor_offset=start(.)-.$cdsstart;.}%>%
# 	addphase(cdsstarts)
# 	shift(rln-1)%>%
# 	trim()%>%
# 	coverage(weight='score')

# fpstartwinds = trspacecds[ribocovtrs%>%head(10000)]%>%
# 	resize(9,'start')%>%
# 	resize(9+6,'end')%>%
# 	shift(-rloffset)%>%
# 	trim
# fpstartwinds = fpstartwinds[width(fpstartwinds)==9+6]
# fpcov[names(fpstartwinds)][fpstartwinds]%>%as.matrix%>%colSums%>%txtplot


# fpstopwinds = trspacecds[ribocovtrs%>%head(10000)]%>%
# 	resize(9,'end')%>%
# 	resize(9+6,'start')%>%
# 	shift(-rloffset)%>%
# 	trim
# fpstopwinds = fpstopwinds[width(fpstopwinds)==9+6]
# fpcov[names(fpstopwinds)][fpstopwinds]%>%as.matrix%>%colSums%>%txtplot


# #suppose rl is 10. 
# tpstartwinds = trspacecds[ribocovtrs%>%head(1000)]%>%
# 	resize(9,'start')%>%
# 	resize(9+6,'end')%>%
# 	shift(+(rln-rloffset-1))%>%
# 	trim
# tpstartwinds = tpstartwinds[width(tpstartwinds)==9+6]
# tpcov[names(tpstartwinds)][tpstartwinds]%>%as.matrix%>%colSums%>%txtplot


# fpcov%>%GRanges%>%subset(score!=0)%>%shift(11)
# fpcov%>%GRanges%>%subset(score!=0)%>%resize(rln)%>%.[!is_out_of_bounds(.)]%>%
# 	resize(1,'end')%>%shift(-(29-11-1))


# tpstopwinds = trspacecds[ribocovtrs%>%head(1000)]%>%
# 	resize(9,'end')%>%
# 	resize(9+6,'start')%>%
# 	shift(+(rln-rloffset-1))%>%
# 	trim
# tpstopwinds = tpstopwinds[width(tpstopwinds)==9+6]
# tpcov[names(tpstopwinds)][tpstopwinds]%>%as.matrix%>%colSums%>%txtplot


# #so e.g. 29 shows 3 pretty clear spots aligning on the start codon, ncluding
# #at all timepoints. It also kind of looks like the phase distribution is screwed though

# ################################################################################
# ########Testing phase specific alighment with PPG motifs
# ################################################################################
# # phs_psitelist->tp_phs_psitelist

# matcolvalues = map_df(.id='peptide',peps,function(ipep){
# 	phs_psitelist%>%
# 	map_df(.id='sample',~.[tr_peptidegr%>%subset(peptide==ipep)%>%resize(9,'start')]%>%
# 		as.matrix%>%
# 		colMeans(na.rm=T)%>%
# 		enframe('position','mean_rdens')
# 	)
# })
# if('1' %in% matcolvalues$peptide) matcolvalues$peptide <- peps[as.numeric(matcolvalues$peptide)]
# #now plot
# for(ipep in peps){
# 	plotfile<- here(paste0('plots/',ipep,'_motif_Psite_Alignment','.pdf'))
# 	pdf(plotfile)
# 	pepchrs = str_split(ipep,'')%>%unlist
# 	p=matcolvalues%>%
# 		filter(peptide==ipep)%>%
# 		separate(sample,c('time','assay','rep'))%>%
# 		ggplot(.,aes(x=position,y=mean_rdens,color=time,group=paste0(time,rep)))+
# 		geom_line()+
# 		# scale_color_manual(values=stagecols)+
# 		scale_x_continuous(paste0('Position'),breaks=1:9,
# 			labels=c('',pepchrs[1],'','',pepchrs[2],'','',pepchrs[3],''))+
# 		scale_y_continuous(paste0('Mean CDSnormed Ribosome Density'))+
# 		geom_vline(xintercept=4,linetype='dashed')+
# 		ggtitle(paste0(ipep,' Motif P-site alignment'))+
# 		theme_bw()
# 	print(p)
# 	dev.off()
# 	message(normalizePath(plotfile))
# }

# ################################################################################
# ########Turns out some reads align only partially ont he transcript
# ########
# ################################################################################
	
# #why out of bounds ranges tho
# oobtr = 'ENSMUST00000021609'
# fpcovlist[[samp1]][[rl]]%>%GRanges%>%subset(seqnames==oobtr)
# trlens[oobtr]
# fpcovlist[[samp1]][[rl]]%>%GRanges%>%subset(seqnames==oobtr)%>%.[28]%>%shift(24)
# trlens[oobtr]==2667
# probpos=2644
# #this can't be, cos it's a 25 bp read length there..

# fpcovlist

# #let's look at the original data
# bamtbl<-"pipeline/deepshapebamdata/E13_ribo_1.bam.reformat"
# bamtbl%>%
# 	str_interp('grep -e $')%>%
# 	fread(select=c(2,6,7))%>%
# 	# filter(V2%>%str_detect(oobtr))
# 	mutate(V2=trimids(V2))%>%
# 	filter(V2==oobtr)
# #yeah this for sure has one there.
# #let's look at the tr bam
# fulloobtr='ENSMUST00000021609.8'


#okay so that out of bounds guy is just a bad alignment, how would one filter