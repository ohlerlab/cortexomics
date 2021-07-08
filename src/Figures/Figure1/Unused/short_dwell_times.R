bamtbl=mainbamtbls[1]
if(!file.exists(here('data/sh_fpcovlist.rds'))){

trlens = exonsgrl%>%width%>%sum
trseqinfo = Seqinfo(seqnames=ribocovtrs,seqlengths=trlens[ribocovtrs])

	
ribobams =c(
	'pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.bam',
	'pipeline/star_transcript/data/E13_ribo_2/E13_ribo_2.bam'
)%>%setNames(.,basename(dirname(.)))

ribobams[1] -> ribobam

if(!file.exists(here('data/sh_fpcovlist.rds'))){
	sh_fpcovlist = ribobams%>%mclapply(mc.cores=10,function(ribobam){
		#
		ribogr <- GenomicAlignments::readGAlignments(ribobam)
		mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
		#
		ribogr %<>% subset(between(readlen,19,24))
		ribogr %<>% as("GenomicRanges")
		ribogr %<>% as.data.frame%>%
			mutate(seqnames=str_extract(seqnames,'\\w+'))%>%
			filter(seqnames %in% ribocovtrs)%>%
			GRanges(seqinfo=trseqinfo)
		ribogr %<>% resize(1,'start')
		split(ribogr,ribogr$readlen)%>%
			lapply(coverage)
	})
	names(sh_fpcovlist)<-names(ribobams)
	saveRDS(sh_fpcovlist,here('data/sh_fpcovlist.rds'))
}else{
	sh_fpcovlist<-readRDS(here('data/sh_fpcovlist.rds'))
	stopifnot(names(sh_fpcovlist)==names(ribobams))
	stopifnot(all(ribocovtrs%in%names(sh_fpcovlist[[1]][['21']])))
	inclusiontable(ribocovtrs,names(sh_fpcovlist[[1]][['21']]))
}

allcodlist<-readRDS(here('data/allcodlist.rds'))
stopifnot(allcodlist@seqinfo@seqnames%>%setequal(ribocovtrs))
stopifnot(allcodlist@seqinfo@seqlengths%>%setequal(sum(width(exonsgrl[ribocovtrs]))))

sh_fpcovlist <- sh_fpcovlist[names(ribobams)]
sh_fpcovlist[names(ribobams[1])]
if(!file.exists(here('data/sh_fpprofilelist.rds'))){
	sh_fpprofilelist <-imap(sh_fpcovlist[names(ribobams)],function(sampfpcov,sampname){
		trsums = sampfpcov%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
		sampfpcov%>%lapply(function(rlfpcov){
			rlfpcov = rlfpcov/(trsums)
			allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
			cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
			message('.')
			out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
			out
		})
	})
	saveRDS(sh_fpprofilelist,here('data/sh_fpprofilelist.rds'))
}else{
	sh_fpprofilelist<-readRDS(here('data/sh_fpprofilelist.rds'))
}


sh_codprof = sh_fpprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
sh_codprof%<>%mutate(position = position - 1 - (FLANKCODS*3))
sh_codprof%<>%mutate(countnonorm=count)
sh_codprof%<>%group_by(sample,readlen,codon)%>%mutate(count= count / median(count))
sh_codprof%<>%filter(!codon %in% c('TAG','TAA','TGA'))

shortoffsets <- tibble(
	offset=c(8),
	compartment='nucl',
	length=21,
)%>%mutate(readlen=paste0('rl',length))

{
# offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_fppos_vs_codon_variance.pdf'
pdf(plotfile,w=12,h=12)
#plotting variance amongst codons at each point.
sh_codprof%>%
# rustprofiledat%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	# filter(count==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	# filter(sample%>%str_detect('ribo_1'))%>%
	filter(!is.nan(count))%>%
	summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(readlen,time,assay,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	filter(position> -numreadlen+3,position < -3)%>%
	# filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}

trna_ab_df_samp = allcodsigmean_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]

codon_data <- trna_ab_df_samp%>%
	filter(fraction=='Total')%>%select(-fraction,-sample)%>%
    select(time,rep,codon,abundance,availability)%>%
    group_by(time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean))

sh_codprof%>%
	filter(position> -as.numeric(readlen)+6,position < -6)%>%
	separate(sample,c('time','assay','rep'))%>%
	left_join(codon_data)%>%
	group_by(position,time,rep)%>%
	nest%>%
	mutate(trnacor=map(data,~tidy(cor.test(.$abundance,.$count))))%>%
	unnest(trnacor)%>%
	arrange(p.value)%>%
	arrange(position)%>%.$estimate%>%txtplot




