

# ################################################################################
# ########Test the quality of the model
# ################################################################################
	


#' ## Testing
#' 
#' We can see if our procedure worked by quickly looking at the fft (fast fourier transform) on the top1k genes
#' binned randomly into 5 groups, before and after the application of our sequence specific cutoffs (to ALL reads),
#' rather than just those around the stop codon

get_periodicity_scores<-function(cdstosamp,reads,topcdsmap,cds,seqshiftisneg=FALSE){
        message('.')

        assert_that(all(has_name(mcols(reads),c('seqshift','cdsshift'))))

        #subset the cds

        cdsmap2use <- topcdsmap[cdstosamp]
        cds2use <- subset(cds,transcript_id%in%cdstosamp)

        cdspsites <- reads%>%
				apply_psite_offset(c('cdsshift'))%>%as("GRanges")%>%
				mapToTranscripts(cds2use%>%split(.,.$transcript_id))%>%
				.[!duplicated(.$xHits)]

        seqcdspsites <- reads%>%
				apply_psite_offset(c('cdsshift','seqshift'))%>%as("GRanges")%>%
				mapToTranscripts(cds2use%>%split(.,.$transcript_id))%>%
				.[!duplicated(.$xHits)]

        #we need psites in transcript space here
        no_seqshift = coverage(cdspsites)%>%
        	{suppressWarnings({unlist(.)})}%>%
        	as.vector%T>%
        	{is3mod=(length(.)%%3 )==0;stopifnot(is3mod)}%>%
        	fft%>%abs%>%`^`(2)%>%{sum(get_frac_inds(.,0.3,0.36))/sum(.) }

        with_seqshift = coverage(seqcdspsites)%>%
        	{suppressWarnings({unlist(.)})}%>%
        	as.vector%T>%
        	{is3mod=(length(.)%%3 )==0;stopifnot(is3mod)}%>%
        	fft%>%abs%>%`^`(2)%>%{sum(get_frac_inds(.,0.3,0.36))/sum(.) }

        data_frame(no_seqshift,with_seqshift)
}

data4readshift <- get_seqforrest_data(topcdsreads,FaFile(REF),nbp=2)

get_predictedshifts(seqshiftmodel_GAonly,data4readshift,prob=0.5)%>%table
mcols(topcdsreads)$seqshift <-  get_predictedshifts(seqshiftmodel_GAonly,data4readshift,prob=0.5)

seqshift_periodicities<- seqnames(topcdsmap)%>%
	unique%>%
	sample%>%		# head(300)%>%
	as.character%>%
	split(seq_along(.)%%5)%>%
	mclapply(F=get_periodicity_scores,topcdsreads,topcdsmap,cds,seqshiftisneg=FALSE)%>%
	bind_rows


t.test(seqshift_periodicities$no_seqshift,seqshift_periodicities$with_seqshift)






# })()



#+ spectral_coefficient_strip_plot, fig.width =4,fig.height=4,out.width=400,out.height=450,dev='pdf',include=TRUE,eval=TRUE

seqshift_periodicities%>%unlist%>%txtplot

stop()

{

cdsread_trmap <- topcdsreads%>%
	# sample(10e3)%>%
	GRanges%>%
	mapToTranscripts(topcdsexons%>%split(.,.$transcript_id))
mcols(cdsread_trmap) <- mcols(topcdsreads)[cdsread_trmap$xHits,]

setpos<- .%>%{strand(.)<-'+';.}
prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')
prestarts <- topcdsmap%>%setpos%>%resize(3,'start')%>%resize(1,'start')


cdsread_trmap$startdist<-start(cdsread_trmap) - start(prestarts[seqnames(cdsread_trmap)])
cdsread_trmap$stopdist<-start(cdsread_trmap) - start(prestops[seqnames(cdsread_trmap)])

#Get riboqc_cutoffs
sample = bam%>%dirname%>%basename 
riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
riboqcdf <- riboqccutoffs%>%fread%>%select(length=read_length,riboqc_shift=cutoff)

starts <- topcdsmap%>%setpos%>%resize(1,'start')%>%setNames(seqnames(.))
cdsread_trmap$startdist <- start(cdsread_trmap) - start(starts[cdsread_trmap@seqnames])

stops <- topcdsmap%>%setpos%>%resize(3,'end')%>%fpend%>%setNames(seqnames(.))
cdsread_trmap$stopdist <- start(cdsread_trmap) - start(stops[cdsread_trmap@seqnames])


fcods=20
lflank=(fcods*3)
rflank=(fcods*3)+2
mpoint=((fcods*2)+1)*3
epoint=((fcods*4)+2)*3

# cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%hist
# cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%between(-lflank,rflank)%>%table
metaplotlabs<-c(paste0('-',lflank/2),'AUG',lflank/2,paste0('mid -',lflank/2),'mid',paste0('mid +',lflank/2),paste0('end -',lflank/2),'stop',paste0('end +',lflank/2))
metaplotbreaks<-c(-lflank/2,0,lflank/2,mpoint-lflank/2,mpoint,mpoint+lflank/2,epoint-(lflank/2),epoint,epoint+(lflank/2))
#

seqshiftfuncs <- list(
	Riboqc = .%>% safe_left_join(riboqcdf,by=c('length'))%>% mutate(startdist=startdist+riboqc_shift,stopdist=stopdist+riboqc_shift),
	CDSmax_shift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift),
	seqshift = .%>% mutate(startdist=startdist+cdsshift+seqshift,stopdist=stopdist+cdsshift+seqshift)
)

#

#' We can also simply look at a metaplot of the frames, (comparing riboqc and cdsmax shift)
plist=lapply(seqshiftfuncs,mymemoise(function(seqshiftfunc){
	#
	p=cdsread_trmap%>%
		# sample(100e3)%>%
		# shift(.$cdsshift)%>%
		mcols()%>%as_tibble	%>%
		select(startdist,stopdist,length,cdsshift,seqshift)%>%
		# filter(between(abs(startdist-stopdist),-lflank,rflank))%>%
		seqshiftfunc%>%
		mutate(
			mdist = startdist - ((floor(((startdist-stopdist+1)/3)/2))*3),
			instart=between(startdist,-lflank,rflank),
			inend=between(stopdist,-lflank,rflank),
			inmiddle=between(mdist,-lflank,rflank)
		)%>%
		filter((instart|inend|inmiddle) & (!(instart&inend&inmiddle)))%>%
		# filter(instart)%>%
		mutate(x = NA,
			x=ifelse(instart,startdist,x),
			x=ifelse(inmiddle,mdist+mpoint,x),
			x=ifelse(inend,stopdist+epoint,x) ) %>%
		mutate(length=as.character(length))%>%
		bind_rows(.,mutate(.,length='all'))%>%
		group_by(x,length)%>%tally%>%
		# filter(length=='all')%>%
		mutate(phase = factor(x%%3))%>%
		# arrange(desc(n))///
		ggplot(aes(x=x,ymax=n,ymin=0,color=phase))+geom_linerange()+
		facet_grid(scale='free',length~.)+
		geom_vline(linetype=2,xintercept=c(0,epoint+2),alpha=I(0.3))+
		scale_x_continuous(labels=metaplotlabs,breaks=metaplotbreaks)+
		# scale_y_continuous(limits=c(0,300))+
		# coord_cartesian(xlim=c(-lflank,rflank))+
		geom_rect(aes(xmin=0,xmax=epoint,ymin=0,ymax=Inf,alpha=I(0.01)),fill='grey',color='white')+
		scale_fill_discrete(guide=F)+
		theme_minimal()	
	p
}))

maxn<-plist%>%map_dbl(~max(.$data$n))%>%max
# for(i in seq_along(plist)) plist[[i]] <- plist[[i]]+scale_y_continuous(limits=c(0,maxn/2))
for(i in seq_along(plist)) plist[[i]] <- plist[[i]]+scale_y_log10(limits=c(100,10e3))

#+ riboshift_comparison plot, fig.width =12,fig.height=12,out.width=1200,out.height=1250,dev='pdf',include=TRUE,eval=TRUE

if(isknitr) {
	ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs))
}else{
	suppressWarnings(dir.create(outfolder,rec=TRUE))
	shiftplotfile <- 'shift_methods_comp.pdf'
	shiftplotfile <- file.path(outfolder,shiftplotfile)
	pdf(shiftplotfile,w=14,h=10)
	print(ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs)))
	dev.off()
	message(normalizePath(shiftplotfile))
}


}


stop()




#







# data4readshift%>%as.data.frame%>%select(matches('^[ft]p'))%>%
# 	apply(2,table)%>%.['1',TRUE]%>%
# 	enframe%>%separate(name,into=c('end','pos','base'),sep='\\.')%>%
# 	group_by(end,pos)%>%mutate(value = value/sum(value))%>%
# 	spread(base,value)





#Now let's check up/down stream of codons





shiftplotfile <- 'shift_methods_comp.pdf'
shiftplotfile <- file.path(outfolder,shiftplotfile)
pdf(shiftplotfile,w=14,h=10)
print(ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs)))
dev.off()
message(normalizePath(shiftplotfile))


library(Biostrings)

cdsseq <- topExonseq[topcdsmap]

for(codon in names(GENETIC_CODE)){
	codonlocs <- vmatchPattern(codon,cdsseq)
	codonmatchgr<-unlist(codonlocs)%>%{GRanges(names(.),.)}%>%setNames(NULL)

	tmp1<-codonmatchgr%>%resize(30,'end')%>%mapFromTranscripts(topcdsmap)%>%head
	seqnames(tmp1)%in%names(topExonseq)
	seqinfo(tmp1)
	topExonseq[tmp1[,NULL]]

	# vmatchPDict(names(GENETIC_CODE)[1:2],cdsseq[1:10])
}

cdsread_trmap%>%.$length%>%table

makeTxDbFromGRanges(gtf_gr)

	
forget(readGAlignments)
bam%>%file.exists
bamfile<-BamFile(bam,yield=)



##Let's make fake, positive control reads

#for each of our 

allstoppsites%>%head%>%as.data.frame



#Simulate Riboseq reads
#length distribution
lengthdist <- c(`26` = 51599L, `27` = 82461L, `28` = 57037L, `29` = 43928L)%>%{./sum(.)}
#initial_frame_dist 
initial_frame_dist <- c(0.6,0.3,0.1)
testexons <- topcdsexons[1:100]
testcds <- topcdsmap[topcdsexons$transcript_id%>%unique]

sample<-base::sample

starts_ends <- testcds%>%head(2)%>%width%>%divide_by(3)%>%map(seq_len)%>%map(sample,1e3,replace=TRUE)%>%
	map(~ .+ sample(c(0,1,2),p=initial_frame_dist,rep=TRUE,size=length(.)))%>%
	map(~data_frame(start=.,end= . + as.numeric(names(lengthdist))[sample(seq_along(lengthdist),rep=TRUE,size=length(.))]))%>%


####Look more modelishly at reads

################################################################################
######### Distribution around codons
################################################################################
	

prestops <- topcdsmap%>%resize(3,'end')%>%resize(1,'start')
stopreads <- subsetByOverlaps(cdsread_trmap,prestops)
leftseg <- fpend(stopreads)%>%downstream_dist_till(prestops)
rightseg <- tpend(stopreads)%>%upstream_dist_till(prestops)

table(stopreads$length)

#so yeah very much left seg influence by right....
table(leftseg,rightseg)%>%apply(1,function(x) x / sum(x))
table(leftseg,rightseg)%>%apply(1,function(x) x / sum(x))%>%.[c(7,8,9),]


codons <- names(GENETIC_CODE)

codondisttabs<-lapply(codons,function(codon){

	codonlocs <- topExonseq[topcdsmap%>%setpos]%>%
		# vmatchPattern(codon,.)%>%
		vmatchPattern(codon,.)%>%
		unlist%>%
		subset((end%%3)==0)%>%
		GRanges(names(.),.)

	codonlocs%<>%mapFromTranscripts(topcdsmap)

	

	lendisttab<-codonlocs%>%mergeByOverlaps(cdsread_trmap)%>%{tibble(dist=start(.$.)-start(.$cdsread_trmap),length=.$length)}

	lendisttab
})

codondisttabs%<>%setNames(codons)
codondisttabs%<>%bind_rows(.id='codon')

codistplot <- here('plots/codondistplots/codondistbycodon.pdf')
codistplot%>%dirname%>%dir.create(showWarnings=FALSE)
pdf(codistplot)
codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
	filter(dist>5)%>%filter(dist<length-5)%>%
	ungroup%>%
	mutate(dist=floor(dist/3))%>%
	group_by(codon,length,dist)%>%summarise(n=sum(n))%>%
	group_by(codon,length)%>%
	mutate(n=n/min(n))%>%
	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
	scale_x_continuous(breaks=1:20)+
	scale_color_discrete(guide=FALSE)+
	theme_bw()
dev.off()
normalizePath(codistplot)%>%message

codistplot <- here('plots/codondistplots/codondist.pdf')
pdf(codistplot)
codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
	filter(dist>5)%>%filter(dist<length-5)%>%
	group_by(codon,length)%>%
	mutate(n=n/min(n))%>%
	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
	scale_x_continuous(breaks=seq(1,30,3),minor_breaks=1:30,name="Codon bp1 - 5' Read End")+
		scale_color_discrete(guide=FALSE)+
	theme_minimal()+
	theme(panel.grid.minor = element_line( size=0.5))
dev.off()
normalizePath(codistplot)%>%message



codonsmagvect<-codondisttabs%>%group_by(codon)%>%summarise(signal=sum(between(dist,6,16)))%>%arrange(desc(signal))%>%.$codon


stopcodons<-codons%>%DNAStringSet%>%translate%>%as.character%>%is_in(c('*'))%>%codons[.]
stopcodonsoratg<-c('ATG',stopcodons)

codistplot <- here('plots/codondistplots/codondist.pdf')
pdf(codistplot)
codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
	filter(codon %in% c(codonsmagvect[1:4],rev(codonsmagvect)[1:4]))%>%
	filter(!codon %in% stopcodonsoratg)%>%
	filter(dist>3)%>%filter(dist<length-3)%>%
	group_by(codon,length)%>%
	mutate(n=n/min(n))%>%
	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
	scale_x_continuous(breaks=seq(1,30,3),minor_breaks=1:30,name="Codon bp1 - 5' Read End")+
		scale_color_discrete()+
	theme_minimal()+
	theme(panel.grid.minor = element_line( size=0.5))
dev.off()
normalizePath(codistplot)%>%message


codondisttabs[[9]]%>%.$dist%>%table%>%
	.[order(as.numeric(names(.)))]%>%txtplot(y=.,x=as.numeric(names(.)))


codondisttabs%>%setNames(codons)%>%bind_rows(.id='codon')
	filter(dist==5,length==29)%>%group_by(codon,dist)%>%tally

codonlocs%>%mergeByOverlaps(cdsread_trmap)%>%{tibble(dist=start(.$.)-start(.$cdsread_trmap),length=.$length)}%>%table%>%
	.[order(as.numeric(names(.)))]%>%txtplot(y=.,x=as.numeric(names(.)))




# TODO - test this
# ' Why am I not getting great periodicity?
# ' What am I not getting 

# length
# total psite count
# psite count at the edges






# TODO - count just outside the CDS here - both to get our edge scores and to allow access to seqs

# and in the middle

# Now derive various states for our ORFs.

