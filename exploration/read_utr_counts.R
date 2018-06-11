root <- ifelse(exists('snakemake'),getwd(),normalizePath('~/projects/cortexomics'))


library(magrittr)
library(assertthat)
select=dplyr::select
slice=dplyr::slice


read_tilecount<-function(x,n_tiles=NULL,fun='mean',exprthresh){

	stopifnot(file.exists(x))

	countdf <- 
		x%>%		#this skips the header, and gets the id, (split by _,1st and 3rd places) the strand, thelength, and the count
		# sprintf('cat %s | tail -n +2  | grep -v "\t0$" | tr "_" "\\\t"   | cut -f 1,3,7,9 ',.)%>%
		sprintf('cat %s | head | tail -n +2  | tr "_" "\\\t"   | cut -f 1,3,7,9 ',.)%>%
		fread(skip=1)%>%
		set_colnames(c('transcript','pos','strand','count'))

	if(is.null(n_tiles)) n_tiles = max(countdf$pos)

	transcripts = countdf[[1]]%>%unique
	countmat <- Matrix::Matrix(data=0,sparse=TRUE,ncol=n_tiles,nrow=length(transcripts))
	countmat[]<-0

	#populate the sparse matrix with the values
	countmat[matrix(c(as.integer(factor(countdf[[1]],levels=transcripts)),countdf$pos),ncol=2)]<-countdf$count
	# countmat[matrix(c(as.integer(factor(countdf[[1]])),countdf$pos),ncol=2)]<-countdf$count
	rownames(countmat)<-transcripts
	countmat

}
message('get size factors so we can scale things right')


#get size factors with gene counts
genecounts <-
	Sys.glob(file.path(root,'data/feature_counts/data/*/*genefeature_counts'))%>%
	setNames(.,basename(.))%>%
	map(fread,skip=1,select=c(1,7))%>%
	Reduce(f=partial(left_join,by=c('Geneid')))
sample_size_factors<-DESeq2::estimateSizeFactorsForMatrix(as.matrix(genecounts[,-1]))
sample_size_factors%<>%setNames(sample_size_factors%>%names%>%dirname%>%basename)
sample_size_factorsdf<-sample_size_factors%>%stack%>%
	set_colnames(c('sizefactor','dataset'))%>%
	mutate(dataset=dataset%>%as.character%>%basename%>%str_replace('.bam',''))%>%
	separate(dataset,c('time','assay','replicate'))


tilefiles=Sys.glob(file.path(root,'data/feature_counts/data/*/*_tile*feature_counts'))

regsamples = tilefiles%>%basename%>%str_replace('_tilesfeature_counts','')
dataset_df<-data_frame(dataset=regsamples,timeassreg=dataset%>%str_replace('\\d+\\.',''))

#read in the tile count data
tilemats <-
	tilefiles%>%
	setNames(regsamples)%>%
	map(read_tilecount)

debug(read_tilecount)

sizefactors = tilemats%>%map(~Matrix::rowSums(.))%>%simplify2array%>%{DESeq2::estimateSizeFactorsForMatrix(.)}


#normalize by library size
tilemats%<>%map2(sizefactors[names(tilemats)],`/`)
#sum replicates
# tilemats <- tilemats%>%split(paste0(dataset_df$timeassreg))
# repnums <- tilemats%>%map(length)
# tilemats <- map(tilemats,Reduce,f=`+`)%>%map2(repnums,`/`)


timecoldict = read_csv('./exploration/tables/stages_colors.tsv',col_names=F)%>%mutate(X1=X1%>%str_replace('\\.',''))%>%{setNames(.$X2,.$X1)}

#now, aggregate over all genes
sumtiledata<-tilemats%>%map(Matrix::colMeans)%>%stack%>%set_colnames(c('count','set'))%>%
	group_by(set)%>%mutate(pos=1:n())
sumtiledata<-	sumtiledata%>%separate(set,c('time','assay','replicate','region'),remove=FALSE)
tilenumfp = sumtiledata%>%filter(region=='fputr')%>%.$pos%>%max
tilenumcds = sumtiledata%>%filter(region=='cds')%>%.$pos%>%max
#offset teh psotiions so we can plot together
sumtiledata%<>%mutate(pos=pos+case_when(region=='fputr'~ 0L,region=='cds'~tilenumfp,region=='tputr'~tilenumfp+tilenumcds))

sumtiledata%>%ggplot(aes(x=as.numeric(pos),y=count,color=assay))+facet_grid(.~time)+geom_line(size=2)+theme_bw()+
	scale_x_continuous(breaks=c(tilenumfp,tilenumfp+(tilenumcds)),labels=c('AUG','Stop'))

sumtiledata%>%filter(pos==1,assay=='ribo',time=='E13')

sumtiledata%>%filter(assay=='ribo')%>%ggplot(aes(x=as.numeric(pos),y=as.numeric(count),color=time,linetype=replicate))+geom_line(size=1)+theme_bw()+
	scale_x_continuous(name = 'Position',breaks=c(tilenumfp,tilenumfp+(tilenumcds)),labels=c('AUG','Stop'))+
	scale_y_continuous(name = 'Mean Read Density')+
	scale_color_manual(values=timecoldict)













#read in values for the utrs and cds as a whole
regfiles<-Sys.glob(file.path(root,'data/feature_counts/data/*/*feature_counts'))%>%
	grep(inv=TRUE,val=TRUE,patt='tile|gene')

datasets = regfiles%>%basename%>%str_replace('feature_counts','')
dataset_df<-data_frame(dataset=datasets,timeassreg=dataset%>%str_replace('\\d+\\.',''))%>%
	separate(dataset,c('time','assay','replicate','region'),remove=FALSE)%>%
	mutate(sample = paste0(time,'_',assay,'_',replicate))
dataset_df%<>%left_join(sample_size_factorsdf)

#read in the tile count data
regcounts <-
	regfiles%>%
	setNames(datasets)%>%
	map(fread,skip=1,select=c(1,7))%>%
	map2(names(.),~ set_colnames(.x,c('transcript_id',.y)))%>%
	Reduce(f=partial(left_join,by=c('transcript_id')))

#Now process - we need to exclude five primer utrs that overlap cds,
#and we 
cds <- import(file.path(root,'data/cds.gtf'))
fputrs <- import(file.path(root,'data/fputrs.gtf'))
tputrs <- import(file.path(root,'data/tputrs.gtf'))

#add gene id to our counts
gtotr<-cds%>%mcols%>%as.data.table%>%distinct(gene_id,transcript_id)
regcounts<-regcounts%>%left_join(gtotr)%>%select(gene_id,transcript_id,everything())

#first exclude things that overlap cds,
strandshift <-function(gr,k){GenomicRanges::shift(gr,ifelse(as.logical(strand(gr)=='-'),-k,k))}
fputrs_ovcds <- fputrs%>%subsetByOverlaps(ignore.strand=TRUE,cds)
regcounts <- regcounts%>%filter(!transcript_id%in%fputrs_ovcds$transcript_id)

tputrs_ovcds <- tputrs%>%subsetByOverlaps(ignore.strand=TRUE,cds)
regcounts <- regcounts%>%filter(!transcript_id%in%tputrs_ovcds$transcript_id)

trcdscounts <- regcounts[,dataset_df$dataset] %>% 
	sweep(2, STATS=dataset_df$sizefactor, FUN='/')%>%
	{.[,str_subset(colnames(.),'cds')]} %>%
	apply(1, sum)
	
regcounts$cdscount <- trcdscounts

#now take the cds with the most counts for each of our genes.
regcounts %<>%
	group_by(gene_id)%>%
	slice(which.max(cdscount))

normregcounts <- regcounts
normregcounts %<>% 	select( gene_id, transcript_id, cdscount, everything())
normregcounts[,-c(1:3)] <- sweep(as.matrix(normregcounts[,-c(1:3)]), 2, STATS=dataset_df$sizefactor, FUN='/')

trr_te<-dataset_df %>%
	split(.,list(.$time,.$region,.$replicate)) %>%
	map('dataset') %>%
	map(~ normregcounts[.[1]] / normregcounts[.[2]] )

reghistdata = normregcounts%>%
	select(-gene_id,-transcript_id,-cdscount)	%>%
	apply(2,function(x){hist(log10(x+1),plot=F)})


histdata<-cbind(
	reghistdata%>%map('mids')%>%stack%>%set_colnames(c('mid','dataset')),
	reghistdata%>%map('counts')%>%stack%>%setNames(c('count','dataset'))%>%select(count)
	)%>%separate(dataset,c('time','assay','replicate','region'),remove=FALSE)

histplot=histdata%>%ggplot(aes(fill=region,x=mid,y=count))+facet_grid(time~assay)+stat_identity(geom='bar',bindwidth=1/100,position='dodge',alpha=I(0.5),width=0.2)+theme_bw()+
	scale_x_continuous(name='Log10(count+1)')
histplot
histplot+coord_cartesian(x=c(1,5))


tehistdata<-cbind(
	trr_te %>% map(function(x){hist(log10(x[[1]]), 30, plot=F)})%>%map('mids')%>%stack%>%set_colnames(c('mid','dataset')),
	count = trr_te %>% map(function(x){hist(log10(x[[1]]), 30, plot=F)})%>%map('counts')%>%stack%>%set_colnames(c('count','dataset'))%>%.[,1]
	)%>%separate(dataset,c('time','region','replicate'),remove=FALSE)

tehistdata%>%head

tehistplot=
	tehistdata%>%
	ggplot(aes(fill=region,x=mid,y=count))+
	facet_grid(time~.)+
	stat_identity(geom='bar',position='identity',alpha=I(0.5),width=0.5)+
	theme_bw()+
	scale_x_continuous(name='RiboSeq Reads / RNAseq Reads')


ggsave(plot=tehistplot,file = file.path(root,'plots/tehist.pdf')%T>%message)



regcountsrepmerge<-dataset_df$dataset%>%split(dataset_df$timeassreg)%>%map(~regcounts[,.])%>%map(rowMeans)%>%simplify2array

#oookay starting ot get tired here.
utrcountprops <-
	regcounts%>%as_data_frame%>%
	gather(dataset,count,-Geneid)%>%
	separate(dataset,c('time','assay','replicate','region'),remove=TRUE)%>%
	group_by(time,assay,replicate,Geneid)%>%
	spread(region,count)

utrcountprops<-
	utrcountprops%>%
	group_by(time,assay,replicate)%>%
	mutate(fpprop = fputrs / (cds+fputrs+tputrs))%>%
	mutate(tpprop = tputrs / (cds+fputrs+tputrs))%>%
	group_by(time,assay,replicate)%>%
	summarise(fpprop = median(fpprop,na.rm=T),tpprop = median(tpprop,na.rm=T),
		fputrs = median(fputrs,na.rm=T),
		cds = median(cds,na.rm=T),
		tputrs = median(tputrs,na.rm=T)
		)
	# summarise(fp_mcountratio = mcount[region=='fputrs']/mcount[region=='cds'],
	# 	tp_mcountratio = mcount[region=='tputrs']/mcount[region=='cds'],
	# 	dataset=unique(dataset))
	# identity

utrcountprops<- libsizes%>%stack%>%set_colnames(c('libsize','ind'))%>%separate(ind,c('time','assay','replicate'))%>%left_join(utrcountprops)


ggplot(utrcountprops,aes(x=log10(libsize),y=fpprop,color=assay))+
	# geom_point()+
	geom_text(aes(label=paste0(time,assay,replicate,sep='_')))+
	scale_x_continuous(limits=c(7,8))+
	scale_y_continuous("Median Proportion of Reads in 5' UTR")+
	theme_bw()+
	ggtitle("Library size vs. Proportion of reads in 5' UTR")%>%
	ggsave(plot=.,file=file.path(root,'plots/libsize_vs_5utr_reads.pdf'))


ggplot(utrcountprops,aes(x=log10(libsize),y=tpprop,color=assay))+
	# geom_point()+
	geom_text(aes(label=paste0(time,assay,replicate,sep='_')))+
	scale_x_continuous(limits=c(7,8))+
	scale_y_continuous("Median Proportion of Reads in 3' UTR")+
	theme_bw()+
	ggtitle("Library size vs. Proportion of reads in 3' UTR")%>%
	ggsave(plot=.,file=file.path(root,'plots/libsize_vs_3utr_reads.pdf'))

cor(libsizes,sample_size_factors$sizefactor)

getOrderedBwSignal <- function(gr,bw,as="NumericList"){

  assert_that(file.exists(bw))
  assert_that(has_extension(bw,'bw'))
  assert_that(is(gr,'GenomicRanges'))
  assert_that(not_empty(gr))

  oldnames = names(gr)
  assert_that(! is.null(oldnames))

  selectionob <- BigWigSelection(gr)
  names_bworder <- map(as.list(selectionob@ranges),names)

  # results <- quietly(import.bw)(con = bw,selection=selectionob,as='NumericList')
  # if(results$warnings=="The universe() getter is deprecated." & is(results$result,'NumericList')) {
  #   results <- results$result
  # }else{
  #   stop('This was deprecated, may need to be changed now')
  # }
  results <- quietly(import.bw)(con = bw,selection=selectionob,as=as)
  if(!is_empty(results$warnings)) stop(results$warnings)
  results <-results$result

  #the results retain the order of the bw selection with chromosome, but they sort the chromosomes
  names_bworder = names_bworder[unique(names(results))]
  names_bworder <- quietly(flatten_chr)(names_bworder)$result

  #all we don't duplicate the old names
  assert_that(length(names_bworder)==length(oldnames))


  #now reorder
  results <- results[match(oldnames,names_bworder)]

  # if(as=='NumericList') results <- sum(results)

  results

}

is_uval <- function(x, k) x == unique(x)[k]
paste_ <- partial(paste,sep='_')


bigwig=sigwigs[[1]]
codons=startcods
normfactors=gnormfacts
#get average signal over a bunch of codons given a bigwig
getsig<-function(bigwig,codons,normfactors,wspan=FIXEDSIZE){

	bsample = bigwig%>%dirname%>%dirname%>%basename
	#get strand from 
	strsym <- switch(str_extract(bigwig,'pos|neg'),
		'pos'='+','neg'='-')

	codons <- codons[strand(codons) == strsym]

	normfacts = normfactors[,bsample,drop=TRUE][match(codons$transcript_id,normfactors$transcript_id)]

	signal <- 
		codons%>%
		resize(wspan+1+wspan,'center')%>%
		{names(.)<-.$transcript_id;.}%>%
		getOrderedBwSignal(bigwig,as='NumericList')


	signal%<>%
		as.list%>%
		map2(normfacts,.f=`/`)%>%
		simplify2array

	signal <- signal[,!apply(signal, 2 ,function(x) any(is.nan(x)))]
	signal <- signal[,!apply(signal, 2 ,function(x) any(! is.finite(x)))]

	signal%>%
		apply(1,mean)%>%
		identity

}
agg_signal <- sigwigs %>% lapply(getsig,startcods,gnormfacts)

#also make fixed points around the tss and stop
FIXEDSIZE <- 60
FIXEDINT <- 5
# FIXEDSIZE <- FIXEDSIZE + 1 + FIXEDSIZE
FIXEDNUM <- (FIXEDSIZE/FIXEDINT) + 1 + (FIXEDSIZE/FIXEDINT)
offsets = seq(-FIXEDSIZE,FIXEDSIZE,by=FIXEDINT)

bigwigs = Sys.glob(file.path(root,'data/bigwigs/*/*/*.chr.bw'))
timecoldict = read_csv(file.path(root,'exploration/tables/stages_colors.tsv'),col_names=F)%>%mutate(X1=X1%>%str_replace('\\.',''))%>%{setNames(.$X2,.$X1)}

#
te_results <- file.path(root,'exploration/tables/riboseqres_P0.txt')
ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols=cols(
  col_character(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double()
)
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)
ribodiffcontrastobs <- te_results%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)
ribodiffgenes <- ribodiffcontrastobs[[1]]%>%select(gene_id=feature_id,everything())%>%filter(adj_p_value < 0.05)




#counts for cds
gnormfacts <- normregcounts%>%filter(gene_id %in%ribodiffgenes)%>%ungroup%>%select(transcript_id, matches('.cds'))%>%{colnames(.)%<>%str_replace('.cds','');.}
gnormfacts <- gnormfacts[apply(gnormfacts,1,median)>2,]
#widths
twidths <- cds%>%GR2DT%>%group_by(transcript_id) %>% summarise(width=sum(width))
#now get density rather than count
gnormfacts[,-1] <-sweep(gnormfacts[,-1] ,1 ,STATS=twidths$width[match(gnormfacts$transcript_id,twidths$transcript_id)],FUN = `/`)

#define codons to look at
startcods<-trstotile %>% 
  filter(transcript_id %in% gnormfacts$transcript_id)%>%
  filter(type=='start_codon')%>%
  DT2GR%>%
  {keepSeqlevels(.,str_subset(.@seqinfo@seqnames,'^chr'),'coarse')}

stopcods<-trstotile %>% 
  filter(transcript_id %in% gnormfacts$transcript_id)%>%
  filter(type=='stop_codon')%>%
  DT2GR%>%
  {keepSeqlevels(.,str_subset(.@seqinfo@seqnames,'^chr'),'coarse')}


#get teh bigwig signal over these
sigwigs <- bigwigs%>% str_subset('pos')
agg_signal <- sigwigs %>% mclapply(mc.cores=20,getsig,startcods,gnormfacts)
tp_agg_signal <-sigwigs %>% mclapply(mc.cores=20,getsig,stopcods,gnormfacts)

agg_signal %<>% setNames(bigwigs%>% str_subset('pos'))  
tp_agg_signal %<>% setNames(bigwigs%>% str_subset('pos'))  

averagogram_data <- list(start=agg_signal,stop=tp_agg_signal)%>%
	modify_depth(2, ~data_frame(density=.)%>%mutate(pos=seq_len(n())-FIXEDSIZE-1))%>%
	map(bind_rows,.id='file')%>%
	bind_rows(.id = 'region')%>%
	mutate(file=basename(file)%>%str_replace('.chr.bw',''))%>%
	separate(file,into=c('time','assay','replicate','strand'),remove=FALSE)%>%
	mutate(sample=paste_(time,assay,replicate))

averagogram_data%<>%
	left_join(sample_size_factorsdf)%>%
	mutate(density = density / sizefactor)


averagogram <- function(x) x %>%	ggplot(aes(x=as.numeric(pos),y=density,color=time,linetype=replicate))+geom_line(size=1)+theme_bw()+
	scale_y_continuous(name = 'Average Relative Density (AU)')+
	scale_color_manual(values=timecoldict)

fp_ribo_averagogram_plot <- 
	averagogram_data%>%filter(region=='start',assay=='ribo')%>%
	averagogram

fp_ribo_averagogram_plot<- fp_ribo_averagogram_plot + 
	scale_x_continuous(name = 'Position Relative to AUG',breaks=offsets[seq(1,length(offsets),len=7)]) +
	ggtitle('RiboSeq Read Density Around Start Codon') 

tp_ribo_averagogram_plot <- 
	averagogram_data%>%filter(region=='stop',assay=='ribo')%>%averagogram

tp_ribo_averagogram_plot <-	
	tp_ribo_averagogram_plot+
	scale_x_continuous(name = 'Position Relative to Stop Codon',breaks=offsets[seq(1,length(offsets),len=7)])+
	ggtitle('RiboSeq Read Density Around Stop Codon')

	
ggsave(fp_ribo_averagogram_plot,file=file.path(root,'plots/fp_ribo_fixed_averageogram.pdf')%T>%message)
ggsave(tp_ribo_averagogram_plot,file=file.path(root,'plots/tp_ribo_fixed_averageogram.pdf')%T>%message)

# plot(1)
# (1:3)%>%debugmap(failif3)

options(error=recover)
options(error=partial(traceback,2))


aug
#testing
vals = augtiles%>%sample(10)%>%rev
assert_that(all(map_dbl(as.list(split(vals,seq_along(vals))),getOrderedBwSignal,bw)==getOrderedBwSignal(vals,bw)))






strandshift <-function(gr,k){GenomicRanges::shift(gr,ifelse(as.logical(strand(gr)=='-'),-k,k))}

augtiles<-startcods%>%
  resize(1,'center')%>%
  rep(each=FIXEDNUM)%>%
  strandshift(offsets)
augtiles<-augtiles[,'transcript_id']
augtiles<-augtiles[str_detect(seqnames(augtiles),'chr')]
names(augtiles)<-paste0(augtiles$transcript_id,'_',offsets)


augsignalmat<-Sys.glob('./data/bigwigs/*/*/*.chr.bw')%>%
  setNames(.,.)%>%
  lapply(getOrderedBwSignal,gr=augtiles)
  mclapply(mc.cores=40,getOrderedBwSignal,gr=augtiles)
augsignalmat<-  simplify2array(augsignalmat)

lapply(seq_along(offsets),function(i){})

colnames(augsignalmat)%<>%basename
rownames(augsignalmat)<-names(augtiles)



#get the data in tall form
augsignaldf<-augsignalmat%>%
  as.data.frame%>%rownames_to_column('tileid')%>%
  gather(dataset,density,-tileid)

#parse file name including strand
augsignaldf=augsignaldf%>%
  mutate(dataset=str_replace(dataset,'.chr.bw',''))%>%
  separate(dataset,into=c('time','assay','replicate','strand'))

#now get only correct strand data
augsignaldf <-
	augsignaldf%>% 
	mutate(strand=ifelse(strand=='neg','-','+'))%>%
	semi_join(data_frame(tileid=names(augtiles),strand=as.vector(strand(augtiles))),by=c('strand','tileid'))

#now we can aggregate across genes for each pos
augsignaldfagg<-augsignaldf<-
  augsignaldf%>%separate(tileid,into=c('transcript','pos'),sep='_')%>%
  group_by(time,assay,replicate,pos)%>%
  summarise(density=mean(density))
  
#and normalize for each sample
augsignaldfagg<-augsignaldf%>%
  left_join(sample_size_factors)%>%
  mutate(count = density / sizefactor)
