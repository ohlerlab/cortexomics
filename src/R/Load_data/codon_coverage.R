conflict_prefer("which", "Matrix")
conflict_prefer("rowSums", "BiocGenerics")
conflict_prefer("colMeans", "BiocGenerics")


ms_id2uprotein_id
best_uprotein_ids

library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
require(txtplot)
#coverageplots

CODOONSOFINTEREST<-DNAStringSet(c('TCA','TCG','TCC','TCT'))
codons2scan <- c(CODOONSOFINTEREST,reverseComplement(CODOONSOFINTEREST))
# codons2scan <- c(CODOONSOFINTEREST,(CODOONSOFINTEREST))
ctrlcodons<-c('TTT','GGG','AAA','CCC')
codons2scan <- c(CODOONSOFINTEREST,ctrlcodons)
#codons2scan <- CODOONSOFINTEREST

FLANKCODS<-15


cdsids2use <- ms_id2protein_id%>%filter(uprotein_id %in% best_uprotein_ids)%>%.$protein_id

cds2use <- cds%>%subset(protein_id %in% cdsids2use)%>%split(.,.$protein_id)

REF <- 'pipeline/my_GRCm38.p5.genome.chr_scaff.fa'
ref<-Rsamtools::FaFile(REF)

cdsseq <- cds2use%>%extractTranscriptSeqs(x=ref)

chrmcdsseqs<-getSeq(x=ref,chrmcds)

allcodons=getGeneticCode()
#get code for mitcondrial coding genes

bamfiles<-Sys.glob(here::here('pipeline/star/data/*/*_*bam'))

bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%
	str_subset(neg=TRUE,'transcript')%>%
	str_subset('_ribo_')

bamfile=bams[1]


cdscodons <- 	lapply(seq_along(allcodons),function(i){
	
	codon=names(allcodons)[[i]]
	message(codon)
	codmatches<-vmatchPattern(pattern=codon,cdsseq)

	nmatches = 	codmatches%>%elementNROWS 

	strands<-strand(cds2use[names(codmatches)])[as(rep(1,length(cds2use)),'IntegerList')]%>%unlist
	strands <- rep(as.vector(strands),nmatches)
	matchgr<-codmatches%>%unlist%>%GRanges(names(.),.,strand=strands)

	matchgr%<>%subset(start %%3 == 1)

	codmatchonchr<-matchgr%>%mapFromTranscripts(cds2use)
	codmatchonchr%<>%subset(width==3)

	codmatchwindows<-codmatchonchr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')

	codmatchwindows

})
cdscodons%<>%setNames(names(allcodons))

codons2scan<-allcodons%>%.[str_detect(names(.),'TTC|GTC|CAC|AAC|ATG')]

#check this worked
stopifnot(all( c('A','T','G') == (cdscodons[['ATG']]%>%getSeq(x=ref)%>%as.matrix%>%apply(2,table)%>%simplify2array%>%head((1+FLANKCODS)*3)%>%tail(3)%>%unlist)%>%names))

bams2scan<-bams%>%.[c(1:2,length(.)-1,length(.))]
psitelist <- mclapply(bams2scan,get_genomic_psites,cds2use%>%unlist)

#Now run through the bam file and get a profile for each codon
profilemats_unproc<-lapply(psitelist,function(psites){
	
	# psites = psitelist[[3]]
	message(bamfile)

	cdscounts_densitys <- countOverlaps(cds2use,psites)
	stopifnot(	all(names(cdscounts_densitys)==names(cds2use)) )
	cdscounts_densitys <- cdscounts_densitys%>%divide_by(sum(width(cds2use)))%>%setNames(names(cds2use))
	
	
	mclapply(mc.cores=10,seq_along(codons2scan),function(i){

		codon=names(codons2scan)[[i]]
		
		# codon='ATG' codon='GTC'

		message(codon)

		codmatchwindows <- cdscodons[[codon]]

		lapply(25:31,function(length_i){
			pospsites<-psites%>%subset(length==length_i)%>%subset(strand=='+')
			negpsites<-psites%>%subset(length==length_i)%>%subset(strand=='-')%>%head(1)

			poswinds<-codmatchwindows%>%subset(strand=='+')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			negwinds<-codmatchwindows%>%subset(strand=='-')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')

			pospsites%>%coverage%>%.[poswinds]%>%length%>%divide_by(1e3)
			negpsites%>%coverage%>%.[negwinds]%>%length%>%divide_by(1e3)

			profiles <- 
			c(
				pospsites%>%coverage%>%.[poswinds],
				negpsites%>%coverage%>%.[negwinds],
				NULL
			)

			nposwinds<-poswinds%>%length
			profiles[-(1:nposwinds)]%<>%revElements

			#get the profiles as a matrix, of nbp cols, and nrow codons
			# profilemat<-profiles%>%setNames(paste0('w',seq_along(.)))%>%as('List')%>%as("DataFrame")%>%as("Matrix")

			# profilemat<-profiles%>%simplify2array%>%t
			windnames<-c(names(poswinds),names(negwinds))
			# profilemat%<>%sweep(1,STATS=cdscounts_densitys[windnames],FUN='/')
			# profilemat%<>%as.matrix
			# profilemat%<>%t
			# profilemat <- profilemat[	!profilemat[,1]%>%map_lgl(is.nan),]
			# profilemat <- profilemat[	profilemat%>%apply(1,function(x) all((is.finite(x))) ) ,]
			np<-length(profiles)
			n_col<-length(profiles[[1]])

			# sapply(i:n_col,function(i) profiles[as(rep(i,np),'IntegerList')]%>%sum%>%add(1)%>%log)
			#now rather than covert to a matrix, we can do this to 
			logsumpvect<-sapply(1:n_col,function(i){
				profiles[as(rep(i,np),'IntegerList')]%>%sum%>%add(1)%>%log%>%sum
			})%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)

			# rustsumvect<-sapply(1:n_col,function(i){
			# 	profiles[as(rep(i,np),'IntegerList')]%>%sum%>%{.>0}%>%sum
			# })%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)

			# sumvect<-sapply(1:n_col,function(i){
			# 	profiles[as(rep(i,np),'IntegerList')]%>%sum%>%sum
			# })%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)


			# invisible(list())
			invisible(list(logsumpvect))
		})
	})
})

profilemats_unproc[[1]][[1]][[1]][[1]]

profilemats_unproc%>%saveRDS('data/codon_profilemats_unproc.rds')

profilemats <- profilemats_unproc%>%setNames(.,bams[1:length(.)])

# stop()

# codonstbl<-read.table('https://raw.githubusercontent.com/zhanxw/anno/master/codon.txt')
# codonstbl%>%write_tsv('ext_data/codons.txt')
codonstbl<-read_tsv('ext_data/codons.txt')
#segregate the 
codonnames <- paste0(codonstbl[[2]],'-',codonstbl[[1]])[match(names(codons2scan),codonstbl[[1]])]

#codoncovdf<-
profilemats%<>%map_df(.id='bam',.%>%map_df(.id='codon',.%>%enframe('position','signal')))

bamsamples<-bams%>%basename%>%str_remove('.bam')

sampleparams<-'src/sample_parameter.csv'%>%fread

profilemats%>%head
profilemats$sample = profilemats$bam%>%match(bams)%>%bamsamples[.]
profilemats$time<-sampleparams$time[match(profilemats$sample,sampleparams$sample_id)]



profilemats$codon%<>%as.numeric
profilemats$codon<-profilemats$codon%>%codonnames[.]

#summarise to get rid of the periodicity signal
codonprofmat<-profilemats%>%mutate(codonpos=ceiling(position/3))%>%group_by(codon,position,bam)

codonprofmat$codonpos%<>%subtract(FLANKCODS+1)
codonprofmat%<>%ungroup


require(plotly)

displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(displaystagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols <- displaystagecols%>%setNames(names(stageconv))
samplestage <- unique(codonprofmat$sample)%>%{setNames(stageconv[str_replace(.,'_.*?_.*?$','')],.)}

samplestage[codonprofmat$sample]%>%unique

highcod = 'AAC'
lowcod = 'GTC'

pdf('tmp.pdf')
codonprofmat%>%
	filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
	# filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
	filter(sample%>%str_detect(c('ribo')))%>%
	{print(ggplot(.,aes(codonpos,signal,group=bam,
		color=samplestage[sample]))+
		scale_color_manual(values=displaystagecols)+
		facet_grid(codon~.)+
		geom_line(aes())+
		theme_bw())}
dev.off()
normalizePath('tmp.pdf')


#' Convert from Rle to one column matrix
#'
Q
setAs("Rle", "Matrix", function(from) {
    rv <- runValue(from)
    nz <- rv != 0
    i <- which(as.vector(from !=0))
    x <- rep(rv[nz], runLength(from)[nz])

	length(i)    

    sparseMatrix(i= i, p=c(0L, length(x)), x=x,
                 dims=c(length(from), 1))
})

#' Convert from DataFrame of Rle to sparse Matrix
#'
setAs("DataFrame", "Matrix", function(from) {
  mat = do.call(cbind, lapply(from, as, "Matrix"))
  colnames(mat) <- colnames(from)
  rownames(mat) <- rownames(from)
  mat
})



