#This script will calculate a set of regions, say all utrs, along with some data about those, and an annotation of which are positive and negate. It will then use the data as well as some caclulated statistics like length and gc content to try and construct a realistic background set for the positive sequences
#can probably jsut generate the sequences directly actually.
message('loading libraries')
library(magrittr)
library(data.table)
library(rtracklayer)
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
library(GenomicRanges)

message('...done')

#load arguments
getwd()
args <- c(
  ribodifffolder = './ribodiff',
  normcountstable='exprdata/allcounts_snorm.tsv',
  annotation='my_gencode.vM12.annotation.gff3',
  genome='my_GRCm38.p5.genome.chr_scaff.fa',
  outputfolder = 'motseqs'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])
dir.create(outputfolder,showWarn=T)

annotation <- 'my_gencode.vM12.annotation.gff3'%>%import
genome <- Rsamtools::FaFile('my_GRCm38.p5.genome.chr_scaff.fa')

#read in ribodiff
ribodiffresfiles <- Sys.glob(file.path(ribodifffolder,'/riboseqres_*.txt'))%>%
	{stopifnot(length(.)>0);.}
ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols<-do.call(cols,c(list(col_character()),rep(list(col_double()),6)))
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)
ribodiffcontrastobs <- ribodiffresfiles%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)
#mark translationally regulated genes
ribodiffcontrastobs%<>%bind_rows(.id='time')


#now define the regulated groups to look at
transregged<-ribodiffcontrastobs%>%group_by(feature_id)%>%filter(sum(adj_p_value<0.05)>1)   
translgeneids <- transregged$feature_id%>%unique
downtranslgeneids <- transregged%>%filter(all(log2fc[adj_p_value<0.05]<0))%>%.$feature_id%>%unique
uptranslgeneids <- transregged%>%filter(all(log2fc[adj_p_value<0.05]>0))%>%.$feature_id%>%unique

#this function chooses a background set given
#a) a logical vector denoting the set and background set
#b) a data frame frame containing relevant variables to use when deciding on the background set
ids<-fread('ids.txt')%>%distinct
normcount <- fread(normcountstable)%>%
	left_join(ids)
normcount%<>%filter(!is.na(gene_id))
normcount%>%colnames

genesetname <- 'transl_reg'
testset <- translgeneids


testsets = list(
	translregged = translgeneids,
	downtranslregged = downtranslgeneids,
	uptranslregged = uptranslgeneids,
	randomset = ribodiffcontrastobs$feature_id%>%sample(1e3)
)

getmatchingset<-function(testset,data){
	data$testvar <- data[[1]] %in% testset
	stopifnot(sum(data[['testvar']]) > 10)
	matchformula <- as.formula(paste0('testvar ~ ',colnames(data)%>%str_subset('total')%>%paste0(collapse='+')) )
	matchobject <- MatchIt::matchit(matchformula,data=data)
	data[[1]][matchobject[[1]]%>%as.numeric]
}

reduce<-GenomicRanges::reduce
for(genesetname in names(testsets)){
	backgroundset <- getmatchingset(
		testset,
		data = normcount%>%select(gene_id,matches('total'))
	)

	#cds set for the genes
	regionname <- 'CDS'
	outseqfile <- file.path(outputfolder,regionname,genesetname,paste0(genesetname,'.fa'))
	outseqfile%>%dirname%>%dir.create(rec=TRUE,showWarnings = FALSE)

	#for gene ids set, get teh non redundant sequence of type regionname from anno and pull
	#the sequence from gen. make sure non overlaps sequence of type subtracttype
	writeseq <- function(regionname,set,outseqfile,anno=annotation,gen=genome,subtracttype=NA){
		subseq <- anno%>%subset(type==subtracttype)
		out = anno %>% 
			subset(type==regionname) %>% 
			subset(gene_id %in% set[99:100]) %>% 
			split(.,.$gene_id)%>% 
			{GenomicRanges::reduce(.)} 

		out%>%setdiff(subseq)

		out
			%>%
			lapply(GenomicRanges::setdiff,subseq)%>%
			GRangesList%>%
			unlist%>%
			{gr=.;getSeq(x=gen,gr)%>%setNames(.,names(gr))} %>%
			split(.,names(.))%>%
			lapply(FUN=Reduce,f=c)%>%
			{Biostrings::DNAStringSet(.)}%>%
			Biostrings::writeXStringSet(outseqfile)
	}
	for(regionname in c('CDS','five_prime_UTR','three_prime_UTR')){
		message(paste0('writing seqs for 'regionname' of 'genesetname))
		subtractseq <- ifelse(regionname=='CDS',NA,'CDS')

		outseqfile <- file.path(outputfolder,regionname,genesetname,paste0(genesetname,'.fa'))
		outseqfile%>%dirname%>%dir.create(rec=TRUE,showWarnings = FALSE)
		writeseq(regionname,set=testset,outseqfile,subtracttype=subtractseq) 

		outseqfile <- outseqfile%>%str_replace('.fa$','_background.fa')
		writeseq(regionname,backgroundset,outseqfile,subtracttype=subtractseq) 
	}
		
}
#methods motif searches

# When searching for motifs, we wished to isolate motifs responsible for translational reuglation
# Since translation and transcriptional regulation often act on the same genes, we constructed
# a background set of genes with similiar expression patterns to our test sets, by using a greeding matching algorithm.
# which considered similiarity in all RNAseq libraries (i.e. not conting mass spec or riboseq data). 
# We then output sequences for 3' UTRs, 5' UTRs, and CDS, with UTR sequences considered only when
# not contained in the coding sequence of any annotated CDS. 










