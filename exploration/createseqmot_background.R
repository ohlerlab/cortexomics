#This script will calculate a set of regions, say all utrs, along with some data about those, and an annotation of which are positive and negate. It will then use the data as well as some caclulated statistics like length and gc content to try and construct a realistic background set for the positive sequences
#can probably jsut generate the sequences directly actually.
message('loading libraries')
library(magrittr)
library(GenomicRanges)
library(data.table)
library(rtracklayer)
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
message('...done')
filter<-dplyr::filter
slice<-dplyr::slice
summarize<-dplyr::summarize

#load arguments
getwd()
sargs <- c(
  ribodifffolder = './ribodiff',
  normcountstable='exprdata/allcounts_snorm.tsv',
  annotation='my_gencode.vM12.annotation.gff3',
  genome='my_GRCm38.p5.genome.chr_scaff.fa',
  outputfolder = 'motif_seqs_backgrounds'
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
transregged <- ribodiffcontrastobs%>%bind_rows(.id='time')%>%group_by(feature_id)%>%filter(sum(adj_p_value<0.05)>1)   
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

allcds <- annotation%>%subset(type=='CDS') 
anno <- annotation%>%subset(type=='three_prime_UTR')%>%head(1000)
subanno <- annotation%>%subset(type=='CDS')%>%head(1000)

subtract_gr<- function(anno,subanno){
	disjoin <- c(anno[,NULL],subanno[,NULL])%>%disjoin(with=TRUE)
	disjoin_sub <- disjoin[max(disjoin$revmap) < length(anno)]
	srevmap <- disjoin_sub$revmap%>%setNames(.,seq_along(.))%>%stack
	outranges <- disjoin_sub[as.numeric(srevmap$name),]
	mcols(outranges) <- mcols(anno)[srevmap$value,]
	outranges
}

tputrs <- import('tputrs.gtf')

genesetname <- 'tranl_reg'
testset <- translgeneids


getmatchingset<-function(testset,data){
	data$testvar <- data[[1]] %in% testset
	stopifnot(sum(data[['testvar']]) > 10)
	matchformula <- as.formula(paste0('testvar ~ ',colnames(data)%>%str_subset('total')%>%paste0(collapse='+')) )
	matchobject <- MatchIt::matchit(matchformula,data=data)
	data[[1]][matchobject[[1]]%>%as.numeric]
}

backgroundset <- getmatchingset(
	testset,
	normcount%>%select(gene_id,matches('total'))
)

regionname <- 'CDS'
outseqfile <- file.path(outputfolder,paste0(regionname,'_',genesetname,'.pdf'))
annotation %>% subset(type==regionname) %>% subset(gene_id %in% testset) %>% reduce %>% getSeq(genome,.) %>% Biostrings::writeXStringSet(outseqfile)

regionname <- 'five_prime_UTR'
outseqfile <- file.path(outputfolder,paste0(regionname,'_',genesetname,'.pdf'))
annotation %>% subset(type==regionname) %>% reduce %>% setdiff(allcds,ignore.strand=T)%>%getSeq(genome,.) %>% Biostrings::writeXStringSet(outseqfile)

regionname <- 'three_prime_UTR'
outseqfile <- file.path(outputfolder,paste0(regionname,'_',genesetname,'.pdf'))
annotation %>% subset(type==regionname) %>% reduce %>% setdiff(allcds,ignore.strand=T)%>%getSeq(genome,.) %>% Biostrings::writeXStringSet(outseqfile)


#methods motif searches

# When searching for motifs, we wished to isolate motifs responsible for translational reuglation
# Since translation and transcriptional regulation often act on the same genes, we constructed
# a background set of genes with similiar expression patterns to our test sets, by using a greeding matching algorithm.
# which considered similiarity in all RNAseq libraries (i.e. not conting mass spec or riboseq data). 
# We then output sequences for 3' UTRs, 5' UTRs, and CDS, with UTR sequences considered only when
# not contained in the coding sequence of any annotated CDS. 


outseqfile <- file.path(outputfolder,paste0(regionname,'_',genesetname,'.pdf'))



#now output each of the sequence types



reggenes <- unique(transregged$feature_id)
controlgenes <- normcount$gene_id[matchobject[[1]]%>%as.numeric]

set <- 
background <- 

#ultimately we want lists of test and non test gene


















