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
#
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

ids<-fread('ids.txt')%>%distinct
normcount <- fread(normcountstable)%>%
	left_join(ids)
normcount%<>%filter(!is.na(gene_id))
normcount%>%colnames


#read in ribodiff
ribodiffresfiles <- Sys.glob(file.path(ribodifffolder,'/riboseqres_*.txt'))%>%
	{stopifnot(length(.)>0);.}
ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols<-do.call(cols,c(list(col_character()),rep(list(col_double()),6)))
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)
ribodiffcontrastobs <- ribodiffresfiles%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)
#mark translationally regulated genes
ribodiffcontrastobs%<>%setNames(ribodiffresfiles%>%str_extract(regex('(?<=_).*(?=.txt)')))
ribodiffcontrastobs%<>%bind_rows(.id='time')


#read in ribodiff
xtailresfiles <- Sys.glob(file.path('xtail','/xtail_*.txt'))%>%
	{stopifnot(length(.)>0);.}
xtailcolsnms=c("feature_id", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", 
"pvalue_v1", "E145_log2TE", "E13_log2TE", "log2FC_TE_v2", "pvalue_v2", 
"log2fc", "p_value", "adj_p_value")

xtailcontrastobs <- xtailresfiles%>%map(read_tsv)
xtailcontrastobs%<>%setNames(xtailresfiles%>%str_extract(regex('(?<=_).*(?=.txt)')))
#mark translationally regulated genes
xtailcontrastobs%<>%bind_rows(.id='time')
#now define the regulated groups to look at

ribodiffcontrastobs%>%group_by(feature_id)%>%mutate(n=sum(adj_p_value<0.05))%>%.$n%>%table

ribodiffcontrastobs%>%filter(time=='P0')%>%filter(adj_p_value<0.05)

ribodiffcontrastobs%>%group_by(feature_id)%>%filter(sum(adj_p_value<0.05)>=1)%>%.$time%>%table

transregged<-ribodiffcontrastobs%>%group_by(feature_id)%>%filter(sum(adj_p_value<0.05)>=1)   
translgeneids <- transregged$feature_id%>%unique
downtranslgeneids <- transregged%>%filter(all(log2fc[adj_p_value<0.05]<0))%>%.$feature_id%>%unique
uptranslgeneids <- transregged%>%filter(all(log2fc[adj_p_value<0.05]>0))%>%.$feature_id%>%unique

xtailtransregged<-xtailcontrastobs%>%group_by(feature_id)%>%filter(sum(adj_p_value<0.05)>=1)   
xtailtranslgeneids <- xtailtransregged$feature_id%>%unique
xtaildowntranslgeneids <- xtailtransregged%>%filter(all(log2fc[adj_p_value<0.05]<0))%>%.$feature_id%>%unique
xtailuptranslgeneids <- xtailtransregged%>%filter(all(log2fc[adj_p_value<0.05]>0))%>%.$feature_id%>%unique

nontransregged <- normcount$gene_id%>%setdiff(translgeneids)%>%setdiff(xtailtranslgeneids)

pdf('../plots/expr_comparison_transl_notransl_hist.pdf'%T>%normalizePath%>%message)
print(
	bind_rows(normcount%>%filter(gene_id %in% translgeneids) %>%select(matches('total'))%>%apply(1,median)%>%data_frame(score=.,set='transl'),
	normcount%>%filter(gene_id %in% nontransregged)%>%select(matches('total'))%>%apply(1,median)%>%
	data_frame(score=.,set='nontransl')) %>%
	ggplot(aes(x=score))+geom_histogram()+facet_grid(scale='free',set~.))
dev.off()
# message(normalizePath('tmp.pdf'))



genesetname <- 'transl_reg'
testset <- translgeneids


testsets = list(
	# translregged = translgeneids,
	# downtranslregged = downtranslgeneids,
	# uptranslregged = uptranslgeneids,
	xtailtranslregged = translgeneids,
	xtaildowntranslregged = downtranslgeneids,
	xtailuptranslregged = uptranslgeneids,
	randomset = ribodiffcontrastobs$feature_id%>%sample(1e3)
)

#also define a set of genes which are non translationally regulated at all


reduce<-GenomicRanges::reduce


#this function chooses a background set given
#a) a logical vector denoting the set and background set
#b) a data frame frame containing relevant variables to use when deciding on the background set

getmatchingset<-function(testset,data){
	data$testvar <- data[[1]] %in% testset
	stopifnot(sum(data[['testvar']]) > 10)
	matchformula <- as.formula(paste0('testvar ~ ',colnames(data)%>%str_subset('total')%>%paste0(collapse='+')) )
	matchobject <- MatchIt::matchit(matchformula,data=data)
	data[[1]][matchobject[[1]]%>%as.numeric]
}

subtract_gr<- function(anno,subanno){
	disjoin <- c(anno[,NULL],subanno[,NULL])%>%disjoin(with=TRUE)
	disjoin_sub <- disjoin[max(disjoin$revmap) < length(anno)]
	srevmap <- disjoin_sub$revmap%>%setNames(.,seq_along(.))%>%stack
	outranges <- disjoin_sub[as.numeric(srevmap$name),]
	mcols(outranges) <- mcols(anno)[srevmap$value,]
	outranges
}

#for gene ids set, get teh non redundant sequence of type regionname from anno and pull
#the sequence from gen. make sure non overlaps sequence of type subtracttype
writeseq <- function(regionname,set,outseqfile,anno=annotation,gen=genome,subtracttype=NA){
	subseq <- anno%>%subset(type==subtracttype)
	out = anno %>% 
		subset(type==regionname) %>% 
		subset(gene_id %in% set) %>%
		subtract_gr(subseq)%>%
		split(.,.$gene_id)%>% 
		{GenomicRanges::reduce(.)}%>% 
		unlist%>%
		{gr=.;getSeq(x=gen,gr)%>%setNames(.,names(gr))}

	out%>%split(.,names(.))%>%
		lapply(FUN=Reduce,f=c)%>%
		{Biostrings::DNAStringSet(.)}%>%
		Biostrings::writeXStringSet(outseqfile)
}

for(regionname in c('CDS','five_prime_UTR','three_prime_UTR')){
	outseqfile <- file.path(outputfolder,regionname,paste0('nontransregged','.fa'))
	outseqfile%>%dirname%>%dir.create(rec=TRUE,showWarnings = FALSE)
	message('.')
	writeseq(regionname,
		set=nontransregged,
		outseqfile,
		subtracttype = ifelse(regionname=='CDS',NA,'CDS')
	)
}

backgroundsets <- mclapply(names(testsets),function(genesetname){
 getmatchingset(
		testset,
		data = normcount%>%
			filter(gene_id %in% c(testset,nontransregged))%>%
			select(gene_id,matches('total'))
		)
	})%>%setNames(names(testsets))

names(testsets)[[1]]->genesetname
mclapply(names(testsets),function(genesetname){
	testset <- testsets[[genesetname]]
	
	backgroundset <- backgroundsets[[genesetname]]

	#cds set for the genes
	regionname <- 'CDS'
	outseqfile <- file.path(outputfolder,regionname,genesetname,paste0(genesetname,'.fa'))
	outseqfile%>%dirname%>%dir.create(rec=TRUE,showWarnings = FALSE)

	
	for(regionname in c('CDS','five_prime_UTR','three_prime_UTR')){
		message(str_interp('writing seqs for ${regionname} of ${genesetname}'))
		subtractseq <- ifelse(regionname=='CDS',NA,'CDS')

		outseqfile <- file.path(outputfolder,regionname,genesetname,paste0(genesetname,'.fa'))
		outseqfile%>%dirname%>%dir.create(rec=TRUE,showWarnings = FALSE)
		writeseq(regionname,set=testset,outseqfile,subtracttype=subtractseq) 

		outseqfile <- outseqfile%>%str_replace('.fa$','_background.fa')
		writeseq(regionname,backgroundset,outseqfile,subtracttype=subtractseq) 
	}
		
})
#methods motif searches

# When searching for motifs, we wished to isolate motifs responsible for translational reuglation
# Since translation and transcriptional regulation often act on the same genes, we constructed
# a background set of genes with similiar expression patterns to our test sets, by using a greeding matching algorithm.
# which considered similiarity in all RNAseq libraries (i.e. not conting mass spec or riboseq data). 
# We then output sequences for 3' UTRs, 5' UTRs, and CDS, with UTR sequences considered only when
# not contained in the coding sequence of any annotated CDS. 










