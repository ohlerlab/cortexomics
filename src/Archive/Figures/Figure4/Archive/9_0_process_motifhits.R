################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
library(rtracklayer)
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}
################################################################################
########This Script produces a shortlisted list of motif hits, integrating various information
################################################################################
#motif isn't just a dumb string of Us.
#Motif's protein seems expressed in our data
#Motif turns up in the cERMIT results for all TE
#Motif turns up in the cERMIT results for only TE up
highcountgnms <- readRDS(here('data/highcountgnms.rds'))
genesofinterest=c('Nes','Tle4','Flna','Satb2','Bcl11b')

cermittomtom = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/teranked/teranked.fa_cermit/tomtom_out/tomtom.tsv'%>%read_tsv

motifnm_df = readLines(here('ext_data/motif_databases/CISBP-RNA/Mus_musculus.meme'))%>%str_subset('MOTIF')%>%str_split(' ')%>%map(tail,2)%>%
	simplify2array%>%t%>%set_colnames(c('motif_id','motif_name'))%>%as_tibble

cermitmotifs = cermittomtom%>%rename('motif_id'='Target_ID')%>%left_join(motifnm_df)
cermitmotifs$cluster = cermitmotifs$Query_ID%>%str_extract('(?<=input_cluster)\\d+')
cermitmotifs$motif_name%<>%str_replace('ENSMUSG00000056951','Rbm4')
cermitmotifs$motif_gene = cermitmotifs$motif_name%>%str_extract('(?<=\\().*?(?=\\))')

cermitmotifs %>% saveRDS(here('data/cermitmotifs.rds'))

cermotmotsum = cermitmotifs%>%group_by(cluster)%>%arrange(`p-value`)%>%filter(motif_gene %in% highcountgnms)%>%summarise(motif_genes = paste0(collapse=',',motif_gene))
cermotmotsum %<>% filter(!is.na(cluster))

{

cermitmemefile = Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/teranked/teranked.fa_cermit/*.meme')
cermitmemefile = Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/terev/terev.fa_cermit/*.meme')
fafile = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/terev/terev.fa'
fimodir = fafile%>%dirname%>%paste0('/fimo')

if(!file.exists(paste0(fimodir,'/fimo.gff'))) system(str_interp('fimo --oc ${fimodir} ${cermitmemefile} ${fafile} '))

fimogr = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/teranked/fimo/fimo.gff'%>%import
fimogr = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/terev/fimo/fimo.gff'%>%import
fimogr$Name%>%str_subset('Satb2')
fimogr$cluster = fimogr$ID%>%str_extract('(?<=cluster)\\d+')
fimogr$gene_name = fimogr$Name%>%str_extract('(?<=_).*?(?=\\|)')%>%str_replace('cluster\\d+_','')


fimogr%>%mcols%>%as.data.frame%>%select(cluster,gene_name)%>%filter(cluster==5)%>%.$gene_name%>%n_distinct

fimogr%>%mcols%>%as.data.frame%>%select(cluster,gene_name)%>%group_by(gene_name,cluster)%>%tally%>%
	group_by(cluster)%>%
	mutate(frac_less_than = rank(n,ties='min')%>%`-`(1)%>%{./length(.)})%>%
	filter(gene_name=="Satb2")%>%
	arrange(-frac_less_than)

}

motifinfo = fimogr%>%mcols%>%as.data.frame%>%select(cluster,gene_name)%>%inner_join(cermotmotsum)%>%
	filter(gene_name%in%genesofinterest)%>%
	rename('cermit_motif':=cluster)%>%
	group_by(cermit_motif,motif_genes)

motifinfo%>%
	summarise(genes_of_interest_bound = paste0(collapse=',',gene_name))%>%
	arrange(cermit_motif)%>%
	as.data.frame%>%
	write_tsv(here('tables/cermit_motif_info.tsv'))

# fimogr%>%mcols%>%as.data.frame%>%select(cluster,gene_name)%>%left_join(cermitmotifs%>%select(cluster,motif_name))%>%

#also just check waht binds satb2
# fimo ../ext_data/motif_databases/CISBP-RNA/Mus_musculus.dna_encoded.meme motseqs/tputr/teranked/teranked.fa
cisbptputrfimo <- import('pipeline/fimo_out/fimo.gff')

cisbptputr_df = cisbptputrfimo%>%mcols%>%as.data.frame%>%select(ID,Name)%>%
	mutate(motif_gene = ID%>%str_extract('(?<=\\().*?(?=\\))'))%>%
	mutate(gene_name = ID%>%str_extract('(?<=\\)\\-\\d{1,3}\\-).*?(?=\\|)'))

cisbptputr_df%>%group_by(gene_name,motif_gene)%>%tally%>%
	group_by(motif_gene)%>%
	mutate(frac_less_than = rank(n,ties='min')%>%`-`(1)%>%{./length(.)})%>%
	filter(gene_name=="Satb2")%>%
	arrange(-frac_less_than)

tputrlendf = long_tputr%>%width%>%sum%>%enframe('gene_name','length')

cisbptputr_df%>%group_by(gene_name,motif_gene)%>%tally%>%
	group_by(motif_gene)%>%
	left_join(tputrlendf)%>%
	mutate(dens = n/length)%>%
	mutate(frac_less_than = rank(dens,ties='min')%>%`-`(1)%>%{./length(.)})%>%
	filter(gene_name=="Satb2")%>%
	arrange(-frac_less_than)%>%
	as.data.frame

################################################################################
########Print motif plots
################################################################################
#we need make_trajplots from 6_	

cermitmotifs <- readRDS(here('data/cermitmotifs.rds'))

motiftrajgenes = (cermitmotifs$motif_gene%>%na.omit%>%unique)%>%c('Samd4b')
for(motifgene in motiftrajgenes){
  # if(motifgene=='Samd4') motfgene = 'Samd1'
  here('plots/motiftraj')%>%dir.create(showWarn=F)
  plotfile<- here(str_interp('plots/motiftraj/${motifgene}.pdf'))
  pdf(plotfile,w=3*4,h=4)
  print(ggarrange(plotlist=list(make_trajplot(motifgene)),nrow=1))
  dev.off()
  message(normalizePath(plotfile))
}






# cermitsumm = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs/tputr/teranked/teranked.fa_cermit/summary.txt'%>%readLines
# cermitsumm%<>%tail(-3)
# cermitlineclusters = cermitsumm%>%str_extract('^\\d+(?= {3})')%>%{for(i in 2:length(.)){if(is.na(.[i])).[i]=.[i-1]};.}
# cermitlinepatterns = cermitsumm%>%str_extract('(?<= {3})\\w+(?= {3})')

# cermitclustpatterndf = tibble(pattern=cermitlinepatterns,cluster=cermitlineclusters)%>%filter(cluster<=10)
# cermitclustpatterndf%<>%filter(!is.na(pattern))

# cermitmotifs%>%select(Query_consensus,motif_name)

# pattern = cermitclustpatterndf$pattern[1] 

# allcermitmatches = lapply(cermitclustpatterndf$pattern%>%setNames(cermitclustpatterndf$cluster),function(pattern){
# 	vmatchPattern(pattern,long_tputrseq,fixed=F)%>%
# 		as.data.frame%>%
# 		transmute(gene=names(long_tputrseq)[group],start,pattern=pattern)
# })%>%bind_rows



allcermitmatches%>%filter(gene=='Satb2')%>%left_join(cermitclustpatterndf)

allcermitmatches%>%left_join(cermitclustpatterndf)%>%group_by(gene,cluster)%>%tally%>%filter(cluster==1)%>%arrange(n)

allcermitmatches%>%left_join(cermitclustpatterndf)%>%filter(cluster==1)%>%.$gene%>%n_distinct

allcermitmatches%>%left_join(cermitclustpatterndf)%>%.$gene%>%n_distinct

################################################################################
########readexcel pum1/2
################################################################################
library(readxl)	
read_xlsx('ext_data/Zhang_2017_pum_sites.xlsx',sheet=1)[[6]]%>%.[nchar(.)==8]%>%na.omit%>%paste0(collapse='\n')%>%cat(file='ext_data/Pum1.txt')
read_xlsx('ext_data/Zhang_2017_pum_sites.xlsx',sheet=2)[[6]]%>%.[nchar(.)==8]%>%na.omit%>%paste0(collapse='\n')%>%cat(file='ext_data/Pum1.txt')
scp projects/cortexomics/ext_data/Pum*  max-login2:~/
mkdir Pum
mv Pum*.txt Pum
sites2meme Pum > Pumbetter.meme


################################################################################
########temp code on max while bih down
################################################################################
TGTANATA
library(magrittr)

seqs = readDNAStringSet('motseqs/tputr/teranked/teranked.fa')
upnames =  readDNAStringSet('motseqs/tputr/terankeduponly/terankeduponly.fa')%>%names
nocnames =  readDNAStringSet('motseqs/tputr/terankednodown/terankednodown.fa')%>%names%>%setdiff(upnames)
downnames =  readDNAStringSet('motseqs/tputr/teranked/teranked.fa')%>%names%>%setdiff(c(upnames,nocnames))

vmatchPattern('TGTANATA',fixed=F,seqs[upnames])%>%as.data.frame%>%nrow%>%divide_by(length(upnames))
vmatchPattern('TGTANATA',fixed=F,seqs[nocnames])%>%as.data.frame%>%nrow%>%divide_by(length(nocnames))
vmatchPattern('TGTANATA',fixed=F,seqs[downnames])%>%as.data.frame%>%nrow%>%divide_by(length(downnames))



	