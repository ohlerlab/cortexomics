library(magrittr)

getwd()

Sys.glob('pipeline//cutadapt_reads/*/*.fastq.gz.cutadaptstats.txt')

system('cat pipeline//cutadapt_reads/*/*.fastq.gz.cutadaptstats.txt | perl -lane \' ~ \\/.*?\\/(.*?)\\/.*\ ([0-9,]+)$/;print $1," ", $2\' ')

system('cat pipeline//cutadapt_reads/E13_ribo_1/*.fastq.gz.cutadaptstats.txt | head -n 10| perl -lane \' ~ \\/.*?\\/(.*?)\\/.*\ ([0-9,]+)$/;print $1," ", $2\' ')


library(tidyverse)
library(data.table)




totalreads <- Sys.glob(paste0('pipeline//cutadapt_reads/*/*.fastq.gz.cutadaptstats.txt'))%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%head(8)%>%tail(1)%>%str_extract('[0-9,]+$')%>%str_replace_all(',','')%>%as.numeric)%>%{setNames(as.vector(.),names(.))}
countedreads <- paste0(pdir,'../Ribo_Lausanne/pipeline/feature_counts/all_feature_counts')%>%fread%>%as.data.frame%>%.[unique(colnames(.))]%>%summarise_at(vars(-feature_id),list(sum))%>%unlist

(round(countedreads[names(totalreads)]/totalreads,3)*100)%>%enframe('sample','usable_reads')%>%arrange(desc(usable_reads))%>%as.data.frame%>%group_by(sample)%>%filter(sample%>%str_detect('OMM'))%>%filter(usable_reads<80)%>%
	slice(1)

setwd('..')

pdirs = c('Ribo_Lausanne','Ovarian_cancer_ribo','cortexomics','Mitosis_Riboseq')
pdir = pdirs[1]
allreadnumstats <- lapply(pdirs,function(pdir){
	samples <- Sys.glob(paste0(pdir,'/pipeline/star/reports/*/*Log.final.out'))%>%dirname%>%basename
	cutadaptfiles <- Sys.glob(paste0(pdir,'/pipeline/cutadapt_reads/*/*.fastq.gz.cutadaptstats.txt'))
	starrepfiles <- Sys.glob(paste0(pdir,'/pipeline/star/reports/*/*Log.final.out'))
	formappednum<-Sys.glob(paste0(pdir,'/pipeline/star/reports/*/Log.final.out'))%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%str_subset('Number of input reads')%>%str_extract('\\d+$')%>%head(1)%>%as.numeric)
	uniquemappednum<-Sys.glob(paste0(pdir,'/pipeline/star/reports/*/Log.final.out'))%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%str_subset('Uniquely mapped reads')%>%str_extract('\\d+$')%>%head(1)%>%as.numeric)
	mappednum<-Sys.glob(paste0(pdir,'/pipeline/star/reports/*/Log.final.out'))%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%str_subset('apped reads')%>%str_extract('\\d+$')%>%head(1)%>%as.numeric)
	totalreads <- cutadaptfiles%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%head(8)%>%tail(1)%>%str_extract('[0-9,]+$')%>%str_replace_all(',','')%>%as.numeric)%>%{setNames(as.vector(.),names(.))}
 	totalreads <- totalreads[samples]%>%setNames(samples)
	countedreads <- paste0(pdir,'/pipeline/feature_counts/all_feature_counts')%>%fread%>%as.data.frame%>%.[unique(colnames(.))]%>%summarise_at(vars(-feature_id),list(sum))%>%unlist
# mapped_reads <- Sys.glob(paste0(pdir,'/pipeline/star/reports/*/*bamstats.txt'))%>%setNames(.,basename(dirname(.)))%>%map_dbl(.%>%readLines%>%str_subset('\treads mapped:')%>%str_extract('\\d+'))
	data.frame(sample=samples,
		totalreadsprefilt = totalreads,
		formap_reads = formappednum,
		processed_read_pc = (round(formappednum[samples]/totalreads[samples],3)*100),
		assignable_reads = countedreads[samples],
		unique_mapped = uniquemappednum
	)%>%mutate(pc_maps_unique = round(unique_mapped/formap_reads,3)*100)
})


allreadnumstats%>%setNames(pdirs)%>%bind_rows(.id='project')%>%group_by(sample)%>%filter(sample%>%
	str_detect('OMM') | (!project %>%str_detect('Laus')))%>%
	mutate(assignable_read_pc=assignable_reads/totalreadsprefilt)%>%
	filter(assignable_read_pc<80 | (is.na(assignable_read_pc)))%>%
	filter(!str_detect(sample,'mappability'))%>%
	filter(!str_detect(sample,'test|AD0|total|80S|Poly|L[57]|chx'))%>%
	filter(!str_detect(project,'Mitosis'))%>%
	ungroup%>%
	# sample_n(10)%>%as.data.frame%>%
	mutate(Million_Input_reads = totalreadsprefilt/1e6)%>%
	mutate(Million_Processed_reads = formap_reads/1e6)%>%
	mutate(Million_Uniquely_Mapped_reads = unique_mapped/1e6)%>%
	mutate(Million_Assigned_Reads = assignable_reads/1e6)%>%
	filter(!is.na(Million_Input_reads))%>%
	identity%>%
	select(project,sample,matches('Million'))%>%write_tsv('project_readuse_stats.tsv')



