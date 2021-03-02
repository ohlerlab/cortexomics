################################################################################
########
################################################################################
'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure1/3_1_codon_profiles_tr.R'

if(!file.exists(here('data/psitecovlist.rds'))){
	psitecovlist = allbamtbls%>%mclapply(mc.cores=1,function(bamtbl){
		bamtbl%>%
			str_interp('grep -e $')%>%
			fread(select=c(2,6,7))%>%
			set_colnames(c('transcript_id','start','length'))%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%ribocovtrs)%>%
			inner_join(offsets%>%select(offset,length))%>%	
			mutate(start = start+offset)%>%
			{GRanges(.$transcript_id,IRanges(.$start,w=1))}%>%
			{seqlevels(.) = ribocovtrs ;.}%>%
			{seqlengths(.) = trlens[ribocovtrs] ;.}%>%
			coverage
	})
	saveRDS(psitecovlist,here('data/psitecovlist.rds'))
}else{
	psitecovlist<-readRDS(here('data/psitecovlist.rds'))
	stopifnot(names(psitecovlist)==names(allbamtbls))
	stopifnot(all(ribocovtrs%in%names(psitecovlist[[1]])))
}




cdsstarts = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1)%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
cdsends = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1,'end')%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
trcds = GRanges(names(cdsstarts),IRanges(cdsstarts,cdsends))

cdsgrl%>%.[order(names(.))]%>%width%>%sum%>%head
trcds%>%.[order(seqnames(.))]%>%width%>%head


neartrs = psitecovlist[[isample]]%>%mean%>%sort%>%tail(1)%>%names


if(!file.exists('data/codposdf.rds')) {
	library(Rsamtools)
	bestcdsseq = cdsgrl[neartrs]%>%sort_grl_st%>%extractTranscriptSeqs(x=FaFile(REF),.)
	vmatchPattern(bestcdsseq,pattern=DNAString('N'))%>%elementNROWS%>%is_in(0)%>%table

	codposdf = lapply(bestcdsseq[neartrs],function(cdsseq){
		codonmat = codons(cdsseq)%>%{cbind(pos = .@ranges@start,as.data.frame(.))}%>%
			identity
	})
	codposdf%<>%bind_rows(.id='tr_id')
	codposdf%>%saveRDS('data/codposdf.rds')
}else{
	codposdf<-readRDS('data/codposdf.rds')	
	stopifnot(all(neartrs%in%codposdf$tr_id))
}


psitecovlist[[1]][trcds]
names(trcds) <- seqnames(trcds)
	
	
isample='E13_ribo_1'


trid2short = ids_nrgname%>%
	group_by(gene_name)%>%distinct(transcript_id)%>%
	transmute(transcript_id,tr_id_short=paste0(substr(gene_name,1,3),'_',rank(transcript_id)))%>%
	{setNames(.$tr_id_short,.$transcript_id)}


nearfolderribo = here('pipeline','near',isample,'ribo-seq')%T>%
	dir.create(.,rec=T,showW=F)
nearfolderribo%>%dir.create
nearfolder = nearfolderribo%>%dirname
#
cdscov <-  psitecovlist[[isample]][trcds[neartrs]]
codon_df <- codposdf%>%subset(tr_id%in%neartrs)
names(cdscov)%<>%trid2short[.]
codon_df$tr_id%<>%trid2short[.]
# cdscov$tr_id = 
#
nearsumdata =cdscov%>%mean%>%multiply_by(3)

#
makenearsumfile<-function(nearsumdata){
	nearsumfile = paste0(nearfolder,'/','genelist','.dat')
	nearsumdata%>%enframe('tr_id','density')%>%
	write_delim(nearsumfile,col_names=F,delim=' ')
	nearsumfile
}
nearsumfile = makenearsumfile(nearsumdata)
#
makeneardata<-function(codon_df,cdscov){
	codcounts = cdscov%>%lapply(as.vector)%>%map(~matrix(.,nrow=3))%>%lapply(colSums)
	neardata = codon_df%>%select(tr_id,codon=x)
	neardata %<>% arrange(tr_id)
	neardata$count = codcounts[unique(codon_df$tr_id)]%>%unlist
	neardata
}
neardata = makeneardata(codon_df,cdscov)
#
makenearfiles<-function(neardata,nearfolderribo){
		i_tr_id = unique(neardata$tr_id)[1]
		map_chr(unique(neardata$tr_id),function(i_tr_id){
			nearfile = paste0(nearfolderribo,'/A-site_',i_tr_id,'.dat')
			neardata%>%
				filter(tr_id==i_tr_id)%>%
				select(codon,count)%>%
				mutate(codon=paste0(codon,' '))%>%
				mutate(count=paste0(' ',count,' '))%>%
				write_tsv(nearfile,col_names=F)
		nearfile
		})
}
nearfiles = makenearfiles(neardata,nearfolderribo)

#install this was a fucking nightmare
#wget https://github.com/stevengj/nlopt/archive/v2.6.2.tar.gz
#stepping int that
#cmake -DCMAKE_INSTALL_PREFIX=$HOME/work/install ; make; make install
#now install NEAR
#git clone https://github.com/jszavits/NEAR
#gfortran -c mt19937.f90
#gfortran NEAR.f90 -I/fast/users/dharnet_m/work/install/include -L/fast/users/dharnet_m/work/install/lib64 -lnlopt -lm -o NEAR mt19937.o
