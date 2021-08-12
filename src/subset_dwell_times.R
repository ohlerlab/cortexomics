{
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source(here("src/Figures/Figure0/0_load_annotation.R"))
}
if(!exists('fpcovlist')) base::source('src/Figures/Figure2/1_load_pos_data.R')
#
limmadf = readxl::read_xlsx("tables/S2.xlsx",3,col_types=c(time='text'))
txnchangegenes = limmadf%>%filter(assay=='all')%>%filter(adj.P.Val<0.05)%>%.$gene_id%>%unique
nochangegenes = limmadf%>%filter(!gene_id %in%txnchangegenes)%>%.$gene_id%>%unique
higexprgenes = limmadf%>%filter(as.numeric(AveExpr) > quantile(as.numeric(AveExpr),0.5))%>%.$gene_id%>%unique
nochange_highexpr = intersect(nochangegenes,higexprgenes)
gid2trid = ids_nrgname%>%distinct(gene_id,transcript_id)%>%{setNames(.$transcript_id,.$gene_id)}
tr_nochange_highexpr =gid2trid[nochange_highexpr]%>%intersect(ribocovtrs)
#
dte_df = readxl::read_xlsx("tables/S2.xlsx",1,col_types=c(time='text'))
dte_df%<>%filter(time=='3')
dte_df$time='E175'
teupgenes = dte_df%>%filter(adj_p_value<0.05,log2fc>0)%>%.$gene_id
tedowngenes = dte_df%>%filter(adj_p_value<0.05,log2fc<0)%>%.$gene_id
dtegenes = c(teupgenes,tedowngenes)
#
#definite chrom and ribo genes.
GTOGO <- 'data/GTOGO.rds'%>%readRDS%>%select(gene_name,go_id,g_id=ensembl_gene_id)
goid_ribo = 'GO:0003735'
goid_chrom = 'GO:0003682'
ribogenes = GTOGO%>%filter(go_id==goid_ribo)%>%.$g_id
chromgenes = GTOGO%>%filter(go_id==goid_chrom)%>%.$g_id
ribchrgenes<-c(ribogenes,chromgenes)
#
STARTBUFF=60
ENDBUFF=60
ribocovtrs <- readRDS(here('data/ribocovtrs.rds'))
exonseq = exonsgrl[ribocovtrs]%>%extractTranscriptSeqs(x=fafile)
allcodons=getGeneticCode()
shift<-GenomicRanges::shift
i=1
innercds = trspacecds%>%
	subset(width>(3+STARTBUFF+ENDBUFF))%>%
	resize(width(.)-STARTBUFF,'end')%>%
	resize(width(.)-ENDBUFF,'start')
FLANKCODS=15
# telleygenegrps = 
# zappgngrps = 
# dtselgenelist = list(
# 	up = ribocovtrs[trid2gid[[ribocovtrs]]%>%is_in(teupgenes)],
# 	down = ribocovtrs[trid2gid[[ribocovtrs]]%>%is_in(tedowngenes)],
# 	nochangehighe = tr_nochange_highexpr%>%intersect(names(innercds)),
# 	all = names(innercds)
# )
#telley 2019 waves
}

{
gid2trid = ids_nrgname%>%filter(transcript_id%in%ribocovtrs)%>%{setNames(.$transcript_id,.$gene_id)}
gnm2trid = ids_nrgname%>%filter(transcript_id%in%ribocovtrs)%>%{setNames(.$transcript_id,.$gene_name)}
telley_ap_early_genes <- readxl::read_xlsx('ext_data/aav2522_Data-S2.xlsx')%>%select(-matches('E14|E15'))%>%
	filter_at(vars(matches('wave')),~.%in%c(1,2))%>%
	.[['Gene symbol']]%>%gnm2trid[.]%>%.[!is.na(.)]
telley_n4_early_genes <- readxl::read_xlsx('ext_data/aav2522_Data-S2.xlsx')%>%select(-matches('E14|E15'))%>%
	filter_at(vars(matches('wave')),~.%in%c(5,6))%>%
	.[['Gene symbol']]%>%gnm2trid[.]%>%.[!is.na(.)]
#telley 2019 waves
telley_ap_late_genes <- readxl::read_xlsx('ext_data/aav2522_Data-S2.xlsx')%>%select(-matches('E12|E13'))%>%
	filter_at(vars(matches('wave')),~.%in%c(1,2))%>%
	.[['Gene symbol']]%>%gnm2trid[.]%>%.[!is.na(.)]
telley_n4_late_genes <- readxl::read_xlsx('ext_data/aav2522_Data-S2.xlsx')%>%select(-matches('E12|E13'))%>%
	filter_at(vars(matches('wave')),~.%in%c(5,6))%>%
	.[['Gene symbol']]%>%gnm2trid[.]%>%.[!is.na(.)]
neurites <- here('ext_data/neurites_zappulo_etal_2017.csv')
neurites%<>%fread(skip=2)
neurite_trids = neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%
	filter(RiboSeq_log2FC_Neurite_Soma>.32)%>%
	.$gene_id%>%gid2trid[.]%>%.[!is.na(.)]
soma_trids = neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%
	filter(RiboSeq_log2FC_Neurite_Soma<.32)%>%
	.$gene_id%>%gid2trid[.]%>%.[!is.na(.)]
dtselgenelist = list(
	all = names(innercds),
	telley_ap_early_genes=telley_ap_early_genes,
	telley_n4_early_genes=telley_n4_early_genes,
	telley_ap_late_genes=telley_ap_late_genes,
	telley_n4_late_genes=telley_n4_late_genes,
	neurite_trids=neurite_trids,
	soma_trids=soma_trids
)
FLANKCODS=15
}

seltrs=dtselgenelist[[4]]
sampfpcov=fpcovlist[[mainsamps[1]]]
rlfpcov=sampfpcov[['29']]
i=1
#get allcodlist granges object descxribing codon positions in the transcripts
if(!file.exists(here('data/allcodlist.rds'))){
	allcodlist <- lapply(seq_along(allcodons)%>%setNames(names(allcodons)),function(i){
		#	
		codon=names(allcodons)[[i]]
		message(codon)
		codmatches<-vmatchPattern(pattern=codon,exonseq[ribocovtrs])#exclude the start ccodon
		#
		nmatches = 	codmatches%>%elementNROWS 
		#
		matchgr<-codmatches%>%unlist%>%GRanges(names(.),.)
		matchgr$cdspos = start(matchgr) - start(trspacecds[as.vector(seqnames(matchgr))])
		matchgr%<>%subset(cdspos %%3 == 0)
		seqlengths(matchgr) = exonsgrl%>%width%>%.[seqlevels(matchgr)]%>%sum
		innercds = trspacecds%>%subset(width>(3+STARTBUFF+ENDBUFF))%>%
			resize(width(.)-STARTBUFF,'end')%>%
			resize(width(.)-ENDBUFF,'start')
		matchgr = matchgr%>%subsetByOverlaps(innercds)
		codmatchwindows<-matchgr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')
		codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
		codmatchwindows%<>%subsetByOverlaps(innercds)
		codmatchwindows
	})
	allcodlist=allcodlist%>%GRangesList%>%unlist
	saveRDS(allcodlist,here('data/allcodlist.rds'))
}else{
	allcodlist<-readRDS(here('data/allcodlist.rds'))
	stopifnot(allcodlist@seqinfo@seqnames%>%setequal(ribocovtrs))
	stopifnot(allcodlist@seqinfo@seqlengths%>%setequal(sum(width(exonsgrl[ribocovtrs]))))
}

if(!exists('fpcovlist')) fpcovlist<-readRDS(here('data/fpcovlist.rds'))
#	
fpcovlist <- fpcovlist[names(allbamtbls)]
sampname=fpcovlist%>%names%>%head(1)
sampfpcov=fpcovlist[[sampname]]
rlfpcov=sampfpcov[['29']]

if(!file.exists(here('data/subfpprofilelist.rds'))){
	fpcodonmats <- 	imap(fpcovlist%>%head(1),function(sampfpcov,sampname){
				trsums = sampfpcov%>%head(1)%>%map(~.[seltrs])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
				sampfpcov['29']%>%lapply(function(rlfpcov){
					rlfpcov = rlfpcov[seltrs]
					rlfpcov = rlfpcov/(trsums)
					# rlfpcov = rlfpcov > mean(rlfpcov)
					allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums>=0])
					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
					('.')
					out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)
					list(seqnames(allcodlistnz))
					out %>% map(~Matrix::Matrix(.,sparse=TRUE))

				})
			})
	object.size(fpcodonmats)/1e6

	fpcodonmats <- 
			imap(fpcovlist,function(sampfpcov,sampname){
				trsums = sampfpcov%>%head(1)%>%map(~.[seltrs])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
				sampfpcov%>%lapply(function(rlfpcov){
					rlfpcov = rlfpcov[seltrs]
					rlfpcov = rlfpcov/(trsums)
					# rlfpcov = rlfpcov > mean(rlfpcov)
					# allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums>=32])
					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
					('.')
					out = rlfpcov[allcodlistnz]%>%
						as.matrix%>%set_rownames(as.vector(seqnames(allcodlistnz)))%>%
						{mat=.;lapply(unique(cods),function(cod)mat[cods==cod,TRUE])}%>%
						lapply(Matrix::Matrix,sparse=TRUE)
					# out %>% map(colMeans)
					out
				})
			})
	
	object.size(subfpprofilelist)/1e6

	saveRDS(subfpprofilelist,here('data/subfpprofilelist.rds'))
}else{
	subfpprofilelist<-readRDS(here('data/subfpprofilelist.rds'))
}



# #RUST
# if(!file.exists(here('data/subfprustprofilelist.rds'))){
# 	subfprustprofilelist <-
# 		mclapply(mc.cores=1,dtselgenelist,function(seltrs){
# 			imap(fpcovlist[mainsamps[1]],function(sampfpcov,sampname){
# 				trsums = sampfpcov['29']%>%map(~.[innercds[seltrs]])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
# 				sampfpcov%>%lapply(function(rlfpcov){
# 					rlfpcov = rlfpcov[seltrs]
# 					rlfpcov = rlfpcov > mean(rlfpcov)
# 					# rlfpcov = rlfpcov > mean(rlfpcov)
# 					allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
# 					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
# 					('.')
# 					out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
# 					out
# 				})
# 			})
# 		})
# 	saveRDS(subfprustprofilelist,here('data/subfprustprofilelist.rds'))
# }else{
# 	subfprustprofilelist<-readRDS(here('data/subfprustprofilelist.rds'))
# }
# fprustprofilelist<-readRDS(here('data/fprustprofilelist.rds'))

# fprustprofilelist%>%saveRDS('data/fprustprofilelist_good.rds')
# rustprofiledat <- subfprustprofilelist%>%map_df(.id='subgroup',.%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon'))))
# codonprofiledat <- subfpprofilelist%>%map_df(.id='subgroup',.%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon'))))
# fprustprofilelist <- subfpprofilelist[[1]]
codonprofiledat <- fprustprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
# codonprofiledat <- codonprofiles%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
codonprofiledat %<>% mutate(position = position - 1 - (FLANKCODS*3))
# codonprofiledat %<>% group_by(subgroup,sample,readlen,codon)%>%mutate(count= count / median(count))
codonprofiledat %<>% group_by(sample,readlen,codon)%>%
	# mutate(count= count / median(count))
	identity
codonprofiledat %<>% filter(!codon %in% c('TAG','TAA','TGA'))
codonprofiledat$readlen%<>%str_replace('^(\\d)','rl\\1')

# subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-11])%>%enframe('codon','directdt')%>%left_join(codondata%>%filter(time=='E13',rep=='1'))%>%
# 	{quicktest(.$directdt,.$dwell_time)}
subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-11])%>%enframe('codon','directdt')%>%left_join(codondata%>%filter(time=='E13',rep=='1'))%>%
	{quicktest(.$directdt,.$availability)}

# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codondata$dwell_time)
# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codonprofiledat$count)
# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codonprofiledat$count)


rustcodon_dts = rustprofiledat%>%
	mutate(length=as.numeric(readlen))%>%
    safe_left_join(offsets%>%select(length,offset))%>%
    filter(position== -offset-3)%>%
	group_by(subgroup,sample,codon)%>%
    summarise(dwell_time = mean (count))

#now plot
plotfile<- here(paste0('plots/','subgroup_dts_density','.pdf'))
pdf(plotfile,w=10,h=10)
rustcodon_dts%>%
	separate(sample,c('time','assay','rep'))%>%
	ggplot(.,aes(x=dwell_time,fill=subgroup))+
	geom_density(alpha=I(0.5))+
	scale_x_continuous(paste0('Codon Dwel Time'))+
	facet_grid(time~rep)+
	ggtitle(paste0('Codon DT for Up vs Down dTE genes'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


tRNA_occ_df_en


