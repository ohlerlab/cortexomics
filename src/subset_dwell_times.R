limmadf = readxl::read_xlsx("tables/S2.xlsx",3,col_types=c(time='text'))
txnchangegenes = limmadf%>%filter(assay=='all')%>%filter(adj.P.Val<0.05)%>%.$gene_id%>%unique
nochangegenes = limmadf%>%filter(!gene_id %in%txnchangegenes)%>%.$gene_id%>%unique
higexprgenes = limmadf%>%filter(as.numeric(AveExpr) > quantile(as.numeric(AveExpr),0.5))%>%.$gene_id%>%unique
nochange_highexpr = intersect(nochangegenes,higexprgenes)
gid2trid = ids_nrgname%>%distinct(gene_id,transcript_id)%>%{setNames(.$transcript_id,.$gene_id)}
tr_nochange_highexpr =gid2trid[nochange_highexpr]%>%intersect(ribocovtrs)

dte_df = readxl::read_xlsx("tables/S2.xlsx",1,col_types=c(time='text'))
dte_df%<>%filter(time=='3')
dte_df$time='E175'
teupgenes = dte_df%>%filter(adj_p_value<0.05,log2fc>0)%>%.$gene_id
tedowngenes = dte_df%>%filter(adj_p_value<0.05,log2fc<0)%>%.$gene_id
dtegenes = c(teupgenes,tedowngenes)


gcodwide<-gcodwide%>%left_join(dte_df%>%select(time,g_id=gene_id,log2fc))


#definite chrom and ribo genes.
GTOGO <- 'data/GTOGO.rds'%>%readRDS%>%select(gene_name,go_id,g_id=ensembl_gene_id)
goid_ribo = 'GO:0003735'
goid_chrom = 'GO:0003682'
ribogenes = GTOGO%>%filter(go_id==goid_ribo)%>%.$g_id
chromgenes = GTOGO%>%filter(go_id==goid_chrom)%>%.$g_id
ribchrgenes<-c(ribogenes,chromgenes)

dtselgenelist = list(
	up = ribocovtrs[trid2gid[[ribocovtrs]]%>%is_in(teupgenes)],
	down = ribocovtrs[trid2gid[[ribocovtrs]]%>%is_in(tedowngenes)],
	nochangehighe = tr_nochange_highexpr
	# all = ribocovtrs
)

seltrs=dtselgenelist[[1]]
sampfpcov=fpcovlist[[mainsamps[1]]]
rlfpcov=sampfpcov[[1]]
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

fpcovlist<-readRDS(here('data/fpcovlist.rds'))


fpcovlist <- fpcovlist[names(allbamtbls)]
if(!file.exists(here('data/subfpprofilelist.rds'))){
	subfpprofilelist <-mclapply(mc.cores=8,dtselgenelist,function(seltrs){
			imap(fpcovlist[mainsamps[c(1)]],function(sampfpcov,sampname){
				trsums = sampfpcov%>%map(sum)%>%map(~.[seltrs])%>%purrr::reduce(.,`+`)#sum over counts for that transcript
				sampfpcov['29']%>%lapply(function(rlfpcov){
					rlfpcov = rlfpcov[seltrs]
					rlfpcov = rlfpcov/(trsums)
					# rlfpcov = rlfpcov > mean(rlfpcov)
					allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
					('.')
					out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
					out
				})
			})
		})
	saveRDS(subfpprofilelist,here('data/subfpprofilelist.rds'))
}else{
	subfpprofilelist<-readRDS(here('data/subfpprofilelist.rds'))
}


if(!file.exists(here('data/subfprustprofilelist.rds'))){
	subfprustprofilelist <-
		mclapply(mc.cores=8,dtselgenelist,function(seltrs){
			imap(fpcovlist[mainsamps[1]],function(sampfpcov,sampname){
				trsums = sampfpcov%>%map(~.[trspacecds[seltrs]])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
				sampfpcov%>%lapply(function(rlfpcov){
					rlfpcov = rlfpcov[seltrs]
					rlfpcov = rlfpcov > mean(rlfpcov)
					# rlfpcov = rlfpcov > mean(rlfpcov)
					allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
					('.')
					out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
					out
				})
			})
		})
	saveRDS(subfprustprofilelist,here('data/subfprustprofilelist.rds'))
}else{
	subfprustprofilelist<-readRDS(here('data/subfprustprofilelist.rds'))
}

# rustprofiledat <- subfprustprofilelist%>%map_df(.id='subgroup',.%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon'))))
# codonprofiledat <- subfpprofilelist%>%map_df(.id='subgroup',.%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon'))))
codonprofiledat <- fpprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
codonprofiledat %<>% mutate(position = position - 1 - (FLANKCODS*3))
# codonprofiledat %<>% group_by(subgroup,sample,readlen,codon)%>%mutate(count= count / median(count))
codonprofiledat %<>% group_by(sample,readlen,codon)%>%mutate(count= count / median(count))
codonprofiledat %<>% filter(!codon %in% c('TAG','TAA','TGA'))

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


