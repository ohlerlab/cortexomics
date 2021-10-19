{
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source(here("src/Figures/Figure0/0_load_annotation.R"))
}
if(!exists('fpcovlist')) base::source('src/Figures/Figure2/1_load_pos_data.R')
#
limmadf = readxl::read_xlsx("tables/S2.xlsx",4,col_types=c(time='text'))
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
exonseq = exonsgrl[ribocovtrs]%>%extractTranscriptSeqs(x=fafileob)
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
	soma_trids=soma_trids,
	nochange=tr_nochange_highexpr
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

# sampname=names(fpcovlist)[[1]]

# #debug
# testtr <- names(cov[['29']])[1]
# maxval = cov[['29']][[testtr]]%>%max
# #first cod window containing that max value
# testcods <- allcodlist%>%subset(seqnames==testtr)%>%resize(3)%>%shift(45)
# mwindind <- rlfpcov[testcods]%>%`==`(18)%>%any%>%which%>%head(1)

# rlfpcov[testcods[mwindind]]

# #there's the max point
# maxpoint <- rlfpcov[1]%>%GRanges%>%subset(score==18)
# #codon that overlaps it
# testcods[mwindind]

# #here's the max point on our new cov object
# cov[['29']][1]%>%GRanges%>%subset(score==18)



# remapmaxpoint <- cov[['29']][1]%>%GRanges%>%subset(score==18)


# #does this max point on new object overlap the same codon remapped??
# mapped_hitcod <- testcods[mwindind]%>%mapToTranscripts(innercds)
# cov[['29']][1]%>%GRanges%>%subset(score==18)%>%
# 	overlapsAny(mapped_hitcod)

# #remapping of the signal makes sense
# start(remapmaxpoint) + start(innercds[testtr]) -1
# start(maxpoint)

# cov[['29']][testcods[mwindind]%>%mapToTranscripts(innercds%>%setstrand)]


#  %>% is_in(cds_codons)

# cds_codons[mwindind]

# cov[['29']][1]

innercdspos <- innercds%>%{strand(.)<-'+';.}
cds_codons <- allcodlist%>%mapToTranscripts(innercdspos)
#
if(!file.exists(here('data/cdsfpcovlist.rds'))){
	cdsfpcovlist <- mcmapply(mc.cores=8,SIMPLIFY=F,fpcovlist,names(fpcovlist),FUN=function(sampfpcov,sampname){
			cov <- fpcovlist[[sampname]]%>%
				map(~.[names(innercds)])%>%
				map(~.[innercds])
			cov
	})
	saveRDS(cdsfpcovlist,here('data/cdsfpcovlist.rds'))
}else{
	cdsfpcovlist<-readRDS(here('data/cdsfpcovlist.rds'))
}

# cds_codons_nz = cds_codons%>%subset(
mcols(cds_codons)=NULL
cds_codons$codon = names(cds_codons)%>%str_split('\\.')%>%map_chr(1)
sampinds <- c(1:10)



if(!file.exists(here('data/cn_norm.rds'))){
	cn_norm <- 	lapply(dtselgenelist['all'],function(seltrs){
			# mcmapply(mc.cores=4,cdsfpcovlist[mainsamps[1:10]],mainsamps[1:10],F=function(cdssampfpcov,sampname){
			mcmapply(mc.cores=4,SIMPLIFY=F,cdsfpcovlist[mainsamps[sampinds]],mainsamps[sampinds],F=function(cdssampfpcov,sampname){
				cdssampfpcov%>%lapply(function(cdsrlfpcov){
					rustcdsrlfpcov <- cdsrlfpcov
					# rustcdsrlfpcov <- rustcdsrlfpcov>mean(rustcdsrlfpcov)
					nz_trs <- any(rustcdsrlfpcov)%>%names(.)[.]
					#
					cds_codons_nz = cds_codons%>%subset(seqnames%in%nz_trs)
					codtrs = seqnames(cds_codons_nz)%>%as.character
					cat('.')
					#calculate ro vals
					ro_cl = rustcdsrlfpcov[cds_codons_nz]%>%
						split(cds_codons_nz$codon)%>%
						lapply(as.matrix)%>%
						map(colMeans)
					#also get evals
					tr_rust_evals <- rustcdsrlfpcov%>%mean
					re_c <- tr_rust_evals[codtrs]%>%
						split(cds_codons_nz$codon)%>%
						map_dbl(mean)
					ro_cl <- ro_cl%>%map_df(.id='codon', enframe, 'position', 'ro_cl')
					re_c <- enframe(re_c, 'codon', 're_c')
					ro_cl%>%left_join(re_c, by='codon')
				})
			})%>%setNames(mainsamps[sampinds])
		})
	saveRDS(cn_norm,here('data/cn_norm.rds'))
}else{
	cn_norm<-readRDS(here('data/cn_norm.rds'))
}

if(!file.exists(here('data/cdsfpcovlist.rds'))){
	rust_roel <- lapply(dtselgenelist['all'],function(seltrs){
		# mcmapply(mc.cores=4,cdsfpcovlist[mainsamps[1:10]],mainsamps[1:10],F=function(cdssampfpcov,sampname){
		mcmapply(mc.cores=4,SIMPLIFY=F,cdsfpcovlist[mainsamps[sampinds]],mainsamps[sampinds],F=function(cdssampfpcov,sampname){
			cdssampfpcov%>%lapply(function(cdsrlfpcov){
					rustcdsrlfpcov <- cdsrlfpcov
					rustcdsrlfpcov <- rustcdsrlfpcov>mean(rustcdsrlfpcov)
					nz_trs <- any(rustcdsrlfpcov)%>%names(.)[.]
					#
					cds_codons_nz = cds_codons%>%subset(seqnames%in%nz_trs)
					codtrs = seqnames(cds_codons_nz)%>%as.character
					cat('.')
					#calculate ro vals
					ro_cl = rustcdsrlfpcov[cds_codons_nz]%>%
						split(cds_codons_nz$codon)%>%
						lapply(as.matrix)%>%
						map(colMeans)
					#also get evals
					tr_rust_evals <- rustcdsrlfpcov%>%mean
					re_c <- tr_rust_evals[codtrs]%>%
						split(cds_codons_nz$codon)%>%
						map_dbl(mean)
					ro_cl <- ro_cl%>%map_df(.id='codon', enframe, 'position', 'ro_cl')
					re_c <- enframe(re_c, 'codon', 're_c')
					ro_cl%>%left_join(re_c, by='codon')
			})
		})%>%setNames(mainsamps[sampinds])
	})
	saveRDS(rust_roel,here('data/rust_roel.rds'))
}else{
	rust_roel<-readRDS(here('data/rust_roel.rds'))
}

#https://www.nature.com/articles/ncomms12915#Sec10 see equation 3
# of RUST paper
frustprofilelist <- rust_roel%>%
# frustprofilelist <- cn_norm%>%
	map_df(.id='set',.%>%map_df(.id='sample',.%>%bind_rows(.id='readlen')))
#
frustprofilelist %<>% mutate(position = position - 1 - (FLANKCODS*3))
frustprofilelist %<>% filter(!codon %in% c('TAG','TAA','TGA'))
frustprofilelist%<>%mutate(count = ro_cl/re_c)
#
frustprofilelist$sample%>%table
frustprofilelist$readlen%>%n_distinct
# frustprofilelist$sample %<>% {mainsamps[as.numeric(.)]}
# frustprofilelist$readlen %<>% as.numeric %>% names(fpcovlist[[1]])[.]
frustprofilelist$nreadlen <-frustprofilelist$readlen%>%as.numeric
frustprofilelist$readlen%<>%str_replace('^(\\d)','rl\\1')
#
kl_df<-frustprofilelist%>%
	filter(set=='all')%>%
	# group_by(readlen)%>%group_slice(4)%>%
	# group_by(position)%>%
	# group_slice(1)%>%
	# filter(position%%3 ==0)%>%
	# .$position
	# filter(position< -6)%>%
	# filter(position> -(nreadlen-6))%>%
	group_by(set,sample,nreadlen,position)%>%
	# summarise(KL=sum(ro_cl * log2(ro_cl/re_c)))
	mutate(ro_cl = ro_cl/sum(ro_cl), re_c = re_c/sum(re_c))%>%
	summarise(KL=sum(ro_cl * log2(ro_cl/re_c)))
#
kl_df%>%
	mutate(phase=position%%3)%>%
	group_by(nreadlen,phase)%>%group_slice(1)%>%{txtplot(.$position,.$KL)}
#get offsets by picking top two KLs per phase, choosing rightmost as psite
kl_offsets <- kl_df%>%
	mutate(phase=position%%3)%>%
	group_by(set,sample,nreadlen,phase)%>%
	filter(position< -6)%>%
	filter(position> -(nreadlen-6))%>%
	mutate(rank=rank(-KL))%>%
	filter(rank%in%c(1:2))%>%
	arrange(sample,nreadlen,phase,position)
kl_offsets <- kl_df%>%
	mutate(phase=position%%3)%>%
	group_by(set,sample,nreadlen,phase)%>%
	arrange(position)%>%
	# mutate(codon_KL = KL+lag(KL)+lag(lag(KL)))%>%
	mutate(codon_KL = KL)%>%
	filter(position< -6)%>%
	filter(position> -(nreadlen-6))%>%
	group_by(sample,nreadlen)%>%
	slice(which.max(codon_KL))%>%
	arrange(sample,nreadlen,phase,position)
#
kl_offsets <- kl_offsets%>%
	group_by(sample,nreadlen,phase)%>%
	# filter(sample=='E13_ribo_1',nreadlen==25)%>%
	slice(which.max(position))%>%
	summarise(offset=-(position+3))%>%
	mutate(readlen=paste0('rl',nreadlen))
#
kl_offsets <- kl_df%>%group_by(nreadlen,position)%>%
	summarise(sumKL = sum(KL))%>%
	filter(position< -6)%>%
	filter(position> -(nreadlen-6))%>%
	slice(which.max(sumKL))%>%
	mutate(offset = (-position)-3)
#
most_freq <- function(x) x%>% table%>%sort%>%names%>%as.numeric%>%tail(1)
#no equalize all timepoints
# kl_offsets <- kl_offsets%>%
	# separate(sample,c('time','assay','rep'),remove=F)%>%
	# filter(assay=='ribo')%>%
	# group_by(nreadlen)%>%
	# group_slice(5)%>%
	# mutate(offset = offset%>%most_freq)
#
clean_fr_sampnames<-function(x) x%>%str_replace('.*_(Poly|80S)(.*)_()','\\2_\\1ribo_\\3')
kl_offsets2plot <- kl_offsets%>%
	# separate(sample,c('time','assay','rep'),remove=F)%>%
	# filter(assay=='ribo')
	mutate(readlen = paste0('rl',nreadlen))

plotfile='plots/rust_fppos_vs_codon_variance.pdf'
pdf(plotfile,w=12,h=3*n_distinct(kl_df$nreadlen))
#plotting variance amongst codons at each point.
# sh_codprof%>%
kl_df%>%
	filter(position< -3)%>%
	filter(position> -(nreadlen-6))%>%
	separate(sample,c('time','assay','rep'),remove=F)%>%
	filter(time%>%is_in(names(stagecols)))%>%
	filter(assay=='ribo')%>%
	{
		qplot(data=.,x=position,y=KL)+
		theme_bw()+
		facet_grid(nreadlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=kl_offsets2plot,aes(xintercept= -(offset)-3),color=I('green'),linetype=2)+
		geom_vline(data=kl_offsets2plot,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		# geom_vline(data=kl_offsets2plot,aes(xintercept= -offset-1),color=I('blue'),linetype=2)+
		# geom_vline(data=kl_offsets2plot,aes(xintercept= -offset-2),color=I('blue'),linetype=2)+
		# geom_vline(data=kl_offsets2plot,aes(xintercept= -offset-3),color=I('green'),linetype=2)+
		# geom_vline(data=kl_offsets2plot,aes(xintercept= -offset-4),color=I('green'),linetype=2)+
		# geom_vline(data=kl_offsets2plot,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
message(normalizePath(plotfile))
#

kl_offsets2plot%>%select(offset,length=nreadlen)%>%write_tsv('ext_data/offsets_rustvar.tsv')