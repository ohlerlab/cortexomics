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

if(!exists('fpcovlist')) fpcovlist<-readRDS(here('data/fpcovlist.rds'))
#	
fpcovlist <- fpcovlist[names(allbamtbls)]
sampname=fpcovlist%>%names%>%head(1)
sampfpcov=fpcovlist[[sampname]]
rlfpcov=sampfpcov[['29']]

width1grs <- function(gr){
	stopifnot(Negate(is.unsorted)(gr))
	isw1 <- width(gr)==1
	broad <- gr[!isw1]
	#vector of integers - 1,2,3 for range 1-3
	narrowstarts <- unlist(as(broad@ranges,'IntegerList'))
	narrow <- {GRanges(
			rep(seqnames(broad),width(broad)),
			IRanges(narrowstarts,w=1)
		)}
	mcols(narrow) <- mcols(broad)[rep(seq_along(broad),width(broad)),,drop=F]
	sort(c(gr[isw1],narrow))
}
regScoreSums<-function(srle,gr){
	scoregr=srle%>%as("GRanges")%>%subset(score!=0)%>%width1grs
	ov = findOverlaps(scoregr,gr)
	sumscores = tibble(reg=ov@to,score=scoregr$score[ov@from])%>%
		group_by(reg)%>%summarise_at('score',sum)
	tibble(reg=1:length(gr))%>%left_join(sumscores, by='reg')%>%
		mutate(score=replace_na(score,0))%>%
		pluck('score')%>%
		setNames(names(gr))
}
#

sampfpcov=fpcovlist[[1]]
if(!file.exists(here('data/trsums.rds'))){
	trsums <- 	mcmapply(mc.cores=8,SIMPLIFY=F,fpcovlist,names(fpcovlist),FUN=function(sampfpcov,sampname){
				sampfpcov = sampfpcov
				sampfpcov%>%map(regScoreSums,innercds)%>%purrr::reduce(.,`+`)
			})
	saveRDS(trsums,here('data/trsums.rds'))
}else{
	trsums<-readRDS(here('data/trsums.rds'))
	stopifnot(names(trsums)==names(fpcovlist))
	stopifnot(names(trsums[[1]])==names(innercds))
}



dtselgenelist[['allhigh']]<-trsums%>%map_df(.id='sample',enframe,'tr_id','count')%>%
	group_by(tr_id)%>%filter(sample%in%mainsamps)%>%
	summarise(high=all(count>32))%>%
	pluck('tr_id')

allcodlistnz = allcodlist%>%
	subset(seqnames%in%names(fpcovlist[[1]][[1]]))%>%
	.[overlapsAny(.,innercds,type='within')]
#	identity
	# head
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
names(allcodlistnz) <- seq_along(allcodlistnz)

zerotrs = which(trsums[[1]]==0)%>%names

# lapply(fpcovlist,function(sampfpcov) sampfpcov%>%map_dbl(function(rlfpcov){rlfpcov[allcodlistnz%>%subset(seqnames%in%zerotrs)]%>%sum%>%sum}))

# #FIND FUCKING INFS AAAAAGH
# infmatrix = fpcodonmats[[1]][[1]][[1]]%>%`==`(Inf)

# (fpcodonmats[[1]][[1]][[1]] * infmatrix)@i%>%head
# (fpcodonmats[[1]][[1]][[1]] * infmatrix)@p%>%head
# (fpcodonmats[[1]][[1]][[1]] * infmatrix)%>%summary%>%head
# fpcodonmats[[1]][[1]][[1]]%>%summary%>%{.[.[,'x']==Inf]}


# fpcodonmats[[1]][[1]][[1]]%>%summary%>%as.data.frame%>%filter(x%in%Inf)
# allcodmat%>%summary%>%as.data.frame%>%filter(x%in%Inf)%>%head
# out[[2]]%>%summary%>%as.data.frame%>%filter(x%in%Inf)%>%head

# fpcodonmats[[1]][[1]][[1]]%>%is.finite%>%not%>%which
# fpcodonmats[[1]][[1]][[1]]

stop()

if(!file.exists(here('data/fpcodonmats.rds'))){
	
	sampfpcov=fpcovlist[[1]]
	sampname=names(fpcovlist)[[1]]
	rlfpcov = sampfpcov[[1]]

	fpcodonmats <- 	mcmapply(mc.cores=1,SIMPLIFY=F,fpcovlist,names(fpcovlist),FUN=function(sampfpcov,sampname){
		trsum = trsums[[sampname]]
		sampfpcov%>%lapply(safely(function(rlfpcov){
			rlfpcov = rlfpcov[names(trsum)]
			rlfpcov = rlfpcov/(trsum)
			cat('.')
			out = rlfpcov[allcodlistnz]
			stopifnot(out%>%max(na.rm=T)%>%max(na.rm=T)%>%is.finite)
			out = out %>%split(cods)%>%lapply(as.matrix)
			out = lapply(out,Matrix::Matrix,sparse=TRUE)
			out
		}))
	})
	saveRDS(fpcodonmats,here('data/fpcodonmats.rds'))
}else{
	fpcodonmats<-readRDS(here('data/fpcodonmats.rds'))
	stopifnot(names(fpcodonmats)==names(fpcovlist))
}
if(fpcodonmats[[1]][[1]]%>%names%>%is_in('result')%>%any){
	fpcodonmats%<>%map_depth(2,'result')
}



# #testing presence of infinities
# fpcodonmats[[1]][[1]][['TTT']]%>%is.finite%>%not%>%rowSums%>%`==`(0)%>%not%>%which
# 5078
# allcodlistnz[5078]
# trsums[[1]]['ENSMUST00000005826']

# allcodmat <- Matrix::Matrix(sparse=T,0,nrow=length(allcodlistnz),ncol=width(allcodlistnz[1]))

# object.size(allcodmat)

# maxl = runLength(rlfpcov)%>%sum%>%max
# lapply(rlfpcov,function(rle)c(rle,Rle(0,maxl-length(rle))))%>%as("RleList")%>%as("GRanges")

###

			# codwindcov = rlfpcov%>%
			# 	as("GRanges")%>%
			# 	subset(score!=0)%>%
			# 	width1grs%>%
			# 	{o=mapToTranscripts(.,allcodlistnz);o$score = .$score[o$xHits];o}%>%
			# 	sort%>%
			# 	width1grs
			# if(any(!is.finite(codwindcov$score)))browser()
			# # seq_int_rle = match(codwindcov@seqnames,seq_along(allcodlistnz))
			# maxseq = length(allcodlistnz)
			# myp = cumsum(c(1,runLength(sort(c(Rle(1:maxseq),seqnames(codwindcov)))))-1)
			# # myp%>%length
			# # ((start(codwindcov)%>%length) == (diff(myp)%>%sum))
			# # myp = cumsum(c(runLength(sort())-1))
			# myspmat <- Matrix::sparseMatrix(
			# 	i=start(codwindcov),
			# 	p=myp,
			# 	x=codwindcov$score,
			# 	dims = c(max(start(codwindcov)),length(myp)-1)
			# )
			# myspmat<-Matrix::t(myspmat)
			# myspmat = lapply(unique(cods)%>%setNames(.,.),function(cod) myspmat[cods==cod,])
			# if(any(!is.finite(codwindcov$score))) browser()
			# myspmat
			
			# # myspmat%>%is.finite%>%`!`%>%rowSums%>%is_in(1)
			# windsize = width(allcodlistnz[1])
			# windsize = 3
			# maxl = runLength(rlfpcov)%>%sum%>%max
			# allcodmat <- Matrix::Matrix(sparse=TRUE,0,nrow=length(allcodlistnz),ncol=windsize)
			# codwindcov <- rlfpcov%>%as("GRanges")%>%subset(score!=0)%>%width1grs
			
			# rlfpcov%>%as("GRanges")%>%subset(score!=0)%>%subsetByOverlaps(allcodlistnz)%>%.$score%>%max


			# stopifnot(all(seqlevels(codwindcov)==names(rlfpcov)))
			# maxseq = length(rlfpcov)
			# intseqnames = seqnames(codwindcov)%>%match(.,unique(.))
			# myp = cumsum(c(1,runLength(sort(c(Rle(1:maxseq),intseqnames))))-1)
			# mysp_covmat <- Matrix::sparseMatrix(
			# 	i=start(codwindcov),
			# 	p=myp,
			# 	x=codwindcov$score,
			# 	dims = c(maxl,length(myp)-1)
			# )
			# allcodseqinds = as.vector(allcodlistnz@seqnames%>%match(.,unique(.)))
			# for(i in 1:windsize){
			# 	indmat <- c(start(allcodlistnz)+i-1,allcodseqinds)
			# 	indmat <- matrix(indmat,ncol=2)
			# 	#get our values
			# 	insx = mysp_covmat[indmat]
			# 	nz = insx!=0
			# 	# dim(allcodmat)
			# 	#enter the nonzero entries
			# 	allcodmat[matrix(c(which(nz),rep(1,sum(nz))),ncol=2)] <- insx[nz]
			# }
			# # allcodmat = lapply(unique(cods)%>%setNames(.,.),function(cod) allcodmat[cods==cod,])
			
			# mym =allcodmat[1:10,1:3]
			# mym[matrix(ncol=2,c(4,2))]<-3
			# mym


			# allcodmat%>%summary%>%as.data.frame%>%mutate(row=1:n())%>%filter(j==1,x%in%Inf)%>%head
			# allcodmat[50735,1]
			# #okay so row 50735 is inf,which is 1725 of the indmat
			# which(insx==Inf)
			# mysp_covmat[indmat[146250,,drop=F]]
			# #okay, so this is infinite - this position, seqnames
			# indmat[146250,,drop=F]
			# #which should correspond to this one			
			# rlfpcov[allcodlistnz[146250]]
			# rlfpcov[allcodlistnz[146250]@seqnames]
			# #which using native subsetting,doesn't have any infinities....
			# mysp_covmat[,allcodseqinds[146250]]%>%head(300)%>%`>`(0)%>%which


			# seqnames(allcodlistnz[323])
			# # all(myspmat[[1]]==allcodmat[[1]])



			# dim(myspmat)
			# length(cods)
###

# dim(allcodmat)
# library(Matrix)
# tn=5^2
# tc = sqrt(tn)
# testmat <- matrix(1:tn,ncol=sqrt(tn))
# testmat[(1:tn)%in%sample(1:tn,floor(tn*0.9))]<-0
# sptestmat <- Matrix(testmat,sparse=TRUE)
# j = rep(seq_along(sptestmat@p[-1]),diff(sptestmat@p))
# j = Rle(j)
# myp = cumsum(c(runLength(sort(c(Rle(0:tc),j)))-1))
# sptestmat@p==myp

# dim(sptestmat)

# cumsum(runLength(j))

# sampfpcov[[5]]%>%lapply()
# Rle2Mat<-function(from) {
#     rv <- runValue(from)
#     nz <- rv != 0
#     i <- as(ranges(from)[nz],'IntegerList')%>%unlist
#     x <- rep(rv[nz], runLength(from)[nz])
#     sparseMatrix(i=i, p=c(0L, length(x)), x=x,
#                  dims=c(length(from), 1))
# }
# out[[1]] %>%Rle2Mat
# out[[1]]%>%runValue

# from=out[[1]]
# # Convert from DataFrame of Rle to sparse Matrix
# #
# setAs("DataFrame", "Matrix", function(from) {
#   mat = do.call(cbind, lapply(from, as, "Matrix"))
#   colnames(mat) <- colnames(from)
#   rownames(mat) <- rownames(from)
#   mat
# })

# names(allcodlistnz) <- seq_along(allcodlistnz)
# codwindcov = rlfpcov%>%
# 	as("GRanges")%>%
# 	subset(score!=0)%>%
# 	width1grs%>%
# 	{o=mapToTranscripts(.,allcodlistnz);o$score = .$score[o$xHits];o}%>%
# 	sort%>%
# 	width1grs
# seq_int_rle = match(codwindcov@seqnames,seq_along(allcodlistnz))
# maxseq = length(allcodlistnz)
# myp = cumsum(c(0,runLength(sort(c(Rle(1:maxseq),seq_int_rle)))-1))
# myspmat <- Matrix::sparseMatrix(
# 	i=start(codwindcov),
# 	p=myp,
# 	x=codwindcov$score
# )
# myspmat%<>%t
# myspmat = lapply(unique(cods)%>%setNames(.,.),function(cod) myspmat[cods==cod,])

# (!(myspmat[['AAA']] == spmat[['AAA']]))

# trsum%>%min

# is.na(spmat[['AAA']])%>%sum
# is.na(myspmat[['AAA']])%>%sum
# out%>%is.na%>%sum

# #why these infinite values


if(!file.exists(here('data/subfpprofilelist.rds'))){

	

	row_gnames <- seqnames(allcodlistnz)%>%split(cods)
	matlist=fpcodonmats[[1]][[1]]
	.x='AAC'
	get_glist_profiles<-function(glist,row_gnames,fpcodonmats){
		matlist = fpcodonmats[[1]][[1]]
		fpcodonmats%>%map_depth(2,function(matlist)map(names(matlist)%>%setNames(.,.),function(.x){
			matlist[[.x]][as.vector(row_gnames[[.x]])%in%unique(glist),]%>%colMeans(na.rm=T)
		}))
	}
	glist=dtselgenelist[['nochange']]
	get_codon_prof_df <- function(glist,row_gnames,fpcodonmats){
		codonprofiledat <- get_glist_profiles(glist,row_gnames,fpcodonmats)	
		codonprofiledat <- codonprofiledat%>%
			map_depth(3,.%>%
				enframe('position','count'))%>%
				map_df(.id='sample',.%>%
					map_df(.id='readlen',.%>%
						bind_rows(.id='codon')))
		# codonprofiledat <- codonprofiles%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
		codonprofiledat %<>% mutate(position = position - 1 - (FLANKCODS*3))
		# codonprofiledat %<>% group_by(subgroup,sample,readlen,codon)%>%mutate(count= count / median(count))
		codonprofiledat %<>% group_by(sample,readlen,codon)%>%
			# mutate(count= count / median(count))
			identity
		codonprofiledat %<>% filter(!codon %in% c('TAG','TAA','TGA'))
		codonprofiledat$readlen%<>%str_replace('^(\\d)','rl\\1')
		codonprofiledat
	}
	subfpprofilelist <- dtselgenelist%>%lapply(get_codon_prof_df,row_gnames,fpcodonmats)
	#subfpprofilelist

	saveRDS(subfpprofilelist,here('data/subfpprofilelist.rds'))
}else{
	subfpprofilelist<-readRDS(here('data/subfpprofilelist.rds'))
}



#RUST
if(!file.exists(here('data/subfprustprofilelist.rds'))){
	subfprustprofilelist <-
		mclapply(mc.cores=1,dtselgenelist['all'],function(seltrs){
			imap(fpcovlist[mainsamps[1]],function(sampfpcov,sampname){
				trsums = sampfpcov%>%map(~.[innercds[seltrs]])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
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
fprustprofilelist<-readRDS(here('data/fprustprofilelist.rds'))


rust_rocl<-readRDS(here('data/subfprustprofilelist.rds'))[['all']]

sampname<-mainsamps[1]
sampfpcov<-fpcovlist[[sampname]]
rlfpcov <- sampfpcov[['29']]
seltrs=TRUE
 

rust_roel <- 	mclapply(mc.cores=1,dtselgenelist['all'],function(seltrs){
			imap(fpcovlist[mainsamps[1]],function(sampfpcov,sampname){
				# trsums = sampfpcov%>%map(~.[innercds[seltrs]])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
				samptrsums = trsums[[sampname]]
				sampfpcov['29']%>%lapply(function(rlfpcov){
					browser()
					rlfpcov = rlfpcov[seltrs]
					innercdsmeans <- samptrsums/width(innercds[names(samptrsums)])
					rlfpcov = rlfpcov > innercdsmeans[names(rlfpcov)]
					
					nz_trs <- names(samptrsums)[samptrsums!=0]
					allcodlistnz = allcodlist%>%
						subset(seqnames%in%nz_trs)
					cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
					('.')
					rustvalmat = rlfpcov[allcodlistnz]%>%split(cods)%>%
						lapply(as.matrix)%>%
						map(colMeans)
					rustvalmat
					
				})
			})
		})



frustprofilelist <- fprustprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
frustprofilelist %<>% mutate(position = position - 1 - (FLANKCODS*3))
frustprofilelist %<>% filter(!codon %in% c('TAG','TAA','TGA'))
frustprofilelist$readlen%<>%str_replace('^(\\d)','rl\\1')



# get_codon_prof_df(names(fpcovlist[[1]][[1]]),row_gnames,fpcodonmats)
# codonprofiledat <- fprustprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
# # codonprofiledat <- codonprofiles%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
# codonprofiledat %<>% mutate(position = position - 1 - (FLANKCODS*3))
# # codonprofiledat %<>% group_by(subgroup,sample,readlen,codon)%>%mutate(count= count / median(count))
# codonprofiledat %<>% group_by(sample,readlen,codon)%>%
# 	# mutate(count= count / median(count))
# 	identity
# codonprofiledat %<>% filter(!codon %in% c('TAG','TAA','TGA'))
# codonprofiledat$readlen%<>%str_replace('^(\\d)','rl\\1')


# subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-11])%>%enframe('codon','directdt')%>%left_join(codondata%>%filter(time=='E13',rep=='1'))%>%
# 	{quicktest(.$directdt,.$dwell_time)}
# subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-11])%>%enframe('codon','directdt')%>%left_join(codondata%>%filter(time=='E13',rep=='1'))%>%
	# {quicktest(.$directdt,.$availability)}

# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codondata$dwell_time)
# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codonprofiledat$count)
# subfpprofilelist[[1]][[1]][['29']]%>%unlist%>%is_in(codonprofiledat$count)


# rustcodon_dts = rustprofiledat%>%
# 	mutate(length=as.numeric(readlen))%>%
#     safe_left_join(offsets%>%select(length,offset))%>%
#     filter(position== -offset-3)%>%
# 	group_by(subgroup,sample,codon)%>%
#     summarise(dwell_time = mean (count))
# #now plot
# plotfile<- here(paste0('plots/','subgroup_dts_density','.pdf'))
# pdf(plotfile,w=10,h=10)
# rustcodon_dts%>%
# 	separate(sample,c('time','assay','rep'))%>%
# 	ggplot(.,aes(x=dwell_time,fill=subgroup))+
# 	geom_density(alpha=I(0.5))+
# 	scale_x_continuous(paste0('Codon Dwel Time'))+
# 	facet_grid(time~rep)+
# 	ggtitle(paste0('Codon DT for Up vs Down dTE genes'))+
# 	theme_bw()
# dev.off()
# message(normalizePath(plotfile))


clean_fr_sampnames<-function(x) x%>%str_replace('.*_(Poly|80S)(.*)_()','\\2_\\1ribo_\\3')
# subfpprofilelist%>%names
geneset='allhigh'
codonvarprofiles <-  subfpprofilelist[[geneset]]%>%
# rustprofiledat%>%
	ungroup%>%
	mutate_at('sample',clean_fr_sampnames)%>%
	group_by(sample,readlen,position)%>%
	filter(!is.nan(count))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(assay,rep,time,readlen,position)%>%
	# mutate(count = count / mean(count))%>%
	summarise(sdsig=sd(count,na.rm=T))%>%
	group_by(assay,readlen,time,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	# filter(position> -6,position < 6)%>%
	filter(position> -numreadlen+6,position< -6)%>%
	filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)
#
codonvaroffsets = codonvarprofiles%>%
	mutate(phase = position %%3)%>%
	group_by(assay,readlen,time)%>%
	# slice(which.max(sdsig+lag(sdsig)+lag(lag(sdsig))))%>%
	slice(which.max(sdsig))%>%
	mutate(offset=-(position+3))%>%
	as.data.frame
codonvaroffsets<-codonvaroffsets%>%group_by(assay,readlen,phase)%>%
	mutate(offset = offset%>%table%>%sort%>%tail(1)%>%names%>%as.numeric)
#the poly codvar offsets look weird
codonvaroffsets%<>%filter(assay=='ribo')
#
codonvaroffsets = subfpprofilelist[[geneset]]%>%
	ungroup%>%
	distinct(sample)%>%
	mutate_at('sample',clean_fr_sampnames)%>%
	separate(sample,c('time','assay','rep'),remove=F)%>%
	left_join(codonvaroffsets,by='time')%>%
	select(sample,readlen,numreadlen,offset,sample)
codonvaroffsets2plot <- codonvaroffsets%>%separate(sample,c('time','assay','rep'))%>%filter(assay=='ribo')
#
{
offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_fppos_vs_codon_variance.pdf'
pdf(plotfile,w=12,h=3*n_distinct(codonvarprofiles$readlen))
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	filter(time%>%is_in(names(stagecols)))%>%
	filter(assay=='ribo')%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-1),color=I('blue'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-2),color=I('blue'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-3),color=I('green'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-4),color=I('green'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}
#

{
codonvaroffsets2plot <- codonvaroffsets%>%separate(sample,c('time','assay','rep'))%>%filter(!assay=='ribo')
offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_fppos_vs_codon_variance_frac.pdf'
pdf(plotfile,w=12,h=3*n_distinct(codonvarprofiles$readlen))
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	filter(assay=='Polyribo')%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}

{
offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_fppos_vs_codon_variance_zoom.pdf'
pdf(plotfile,w=6,h=3*1)
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	filter(time%>%is_in(names(stagecols)))%>%
	filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29'))%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot%>%filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29')),aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=codonvaroffsets2plot%>%filter(time%>%is_in(c('E13','E175')),readlen%>%is_in('rl29')),aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}
