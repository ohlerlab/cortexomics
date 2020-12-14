source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
iso_tx_countdata <-readRDS(here('data/iso_tx_countdata.rds'))




#plot for TE up, TE down, extrTEup, extrTEdown, and the clusters.
#We should do the TE up/down as ranked lists as well.

highcountgnms <- readRDS(here('data/highcountgnms.rds'))
allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
allxtail = allxtail%>%filter(gene_name %in% highcountgnms)
techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
extrtechangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(2)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
techangegenes = techangedf%>%filter(up==1|down==1)%>%.$gene_name
teupgenes = techangedf%>%filter(up==1)%>%.$gene_name
tedowngenes = techangedf%>%filter(down==1)%>%.$gene_name
extrtechangegenes = extrtechangedf%>%filter(up==1|down==1)%>%.$gene_name
extrteupgenes = extrtechangedf%>%filter(up==1)%>%.$gene_name
extrtedowngenes = extrtechangedf%>%filter(down==1)%>%.$gene_name
notechangegenes = techangedf$gene_name%>%setdiff(techangegenes)
pcadf = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/data/pca_t0.tsv'%>%read_tsv
terankedgnms = pcadf%>%select(gene_name=gene.name,PC1_TE)%>%arrange(-PC1_TE)%>%.$gene_name
techangescores = pcadf%>%select(gene_name=gene.name,PC1_TE)%>%arrange(-PC1_TE)%>%{setNames(.$PC1_TE,.$gene_name)}
techangescores%<>%{rank(.)/(length(.)+1)}%>%{qnorm(.)}

long_exons = exons%>%split(.,.$transcript_id)%>%width%>%sum%>%enframe%>%filter(!is.na(name))%>%mutate(g_id=trid2gid[[name]])%>%
	group_by(g_id)%>%
	slice(which.max(value))%>%.$name
long_exons = (exons%>%split(.,.$transcript_id))[long_exons]
names(long_exons)%<>%trid2gnm[[.]]
long_cds = cds%>%split(.,.$transcript_id)%>%width%>%sum%>%enframe%>%filter(!is.na(name))%>%mutate(g_id=trid2gid[[name]])%>%
	group_by(g_id)%>%
	slice(which.max(value))%>%.$name
long_cds = (cds%>%split(.,.$transcript_id))[long_cds]
names(long_cds)%<>%trid2gnm[[.]]
long_fputr = fputrs%>%width%>%sum%>%enframe%>%filter(!is.na(name))%>%mutate(g_id=trid2gid[[name]])%>%
	group_by(g_id)%>%
	slice(which.max(value))%>%.$name
long_fputr = fputrs[long_fputr]
names(long_fputr)%<>%trid2gnm[[.]]
long_tputr = tputrs%>%width%>%sum%>%enframe%>%filter(!is.na(name))%>%mutate(g_id=trid2gid[[name]])%>%
	group_by(g_id)%>%
	slice(which.max(value))%>%.$name
long_tputr = tputrs[long_tputr]
names(long_tputr)%<>%trid2gnm[[.]]

library(Rsamtools)

long_cdsseq = GenomicFeatures::extractTranscriptSeqs(long_cds,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))
long_fputrseq = GenomicFeatures::extractTranscriptSeqs(long_fputr,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))
long_tputrseq = GenomicFeatures::extractTranscriptSeqs(long_tputr,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))
long_trsseq = GenomicFeatures::extractTranscriptSeqs(long_exons,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))

regseqlist = list(
	cds = long_cdsseq,
	fputr = long_fputrseq,
	tputr = long_tputrseq,
	tr = long_trsseq
)

setlist = list(
	random=sample(terankedgnms,1e3),
	teranked=terankedgnms,
	terev=rev(terankedgnms),
	terankeduponly=terankedgnms%>%intersect(teupgenes),
	terankednodown=terankedgnms%>%intersect(c(teupgenes,notechangegenes)),
	techangegenes=techangegenes,
	teupgenes=teupgenes,
	tedowngenes=tedowngenes,
	extrtechangegenes=extrtechangegenes,
	extrteupgenes=extrteupgenes,
	extrtedowngenes=extrtedowngenes,
	teupgeneslmatch=teupgenes,
	tedowngeneslmatch=tedowngenes
)

testset = teupgenes%>%head(200)
# mdata = mdata
getmatchingset<-function(testset,mdata,ratio=1){
	mdata$testvar <- mdata[[1]] %in% testset
	stopifnot(sum(mdata[['testvar']]) > 10)
	matchformula <- as.formula(paste0('testvar ~ ',colnames(mdata)%>%str_subset(neg=T,'testvar|gene_name|gene_id')%>%paste0(collapse='+')) )
	matchobject <- MatchIt::matchit(matchformula,data=mdata[,-1],ratio=ratio)
	mdata[[1]][matchobject[[1]]%>%as.numeric]
}

mdata = enframe(sum(width(long_tputr)),'gene_name','length')
fmdata = enframe(sum(width(long_fputr)),'gene_name','length')


ctlsetlist = list(
	random=sample(terankedgnms,1e3),
	teranked=NULL,
	terankeduponly=NULL,
	terankednodown=NULL,
	techangegenes=notechangegenes,
	teupgenes=notechangegenes,
	teupgeneslmatch=getmatchingset(teupgenes,mdata%>%filter(gene_name%in%c(teupgenes,notechangegenes,tedowngenes)),ratio=3),
	tedowngenes=notechangegenes,
	tedowngeneslmatch=getmatchingset(tedowngenes,fmdata%>%filter(gene_name%in%c(tedowngenes,notechangegenes,tedowngenes)),ratio=3),
	extrtechangegenes=notechangegenes,
	extrteupgenes=notechangegenes,
	extrtedowngenes=notechangegenes
)

# stop()

regnm='tputr'
setnm='terev'
fafile = str_interp('pipeline/motseqs/${regnm}/${setnm}/${setnm}.fa')
if(!file.exists(fafile)){
	library(Biostrings)
	for(regnm in names(regseqlist)){
		for(setnm in names(setlist)){
			set = setlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
			seqs = regseqlist[[regnm]][set]
			if(str_detect(setnm,'teranked')){
				names(seqs)%<>%{paste0(.,'|',round(techangescores[set][.],3))}
			}
			if(str_detect(setnm,'terev')){
				names(seqs)%<>%{paste0(.,'|',round(-techangescores[set][.],3))}
			}
			fafile = str_interp('pipeline/motseqs/${regnm}/${setnm}/${setnm}.fa')
			fafile %>% dirname %>%dir.create(showWarn=F,rec=TRUE)
			seqs%>%{Biostrings::writeXStringSet(.,fafile)}
			controlset = ctlsetlist[[setnm]]%>%intersect(names(regseqlist[[regnm]]))
			controlseq = regseqlist[[regnm]][controlset]
			ctlfafile = str_interp('pipeline/motseqs/${regnm}/${setnm}/${setnm}.background.fa')
			controlseq%>%{Biostrings::writeXStringSet(.,ctlfafile)}
			#	
			message(normalizePath(fafile))
		}
	}
}

motifmatches = vmatchPattern('TGTA',long_tputrseq)%>%lapply(as.data.frame)%>%bind_rows(.id='gene')

long_tputrseq['Satb2']

readDNAStringSet('pipeline/motseqs/fputr/tedowngeneslmatch/tedowngeneslmatch.fa')%>%nchar%>%log2%>%txtdensity
readDNAStringSet('pipeline/motseqs/fputr/tedowngeneslmatch/tedowngeneslmatch.background.fa')%>%nchar%>%log2%>%txtdensity

readDNAStringSet('pipeline/motseqs/fputr/teupgeneslmatch/teupgeneslmatch.fa')%>%nchar%>%log2%>%.[between(.,4,11)]%>%txtdensity
readDNAStringSet('pipeline/motseqs/fputr/teupgeneslmatch/teupgeneslmatch.background.fa')%>%nchar%>%log2%>%.[between(.,4,11)]%>%txtdensity

readDNAStringSet('pipeline/motseqs/tputr/tedowngeneslmatch/tedowngeneslmatch.fa')%>%nchar%>%log2%>%.[between(.,4,15)]%>%txtdensity
readDNAStringSet('pipeline/motseqs/tputr/tedowngeneslmatch/tedowngeneslmatch.background.fa')%>%nchar%>%log2%>%.[between(.,4,15)]%>%txtdensity

readDNAStringSet('pipeline/motseqs/tputr/teupgeneslmatch/teupgeneslmatch.fa')%>%nchar%>%log2%>%.[between(.,4,15)]%>%txtdensity
readDNAStringSet('pipeline/motseqs/tputr/teupgeneslmatch/teupgeneslmatch.background.fa')%>%nchar%>%log2%>%.[between(.,4,15)]%>%txtdensity


#get the files on the max cluster
#conda activate biopython
#rsync  --prune-empty-dirs --progress -avzh --copy-links  --append-verify --include '*/' --include '*terev*' --include '*ranked*.fa' --exclude '*'  dharnet_m@login-1.research.hpc.bihealth.org:/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs  .
#Command to run on max cluster
#for fa in $( find  motseqs/ ! -iname 'teranked.fa' ! -iname '*background*' -iname '*ranked*' -iname '*.fa' ) ; do echo $fa "\n"; python /fast/AG_Akalin/wkopp/alison_sciatac_2018/scipipe_v4/scripts/cermit.py  -fasta ${fa} -output ${fa}_cermit -motifdb Mus_musculus.dna.meme  ;done
#send results back
#rsync   --prune-empty-dirs --progress -avzh --copy-links  --append-verify --include '*/' --exclude '*.fa' motseqs dharnet_m@login-1.research.hpc.bihealth.org:/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/motseqs




# #now add sccores for smart
# set = 'teranked'
# reg = 'fputr'



# python ssmart.py -f teranked.fa -o teranked.ssmart.out-c 2 -p $SMARTBIN
# python run_rnaprofiling.py -i teranked.fa -o rnaproffile -n {threads} 
                 