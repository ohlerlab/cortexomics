library(data.table)

source(here('src/R/Rprofile.R'))
trid2gid = metainfo%>%
	distinct(transcript_id,gene_id)%>%
	{safe_hashmap(.[['transcript_id']],.[['gene_id']])}

#get exons
gtf = here('pipeline/my_gencode.vM12.annotation.gtf')
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
if(is.function('exons'))exons <- gtf_gr%>%subset(type=='exon')
if(is.function('cds'))cds <- gtf_gr%>%subset(type=='CDS')

trlens = exons%>%split(.,.$transcript_id)%>%width%>%sum
cdslens = cds%>%split(.,.$transcript_id)%>%width%>%sum
deepshapedata = fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/deepshape/tr_pos_ribosig')
deepshapedata = deepshapedata%>%mutate(V1=str_replace(V1,'\\.\\d+',''))%>%set_colnames(c('seqnames','start','score'))%>%
	mutate(end=start)%>%GRanges
seqlengths(deepshapedata) = cdslens[names(seqlengths(deepshapedata))]
deepshapecov = deepshapedata%>%coverage(weight='score')
deepshapecovsums = sum(deepshapecov)

dpfiles = Sys.glob(here('pipeline/deepshapeprime/*/run*'))%>%
	data.frame(file=.)%>%
	mutate(sample=basename(dirname(file)))%>%
	mutate(name=basename(file))%>%
	mutate(number=str_extract(name,'\\d+'))%>%
	arrange(desc(as.numeric(number)))%>%
	group_by(sample)%>%slice(1)%>%{setNames(.$file,.$sample)}

dpexprdata = dpfiles%>%map_df(.id='sample',fread)
dpexprdata%<>%set_colnames(c('sample','transcript_id','TPM','diff','count'))
dpexprdata%<>%mutate(transcript_id=str_replace(transcript_id,'\\.\\d+',''),gene_id = trid2gid[[transcript_id]])

#testing aggreement with previous highcount genes
dpexprdata%>%group_by(gene_id %in% highcountgenes,gene_id)%>%filter(sample%in%mainribosamps)%>%summarise(highcount=any(count>=32))%>%summarise(as.data.frame(table(highcount)))


ms_id2protein_id



dptpms = mean(deepshapecov)%>%{./sum(.)}%>%multiply_by(1e6)%>%enframe('transcript_id','TPM')

#load the deepshape results
deepshapedf = fread(here('pipeline/deepshape/DeepShapePrimeOutputs/runlog_199.txt'))%>%head(-1)%>%
	select(transcript_id=1,TPM=2)%>%
	mutate(transcript_id=str_replace(transcript_id,'\\.\\d+',''))


# #now calculate spec. coefficients.
# deepshapespevals = deepshapecovinorfquant%>%mclapply(.%>%as.vector%>%ftestvect)
# deepshapespecs = deepshapespevals%>%map_dbl(1)%>%sqrt%>%enframe('transcript_id','spec_coef')
# deepshapespec_pM = deepshapespecs%>%mutate(spec_pM = spec_coef / (sum(spec_coef,na.rm=T)) * 1e6)


# deepshapespewave = deepshapecovinorfquant%>%
# 		mclapply(mc.cores=10,function(v)DWPT(as.vector(v)))
# deepshapespewavemn = deepshapespewave%>%map_dbl(mean)

# txtplot(log1p(deepshapespewavemn),log1p(deepshapecovinorfquant%>%mean))

################################################################################
########Just playing with the stats here
################################################################################
	


txtplot(deepshapecovinorfquant%>%runLength%>%sum,(deepshapespevals%>%map_dbl(1)))


txtplot(deepshapecovinorfquant%>%mean,sqrt((deepshapespevals%>%map_dbl(1))))


txtdensity(lm(deepshapecovinorfquant%>%mean~sqrt((deepshapespevals%>%map_dbl(1))))$residuals)

library(outliers)

cooksd <- cooks.distance(lm(deepshapecovinorfquant%>%mean~sqrt((deepshapespevals%>%map_dbl(1)))))

isoutlier = cooksd > 4/length(cooksd)
