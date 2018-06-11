#load up the counts for utr regions, cds


#first check if we can find differential expression in UTRs at all
#then check for it in 
featurecountfiles <- Sys.glob(file.path(root,'data/feature_counts/data/*/*feature_counts'))
timepoints	<- c("E13", "E145", "E16", "E175", "P0")
# comptimepoints = c(timepoints[c(1,5)])
comptimepoints = c(timepoints[c(1,4)])

#read in all the counts for different region files
regioncounts <-
	featurecountfiles %>%
	.[!str_detect(.,'tile')] %>%
	setNames(.,basename(.)) %>%
	map(fread,skip=1,select=c(1,7)) %>%
	map2(names(.),~ set_colnames(.x,c('Geneid',.y)))%>%
	Reduce(f=partial(left_join,by=c('Geneid')))

regioncounts %<>% as_data_frame
colnames(regioncounts)%<>%str_replace('feature_counts','')

assert_that( all((regioncounts%>%colnames%>%head)==c("Geneid", "E13_ribo_1.cds", "E13_ribo_1.fputrs", 
"E13_ribo_1.gene", "E13_ribo_1.tputrs", 
"E13_ribo_2.cds") ))

#
coldata=colnames(regioncounts)%>%data_frame(set=.)%>%slice(-1)%>%separate(set,c('time','assay','replicate','region'),remove=FALSE)
coldata$assay%<>%factor(c('total','ribo'))
cdscoldata = coldata%>%filter(region=='cds',time %in% comptimepoints)%>%{stopifnot(nrow(.)>3);.}
cdsdds <- DESeqDataSetFromMatrix(regioncounts[,dcols],colData=cdscoldata,design = ~ time*assay)

fputrcoldata = coldata%>%filter(region=='fputrs',time %in% comptimepoints)
dcols = c(fputrcoldata%>%.$set)
fputrdds <- DESeqDataSetFromMatrix(regioncounts[,dcols],colData=fputrcoldata,design = ~ time*assay)

tputrcoldata = coldata%>%filter(region=='tputrs',time %in% comptimepoints)
dcols = c(fputrcoldata%>%.$set)
tputrdds <- DESeqDataSetFromMatrix(regioncounts[,dcols],colData=tputrcoldata,design = ~ time*assay)

sizefacts <- regioncounts[,cdscoldata$set]%>%as.matrix%>%estimateSizeFactorsForMatrix



####Now do analysis for fp utrs
sizeFactors(fputrdds) <- sizefacts
sizeFactors(tputrdds) <- sizefacts

has_fp_counts <- ! counts(fputrdds)%>%apply(1,function(x) any(x==0))
has_cds_counts <- ! counts(cdsdds)%>%apply(1,function(x) any(x==0))
genes_to_use <- has_fp_counts & has_cds_counts

fputrdds <- fputrdds[genes_to_use,]
cdsdds <- cdsdds[genes_to_use,]

normalizationFactors(fputrdds) <- counts(cdsdds)
normalizationFactors(fputrdds) <- counts(fputrdds)

fputrdds <- estimateDispersionsGeneEst(fputrdds)
dispersions(fputrdds) <- mcols(fputrdds)$dispGeneEst
fputrdds <- nbinomWaldTest(fputrdds)

# assert_that(all(resultsNames(fputrdds) == c("Intercept", "as.numeric.time.", "assay_ribo_vs_total", "as.numeric.time..assayribo")))
intres <- results(fputrdds,paste0('time_',comptimepoints[2],'_vs_',comptimepoints[1])))



intres <- results(fputrdds,paste0('time_',comptimepoints[2],'_vs_',comptimepoints[1])))

assert_that(intres%>%as.data.frame%>%filter(padj<0.05)%>%nrow%>%`>`(40))

te_results <- file.path(root,'exploration/tables/riboseqres_P0.txt')

ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols=cols(
  col_character(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double()
)
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)

ribodiffcontrastobs <- te_results%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols

cor.test(ribodiffcontrastobs[[1]]$log2fc[match(regioncounts$Geneid, ribodiffcontrastobs[[1]]$feature_id)],intres$log2FoldChange)


####and now for tp utsrs
tputrcoldata = coldata%>%filter(region=='tputrs',time %in% c('E13','P0'))
dcols = c(fputrcoldata%>%.$set)
tputrcoldata
tputrdds <- DESeqDataSetFromMatrix(regioncounts[,dcols],colData=tputrcoldata,design = ~ as.numeric(time)*assay)

tputrdds <- DESeq(tputrdds)

assert_that(resultsNames(tputrdds) == c("Intercept", "as.numeric.time.", "assay_ribo_vs_total", "as.numeric.time..assayribo"))

intres <- results(fputrdds,list('as.numeric.time..assaytotal'))

assert_that(intres%>%as.data.frame%>%filter(padj<0.05)%>%nrow%>%`>`(40))

te_results <- file.path(root,'exploration/tables/riboseqres_P0.txt')
te_results%>%read_tsv

ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols=cols(
  col_character(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double()
)
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)

ribodiffcontrastobs <- te_results%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols

cor.test(ribodiffcontrastobs[[1]]$log2fc[match(regioncounts$Geneid, ribodiffcontrastobs[[1]]$feature_id)],intres$log2FoldChange)

