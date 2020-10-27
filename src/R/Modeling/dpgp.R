


dpgpdata<-timeff_ciddf%>%filter(time!='E13')%>%select(uprotein_id,ntime,logFC)%>%
	spread(ntime,logFC)

dpgpdata%<>%dplyr::rename('gene'=uprotein_id)

dpgpinput%>%dirname%>%dir.create
dpgpdata%>%write_tsv(dpgpinput)


dpgpbin <- '~/work/Applications/DP_GP_cluster/bin/DP_GP_cluster.py'%T>%{stopifnot(file.exists(.))}
dpgpinput <- 'pipeline/dpgp/ttest.txt'
dpgpoutput <- 'pipeline/dpgp/ttest'
condabin <-'~/work/miniconda3/envs/dpgp/bin/conda activate'
dpgpiter=1e3

dpgpcmd<-str_interp('source activate dpgp; ${dpgpbin} -i ${dpgpinput}  -o ${dpgpoutput} -p png -n ${dpgpiter}  --fast --do_not_mean_center --unscaled --save_cluster_GPs')
dpgpcmd%>%cat(file='tmp.sh')
system('bash tmp.sh')

dpgpoutfiles<-Sys.glob(str_interp('${dpgpoutput}*'))

 # BiocManager::install(c('ComplexHeatmap','circlize','colorspace','GetoptLong'))
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

dpgpoptclust <- fread(paste0(dpgpoutput,'_optimal_clustering.txt'))
dpgpsimmatdf <- fread(paste0(dpgpoutput,'_posterior_similarity_matrix.txt'))

dpgpsimmat<-dpgpsimmatdf[,-1]%>%as.matrix%>%set_rownames(dpgpsimmatdf[[1]])
dpgpsimmat%<>%{.[match(dpgpoptclust$gene,rownames(.),),]}
dpgpsimmat%<>%{.[,match(dpgpoptclust$gene,colnames(.))]}


pdf('tmp.pdf')
Heatmap(dpgpsimmat[1:1000,1:1000],cluster_columns=F,cluster_rows=F,show_row_names=F,show_column_names=F)
dev.off()

dpgpsimmat%>%table