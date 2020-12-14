



# some dots seem to have been replaced with commas, need to match the supp table gene names
t_timecoretbl$gene_name%<>%str_replace(',','.')


#knnMap$M           : is a matrix 16441x2061 containing the RPM values of all genes (16441) for all cells (2061)
#knnMap$grid        : is a precomputed grid of size 512x512 where genes expression are infered using 15 nearest neighboring cells.
#knnMap$grid$(xi,yi): are the coordinates of the elements in the grid
#knnMap$grid$knn    : are precomputed indices of the 15 nearest neighbooring cells of each pixel of the grid


timecellscores = t_timecoretbl$Specificity %*% knnMap$M[t_timecoretbl$gene_name,]
diffcellscores = t_diffcoretbl$Specificity %*% knnMap$M[t_diffcoretbl$gene_name,]

ss_emat <- projmemoise(fread)(Sys.glob(here('ext_data/GSE11*')))
rpmmat = ss_emat[,-1]%>%sweep(2,STAT=colSums(.),FUN='/')%>%multiply_by(1e6)

all(ss_emat[,1]==rownames(knnMap$M))

cor(rpmmat[,1],knnMap$M[,1])

apply(rpmmat,2,F=cor,knnMap$M[,1])%>%sort%>%tail

bestcountexprdata = countexprdata[rownames(countexprdata)%in%(metainfo%>%filter(isbest)%>%.$protein_id)]
rownames(bestcountexprdata) = tibble(protein_id=rownames(bestcountexprdata))%>%left_join(metainfo%>%filter(isbest)%>%distinct(gene_name,protein_id,transcript_id,uprotein_id))%>%
  .$gene_name
bulkdata=assayData(bestcountexprdata)$exprs
bulkdata=bulkdata%>%sweep(2,STAT=colSums(.),FUN='/')%>%multiply_by(1e6)

bt_timecoretbl = t_timecoretbl%>%filter(gene_name%in%rownames(bulkdata))
bt_diffcoretbl = t_diffcoretbl%>%filter(gene_name%in%rownames(bulkdata))
bulktimecellscores = bt_timecoretbl$Specificity %*% bulkdata[bt_timecoretbl$gene_name,]
bulkdiffcellscores = bt_diffcoretbl$Specificity %*% bulkdata[bt_diffcoretbl$gene_name,]

#now plot
plotfile<- here(paste0('plots/','hm_cellrpm_tdscores','.pdf'))
grDevices::pdf(plotfile)
qplot(y=diffcellscores,x=timecellscores)+
	scale_x_continuous(paste0('time score'))+
	scale_y_continuous(paste0('differentiation score'))+
	ggtitle(paste0('Time/DiffScore based on heatmap cell RPMs'))+
	theme_bw()
dev.off()
normalizePath(plotfile)


bt_timecoretbl = t_timecoretbl%>%filter(gene_name%in%rownames(bulkdata))
bt_diffcoretbl = t_diffcoretbl%>%filter(gene_name%in%rownames(bulkdata))
bulktimecellscores = bt_timecoretbl$Specificity %*% bulkdata[bt_timecoretbl$gene_name,]
bulkdiffcellscores = bt_diffcoretbl$Specificity %*% bulkdata[bt_diffcoretbl$gene_name,]

timecellscores = bt_timecoretbl$Specificity %*% knnMap$M[bt_timecoretbl$gene_name,]
diffcellscores = bt_diffcoretbl$Specificity %*% knnMap$M[bt_diffcoretbl$gene_name,]

#now plot
plotfile<- here(paste0('plots/','hm_cellrpm_tdscores_withbulk','.pdf'))
grDevices::pdf(plotfile)
qplot(y=diffcellscores,x=timecellscores)+
	scale_x_continuous(paste0('time score'))+
	scale_y_continuous(paste0('differentiation score'))+
	ggtitle(paste0('Time/DiffScore based on heatmap cell RPMs'))+
	theme_bw()
dev.off()
normalizePath(plotfile)



foldchangecatdf<-readRDS(here('data/foldchangecatdf.rds'))

fcdf <- foldchangecatdf%>%group_by(gene_id)%>%filter(time=='P0')%>%select(gene_id,translational_logFC)%>%
 left_join(metainfo%>%distinct(gene_id,gene_name))



# Example to display the mean of several genes
allTEchangedf%>%filter(up==1)%>%
	left_join(fcdf)%>%filter(abs(translational_logFC)>1)%>%.$gene_name%>%intersect(rownames(knnMap$M))%>%
	{plotKnnMap(knnMap,colMeans(knnMap$M[.,]))}
title('All TE Up genes, > 2-fold change ')

allTEchangedf%>%filter(down==1)%>%left_join(fcdf)%>%filter(abs(translational_logFC)>1)%>%.$gene_name%>%intersect(rownames(knnMap$M))%>%
	{plotKnnMap(knnMap,colMeans(knnMap$M[.,]))}
title('All TE Down genes, > 2fold Change')

allTEchangedf%>%filter(up==0,down==0)%>%.$gene_name%>%intersect(rownames(knnMap$M))%>%
	{plotKnnMap(knnMap,colMeans(knnMap$M[.,]))}
title('All TE stable genes')



