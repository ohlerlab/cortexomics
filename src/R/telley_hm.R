

# Load the data
knnMap <- readRDS(here("ext_data/telley_sentdata/knnMap.rds"))
#knnMap$M           : is a matrix 16441x2061 containing the RPM values of all genes (16441) for all cells (2061)
#knnMap$grid        : is a precomputed grid of size 512x512 where genes expression are infered using 15 nearest neighboring cells.
#knnMap$grid$(xi,yi): are the coordinates of the elements in the grid
#knnMap$grid$knn    : are precomputed indices of the 15 nearest neighbooring cells of each pixel of the grid





# knnMap: the knnMap object containing the data
#      z: a numeric vector of expression of the gene to display for all cells (should be one row of knnMap$M)
plotKnnMap <- function(knnMap,z=knnMap$M[,"Dcx"]) {
  # The color palette we use
  pal10 <- c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")
  
 # for each pixel in the grid, compute the average expression of the 15 nearest neighbors
 knnMap$grid$knn.avg <- rowMeans(array(z[knnMap$grid$knn],dim(knnMap$grid$knn)))
 
 # then, break the values in 11 levels to match the size of our color palette
 knnMap$grid$knn.avg.lev <- cut(knnMap$grid$knn.avg,breaks = seq(-0.1,max(knnMap$grid$knn.avg,1),length.out = 11))

 # transform the flat grid format into a raster.image to display
 m <- as.raster(matrix(NA_character_,512,512))
 m[cbind(knnMap$grid$yi,knnMap$grid$xi)] <- pal10[as.integer(knnMap$grid$knn.avg.lev) - min(as.integer(knnMap$grid$knn.avg.lev)) + 1]
 m <- m[rev(seq(nrow(m))),]

 # display the image
 plot(NA,NA,xlim=c(0,1.15),ylim=c(0,1),ylab="Differentiation score",xlab="Birth score",axes = FALSE,asp=1)
 axis(1,c(0,1),c("E12","E15"))
 axis(2,c(0,1),c("AP","N4d"))
 rasterImage(m,0,0,1,1)
 
 # display colorlegend
 P <- as.raster(matrix(rev(pal10),ncol=1))
 rasterImage(P,1.05,0.7,1.1,1)
 rect(1.05,0.7,1.1,1)
 text(1.11,1,"high",c(0,1))
 text(1.11,0.7,"low",c(0,0))
}



# Example to display individual genes
plotKnnMap(knnMap,knnMap$M["Pax6",])
plotKnnMap(knnMap,knnMap$M["Eomes",])
plotKnnMap(knnMap,knnMap$M["Dcx",])


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



