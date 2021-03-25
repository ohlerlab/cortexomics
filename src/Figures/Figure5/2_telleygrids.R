#src/R/Load_data/load_telley_etal.R
source('src/R/Rprofile.R')

hclustob <- readRDS(here('data/hclustob.rds'))
# Load the data
knnMap <- readRDS(here("ext_data/telley_sentdata/knnMap.rds"))


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

kngenes = knnMap$M%>%rownames
igene_grps = hclustob$cluster%>%split(names(.),.)
igene_grps = c(list(te_up=teupgenes,te_down=tedowngenes),igene_grps)
igene_grps = c(list(Satb2='Satb2',Pum2='Pum2',Flna='Flna',Tle4='Tle4',Nes='Nes',Bcl11b='Bcl11b'),igene_grps)


imap(igene_grps,function(igenes,grpname){
igenes = igene_grps[[grpname]]
message(igenes%>%is_in(kngenes)%>%table%>%print%>%capture.output%>%paste0(collapse='\n'))
igenes = intersect(kngenes,igenes)
#now plot
plotfile<- here(paste0('plots/','/Figures/Figure5/telley_timediff_grids/',grpname,'.pdf'))
# dir.create(dirname(plotfile))
pdf(plotfile)
plotKnnMap(knnMap,knnMap$M[igenes,])
dev.off()
message(normalizePath(plotfile))
})





# Example to display individual genes

plotKnnMap(knnMap,knnMap$M["Eomes",])
plotKnnMap(knnMap,knnMap$M["Dcx",])
