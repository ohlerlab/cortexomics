library(GenomicAlignments)


rnabam <- 'star/data/E13_total_1/E13_total_1.bam'
annofile<-'my_gencode.vM12.annotation.gtf'


anno<-rtracklayer::import(annofile)

genes <- subset(anno,type=='gene')

counts <-  fread('feature_counts/all_feature_counts')

topgenes<-counts[,c('E13_total_1','feature_id')]%>%arrange(desc(E13_total_1))%>%head(1000)%>%.$feature_id


greads<-rtracklayer::import(rnabam,which=genes%>%subset(gene_id%in%topgenes[1:10]))


#and the magnitudes of the FFT at the closest frequency to 0.33
take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
	require(multitaper)
     if(length(x)<25){
          remain<-50-length(x)
          x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
     }
     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
          
     closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
     
     freq_max_3nt<-resSpec1$freq[closestfreqind]
     Fmax_3nt<-resSpec1$mtm$Ftest[closestfreqind]
     spect_3nt<-resSpec1$spec[closestfreqind]
     return(c(Fmax_3nt,spect_3nt))
     
}
#this takes in a numeric vector of psite values (i.e. the number of psites at a given bp),
#and possibly the number of slepians and bandwidth for the multitaper test, and returnst the
#pvalue of the FFT.
ftestvect<-function(psit,k=24,bw=12){
	psit <- as.vector(psit)
	if(sum(psit!=0) < 3) return(NA)
	if(sum(psit) < 10) return(NA)
		
	slepians_values<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
	vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = slepians_values)
	pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
	return(c(pval))
}


exons<-subset(anno,type=='exon')%>%split(.,.$transcript_id)


library(GenomicFeatures)

transcriptvects<-mapToTranscripts(resize(as(greads,'GRanges'),1),exons)%>%coverage%>%.[sum(.)>0]
transcriptvects%<>%lapply(as.vector)

transcriptvects%>%lapply(ftestvect)