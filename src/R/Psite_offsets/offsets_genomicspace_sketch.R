#simulatae some data

#frame preference vector
#read length distribution


#see if the test bam can be corrected by my script

toptrs <- topcds$transcript_id%>%unique%>%head(10)#top 100 genes
#fakeread<-topcdsreads[qwidth(topcdsreads)==30][1]
#fakeread%>%saveRDS(here('data/fakeread.rds'))
fakeread<-readRDS(here('data/fakeread.rds'))
fakeread%<>%GRanges

#get their cds
cdstosim<-topcdsmap[toptrs]

i=1
p1=seq(start(cdstosim[i])+0,end(cdstosim[i]),by=3)
p2=seq(start(cdstosim[i])+1,end(cdstosim[i]),by=3)
p3=seq(start(cdstosim[i])+2,end(cdstosim[i]),by=3)
codonstrengths = runif(length(p1),10,100)
pstrengths = c(0.6,0.3,0.1)

psites <- floor(unlist(t((codonstrengths)%*%t(pstrengths)))%>%as.vector)
shifts <- rep(0:(length(psites)-1),psites)

fakereads <- rep(fakeread,length(shifts))
start(fakereads)<-start(cdstosim[i])
fakereads%<>%shift(shifts)

#Now we can shift these as we 
lenfreqs <- c(`26` = 0.191948671893167, `27` = 0.297173585890137,
`28` = 0.243786738810325, `29` = 0.182682147043328, `30` = 0.084408856363043)
gensizes <- function(n) sample(names(lenfreqs),n,prob=lenfreqs,rep=TRUE)%>%as.numeric

fakereads%>%resize(gensizes(length(fakereads)))

#

#make fake reads shifted a bit

#get some fake reads of each size

#

#I want a functiont that outputs a phase Rle given some cds
mycds <- GRanges(c('1:1-3','1:7-10','14-16'))

