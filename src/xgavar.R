





#make fake dataset for this.


nreads = 1e5
# nseqcols=4
randomseqmat <- function(nreads){
	seqmat <- matrix(0,nrow=nreads,ncol=4)
	seqmat[matrix(c(1:nreads,sample(1:4,nreads,replace=T)),ncol=2)]<-1	
	seqmat
}

seqcols = replicate(4,randomseqmat(nreads),simp=F)%>%
	lapply(as.data.frame)%>%
	bind_cols

codons = GENETIC_CODE%>%names
codon_dwells = (10+(1:64))%>%setNames(codons)

asitecodons = sample(codons,nreads,replace=T,prob=codon_dwells)
asitecodons%>%table%>%sort

othercodons = matrix(sample(codons,nreads*7,replace=T),ncol=7)
othercodons[,4] <- asitecodons


seqs2shift = seqcols%>%unique%>%sample_n(floor(nrow(.)/10))

#get original sequences to shift
rds2shift <- row.match(seqcols,seqs2shift)%>%is.na%>%`!`

centercodons = othercodons[,c(3,4,5)]
centercodons[rds2shift] = othercodons[rds2shift,c(2,3,4)]

v=centercodons[,3]

entropy <- function(v) {
	p = table(v)
	p = p/sum(p)
	p = p[p!=0]
	p = sum(p*log(p))
	- p
}

centercodons[!rds2shift,]%>%apply(2,entropy)%>%txtplot
centercodons[rds2shift,]%>%apply(2,entropy)%>%txtplot

install.packages("drat", repos="https://cran.rstudio.com")
drat:::addRepo("dmlc")
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")

require(xgboost)

require(magrittr)
x <- model.matrix(Species ~ .^2, iris)[,-1]
bstSparse <- xgboost(data = train$data, label = train$label, max.depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")
dtrain <- xgb.DMatrix(scale(x), label = as.numeric(iris$Species) - 1)
dtrain <- xgb.DMatrix(scale(x), label = as.numeric(iris$Species) - 1)
param <- list(booster = "gblinear", objective = "multi:softprob", num_class = 3,
lambda = 0.0003, alpha = 0.0003, nthread = 2)
# For the default linear updater 'shotgun' it sometimes is helpful
# to use smaller eta to reduce instability
bst <- xgb.train(param, dtrain, list(tr=dtrain), nrounds = 70, eta = 0.5,
callbacks = list(cb.gblinear.history()),obj = function(...){browser()})


# sample(codons,10000,replace=T,prob=codon_dwells)%>%entropy
# sample(codons,10000,replace=T,prob=c(rep(1,63),20)%>%{./sum(.)})%>%entropy
