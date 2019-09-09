################################################################################
########This script takes in a matrix of ms values and then outputs posterior
########Values 
################################################################################
	
# devtools::install_github("const-ae/proDD")
library(proDD)

# Packages for plotting
library(ggplot2)
library(pheatmap)
library(viridisLite)
library(tidyverse,magrittr)
set.seed(1)

# The experimental_design is a vector that assignes each sample to one condition
experimental_design <- matchedms_mat%>%colnames%>%str_extract('[^_]+')

# The generate_synthetic_data can be customized a lot, but here we will 
# use it in its most basic form
#syn_data <- generate_synthetic_data(n_rows=1234, experimental_design=experimental_design,
#                                   frac_changed = 0.1)
#syn_data
# The data matrix, where non-observed values are coded as NA
#X <- syn_data$X
X <- matchedms_mat[,TRUE]%>%
	{.[apply(.,1,is.na)%>%colMeans%>%`!=`(1),]}
ms_id2row <- data_frame(ms_id=rownames(X),id=1:nrow(X))
X%<>%set_rownames(1:nrow(.))

# Median normalization to remove sample effects 
X <- median_normalization(X)

# The columns are the samples and each row is a protein
params <-  fit_parameters(X[TRUE,TRUE],experimental_design,n_subsample=2000,verbose=TRUE)



posteriors <- sample_protein_means(X, params, verbose=TRUE)

proddtest <- proDD::test_diff(posteriors[['E13']],posteriors[['P0']])

proddtest$adj_pval%>%`<`(0.05)%>%table

# save.image('proDD_fit.RData')

# load('proDD_fit.RData')




rowmeans <- posteriors$A%>%rowMeans
rowcoefvars<-posteriors$A%>%apply(1,sd)%>%divide_by(rowmeans)

qplot(rowmeans,rowcoefvars)+geom_smooth()


weirdinds <- which(( rowmeans < 26 ) &(rowcoefvars>0.05))
X[weirdinds,]%>%head
posteriors$A[weirdinds,]%>%head%>%rowMeans
posteriors$A[sample(weirdinds,1),]%>%hist(50)


#Now, for limma we can simply replace with three imputed values, and for the kinetic moddeling we can directly use our
#prot values.
library(testthat)

posteriorsumsplit <- posteriors%>%map(~ cbind(data_frame( 
	mean = apply(.,1,mean),
	var = apply(.,1,var),
	ymax=apply(.,1,quantile,0.975),
	ymin=apply(.,1,quantile,0.125),
)))

posteriorsum <- posteriorsumsplit%>%as.list%>%
	bind_rows(.id='time')%>%
	#mutate(ymin=mean-(sd*1.96),ymax=mean+(sd*1.96))
	identity
posteriorsum%<>% mutate(precision = 1 / var)
posteriorsum%<>%mutate(ms_id=rep(ms_id2row$ms_id,5))
posteriorsum%>%head
posteriorsum%<>%mutate_at(vars(mean,sd,ymin,ymax),~ 2^.)
posteriorsum%>%head

test_that("Imputation Looks reasonable, not just ",{
	#first let's find a gene with NAs at the first tp and upward slope
	allnatp1 <- apply(matchedms_mat[,c(1:3)],1,is.na)%>%colMeans%>%`==`(1)
	notallnatp5 <- apply(matchedms_mat[,c(13:15)],1,is.na)%>%colMeans%>%`==`(0)

	slopes <- matchedms_mat[,]%>%as.data.frame%>%rownames_to_column('id')%>%mutate(id=as.numeric(as.factor(id)))%>%gather(dataset,signal,-id)%>%mutate(tp=as.numeric(as.factor(str_extract(dataset,'[^_]+'))))%>%
		split(.,.$id)%>%
		map_dbl(~lm(.$signal ~ .$tp)$coef[2])

	testrows <- which(allnatp1 & notallnatp5 & (slopes > quantile(slopes,0.9,na.rm=TRUE)))
	testrow<-testrows%>%sample(1)

	notallnatp1 <- apply(matchedms_mat[,c(1:3)],1,is.na)%>%colMeans%>%`==`(0)
	allnatp5 <- apply(matchedms_mat[,c(13:15)],1,is.na)%>%colMeans%>%`==`(1)

	testrowsneg <- which(notallnatp1 & allnatp5 & (slopes < quantile(slopes,0.1,na.rm=TRUE)))
	testrowneg<-testrowsneg%>%sample(1)


})

if(!exists('plotting') plotting = TRUE)

if(plotting){


	plotfile=here('plots/protein_impute/testimputtraj.pdf')
	plotfile%>%dirname%>%dir.create
	pdf(plotfile)
	p=matched_ms%>%filter(ms_id==names(testrow))%>%
			{
				ggplot(data=.,)+
				geom_linerange(
					data=posteriorsum%>%filter(ms_id==names(testrow)),
					aes(ymin=ymin,ymax=ymax,x=as.numeric(as.factor(time))))+
				theme_bw()+
				geom_point(data=.,aes(x=as.numeric(as.factor(time)),y=signal))
		}
	p
	dev.off()
	message(normalizePath(plotfile))

	plotfile=here('plots/protein_impute/testimputtraj_neg.pdf')
	plotfile%>%dirname%>%dir.create
	pdf(plotfile)
	for(ind in sample(seq_along(testrowsneg),20)){
	testrowneg<-testrowsneg%>%.[ind]
	p=matched_ms%>%filter(ms_id==names(testrowneg))%>%
			{
				ggplot(data=.,)+
				geom_linerange(
					data=posteriorsum%>%filter(ms_id==names(testrowneg)),
					aes(ymin=ymin,ymax=ymax,x=as.numeric(as.factor(time))))+
				theme_bw()+
				geom_point(data=.,aes(x=as.numeric(as.factor(time)),y=signal))+
				# scale_y_log10()+
				ggtitle(matched_ms%>%filter(ms_id==names(testrowneg))%>%.$gene_name%>%.[1])
		}
	print(p)
	}
	dev.off()
	message(normalizePath(plotfile))


	plotfile=here('plots/protein_impute/testimputtraj_random.pdf')
	plotfile%>%dirname%>%dir.create
	pdf(plotfile)
	for(ms_id_i in sample(matched_ms$ms_id,20)){
	p=matched_ms%>%filter(ms_id==ms_id_i)%>%
			{
				ggplot(data=.,)+
				geom_linerange(
					data=posteriorsum%>%filter(ms_id==ms_id_i),
					aes(ymin=ymin,ymax=ymax,x=as.numeric(as.factor(time))))+
				theme_bw()+
				geom_point(data=.,aes(x=as.numeric(as.factor(time)),y=signal))+
				# scale_y_log10()+
				ggtitle(matched_ms%>%filter(ms_id==ms_id_i)%>%.$gene_name%>%.[1])
		}
	print(p)
	}
	dev.off()
	message(normalizePath(plotfile))

}
save(posteriors,posteriorsum,X,params,file='data/proDD.data')

#' Okay so, indeed this method doesn't take time point variation into account, so we end up with strange results

