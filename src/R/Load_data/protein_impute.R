# devtools::install_github("const-ae/proDD")
library(proDD)

# My package
library(proDD)

# Packages for plotting
library(ggplot2)
library(pheatmap)
library(viridisLite)
library(tidyverse,magrittr)
set.seed(1)

# The experimental_design is a vector that assignes each sample to one condition
experimental_design <- c("A", "A", "A", "B", "B", "B")
experimental_design <- matchedms_mat%>%colnames%>%str_extract('[^_]+')

# The generate_synthetic_data can be customized a lot, but here we will 
# use it in its most basic form
#syn_data <- generate_synthetic_data(n_rows=1234, experimental_design=experimental_design,
#                                   frac_changed = 0.1)
#syn_data
# The data matrix, where non-observed values are coded as NA
#X <- syn_data$X
X <- matchedms_mat[,TRUE]%>%
	{.[apply(.,1,is.na)%>%colMeans%>%`!=`(1),]}%>%
	set_rownames(1:nrow(.))

# Median normalization to remove sample effects 
X <- median_normalization(X)

# The columns are the samples and each row is a protein
params <-  fit_parameters(X[TRUE,TRUE],experimental_design,n_subsample=2000,verbose=TRUE)



posteriors <- sample_protein_means(X, params, verbose=TRUE)


save.image('proDD_fit.RData')



rowmeans <- posteriors$A%>%rowMeans
rowcoefvars<-posteriors$A%>%apply(1,sd)%>%divide_by(rowmeans)

qplot(rowmeans,rowcoefvars)+geom_smooth()


weirdinds <- which(( rowmeans < 26 ) &(rowcoefvars>0.05))
X[weirdinds,]%>%head
posteriors$A[weirdinds,]%>%head%>%rowMeans
posteriors$A[sample(weirdinds,1),]%>%hist(50)


#Now, for limma we can simply replace with three imputed values, and for the kinetic moddeling we can directly use our
#prot values.



