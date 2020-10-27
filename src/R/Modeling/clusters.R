library(tidyverse)
library(magrittr)
 dt <- 'allCounts.txt'
dt%>%read_tsv(col_names=FALSE)->dt



add_noise<-function(.){
	for (i in 5:11){
		.[[i]]= .[[i]]+((.[[i]]/100)*rnorm(length(.[[1]])))
	};.
}

%>%

getwd()
dt%>%
	list%>%rep(100)%>%
	map(.%>%mutate(X2=X2+1e6,X3=X3+1e6)%>%add_noise)%>%
	bind_rows%>%
	write_tsv('biggercounts.txt')

library(bnlearn)
data(learning.test)
data(gaussian.test)
learning.test%>%str

dag = model2network("[A][B]")

bn.fit(dag,gaussian.test)