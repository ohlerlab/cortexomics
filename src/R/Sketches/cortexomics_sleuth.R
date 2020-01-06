library(parallel)
library(stringr)
library(tidyverse)
library(sleuth)

kallistodir = "~/projects/cortexomics/data/kallisto_out/"

sample_id = dir(kallistodir)
metacols = c('time','assay','replicate')
s2c <-  
	data.frame(sample = sample_id )%>%
	mutate(path = file.path(kallistodir,sample_id))%>%
	mutate(tmp = sample_id%>%str_replace('_28_30','-28-30'))%>%
	separate(tmp,metacols,sep='_')

so <- sleuth_prep(s2c , extra_bootstrap_summary=TRUE, num_cores = 4)



h5files <- Sys.glob(file.path(kallistodir,'*/*tsv'))


h5file <- "~/projects/cortexomics/kallisto_out/E13_total_1/abundance.h5"
file.exists(h5file)

h5files <- Sys.glob(file.path(kallistodir,'*/*h5'))

tryread <- h5files%>%mclapply(safely(read_kallisto_h5))
