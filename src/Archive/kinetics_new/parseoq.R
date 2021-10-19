    library(R.utils)

    library(GenomicRanges)
    library(tidyverse)
    args = R.utils::commandArgs(trailingOnly=TRUE,asValues=TRUE,defaults =list(
      # oqfile='/fast/AG_Ohler/florian/OrfQuant_induced_neurons/ribo_0h_A.bam_for_SaTAnn_final_ORFquant_results',
      oqfile='/fast/AG_Ohler/florian/OrfQuant_induced_neurons/ribo_48h_A.bam_for_SaTAnn_final_ORFquant_results',
      outfile='oqout.gtf'
    ))
    if(interactive())args=args[-1]
    # Turn arguments into R variables
    keys <- attachLocally(args)
    cat("Command-line arguments attached to global environment:\n");
    print(keys);
    str(mget(keys, envir=globalenv()))

    #load the file
    oqobject <- get(load(oqfile))

    #get the transcript space granges object with the orf data
    oqmcols = mcols(oqobject$ORFs_tx)
    #
    nonlistcols = oqmcols%>%
        sapply(mode)%>%
        `!=`('S4')

    mcols(oqobject$ORFs_gen ) <- oqobject$ORFs_tx%>%
        mcols%>%.[oqobject$ORFs_gen%>%names,]

    stopifnot(setequal(oqobject$ORFs_tx$ORF_id_tr,
        oqobject$ORFs_gen$ORF_id_tr))

    #now output as a gtf
    oqobject$ORFs_gen[,nonlistcols]%>%{rtracklayer::export(.,outfile)}