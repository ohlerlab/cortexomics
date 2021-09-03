library(readxl)
#which housekeepers to normalize with
HKGENES2USE <- c("5S rRNA", "18S rRNA")


################################################################################
################################################################################
if(!exists('safe_left_join'))base::source(here::here('src/R/Rprofile.R'))
if(!exists('iso_tx_countdata')) load(here('data/1_integrate_countdata.R'))
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
library(GenomicFeatures)


setwd(here())
cdsseq = extractTranscriptSeqs(x=fafile,cdsgrl)

################################################################################
######## Codon optimality scores based on codon frequencies, gene expression levels
################################################################################source(here("src/R/Rprofile.R"))
stagecolsdisplay <- c(E12.5 = "#214098", E14 = "#2AA9DF", E15.5 = "#F17E22", E17 = "#D14E28", P0 = "#ED3124")
stagecols <- stagecolsdisplay %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stageconv <- names(stagecols) %>% setNames(c("E13", "E145", "E16", "E175", "P0"))

codonfreqs <- memoise(oligonucleotideFrequency)(cdsseq, 3, step = 3)
rownames(codonfreqs) <- names(cdsseq)
overallcodonfreqs <- codonfreqs %>% colSums()

codon_is_optimal <- overallcodonfreqs %>%
    enframe("codon", "freq") %>%
    mutate(AA = as.character(translate(DNAStringSet(codon)))) %>%
    group_by(AA) %>%
    mutate(optimal = freq == max(freq)) %>%
    {
        setNames(.$optimal, .$codon)
    }
#
optimalcodons <- names(codon_is_optimal)[codon_is_optimal]
pc_optcodons <- ((codonfreqs[, optimalcodons] %>% rowSums()) / (codonfreqs %>% rowSums())) %>%
    setNames(rownames(codonfreqs)) %>%
    enframe("transcript_id", "pc_opt")
#
cdsexprvals <- iso_tx_countdata$abundance
#
cdsexprdf <- cdsexprvals%>%
    as.data.frame() %>%
    rownames_to_column('transcript_id') %>%
    gather(sample, abundance, -transcript_id) %>%
    separate(sample, c("time", "assay", "rep")) %>%
    group_by(transcript_id, time, assay) %>%
    summarise(abundance = mean(abundance)) %>%
    filter(assay == "total") %>%
    spread(time, abundance) %>%
    arrange(match(transcript_id, rownames(codonfreqs)))
#

times=colnames(cdsexprdf)[-c(1:2)]%>%setNames(.,.)
itime=times[1]
weighted_codon_usage <- lapply(times, function(itime) {
    # stopifnot(rownames(codonfreqs)==cdsexprdf$transcript_id)
    (codonfreqs[cdsexprdf$transcript_id, ] * cdsexprdf[[itime]]) %>%
        colSums() %>%
        {
            .
        } %>%
        enframe("codon", "weightedusage")
}) %>% bind_rows(.id = "time")
weighted_codon_usage%<>%mutate(AA = GENETIC_CODE[str_extract(codon, "[^\\-]+$")])
weighted_codon_usage%<>%group_by(time,AA)%>%mutate(w_cAI = weightedusage/sum(weightedusage))%>%ungroup

codon_usage <-  (codonfreqs[cdsexprdf$transcript_id, ] * 1) %>%
        colSums() %>%
        {
            .
        } %>%
        enframe("codon", "freq")
codon_usage%<>%mutate(AA = GENETIC_CODE[str_extract(codon, "[^\\-]+$")])
codon_usage%<>%group_by(AA)%>%mutate(cAI = freq/sum(freq))%>%ungroup
# test_that('freq and weighted usage are very similiar',{
#     txtplot(codonfreqs%>%colSums,weighted_codon_usage%>%filter(time=='P0')%>%.$weightedusage)
#     expect_true(cor(codonfreqs%>%colSums,weighted_codon_usage%>%filter(time=='P0')%>%.$weightedusage) > .95)

# })


################################################################################
######## Reading the tRNA data
################################################################################


#which reports to take data from
alltRNAreps <- Sys.glob(here("ext_data/tRNA_data/*Rep*/*")) %>%
    str_subset("/(80S|Poly|Total)[^/]+xlsx$")
stopifnot(alltRNAreps%>%length%>%`>`(3))
# data is pulled from a combination of report and 'test/'ctrl' status
datasources <- tibble(
    file = c(alltRNAreps, alltRNAreps[1]),
    case = c(rep("Test", length(alltRNAreps)), "Control")
)
datasources %<>% mutate(
    sample =
        ifelse(case == "Test",
            str_extract(basename(file), regex(".*(?= VS )")),
            str_extract(file, regex("(?<= VS ).*(?=.xlsx)"))
        )
)
# Now iterate over these
trnaexprdf <- map2(datasources$file, datasources$case, .f = function(file, case) {
    excel_sheets(file)
    hksheet <- readxl::read_xlsx(file, sheet = "Choose Housekeeping Genes")
    hks2use_rows <- which(hksheet[[1]][1:4] %in% HKGENES2USE)
    casecol <- if (case == "Test") 4 else 27
    # hknormfacts <- list(
    #     rep1 = mean(as.numeric(hksheet[[casecol]][hks2use_rows])),
    #     rep2 = mean(as.numeric(hksheet[[casecol + 1]][hks2use_rows]))
    # )
    # # calibfactrow <- which(hksheet[[1]]=="Calibration factor")
    # calibfactrow <- which(hksheet[[1]] == "PPC")
    # ppcnormfacts <- list(
    #     rep1 = mean(as.numeric(hksheet[[casecol]][calibfactrow])),
    #     rep2 = mean(as.numeric(hksheet[[casecol + 1]][calibfactrow]))
    # )
    exprsheet <- if (case == "Test") "Test Sample Data" else "Control Sample Data"
    #
    message(file)
    suppressMessages({
        exprsheet <- read_xlsx(file, sheet = exprsheet)
    })
    #
    exprdf <- exprsheet[c(-1), 1:4] %>%
        as.data.frame() %>%
        head(-2) %>%
        set_colnames(c("decoder", "well", "rep1", "rep2"))
    stopifnot('18S rRNA' %in% exprdf$decoder)
    stopifnot('Undetermined' %in% exprdf$rep1)
    stopifnot('Undetermined' %in% exprdf$rep1)
    suppressWarnings({
    exprdf$rep1 %<>% as.numeric
    })
    suppressWarnings({
        exprdf$rep2 %<>% as.numeric
    })
    exprdf
})

allqpcrsig <- trnaexprdf %>%
    setNames(datasources$sample) %>%
    bind_rows(.id = "sample")
#
codonfromdecoder <- function(x) x %>%
    str_replace("-\\d+$", "") %>%
    str_extract("")
allqpcrsig %<>% mutate(anticodon = decoder %>% str_replace("-\\d+$", ""))
allqpcrsig %<>% mutate(iscodon = str_detect(decoder, "\\w+-\\w+$"))
allqpcrsig %<>% mutate(iscodon = ifelse(str_detect(decoder, "Spike-In"),FALSE,iscodon))
allqpcrsig %<>% mutate(time = sample %>% str_extract("(?<= ).*$"))
allqpcrsig <- gather(allqpcrsig, rep, Ct, rep1:rep2)
allqpcrsig %<>% mutate(fraction = str_extract(sample, "\\w+"))
stopifnot(all(allqpcrsig$rep%in%c('rep1','rep2')))
#normalize to hkgenes
# allqpcrsig %<>% mutate(Ct = replace_na(Ct,Inf))
allqpcrsig <- allqpcrsig%>% group_by(sample, rep) %>% mutate(Ct = Ct,
    hk_expr = mean(Ct[decoder %in% HKGENES2USE]),
    dCt = Ct - hk_expr
)



#function to add codons
addcodon <- function(x) x %>%
        mutate(codon = ifelse(iscodon,str_extract(anticodon, "[^-]+$") %>% DNAStringSet() %>%
        reverseComplement() %>%
        as.character() %>% qs("S+")))
#remove mitochondrial ones, and spike ins, norm genes etc.
all_deco_sig <- allqpcrsig%>%filter(iscodon)%>%filter(!anticodon%>%str_detect('^mt'))
all_deco_sig <- all_deco_sig %>% filter(!anticodon %>% str_detect("iMet"))
all_deco_sig%<>%addcodon
#now add weighted codon usage
all_deco_sig%<>%select(-matches('weightedusage'))
all_deco_sig %<>% safe_left_join(allow_dup=F,weighted_codon_usage,by=c('codon','time')) 
#join up with freqs as well
all_deco_sig %<>% safe_left_join(overallcodonfreqs%>%enframe('codon','freq'))
all_deco_sig %<>% mutate(fraction = str_extract(sample, "\\w+"))
all_deco_sig%>%group_by(time,fraction,codon)%>%group_by()
logSumExp2 <- function(x) logSumExp(log(2)*x)/log(2)
 
decodtabledat <- all_deco_sig%>%
    filter(sample%>%str_detect('Total'))%>%select(sample,decoder,well,Ct,dCt)%>%
    mutate(sample=paste0(sample,'_',str_replace(rep,'rep','')))%>%
    mutate(sample = str_replace(sample,' ','_'))%>%
    ungroup%>%select(-rep)

decodtabledat%>%select(-dCt)%>%spread(sample,Ct)%>%write_tsv('tables/tRNA_decoder_data.tsv')
decodtabledat%>%select(-Ct)%>%spread(sample,dCt)%>%write_tsv('tables/tRNA_decoder_data_norm.tsv')


allcodsig_isomerge <- all_deco_sig %>%
    group_by(fraction, iscodon,time, sample, anticodon,rep,weightedusage) %>%
    # filter(anticodon=='Ala-CGC')%>%
    #note here we flip the dCt so the infinities don't effect the summation
    summarise(Ct = -logSumExp2(-Ct),dCt = -logSumExp2(-dCt),abundance= - dCt)%>%
    mutate(abundance=ifelse(abundance== -Inf,NA,abundance))
#get availability - the residuals of the linear fit between abundance and weighted usage
allcodsig_isomerge%<>%group_by(fraction,time)%>%nest%>%
    mutate(availability = map(data,~ {
        l=nrow(.);
        lm(data=.,abundance ~ weightedusage)$residuals%>%
        {.[as.character(seq_len(l))]}
    }))%>%
    unnest
#now get the mean across replicates
allcodsigmean_isomerge <- allcodsig_isomerge %>%
    group_by(fraction, time, iscodon, sample, anticodon, rep,weightedusage)%>%
    summarise_at(vars(one_of(c('Ct','dCt','abundance','availability'))),list(mean),na.rm=T)
allcodsigmean <- all_deco_sig%>%
    group_by(fraction, time, iscodon, sample, anticodon, rep,weightedusage)%>%
    summarise_at(vars(one_of(c('Ct','dCt','abundance','availability'))),list(mean),na.rm=T)
# add codon info
allcodsigmean_isomerge %<>% addcodon
allcodsigmean %<>% addcodon
#add AA
# allcodsigmean %<>% mutate(AA = GENETIC_CODE[str_extract(codon, "[^\\-]+$")] %>% qs("S+"))
allcodsigmean_isomerge %<>% mutate(AA = GENETIC_CODE[str_extract(codon, "[^\\-]+$")] %>% qs("S+"))
#
allcodsigmean_isomerge%>%filter(sample%>%str_detect('Total'))%>%select(sample,codon,availability,Ct,abundance,availability)%>%
    write_tsv('tables/tRNA_decoder_data.tsv')

 all_deco_sig%>%
    filter(fraction=='Total')%>%
    select(-hk_expr,-fraction,-weightedusage,-freq)%>%
    write_tsv('tables/isodecoder_data.tsv')

allcodsigmean_isomerge%>%group_by(fraction,time)%>%group_slice(1)%>%filter(!is.na(Ct))%>%select(AA,codon,abundance)%>%arrange(desc(abundance))
allcodsigmean_isomerge%>%group_by(fraction,time)%>%group_slice(1)%>%filter(!is.na(Ct))%>%select(AA,codon,abundance)%>%arrange(abundance)

if(FALSE){
    
# get normfacts
# normfacts <- allcodonsig%>%group_by(sample,rep)%>%summarise(normfact = mean(abundance[decoder%in%HKGENES2USE]))
# Now normalize
test_that("HTqPCR agrees with me",{

    testvals<-all_deco_sig%>%group_by(decoder)%>%group_slice(1)%>%
        filter(sample=='Total E16'|sample=='Total E13')%>%as.data.frame%>%
        select(decoder,sample,dCt,Ct,rep)%>%group_by(sample)%>%summarise(dCt = mean(dCt))
    testvals%>%    {print(.);.}%>%
        summarise(
            logfoldchange =  -diff(dCt),
            linear_foldchange = 2^logfoldchange
        )

    qpcrsamples <- readLines('ext_data/tRNA_data/Raw data/alldata.csv',2)%>%tail(1)%>%str_split(';')%>%.[[1]]%>%tail(-2)
    qpcrdata <- fread('ext_data/tRNA_data/Raw data/alldata.csv',skip=2)%>%set_colnames(c(colnames(.)[1:2],qpcrsamples))

    ctdatfiles<-qpcrdata%>%
        gather(sample,Ct,-`Transcript name`,-Well)%>%
        mutate(flag=NA,type=NA)%>%
        mutate(Ct = str_replace(Ct,',','.'))%>%
        mutate(Ct = as.numeric(Ct))%>%
        mutate(Ct = replace_na(Ct,NA))%>%
        split(.$sample)%>%
        map(.%>%select(`Transcript name`,Well,Ct,flag,type))%>%
        imap(function(x,nm){file=paste0('ext_data/tRNA_data/',nm,'.txt');write_tsv(x,file,col_names=F);file})%>%
        unlist
    warnings()

    # totdatfiles <- ctdatfiles%>%.[str_detect(.,'Total')]

    # ctobject <- readCtData(totdatfiles,column.info=list(flag=4, feature=1, type=5, position=2,Ct=3),header=F,n.features=192)
    # mygroups <- totdatfiles%>%names%>%str_replace('\\-\\d','')

    # d.norm <- normalizeCtData(ctobject, norm = "deltaCt",deltaCt.genes = HKGENES2USE)


    # stopifnot(`<`(testvals$dCt[1]- mean(exprs(d.norm)['Ala-AGC-1',c(1,2,5,6)][1:2]),0.001))
    # stopifnot(`<`(testvals$dCt[2]- mean(exprs(d.norm)['Ala-AGC-1',c(1,2,5,6)][3:4]),0.001))

})


allcodsig_isomerge%>%left_join(enframe('codon','freq'),by='codon')

test_that("These values look right, abundance correlates with usage on log scale before and after",{
    
    expect_true(allqpcrsig%>%filter(!iscodon)%>%.$decoder%>%unique%>%identical(c("U6", "5S rRNA", "18S rRNA", "Spike-In", "PPC", "GDC", "Blank")))
    expect_true(overallcodonfreqs%>%enframe('codon','freq')%>%safe_left_join(weighted_codon_usage%>%filter(time=='P0'))%>%{cor(.$freq,.$weightedusage)}%>%between(0.95,0.99))

    ctdatfiles%<>%mutate(isMT = str_detect(anticodon,'^mt'))

    ct_usage_correlation <- ctdatfiles %>% 
        mutate(isMT = str_detect(anticodon,'^mt'))%>%
        filter(iscodon,!isMT)%>%    
        filter(sample=='Poly E13')%>% mutate(codon = ifelse(iscodon,str_extract(anticodon, "[^-]+$") %>% DNAStringSet() %>%
        reverseComplement() %>%
        as.character() %>% qs("S+"))
    ,NA)%>%
    group_by(codon,decoder,sample)%>%
    # filter(rep=='rep2')%>%
    summarise(Ct = log2(mean(2^Ct)))%>%
    group_by(codon,sample)%>%
    summarise(Ct=-logSumExp(-Ct,na.rm=T),na.rm=T)%>%
    safe_left_join(overallcodonfreqs%>%enframe('codon','freq'),by='codon')%>%
    filter(is.finite(Ct))%>%
    safe_left_join(weighted_codon_usage%>%filter(time=='E13'),by=c('codon'))%>%
    {cor.test(.$Ct,.$weightedusage)}%>%tidy%>%.$estimate
    expect_true(between(ct_usage_correlation,-.7,-.5))

    allcodsigmean %>%
    mutate(as.character(translate(DNAStringSet(codon)))) %>%
    head(1) %>%
    {
        stopifnot(.$anticodon == "Ala-AGC" & (.$codon == "GCT"))
    }

    allcodsigmean_isomerge %>%
        mutate(as.character(translate(DNAStringSet(codon)))) %>%
        head(1) %>%
        {
            stopifnot(.$anticodon == "Ala-AGC" & (.$codon == "GCT"))
        }
    #
    glutRNAabund <- allcodsigmean_isomerge %>%
        filter(codon %>% str_detect("CTC")) %>%
        .$abundance
    valtRNAabund <- allcodsigmean_isomerge %>%
        filter(codon %>% str_detect("AAC")) %>%
        .$abundance
    stopifnot(valtRNAabund < glutRNAabund)

})

}
# https://github.com/dbgoodman/ecre_cds_analysis/blob/master/codonR/tAI.R
#

get_ws <- function(abundance, codons, # tRNA gene copy number
                   s = c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89) # sharing rules
) {
    ordered_codons <- names(GENETIC_CODE)
    tRNA_in <- abundance %>% setNames(codons)
    stopifnot(all(names(tRNA_in) %in% ordered_codons))
    stopifnot(!anyDuplicated(names(tRNA_in)))
    # Make sure we have them all
    tRNA <- tRNA_in[ordered_codons] %>% setNames(ordered_codons)
    stopifnot(names(tRNA) == ordered_codons)
    tRNA <- 2^tRNA # put on a linear scale
    tRNA[is.na(tRNA)] <- 0
    p <- 1 - s
    #
    # initialise w vector
    W <- NULL # don't confuse w (lowercase) and W (uppercase)
    #
    # obtain absolute adaptiveness values (Ws)
    for (i in seq(1, 61, by = 4)) {
          W <- c(
              W,
              p[1] * tRNA[i] + p[5] * tRNA[i + 1], # INN -> NNT, NNC, NNA
              p[2] * tRNA[i + 1] + p[6] * tRNA[i], # GNN -> NNT, NNC
              p[3] * tRNA[i + 2] + p[7] * tRNA[i], # TNN -> NNA, NNG
              p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2]
          )
      } # CNN -> NNG
    #
    # check methionine
    # W[36] = p[4]*tRNA[36]
    #
    # get rid of stop codons (11, 12, 15) and methionine (36)
    # W = W[-c(11,12,15,36)]
    W <- W[-c(11, 12, 15)]
    #
    # get ws
    w <- W / max(W, na.rm = T)
    #
    if (sum(w == 0) > 0) {
        ws <- w[w != 0] # zero-less ws
        gm <- exp(sum(log(ws)) / length(ws)) # geometric mean
        w[w == 0] <- gm # substitute 0-ws by gm
    }
    # p
    return(log2(w[codons]) %>% setNames(codons))
}




weighted_aa_usage <- lapply(times, function(itime) {
    aas <- DNAStringSet(colnames(codonfreqs))%>%translate%>%as.character
    split(1:ncol(codonfreqs),aas)%>%lapply(function(cols)codonfreqs[,cols,drop=F]%>%rowSums)%>%bind_cols%>%
    as.matrix%>%    
     multiply_by( cdsexprdf[[itime]])%>%
     colSums%>%{./sum(.)}%>%
     enframe("AA", "aa_weightedusage")

}) %>% bind_rows(.id = "time")


tRNA_abchange <- allcodsigmean_isomerge %>%
    group_by(fraction, codon) %>%
    filter(!is.na(abundance)) %>%
    filter(is.finite(abundance), !is.na(abundance)) %>%
    nest() %>%
    mutate(abundancechange = map_dbl(data, ~ {
        lm(data = ., abundance ~ seq_along(time))$coef[2]
    }))

usage_v_abundance_df <- weighted_codon_usage %>%
    filter(time == "E13") %>%
    ungroup() %>%
    select(-time) %>%
    left_join(tRNA_abchange %>% filter(is.finite(abundancechange)) %>% select(fraction, codon, abundancechange)) %>%
    filter(!is.na(abundancechange), !is.na(weightedusage))

wus_tab_cors <- usage_v_abundance_df %>%
    filter(!codon %in% c("TAG", "TAA", "TGA")) %>%
    group_by(fraction) %>%
    nest() %>%
    mutate(
        cor = map_dbl(data, ~ cor(.$abundancechange, .$weightedusage)),
        pval = map_dbl(data, ~ cor.test(.$abundancechange, .$weightedusage)$p.value)
    )
  
################################################################################
########Model differences between fractions
################################################################################
    

allcodsigmean_isomerge$ntime<-allcodsigmean_isomerge$time%>%as_factor%>%as.numeric
sigcol=sym('abundance')
tRNAlmlist <- lapply(list(sym('abundance'), sym('availability')),function(sigcol){
# tRNAlmlist <- lapply(list(sym('abundance')),function(sigcol){
    #prep data for the lm
    tRNAlmdf<-allcodsigmean_isomerge%>%
        ungroup%>%
        mutate(fraction=factor(fraction,c('Total','Poly','80S')))%>%
        # filter(codon%in%unique(codon)[1:2])%>%
        mutate(!!sigcol := !!sigcol%>%ifelse(is.nan(.),-Inf,.))%>%
        group_by(sample,time,ntime,fraction)%>%
        # filter(!!sigcol%>%{!is.finite(.)})
        mutate(!!sigcol := !!sigcol%>%{ifelse(!is.finite(.),min(vals(.))-3,.)})%>%
        group_by(sample,time,ntime,fraction,codon)%>%
        summarise(sigcol := matrixStats::logSumExp(!!sigcol))
    #now do the lm
    tRNAlm <- lm(data=tRNAlmdf,sigcol ~ time*fraction*codon)
    #parse teh coefficients
    # tRNAeffectdf<-tRNAlm$coef%>%enframe%>%filter(!name%>%str_detect('ntercept'))%>%
    #     # mutate(timedep=str_detect(name,'ntime'))%>%
    #     mutate(time=str_extract(name,'(?<=time)(\\w+)'))%>%
    #     mutate(fraction=str_extract(name,'(?<=fraction)(\\w+)'))%>%
    #     mutate(codon=str_extract(name,'(?<=codon)(\\w+)'))%>%
    tRNAlmdf_vm = tRNAlmdf
    tRNAlmdf_vm$fraction%<>%factor(levels=c('80S','Poly','Total'))
    tRNAlm_vmono <- lm(data=tRNAlmdf_vm,sigcol ~ time*fraction*codon)
    names(tRNAlm_vmono$coef)%<>%str_replace('Total','Total_vmono')
    names(tRNAlm_vmono$coef)%<>%str_replace('Poly','Poly_vmono')
    coefs = c(tRNAlm$coef,tRNAlm_vmono$coef%>%.[names(.)%>%str_detect('vmono')])
    tRNAeffectdf<-coefs%>%enframe%>%filter(!name%>%str_detect('ntercept'))%>%
        # mutate(timedep=str_detect(name,'ntime'))%>%
        mutate(time=str_extract(name,'(?<=time)(\\w+)'))%>%
        mutate(fraction=str_extract(name,'(?<=fraction)(\\w+)'))%>%
        mutate(codon=str_extract(name,'(?<=codon)(\\w+)'))%>%
        select(name,time,fraction,codon,value)
    tRNAeffectdf$time[((is.na(tRNAeffectdf$time))&(!is.na(tRNAeffectdf$fraction))&(!is.na(tRNAeffectdf$codon)))]<-'E13' 
    tRNAeffectdf%>%filter(!is.na(time),!is.na(fraction),!is.na(codon),!is.na(value))%>%select(-name)
})
#
tRNAenrichdf <- safe_left_join(tRNAlmlist[[1]],tRNAlmlist[[2]],by=c('time','fraction','codon'))%>%
    dplyr::rename('abundance_enrich':=value.x,'availability_enrich':=value.y)


# pdf <- cairo_pdf
# fractions <- c("Poly", "Total")
# for (ifraction in fractions) {
#     plotfile <- str_interp("plots/figures/figure2/trna_codons/trna_abchange_exprusage_${ifraction}.pdf")
#     pdf(plotfile)
#     usage_v_abundance_df %>%
#         filter(!codon %in% c("TAG", "TAA", "TGA")) %>%
#         filter(fraction == ifraction) %T>%
#         {
#             message(nrow(.))
#         } %>%
#         nest() %>%
#         mutate(labl = paste0(
#             "r = ", map_dbl(data, ~ cor(.$abundancechange, .$weightedusage)) %>% round(3), "\n",
#             "p = ", map_dbl(data, ~ cor.test(.$abundancechange, .$weightedusage)$p.value %>% round(3))
#         )) %>%
#         unnest() %>%
#         {
#             print(
#                 qplot(data = ., x = .$abundancechange, y = .$weightedusage, label = codon, geom = "blank") +
#                     scale_x_continuous("tRNA Abundance abundance Slope E13 - P0") +
#                     scale_y_continuous("RiboExpression Weighted E13 ") +
#                     geom_text(data = distinct(., labl), hjust = 0, vjust = 1, x = -Inf, y = Inf, aes(label = labl)) +
#                     geom_smooth(method = "lm") +
#                     # facet_grid(tRNA_time~.)+
#                     geom_text() +
#                     geom_text(data = wus_tab_cors, aes(x = 0, y = 0, label = paste0("r = ", round(cor, 2)))) +
#                     theme_bw()
#             )
#         }
#     dev.off()
#     message(normalizePath(plotfile))
# }


# # are and global frequency correlated? - AA level?
# pdf("plots/figures/figure2/trna_codons/AAusage_vs_summed_tRNAab.pdf")
# allcodsigmean_isomerge %>%
#     left_join(enframe(codonfreqs %>% colSums(), "codon", "freq")) %>%
#     group_by(fraction, time, AA) %>%
#     # slice(which.max(abundance))%>%
#     # slice(which.max(freq))%>%
#     summarise(abundance = log2(sum(2^abundance)), freq = sum(freq)) %>%
#     # summarise(abundance = mean(abundance),freq=sum(freq))%>%
#     group_by(fraction, time) %>%
#     nest() %>%
#     mutate(labl = paste0(
#         "r = ", map_dbl(data, ~ cor(use = "complete", .$freq, .$abundance)) %>% round(3), "\n",
#         "p = ", map_dbl(data, ~ cor.test(use = "complete", .$freq, .$abundance)$p.value %>% round(3))
#     )) %>%
#     unnest() %>%
#     {
#         ggplot(., aes(x = abundance, y = freq)) + geom_point() + facet_grid(time ~ fraction) +
#             theme_bw() + geom_text(data = distinct(., time, fraction, labl), hjust = 0, vjust = 1, x = -Inf, y = Inf, aes(label = labl))
#     } +
#     ggtitle("Amino acid usage frequency vs summed tRNA abundance")
# dev.off()
# normalizePath("plots/figures/figure2/trna_codons/AAusage_vs_summed_tRNAab.pdf")

# # Are occupancy and global frequency correlated? - AA level??
# pdf("plots/figures/figure2/trna_codons/AAweightedusage_vs_summed_tRNAab.pdf")
# allcodsigmean_isomerge %>%
#     left_join(weighted_codon_usage) %>%
#     group_by(fraction, time, AA) %>%
#     summarise(abundance = log2(sum(2^abundance)), weightedusage = sum(weightedusage)) %>%
#     filter(is.finite(abundance)) %>%
#     group_by(fraction, time) %>%
#     nest() %>%
#     mutate(labl = paste0(
#         "r = ", map_dbl(data, ~ cor(use = "complete", .$weightedusage, .$abundance)) %>% round(3), "\n",
#         "p = ", map_dbl(data, ~ cor.test(use = "complete", .$weightedusage, .$abundance)$p.value %>% round(3))
#     )) %>%
#     unnest() %>%
#     {
#         ggplot(., aes(x = abundance, y = weightedusage)) + geom_point() + facet_grid(time ~ fraction) +
#             theme_bw() + geom_text(data = distinct(., time, fraction, labl), hjust = 0, vjust = 1, x = -Inf, y = Inf, aes(label = labl))
#     } +
#     ggtitle("Amino acid expr weighted usage frequency vs summed tRNA abundance")
# dev.off()
# normalizePath("plots/figures/figure2/trna_codons/AAweightedusage_vs_summed_tRNAab.pdf")


# # Are occupancy and global abundance correlated?
# pdf("plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf")
# allcodsigmean_isomerge %>%
#     mutate(abundance = abundance) %>%
#     left_join(weighted_codon_usage) %>%
#     filter(is.finite(abundance)) %>%
#     group_by(fraction, time) %>%
#     nest() %>%
#     mutate(labl = paste0(
#         "r = ", map_dbl(data, ~ cor(use = "complete", .$weightedusage, .$abundance)) %>% round(3), "\n",
#         "p = ", map_dbl(data, ~ cor.test(use = "complete", .$weightedusage, .$abundance)$p.value %>% round(3))
#     )) %>%
#     unnest() %>%
#     {
#         ggplot(., aes(x = abundance, y = weightedusage)) + geom_point() + facet_grid(time ~ fraction) +
#             theme_bw() + geom_text(data = distinct(., time, fraction, labl), hjust = 0, vjust = 1, x = -Inf, y = Inf, aes(label = labl))
#     } +
#     ggtitle("Codon usage frequency vs summed tRNA abundance")
# dev.off()
# normalizePath("plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf")

