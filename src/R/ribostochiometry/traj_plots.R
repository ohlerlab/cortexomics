library(here)
library(assertthat)
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
vlookup <- function(query,dicttable,key,vals){
  dict = dicttable%>%ungroup%>%distinct_(key,vals)
  stopifnot(!anyDuplicated(dict))
  data_frame(tmp=query)%>%
    left_join(dict,by=c('tmp'=key))%>%
    .[[vals]]
}

str_split_fast = function(x,sep=';') x %>% {str_split(.,sep,n = str_count(.,';')%>%max%>%add(1))}
#' setup, eval = TRUE
sep_element_in<-function(colonlist,ridssplit,sep=';'){
  assert_that(is.character(colonlist))
  assert_that(is.character(ridssplit))
  values<-colonlist%>%str_split_fast
  
  inds <- rep(seq_along(colonlist),lengths(values))
  
  values<-flatten_chr(values)
  
  data_frame(inds,values)%>%
    mutate(match = values %in% ridssplit)%>%
    group_by(inds)%>%
    summarize(match=any(match))%>%
    .$match
  
}



#load fractionated ms data
ms_tall <- map(c('./pipeline/ms_tables/ms_LFQ_total_ms_tall.tsv',
       './pipeline/ms_tables/ms_LFQ_poly_ms_tall.tsv',
       './pipeline/ms_tables/ms_LFQ_80S_ms_tall.tsv'),
     fread)
ms_tall <- ms_tall%>%bind_rows
message('Labeling protein IDs with annotations')
pids = ms_tall%>%ungroup%>%distinct(Protein_IDs)

stopifnot(ms_tall$fraction%>%unique%>%length == 3)


#do MEDIAN NORMALIZATION on our signal
ms_tall%<>%group_by(time,fraction,replicate)%>%mutate(signal = signal / median(signal,na.rm=T))

#define our catagories
message('looking up protein ID annotations')
#'get info on the ribosomal subu,nts from mats table
rids <- read_tsv(here('./ext_data/riboprotids.tsv'))
lridssplit<-rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
sridssplit<-rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
#get info on proteins so we can catagorize them
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
transreggoterm <- "GO:0006417"
transregprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
                                  filters=('go'),values=transreggoterm,
                                  mart = mart)
translationgoterm <- "GO:0006417"
translationprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
                                     filters=('go'),values=translationgoterm,
                                     mart = mart)[[1]]
ebp1pid = ms_tall$Protein_IDs[match('Pa2g4',ms_tall$gene_name)]
#intersect iwth pids
pids%<>%mutate(pcat = case_when(
  (Protein_IDs==ebp1pid) ~ "Ebp1",
  sep_element_in(Protein_IDs,sridssplit) ~ "Rps",
  sep_element_in(Protein_IDs,lridssplit) ~ "Rpl",
  sep_element_in(Protein_IDs,translationprotids) ~ "translation-associated",
  TRUE ~ "other"
))

#get only our select groups
selms_tall <- ms_tall %>% inner_join(filter(pids,pcat!='other'))


#now summarise 
selms_tall$fraction%>%unique
selms_tall_summ <- selms_tall %>% 
    # mutate(signal=log10(signal))%>%
  group_by(pcat,time,fraction)%>%
  summarize(
      n_distinct(signal),
      upper=quantile(signal,0.75,na.rm=T),
      lower=quantile(signal,0.25,na.rm=T),
      signal=median(signal,na.rm=T)
          )

pcats <- c('Ebp1','Rpl','Rps','translation-associated')
plotcats <- c('Ebp1','Rpl median','Rps median','Translation-associated median')
plotcols <- c('#EA1D3D','darkblue','#FADE4E','#010702')

selms_tall_summ$pcat <- list(selms_tall_summ$pcat) %>% c(setNames(plotcats,pcats)) %>% do.call(recode,args = .)

stopifnot(selms_tall_summ$signal < 100)

# signame <- paste0('Log10 ', unique(selms_tall$sigtype))
signame <- paste0(unique(selms_tall$sigtype), 'expression\nvs.\n proteome median=1')

selms_tall_summ$fraction%<>%recode(poly='Poly',total='Total')
selms_tall_summ$time %<>% str_replace('p','.')
dev.new(w=12,h=4)

selms_tall_summ %>%
  ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,ymax=upper,ymin=lower,color=pcat))+
  geom_point(size=4)+
  geom_line()+
  geom_linerange()+
  scale_x_continuous(name='',labels = unique(selms_tall_summ$time))+
  facet_grid(cols = vars(fraction))+
  scale_color_manual(name='',values = setNames(plotcols,plotcats))+
  # scale_y_continuous(name=signame)+
  scale_y_log10(name=signame)+
  theme_minimal()+
  theme(panel.grid = element_blank(),text = element_text(size=3))

# 
# library('ggpubr')
# plist <- lapply(1:3,function(i)qplot(1:10,1:10))
# do.call(ggpubr::ggarrange,c(plist,list(args = (nrow=1))))
# 
# ggarrange(plist[[1]],plist[[2]],plist[[3]],nrow=1)
