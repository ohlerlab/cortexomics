## ------------------------------------------------------------------------
matches = dplyr::matches
res_file <-
  file.path(modeldir,"intermediate_results",
            paste0(my_contrast_name, "_results.rds"))

isdeseq<-is.null(parameters$alt_res_file)

res <- readRDS(file = res_file)
res %<>% select(-matches('gene_name'))
get_de_summary(res, my_alpha, my_log2fc_th) %>%
  dplyr::select(category = category_nice, number) %>%
  knitr::kable()

#show the beta names and model coeefs for the contrast as well
#nb this may not work with betaPrior= FALSE
#TODO - make this nicer for external programs
if(isdeseq){

    modelcoeffs <- res@elementMetadata@listData$description[2]%>%str_extract('[^:]*$')%>%str_split(',')%>%.[[1]]

  knitr::kable(
    data_frame(
      modelterms=resultsNames(dds),
      modelcoeffs
  ))
}



## ------------------------------------------------------------------------
plot_DiffMA(res, alpha = my_alpha, ylim = c(z_min,z_max))


## ------------------------------------------------------------------------
up_regul_features_nl <-
  get_sig_log2fc(res, red_feature_annot, alpha =  my_alpha, 
                 up_log2fc_th = my_log2fc_th, down_log2fc_th = -Inf, links = F) 

up_file <- file.path(modeldir,"tables/", paste0(my_contrast_name, "_up_regul_features.csv"))
write_csv(up_regul_features_nl, path = up_file)

rds_up_file <- 
  file.path(modeldir,"intermediate_results",
            paste0(my_contrast_name, "_up_regul_features.rds"))
saveRDS(up_regul_features_nl, rds_up_file)

up_regul_features <-
  get_sig_log2fc(res, red_feature_annot, alpha =  my_alpha, 
                 up_log2fc_th = my_log2fc_th, down_log2fc_th = -Inf, links = T)


## ------------------------------------------------------------------------
library(DT)
up_regul_features %>%
  tablefunc(rownames = F, 
                escape = F,
                extensions = 'Buttons', 
                options = list(
                  dom='Bfrtip',
                  buttons = list(
                    list(
                      extend='csv',
                      buttons=c('csv'),
                      text='download')
                  )
               ,.test=is_markdowntest )
)


## ----plot_up_features, fig.height=4.5------------------------------------

if (dim(up_regul_features_nl)[1] > 0) {

  goi_expr <-
    dplyr::left_join(up_regul_features_nl %>% 
                       dplyr::arrange(adj_p_value) %>%
                       head(n=12) %>%
                       select(-matches('gene_name') )%>%
                       dplyr::select(feature_id,gene_name), 
                     reg_log2_cpk,
                     by='feature_id') %>%
    dplyr::left_join(dplyr::select(sample_annot, sample_id, group),
                     by='sample_id')
  
  goi_expr %>%
    ggplot(aes(x = group, y = reg_log2_cpk, color = group)) +
    geom_boxplot(outlier.shape=NA,show.legend=FALSE) +
    geom_jitter(height = 0, width = 0.1, color='black', size=1) +
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    scale_x_discrete(limits = unique(goi_expr$group)) +
    facet_wrap(~gene_name,scales='free_y',ncol=4)  
}


## ------------------------------------------------------------------------
down_regul_features_nl <-
  get_sig_log2fc(res, red_feature_annot, alpha =  my_alpha, 
                 up_log2fc_th = +Inf, down_log2fc_th = -my_log2fc_th, links = F) 

down_file <- file.path(modeldir,"tables", paste0(my_contrast_name, "_down_regul_features.csv"))
write_csv(down_regul_features_nl, path = down_file)

rds_down_file <-
  file.path(modeldir,"intermediate_results",
            paste0(my_contrast_name, "_down_regul_features.rds"))
saveRDS(down_regul_features_nl, rds_down_file)

down_regul_features <-
  get_sig_log2fc(res, red_feature_annot, alpha =  my_alpha, 
                 up_log2fc_th = +Inf, down_log2fc_th = -my_log2fc_th, links = T)


## ------------------------------------------------------------------------
down_regul_features %>%
  tablefunc(rownames = F, 
                escape = F,
                extensions = 'Buttons', 
                options = list(
                  dom='Bfrtip',
                  buttons = list(
                    list(
                      extend='csv',
                      buttons=c('csv'),
                      text='download')
                  )
                ),.test=is_markdowntest
  )


## ----plot_down_features, fig.height=4.5----------------------------------

if (dim(down_regul_features_nl)[1] > 0) {

  goi_expr <-
    dplyr::left_join(down_regul_features_nl %>% 
                       dplyr::arrange(adj_p_value) %>%
                       head(n=12) %>%
                      select(-matches('gene_name'))%>%
                       dplyr::select(feature_id,gene_name), 
                     reg_log2_cpk,
                     by='feature_id') %>%
    dplyr::left_join(dplyr::select(sample_annot, sample_id, group),
                     by='sample_id')
  
  goi_expr %>%
    ggplot(aes(x = group, y = reg_log2_cpk, color = group)) +
    geom_boxplot(outlier.shape=NA,show.legend=FALSE) +
    geom_jitter(height = 0, width = 0.1, color='black', size=1) +
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    scale_x_discrete(limits = unique(goi_expr$group)) +
    facet_wrap(~gene_name,scales='free_y',ncol=4)  
}

