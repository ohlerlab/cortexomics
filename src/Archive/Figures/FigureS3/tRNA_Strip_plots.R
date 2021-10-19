{
dir.create('plots/Figures/FigureS3/',showWarn=F,rec=T)
################################################################################
######## Now plot it
################################################################################
# how each decoder
plotfile <- "plots/Figures/FigureS3/trna_sig_total_alldecod.pdf"
pdf(plotfile, w = 9, h = 16)
allcodonsig %>%
    group_by(decoder) %>%
    mutate(abundance = -dCt)%>%
    mutate(dmean = mean(abundance)) %>%
    group_by(anticodon) %>%
    arrange(dmean) %>%
    mutate(decoder = as_factor(decoder)) %>%
    filter(sample %>% str_detect("Total")) %>%
    group_by(anticodon) %>%
    # group_by(decoder,)
    ggplot(data = ., aes(x = decoder, color = time, y = abundance)) +
    geom_point() +
    scale_color_manual(values = stagecols) +
    scale_y_continuous(name = "- deltaCt", breaks = seq(0, -30, by = -5)) +
    coord_flip() +
    theme_bw()
# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)
#
# how each decoder
plotfile <- "plots/Figures/FigureS3/trna_sig_poly_alldecod.pdf"
pdf(plotfile, w = 9, h = 16)
allcodonsig %>%
    group_by(decoder) %>%
    mutate(abundance = -dCt)%>%
    mutate(dmean = mean(abundance)) %>%
    group_by(anticodon) %>%
    arrange(dmean) %>%
    mutate(decoder = as_factor(decoder)) %>%
    filter(sample %>% str_detect("Poly")) %>%
    group_by(anticodon) %>%
    # group_by(decoder,)
    ggplot(data = ., aes(x = decoder, color = time, y = abundance)) +
    geom_point() +
    scale_color_manual(values = stagecols) +
    scale_y_continuous(name = "- deltaCt", breaks = seq(0, -30, by = -5)) +
    coord_flip() +
    theme_bw()
dev.off()
normalizePath(plotfile)


plotfile <- "plots/Figures/FigureS3/trna_sig_allfrac_mergedecod.pdf"
pdf(plotfile, w = 9*2, h = 16)
allcodsigmean_isomerge %>%
    # filter(sample %>% str_detect("Total")) %>%
    group_by(anticodon) %>%
    mutate(cmean = mean(abundance)) %>%
    ungroup() %>%
    arrange(cmean) %>%
    mutate(anticodon = as_factor(anticodon)) %>%
    # group_by(decoder,)
    mutate(fraction=fraction%>%factor(.,c('Total','Poly','80S')))%>%
    ggplot(data = ., aes(x = anticodon, color = stageconv[time], y = abundance)) +
    facet_grid( . ~ fraction,scale='free_x')+
    geom_point() +
    scale_color_manual(values = stagecols) +
    scale_y_continuous(name = "- deltaCt", breaks = -seq(0, 20, by = 2.5)) +
    coord_flip() +
    theme_bw()
# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)

#
plotfile <- "plots/Figures/FigureS3/trna_sig_poly_mergedecod.pdf"
pdf(plotfile, w = 9, h = 16)
allcodsigmean_isomerge %>%
    filter(sample %>% str_detect("Poly")) %>%
    group_by(anticodon) %>%
    mutate(cmean = mean(abundance)) %>%
    ungroup() %>%
    arrange(cmean) %>%
    mutate(anticodon = as_factor(anticodon)) %>%
    # group_by(decoder,)
    ggplot(data = ., aes(x = anticodon, color = stageconv[time], y = abundance)) +
    geom_point() +
    scale_color_manual(values = stagecols) +
    scale_y_continuous(name = "- deltaCt", breaks = -seq(0, 20, by = 2.5)) +
    coord_flip() +
    theme_bw()
dev.off()
normalizePath(plotfile)


# group them on same row
plotfile <- "plots/Figures/FigureS3/isomerge_tRNA_sig_strip_poly.pdf"
pdf(plotfile, w = 9, h = 16)
allcodsigmean_isomerge %>%
    filter(sample %>% str_detect("Total")) %>%
    group_by(anticodon) %>%
    mutate(cmean = mean(abundance)) %>%
    ungroup() %>%
    arrange(cmean) %>%
    mutate(anticodon = as_factor(anticodon)) %>%
    # group_by(decoder,)
    ggplot(data = ., aes(x = anticodon, color = time, y = abundance)) +
    geom_point() +
    scale_color_manual(values = stagecols) +
    scale_y_continuous(breaks = seq(0, 20, by = 2.5)) +
    coord_flip() +
    theme_bw()
dev.off()
normalizePath(plotfile)
# allcodsigmean$ntime<-allcodsigmean$time%>%as_factor%>%as.numeric
# allcodsigmean$fraction%<>%factor(c('Total','Poly','80S'))
}