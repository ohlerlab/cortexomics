

heatscatter = LSD::heatscatter
debugcompdf%>%{heatscatter(cor=TRUE,xlab='RPF Density',.$ldens, ylab='iBAQ',y=.$ibaq_E13_1 ,cex.lab=1.5)}
debugcompdf%>%{heatscatter(cor=TRUE,xlab='log2(sqrt(TrP))',log2(.$spec), ylab = 'iBAQ',y=.$ibaq_E13_1 ,cex.lab=2)}

par(mfrow=c(1,2))
deepshape_orfquant_compfilt%>%{heatscatter(cor=TRUE,xlab='log2(DpPrime TPMs)',.$l2TPM, ylab='iBAQ',y=.$ibaq_E13_1 ,cex.lab=1.5)}
deepshape_orfquant_compfilt%>%{heatscatter(cor=TRUE,xlab='log2(ORF_pM)',(.$l2ORFs_pM), ylab = 'iBAQ',y=.$ibaq_E13_1 ,cex.lab=2)}

#ORFquant sucks
inclusiontable(orfquantdf$gene_name,deepshapedf$gene_name)%>%reshape2::melt()%>%set_colnames(c('In ORFquant','In DpPrime','Number'))%>%grid.table

par(mfrow=c(1,2))
deepshape_orfquant_compfilt%>%{heatscatter(cor=TRUE,xlab='log2(DpPrime TPMs)',.$l2TPM, ylab='iBAQ',y=.$ibaq_E13_1 ,cex.lab=1.5)}
deepshape_orfquant_compfilt%>%{heatscatter(cor=TRUE,xlab='log2(ORF_pM)',(.$l2ORFs_pM), ylab = 'iBAQ',y=.$ibaq_E13_1 ,cex.lab=2)}

#Isoform specific expression rocks
par(mfrow=c(1,2))
tpmaggdf%>%filter(TPM>1)%>%filter(ms_id%in%ids2comp)%>%{heatscatter(cor=TRUE,xlab='log2(DpPrime TPMs)',log2(.$TPM),ylab='iBAQ',y=.$ibaq_E13_1 ,cex.lab=1.5)}
tpm_kallistodf%>%filter(TPM>1)%>%filter(ms_id%in%ids2comp)%>%{heatscatter(cor=TRUE,xlab='log2(RNseq TPMs) - Kallisto',log2(.$TPM),ylab='iBAQ',y=.$ibaq_E13_1 ,cex.lab=1.5)}


