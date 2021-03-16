
cairo_pdf(h=5,w=5,'plots/figures/figure2/cumdist_TEchange_fplen.pdf')
allTEchangedf%>%
	left_join(transcript_fp_utrlen%>%left_join(tr_gid_df)%>%group_by(gene_name)%>%slice(which.max(fputr_length)))%>%
	arrange(fputr_length)%>%
	mutate(TE_change = case_when(up==1 ~ 'up',down==1~ 'down',TRUE ~ 'None'))%>%
	group_by(TE_change)%>%
	mutate(`Cumulative Fraction`= (1:n() / n()) )%>%
	ggplot(data=.,aes(y=`Cumulative Fraction`,color=TE_change,x=log10(fputr_length)))+geom_line()+
	theme_bw()
dev.off()
normalizePath('plots/figures/figure2/cumdist_TEchange_fplen.pdf')

cairo_pdf(h=5,w=5,'plots/figures/figure2/cumdist_TEchange_tplen.pdf')
allTEchangedf%>%
	left_join(transcript_tp_utrlen%>%left_join(tr_gid_df)%>%group_by(gene_name)%>%slice(which.max(tputr_length)))%>%
	arrange(tputr_length)%>%
	mutate(TE_change = case_when(up==1 ~ 'up',down==1~ 'down',TRUE ~ 'None'))%>%
	group_by(TE_change)%>%
	mutate(`Cumulative Fraction`= (1:n() / n()) )%>%
	ggplot(data=.,aes(y=`Cumulative Fraction`,color=TE_change,x=log10(tputr_length)))+geom_line()+
	theme_bw()
dev.off()
normalizePath('plots/figures/figure2/cumdist_TEchange_tplen.pdf')