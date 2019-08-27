

#Plotting different segment parts next to one another
pdf(here('plots/integratae_exprdata2/endcount_vs_center.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
allsegcounts%>%
	filter(sample=='E13_ribo_1')%>%group_by(gene_id)%>%
	filter(protein_id==sample(protein_id,1))%>% 
	{qplot(data=.,x=log10(1+.$centercount),y=log10(1+.$enddcount))}
allsegcounts%>%
	filter(sample=='P0_ribo_2')%>%group_by(gene_id)%>%
	filter(protein_id==sample(protein_id,1))%>% 
	{qplot(data=.,x=log10(1+.$centercount),y=log10(1+.$enddcount))}
dev.off()


#Plotting different segment parts next to one another
pdf(here('plots/integratae_exprdata2/centercount_vs_spec.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
allsegcounts%>%
	filter(sample=='E13_ribo_1')%>%group_by(gene_id)%>%
	filter(protein_id==sample(protein_id,1))%>% 
	{qplot(data=.,x=log10(1+.$centercount),y=log10(1+.$spec_coef))}
allsegcounts%>%
	filter(sample=='P0_ribo_2')%>%group_by(gene_id)%>%
	filter(protein_id==sample(protein_id,1))%>% 
	{qplot(data=.,x=log10(1+.$centercount),y=log10(1+.$spec_coef))}
dev.off()

pdf(here('plots/integratae_exprdata2/ribo_ms_cor.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
qplot(x=top_protein_ids$ms_cor,fill=I('blue'),color=I('black'))+theme_bw()+theme(title=element_text(size=20))+ggtitle('distribution of Riboseq - MS correlation')+scale_x_continuous(name='correlation Riboseq MS')
dev.off()


pdf(here('plots/integratae_exprdata2/spec_count_comptable.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
spec_count_comp_table%>%.$specbetter%>%table
dev.off()

pdf(here('plots/integratae_exprdata2/spec_count_ddiff_hist.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
# qplot(data=spec_count_comp_table,x=ms_cor.x,y=ms_cor.y,geom='point',color=I('black'),fill=I("blue"))+theme_bw()
spec_count_comp_table%>%filter(!is.na(ms_cor.x),!is.na(ms_cor.y))%>%select(spectral_correlation=ms_cor.x,count_correlation=ms_cor.y)%>%gather(cortype,cor,-ms_id)%>%	
	qplot(data=.,x=cor,color=cortype,fill=cortype,alpha=I(0.5),geom='density')+scale_fill_manual(values=c('blue','red'))+
	theme_bw()
dev.off()

