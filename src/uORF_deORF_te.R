

#this qill quickly check if we see regulation in opposing directions in our 
defile <- '/fast/groups/ag_ohler/dharnet_m/cortexomics/pipeline/xtail/xtail_P0.txt'%>%fread

xtailres <- 'xtail/xtail_E175.txt'%>%fread

dev.off()
plot(1)


xtailres$feature_id%<>%str_replace('uORF_ ','uORF_')

uorfgene_te_tbl<-inner_join(
	xtailres%>%filter(str_detect(feature_id,'uORF'))%>%select(feature_id,log2fc)%>%
		mutate(feature_id=str_split_fixed(feature_id,'_',2)%>%.[,2] )
	,
	xtailres%>%filter(!str_detect(feature_id,'uORF'))%>%select(feature_id,log2fc),
		by='feature_id'
)%>%set_colnames(c('feature_id','uORF_te_fc','gene_te_fc'))

svglite(h=5,w=6,'../plots/techange_uORFs_vs_CDS.svg'%>%normalizePath%T>%message)
uorfgene_te_tbl%>%ggplot(aes(x=uORF_te_fc,y=gene_te_fc))+geom_point()+
scale_x_continuous(limits=c(-5,5))+
scale_y_continuous(limits=c(-5,5))+
theme_minimal()+
ggtitle('Changes in TE in uORF vs CDS ')
dev.off()
