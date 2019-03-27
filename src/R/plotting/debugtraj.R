


ids <- fread('ids.txt')
satb2id <- ids%>%filter(gene_name=='Satb2')%>%.$gene_id

cdata <- fread('feature_counts/all_feature_counts')


sb2cdata <- cdata%>%filter(feature_id==satb2id)


e13ribo <- (sb2cdata$E13_ribo_1+sb2cdata$E13_ribo_2)/2
e13tot <- (sb2cdata$E13_total_1+sb2cdata$E13_total_2)/2


P0ribo <- (sb2cdata$P0_ribo_1+sb2cdata$P0_ribo_2)/2
P0tot <- (sb2cdata$P0_total_1+sb2cdata$P0_total_2)/2


log2(P0ribo/e13ribo)

(P0tot/e13tot)




trdata <- fread('exprdata/transformed_data.txt')

fread('ms_tables/ms_LFQ_total_ms_tall.tsv')%>%filter(gene_name=='Satb2')%>%filter(!fraction_missing)%>%
	group_by(time)%>%summarise(sig=mean(signal))%>%
	{.$sig[5] / .$sig[1]}


fread('exprdata/transformed_data.txt')%>%filter(gene_name=="Satb2")%>%
	gather%>%
	separate(key,c('time','assay','rep'))%>%
	group_by(time,assay)%>%summarise(sig=mean(as.numeric(value)))%>%
	filter(assay%in%c('ribo','MS'),time %in% c('E13','P0'))%>%arrange(assay)
