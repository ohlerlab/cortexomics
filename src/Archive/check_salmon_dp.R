salmonribo = fread('pipeline/salmon/data/E13_ribo_1/quant.sf')
salmonribo$Name%<>%str_extract('ENSMUST\\d+')


salmonribo%>%left_join(enframe(iso_tx_countdata$abundance[,'E13_ribo_1'],'Name','dpTPM'))%>%
	filter(dpTPM>1,TPM>1)%>%
	{quicktest(log2(.$TPM),log2(.$dpTPM))}


salmonribo%>%left_join(enframe(iso_tx_countdata$counts[,'E13_ribo_1'],'Name','dpcount'))%>%
	filter(dpcount>1,NumReads>1)%>%
	{quicktest(log2(.$dpcount),log2(.$NumReads))}

	