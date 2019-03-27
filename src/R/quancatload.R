

quancatdf <- fread('../ext_data/quancat_proteinGroups.txt')

quancatdf%>%colnames%>%str_subset('soma')

quancatdf%>%
	filter(!is.na(log2pro_somaH.soma_proM.4.6))%>%
	filter(!is.na(log2pro_somaM.soma_proH.8.2))%>%
	nrow