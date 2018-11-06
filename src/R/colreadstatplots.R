readlengths <- 
	Sys.glob('collapse_reads/*/*collreadstats.txt')%>%
	grep(v=T,inv=T,patt='test')%>%
	# str_subset('Poly|80S')%>%
	map(.%>%{
		message(.)
		lines <- readLines(.)
		tabstart = lines%>%str_detect('>read_lengths_unique')%>%which
		tab <- fread(.,skip=tabstart+2)
		tab$sample <- dirname(.)%>%basename
		tab%<>%set_colnames(c('readlength','count','sample'))
		tab
	})%>%bind_rows

readlengths$readlength = readlengths$readlength - 8

pdf(w=14,h=3,'../plots/readlendists_Poly_80S.pdf')
readlengths%>%
	# filter(sample%>%str_detect('ribo'))%>%
	filter(sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_Poly_80S.pdf')


pdf(w=14,h=3,'../plots/readlendists_totalribo.pdf')
readlengths%>%
	filter(sample%>%str_detect('ribo'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_totalribo.pdf')


pdf(w=14,h=3,'../plots/readlendists_all.pdf')
readlengths%>%
	filter(sample%>%str_detect('ribo|80S|Poly'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_all.pdf')


readlengths%>%group_by(sample)%>%summarize(max_rl = readlength[which.max(count)])%>%filter(!str_detect(sample,'total'))%>%arrange(desc(max_rl))


pdf(w=14,h=3,'../plots/readlendists_80SE13.pdf')
readlengths%>%
	filter(sample%>%str_detect('80SE13'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_80SE13.pdf')


readlengths%>%group_by(sample)%>%summarize(max_rl = readlength[which.max(count)])%>%filter(!str_detect(sample,'total'))%>%arrange(desc(max_rl))