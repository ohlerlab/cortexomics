psitereformat = fread('pipeline/deepshapebamdata/E13_ribo_1.bam.reformat')
psitereformat %<>% set_colnames(c('readname','transcriptname','genename','CDSstart','CDSstop','mapstart','readlen'))

psitereformat%<>%mutate(pos = mapstart - CDSstart)

psitereformatgr = psitereformat%>%mutate(pos = mapstart - CDSstart)%>%{GRanges(.$transcriptname,IRanges(start=.$pos,width=1))}


left_join(data.frame(pos = -50:50),psitereformat,by='pos')%>%.$pos%>%table%>%.[as.character(-50:50)]%>%txtplot

left_join(data.frame(pos = -50:50),psitereformat,by='pos')%>%.$pos%>%table%>%sort



left_join(data.frame(pos = -20:20),psitereformat%>%filter(readlen==31),by='pos')%>%.$pos%>%table%>%sort




