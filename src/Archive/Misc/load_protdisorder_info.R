
disorderdf = readRDS('Downloads/df.rds')



disorderdf%>%mutate(parsedraw = RAW%>%str_extract_all('\\d+\\s+\\w+\\s+[0-9\\.]+(?=\r\n)'))%>%
	select(-RAW)%>%
	unnest(parsedraw)%>%
	# filter(parsedraw%>%str_count(' ')%>%`!=`(2))
	separate(parsedraw,c('pos','residue','score'),sep=' ')






1:10 + 3
add = `+`

(1:10)%>%add(3)

(1:10)%>%{
	x = .
	y = x +5
	add(y,3)
}


add(1:10,3)


disorderdf$RAW%>%str_extract_all('\\d+\\s+\\w+\\s+[0-9\\.]+(?=\\\\r\\\\n)')

mutate(raw = str_extract_all(raw,'\\d+\\s+\\w+\\s+[0-9\\.]+(?=\\\\r\\\\n)'))