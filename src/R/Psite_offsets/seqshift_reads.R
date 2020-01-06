


# fcountsdf<-'pipeline/feature_counts/all_feature_counts'%>%fread%>%select(gene_id=feature_id,fcount=thissamp)
# p2gdf<-with(cds,data_frame(protein_id,gene_id))%>%distinct

# fcountcompdf<-segcountdf%>%safe_left_join(.,p2gdf,by=c('protein_id'))%>%
# 	left_join(fcountsdf)

# fcountcompdf%>%head

# fcountcompdf%>%{cor(.$fcount,.$total,use='complete')}

# fcountcompdf%>%{cor(log10(1+.$fcount),log10(1+.$total),use='complete')}

# fcountcompdf%>%{txtplot(log10(1+.$fcount),log10(1+.$total))}

# ns_all_genes%>%map('result')%>%map_lgl(is.null)

# segcountdf%>%head



0
# mapped<-mapped%>%head(5)
# reads<-reads%>%head(5)

# 	strand(mapped)<-'+'

# mapped$seqoffset <- get_seqforrest_traindata(mapped,exonseq,trim=FALSE)%>%
# 	predict(psite_model$seqshiftmodel,.)%>%
# 	.$prediction%>%
# 	as.character%>%
# 	as.numeric



# mapFromTranscripts(mapped,exons2use)
# reads




####
#alternative approaches - try looking at independance of end distances
#Try looking at artifical positive data
#look at read length distribution around pos and neg







