
msid2pid = metainfo%>%
	distinct(ms_id,protein_id)%>%
	{safe_hashmap(.[['ms_id']],.[['protein_id']])}
msid2upid = metainfo%>%
	distinct(ms_id,uprotein_id)%>%
	{safe_hashmap(.[['ms_id']],.[['uprotein_id']])}
upid2msid = metainfo%>%
	distinct(ms_id,uprotein_id)%>%
	{safe_hashmap(.[['uprotein_id']],.[['ms_id']])}


#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))
#get our matrix and design df
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)

sample_info_df <- data.frame(name = colnames(matchedms_mat),
                             stringsAsFactors = FALSE)
sample_info_df$time <- sample_info_df$name%>%str_extract('[^_]+')
sample_info_df$assay <- sample_info_df$name%>%str_extract('(?<=_)[^_]+')

bestmsmat=matchedms_mat[upid2msid[[best_uprotein_ids]],]

proDAfitms <- proDA::proDA(bestmsmat, design = ~ time, 
             col_data = sample_info_df%>%filter(assay=='MS'), reference_level = "E13")
prodalfcs <- proDA::test_diff(proDAfitms,contrast='timeP0')

#make combined matrix
proDAdatmat <- cbind(allcountmat[bestprotids,]%>%add(.1),bestmsmat)%>%set_rownames(rownames(bestmsmat))

proda_ms_pred <- proDAfitms$coefficients
proda_ms_pred[,2:5]=(proda_ms_pred[,2:5] + proda_ms_pred[,1])
colnames(proda_ms_pred)[1]='timeE13'
colnames(proda_ms_pred)%<>%str_replace('time(.*)','\\1_MS')


################################################################################
########Do fold changes at the 
################################################################################
