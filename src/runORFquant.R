
orfquant_anno_file <-prepare_annotation_files(annotation_directory = "annotation_directory/" ,genome_seq = my_fasta_file,gtf_file =  my_gtf_file,scientific_name = "Human.test",annotation_name = "my_annotation",export_bed_tables_TxDb = F,create_TxDb = T)
stopifnot(file.exists(orfquant_anno_file))
#run orfquant
ORFquant_results <- run_ORFquant(for_ORFquant_file = my_orfquant_psites_file,  annotation_file = orfquant_anno_file, n_cores = 4)
