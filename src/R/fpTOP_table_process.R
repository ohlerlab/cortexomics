toptable <- fread('../ext_data/tops_yamashita_etal_2008/supp_data_table_2.csv')
colnames(toptable)[5]<-'sequences'
toptable%<>% mutate(seq = strsplit(sequences,','))%>%unnest(.id='tss_num')

%>%write_tsv('topsseqs_unnest.tsv')

toptable %>% head%>% {paste0('>',.$Refseq )} 

dharnet_m@med-login2:pipeline$ cat ../ext_data/tops_yamashita_etal_2008/supp_data_tab1.csv | cut -f 1,8 -d ',' | sed 's/"/\n>/' | sed 's/,/\n/g' | tail -n +4 > ../ext_data/tops_yamashita_etal_2008/supp1seqs.fa