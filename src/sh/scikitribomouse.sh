#ata and Statistical Analysis in This Paper
#For the wild-type S. cerevisiae analysis and validation, the Riboseq (flash-freeze protocol) 
#and RNA-seq (Ribo-zero protocol) data were from Weinberg et al(Weinberg et al., 2016). 
#The accession numbers are GSM1289257, GSM1289256. For the CHX comparison, the CHX-treated data is SRR948553 and the RNA-seq data is SRR948551, from McManus et al(McManus et al., 2014). The reference genome of S. cerevisiae used is S288C R64-2-1. The gene annotation file was the SGD annotation downloaded from UCSC. For the E. coli analysis, the Riboseq (RelE protocol) and RNA-seq data were from Hwang et al(Hwang and Buskirk, 2017). The accession number is GSE85540. The reference genome of E. coli used is the MG1655 genome. For more details of how these data were generated, please refer to the original papers. All the figures in the paper were plotted using matplotlib(Hunter, 2007) (v2.0.0) and seaborn(Waskom and Wagner, 2017) (v0.7.1). The Pearson correlation and Spearman correlation are denoted as r and ρ, respectively.
# processing the yeast data https://www.nature.com/articles/s41594-018-0080-2?WT.feed_name=subjects_translation
# nmismatch=1
# adapter="TCGTATGCCGTCTTCTGCTTG"
# STAR --runThreadN 8 --runMode genomeGenerate --genomeDir index --genomeFastaFiles Yeast.saccer3.fa --genomeSAindexNbases 11 --genomeChrBinNbits 12

# SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM "
# align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical --sjdbOverhang 35"
# STAR --runThreadN 8 --genomeDir index --readFilesIn <(zcat input/srr104951/SRR1049521.trim.fastq.gz) --outFileNamePrefix 'weinberg_yeast_ribo' ${SAM_params} ${align_params}
# STAR --runThreadN 8 --genomeDir index --readFilesIn <(zcat input/srr1049520/SRR1049520.trim.fastq.gz) --outFileNamePrefix 'weinberg_yeast_rna' ${SAM_params} ${align_params}

# if [ ! -f decoys.txt ]; then
#   grep "^>" Yeast.saccer3.fa | cut -d " " -f 1 > decoys.txt
#   sed -i.bak -e 's/>//g' decoys.txt
# fi

# salmonbin=/fast/work/users/dharnet_m/Applications/salmon-latest_linux_x86_64/bin/salmon
# if [ ! -f salmon_index ]; then
#    $salmonbin index -t <(cat yeast_transcripts.fa Yeast.saccer3.fa) -k21 -d decoys.txt -p 4 -i salmon_index
# fi

# rm -rf weinberg_yeast_rna_salm
# $salmonbin quant -i salmon_index -l SF -r input/srr1049520/SRR1049520.fastq.gz  --output weinberg_yeast_rna_salm -p 4
# # $salmonbin quant -t yeast_transcripts.fa -l SF -a weinberg_yeast_rnaAligned.out.sort.bam -o weinberg_yeast_rna_salm -p 4 --seqBias


# #create a lexographic sam header
# samtools view -H $bam | head -n 1 > sqheadersorted
# samtools view -H $bam | grep SQ | perl -lane '/SN:(chr\w+)/;print $1,"\t",$_' | sort -k 1,1 | cut -f 2,3,4,5,6 >> sqheadersorted
# samtools view -H weinberg_yeast_riboAligned.out.sort.bam | grep -v -Pe '@HD|@SQ' >> sqheadersorted
# cat sqheadersorted

# samtools reheader sqheadersorted weinberg_yeast_riboAligned.out.bam | samtools sort - -o weinberg_yeast_riboAligned.out.sort.bam

# samtools view -H weinberg_yeast_riboAligned.out.sort.bam

sktest_rnafold.txt
sktest.gtf

pipeline/sktesttrs.fa
#trlens[sktesttrs]%>%lapply(function(n)rep(0,n))%>%imap_chr(~paste0(.y,'\t',paste(.x,collapse=' ')))%>%paste0(collapse='\n')%>%cat(file='pipeline/sktesttrs.rnafold.txt')
#exonsgrl[sktesttrs] %>% unlist %>%resize(width(.)+1e3)%>%gaps%>%export('pipeline/sklearntest_unmap.bed')
# gtf=../annotation/gencode.vM12.annotation.gtf
#grep -f sktesttrs.txt ../annotation/gencode.vM12.annotation.gtf  > sktest.gtf

gtf=sktest.gtf
fasta=my_GRCm38.p5.genome.chr_scaff.fa
prefix=E13_ribo_1_sk
rnafold=sktesttrs.rnafold.txt
RNA=salmon/data/E13_ribo_1/quant.sf
index_folder=sk_index
unmap=sklearntest_unmap.bed
bam=star/data/E13_ribo_1/E13_ribo_1.bam
output=E13_ribo_1__skout

python /fast/work/groups/ag_ohler/dharnet_m/cortexomics/Applications/scikit-ribo/scikit_ribo/scikit-ribo-build.py \
-g $gtf \
-f $fasta \
-p $prefix \
-r $rnafold \
-t $RNA \
-o $index_folder


gtf=Yeast.sacCer3.sgdGene.sort.gtf
fasta=Yeast.saccer3.fa
prefix=E13_ribo_1__sk
rnafold=yeast_sacCer3_rnafold.txt
RNA=E13_ribo_1_salm/quant.sf
index_folder=sk_index
unmap=2013-01-23.sc_o12_26_2_crossmap.nochrm.bed
bam=E13_ribo_1_Aligned.out.sort.sample.bam
output=E13_ribo_1_skout

scikit-ribo-run.py \
-i $bam \
-f $index_folder \
-p $prefix \
-o $output \
-u $unmap


#
#create a lexographic sam header
bam
samtools view -H $bam | head -n 1 > sqheadersorted
samtools view -H $bam | grep SQ | perl -lane '/(.*?)\|/;print $1' >> sqheadersorted
samtools view -H $bam | grep -v -Pe '@HD|@SQ' >> sqheadersorted
cat sqheadersorted
samtools reheader sqheadersorted $bam | samtools sort - -o weinberg_yeast_riboAligned.out.sort.bam





