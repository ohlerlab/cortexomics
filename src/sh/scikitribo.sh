################################################################################
########Part1
################################################################################
	
#ata and Statistical Analysis in This Paper
#For the wild-type S. cerevisiae analysis and validation, the Riboseq (flash-freeze protocol) 
#and RNA-seq (Ribo-zero protocol) data were from Weinberg et al(Weinberg et al., 2016). 
#The accession numbers are GSM1289257, GSM1289256. For the CHX comparison, the CHX-treated data is SRR948553 and the RNA-seq data is SRR948551, from McManus et al(McManus et al., 2014). The reference genome of S. cerevisiae used is S288C R64-2-1. The gene annotation file was the SGD annotation downloaded from UCSC. For the E. coli analysis, the Riboseq (RelE protocol) and RNA-seq data were from Hwang et al(Hwang and Buskirk, 2017). The accession number is GSE85540. The reference genome of E. coli used is the MG1655 genome. For more details of how these data were generated, please refer to the original papers. All the figures in the paper were plotted using matplotlib(Hunter, 2007) (v2.0.0) and seaborn(Waskom and Wagner, 2017) (v0.7.1). The Pearson correlation and Spearman correlation are denoted as r and ρ, respectively.
# processing the yeast data https://www.nature.com/articles/s41594-018-0080-2?WT.feed_name=subjects_translation
#fasterq-dump -O input/srr1049520 -f -x -p -e 4 GSE118953
#fasterq-dump SRR1049521  -O input/srr104951/ -x -p -e 4

#download yeast genome here, with gtf files etc
#https://github.com/schatzlab/scikit-ribo/tree/master/example
nmismatch=1
adapter="TCGTATGCCGTCTTCTGCTTG"
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir index --genomeFastaFiles Yeast.saccer3.fa --genomeSAindexNbases 11 --genomeChrBinNbits 12


#trim 8 basepairs (The barcode, and an extra one) fromt he riboseq
zcat input/srr104951/SRR1049521.fastq.gz  | ~/work/bin/collapse_reads.pl  SRR1049521 | cutadapt  -j 8 -g '^N{9}'  --match-read-wildcards - | cutadapt -j 8 -a TCGTATGCCGTCTTCTGCTTG -m 10 - > input/srr104951/SRR1049521.trim.fastq.gz
zcat input/srr1049520/SRR1049520.fastq.gz  | ~/work/bin/collapse_reads.pl  SRR1049520 | cutadapt  -j 8 -g '^N{9}'  --match-read-wildcards - | cutadapt -j 8 -a TCGTATGCCGTCTTCTGCTTG -m 10 - > input/srr1049520/SRR1049520.trim.fastq.gz

#align
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM "
align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical --sjdbOverhang 35"
STAR --runThreadN 8 --genomeDir index --readFilesIn <(zcat input/srr104951/SRR1049521.trim.fastq.gz) --outFileNamePrefix 'weinberg_yeast_ribo' ${SAM_params} ${align_params}
STAR --runThreadN 8 --genomeDir index --readFilesIn <(zcat input/srr1049520/SRR1049520.trim.fastq.gz) --outFileNamePrefix 'weinberg_yeast_rna' ${SAM_params} ${align_params}

if [ ! -f decoys.txt ]; then
  grep "^>" Yeast.saccer3.fa | cut -d " " -f 1 > decoys.txt
  sed -i.bak -e 's/>//g' decoys.txt
fi

salmonbin=/fast/work/users/dharnet_m/Applications/salmon-latest_linux_x86_64/bin/salmon
if [ ! -f salmon_index ]; then
   $salmonbin index -t <(cat yeast_transcripts.fa Yeast.saccer3.fa) -k21 -d decoys.txt -p 4 -i salmon_index
fi

rm -rf weinberg_yeast_rna_salm
$salmonbin quant -i salmon_index -l SF -r input/srr1049520/SRR1049520.fastq.gz  --output weinberg_yeast_rna_salm -p 4
$salmonbin quant -i salmon_index -l U -r input/srr104951/SRR1049521.fastq.gz  --output weinberg_yeast_ribo_salm -p 4
# $salmonbin quant -t yeast_transcripts.fa -l SF -a weinberg_yeast_rnaAligned.out.sort.bam -o weinberg_yeast_rna_salm -p 4 --seqBias


#create a lexographic sam header
# samtools view -H weinberg_yeast_riboAligned.out.sort.sample.bam | head -n 1 > sqheadersorted
# samtools view -H weinberg_yeast_riboAligned.out.sort.sample.bam | grep SQ | perl -lane '/SN:(chr\w+)/;print $1,"\t",$_' | sort -k 1,1 | cut -f 2,3,4,5,6 >> sqheadersorted
# samtools view -H weinberg_yeast_riboAligned.out.sort.bam | grep -v -Pe '@HD|@SQ' >> sqheadersorted
# cat sqheadersorted
# samtools reheader sqheadersorted weinberg_yeast_riboAligned.out.bam | samtools sort - -o weinberg_yeast_riboAligned.out.sort.bam
# samtools view -H weinberg_yeast_riboAligned.out.sort.bam

################################################################################
########Part 2
################################################################################
	


gtf=Yeast.sacCer3.sgdGene.sort.gtf
fasta=Yeast.saccer3.fa
prefix=weinberg_yeast_ribo_sk
rnafold=yeast_sacCer3_rnafold.txt
RNA=weinberg_yeast_rna_salm/quant.sf
index_folder=sk_index
unmap=2013-01-23.sc_o12_26_2_crossmap.nochrm.bed
bam=weinberg_yeast_riboAligned.out.sort.bam
output=weinberg_yeast_ribo_skout
/fast/work/groups/ag_ohler/dharnet_m/cortexomics/Applications/scikit-ribo/scikit_ribo/scikit-ribo-build.py \
-g $gtf \
-f $fasta \
-p $prefix \
-r $rnafold \
-t $RNA \
-o $index_folder


gtf=Yeast.sacCer3.sgdGene.sort.gtf
fasta=Yeast.saccer3.fa
prefix=weinberg_yeast_ribo_sk
rnafold=yeast_sacCer3_rnafold.txt
RNA=weinberg_yeast_rna_salm/quant.sf
index_folder=sk_index
unmap=2013-01-23.sc_o12_26_2_crossmap.nochrm.bed
# bam=weinberg_yeast_riboAligned.out.sort.bam
bam=weinberg_yeast_riboAligned.out.sort.sample.bam
output=weinberg_yeast_ribo_skout

scikit-ribo-run.py \
-i $bam \
-f $index_folder \
-p $prefix \
-o $output \
-u $unmap

