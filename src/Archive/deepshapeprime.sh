set -ex
inputfastq=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/trim_reads/E13_ribo_1/E13_ribo_1.fastq.gz
# inputfastq=test.fastq.gz
prefix=$(basename $inputfastq)
prefix=${prefix%.fastq.gz}
dpdir=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/Applications/DeepShape/DeepShape-prime/

trfasta=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa
ln -fs $trfasta
trfasta=$(basename $trfasta)
GTF=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.annotation.gtf
ln -fs $GTF
GTF=$(basename $GTF)

transcriptbam=align/transcript/${prefix}transcript_Aligned.out.bam

# samtools view \
# 	$transcriptbam | \
# 	awk '{if(length($10)>=25 && length($10)<=36){print $1"\t"$3"\t"$4"\t"length($10);}}' | \
# 	awk -F"|" '{printf $1"\t"$2; for(i=8;i<=NF-1;i++){if($i~/CDS:/) {gsub("CDS:","",$i); gsub("-","\t",$i); printf "\t"$i;}} printf $NF"\n"}' \
# 	> ${prefix}.bam.reformat

# python ${dpdir}/ParseBam2RiboPos.py \
# 	$GTF \
# 	$trfasta \
# 	${prefix}.bam.reformat \
# 	${prefix}.bam.reformat.absolute

mkdir -p DeepShapePrimeOutputs

python ${dpdir}/DeepShape-prime.py \
	${prefix}.bam.reformat.absolute \
	$trfasta \
	./DeepShapePrimeOutputs \
	200

