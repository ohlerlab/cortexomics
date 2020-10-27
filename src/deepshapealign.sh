set -ex
inputfastq=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/trim_reads/E13_ribo_1/E13_ribo_1.fastq.gz
# inputfastq=test.fastq.gz
prefix=$(basename $inputfastq)
prefix=${prefix%.fastq.gz}


trfasta=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa
ln -fs $trfasta
trfasta=$(basename $trfasta)
GTF=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.annotation.gtf
ln -fs $GTF
GTF=$(basename $GTF)
mkdir -p align

mkdir -p align/contaminant/

filter_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"

SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM"

STAR --runThreadN 15 --genomeDir StarIndex/contaminant/ --readFilesIn <( zcat $inputfastq) --outFileNamePrefix align/contaminant/${prefix}. --outStd BAM_Unsorted --outReadsUnmapped Fastx --outSAMmode NoQS ${filter_params} > ${prefix}.contam.bam

mkdir -p align/transcript

filter_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM"
STAR --runThreadN 15 --genomeDir StarIndex/transcript/ --readFilesIn align/contaminant/${prefix}.Unmapped.out.mate1 --outFileNamePrefix align/transcript/${prefix}transcript_ ${SAM_params} ${filter_params}

samtools view align/transcript/${prefix}transcript_Aligned.out.bam | awk '{if(length($10)>=25 && length($10)<=36){print $1"\t"$3"\t"$4"\t"length($10);}}' | awk -F"|" '{printf $1"\t"$2; for(i=8;i<=NF-1;i++){if($i~/CDS:/) {gsub("CDS:","",$i); gsub("-","\t",$i); printf "\t"$i;}} printf $NF"\n"}' > ${prefix}.bam.reformat

python /fast/work/groups/ag_ohler/dharnet_m/cortexomics/Applications/DeepShape/DeepShape-prime/ParseBam2RiboPos.py \
	${GTF} \
	${trfasta%.fa}_filter.fa \
	${prefix}.bam.reformat \
	${prefix}.bam.reformat.absolute