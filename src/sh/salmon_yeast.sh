set -x 
#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz

if [ ! -f decoys.txt ]; then
  grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
  sed -i.bak -e 's/>//g' decoys.txt
fi

#cat gencode.vM25.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz> gentrome.fa.gz
#conda create -c bioconda  -n salmon salmon
#conda activate salmon
#if [ ! -f salmon_index ]; then
#    salmon index -t gentrome.fa.gz -d decoys.txt -p 4 -i salmon_index --gencode
#fi
#salmoi salmon_index -l SR -r Glial_synaptogenesis/GABAergic_noGlia_rep1_S1_R1_001.fastq.gz --validateMappings -o transcripts_quant -p 4 --seqBias
#salmon quant -i salmon_index -l SR -r Glial_synaptogenesis/GABAergic_noGlia_rep1_S1_R1_001.fastq.gz --validateMappings -o transcripts_quant -p 4 --seqBias
#salmon quant -i salmon_index -l SR -r Glial_synaptogenesis/GABAergic_noGlia_rep1_S1_R1_001.fastq.gz --validateMappings -o transcripts_quant_nb -p 4

fastqs=Glial_synaptogenesis/*R1_001.fastq.gz

for i in $fastqs; do 
    samplename=${i##*/}
    samplename=${samplename%_S1_R1_001.fastq.gz}
    salmon quant -i salmon_index -l SR -r $i --validateMappings -o $samplename -p 4 --seqBias
done

#make table of transcript ids to gene ids
zgrep -e '>' gencode.vM25.transcripts.fa.gz | cut -d '|' -f 1,2 | tr -d '>' | tr '|' '\t' > trid_gid.tsv

Rscript makecounts.R

