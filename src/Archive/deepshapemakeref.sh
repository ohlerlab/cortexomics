set -ex
GTF=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.annotation.gtf
trfasta=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa
translations=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_translations.fa
dpdir=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/Applications/DeepShape/PreprocessRefs/
rrna_fasta=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.transcripts.fa
tRNAfasta=/fast/work/groups/ag_ohler/dharnet_m/cortexomics/contaminants/mm10-tRNAs.fa

ln -fs $GTF
ln -fs $trfasta
ln -fs $translations
ln -fs $dpdir
ln -fs $rrna_fasta
ln -fs $tRNAfasta

GTF=$(basename $GTF)
trfasta=$(basename $trfasta)
translations=$(basename $translations)
dpdir=$(basename $dpdir)
rrna_fasta=$(basename $rrna_fasta)
tRNAfasta=$(basename $tRNAfasta)

python ${dpdir}/filter_gencode_transcript.py $GTF $trfasta $translations

python ${dpdir}/select_rrna_from_ncrna.py $rrna_fasta ${GTF%.gtf}.rrna.fa

python ${dpdir}/build_contaminant.py ${GTF%.gtf}.rrna.fa $tRNAfasta Mus_musculus_tRNA mm10_contaminant.fa

python ${dpdir}/find_transcript_synonyms.py ${trfasta%.fa}_filter.fa synonym.txt

mkdir -p StarIndex

mkdir -p StarIndex/contaminant

STAR --runThreadN 15 --runMode genomeGenerate --genomeDir StarIndex/contaminant --genomeFastaFiles mm10_contaminant.fa --genomeSAindexNbases 8 --genomeChrBinNbits 11
mkdir -p StarIndex/transcript
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir StarIndex/transcript/ --genomeFastaFiles ${trfasta%.fa}_filter.fa --genomeSAindexNbases 11 --genomeChrBinNbits 12
