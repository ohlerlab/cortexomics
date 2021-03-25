# mv transcript_star/data/ribo_0h/ribo_0h.bam{,.bak}
# mv transcript_star/data/rnaseq_0h/rnaseq_0h.bam{,.bak}
# mv gencode.v24lift37.annotation.cdsrange.txt{,.bak}
#samtools index transcript_star/data/ribo_0h/ribo_0h.bam.bak 
#samtools index transcript_star/data/rnaseq_0h/rnaseq_0h.bam.bak
#mv gencode.v24lift37.annotation.transcript.fa{.bak}

testtr=$(samtools view transcript_star/data/ribo_0h/ribo_0h.bam.bak | cut -f 3 | uniq -c | head -n 1000 | awk '$1 > 1000' | head -n 2 | tail -n 1 | awk '{print $2}'  )
echo $testtr
head gencode.v24lift37.annotation.cdsrange.txt.bak

grep -e "$testtr" gencode.v24lift37.annotation.cdsrange.txt.bak

grep -e "$testtr" gencode.v24lift37.annotation.cdsrange.txt.bak > gencode.v24lift37.annotation.cdsrange.txt
head gencode.v24lift37.annotation.cdsrange.txt


samtools view -c  transcript_star/data/ribo_0h/ribo_0h.bam.bak $testtr 
samtools view transcript_star/data/ribo_0h/ribo_0h.bam.bak $testtr -bh > transcript_star/data/ribo_0h/ribo_0h.bam
samtools index transcript_star/data/ribo_0h/ribo_0h.bam
samtools view -c  transcript_star/data/ribo_0h/ribo_0h.bam $testtr

samtools view -c   transcript_star/data/ribo_0h/ribo_0h.bam.bak $testtr
samtools view transcript_star/data/rnaseq_0h/rnaseq_0h.bam.bak $(head  gencode.v24lift37.annotation.cdsrange.txt -n 1 | cut -f 1) -bh > transcript_star/data/rnaseq_0h/rnaseq_0h.bam
samtools index transcript_star/data/rnaseq_0h/rnaseq_0h.bam

sed -e 's/\(^>.*$\)/#\1#/' gencode.v24lift37.annotation.transcript.fa.bak  | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep $testtr -A 1 > gencode.v24lift37.annotation.transcript.fa

grep -e "$testtr" gencode.v24lift37.annotation.transcript.fa 

snakemake -f salmon_quant/rnaseq_0h/quant.sf -p
grep -e $testtr salmon_quant/rnaseq_0h/quant.sf

snakemake -f  ribomap/ribo_0h/ribo_0h.ribomap

samtools view -c transcript_star/data/ribo_0h/ribo_0h.bam $testtr


mv transcript_star/data/ribo_0h/ribo_0h.bam{.bak,}
mv transcript_star/data/rnaseq_0h/rnaseq_0h.bam{.bak,}
mv gencode.v24lift37.annotation.cdsrange.txt{.bak,}
samtools index transcript_star/data/ribo_0h/ribo_0h.bam 
samtools index transcript_star/data/rnaseq_0h/rnaseq_0h.bam
mv gencode.v24lift37.annotation.transcript.fa{.bak,}


snakemake -f  ribomap/ribo_0h/ribo_0h.ribomap 2> err > out

ribobam=transcript_star/data/ribo_0h/ribo_0h.bam
rnabam=transcript_star/data/rnaseq_0h/rnaseq_0h.bam

testseq=CCTTTCCTGACTCAACAGAGGCCTTTGACC
samtools view $ribobam | grep -m1 -e $testseq > testsamline
head testsamline
testtr=ENST00000628980

grep -e CCTTTCCTGACTCAACAGAGGCCTTTGACC -m1 out

cdsrangefl=gencode.v24lift37.annotation.cdsrange.txt
grep -e $testtr -m1 out $cdsrangefl

#out of bound line
 grep -e $testseq -m1 out
 #id of the sequence it'smapping to

 73870

#Which seems to point to this
#ENST00000462316

which has 571 characters in the RNA fasta