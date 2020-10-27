

STAR \
	--runMode genomeGenerate \
	--runThreadN 8 \
	--genomeFastaFiles $1 \
	--sjdbGTFfile $2 \
	--sjdbOverhang 20. \
	--genomeDir $3 \