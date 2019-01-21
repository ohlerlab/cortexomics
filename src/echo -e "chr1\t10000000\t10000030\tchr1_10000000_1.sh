echo -e "chr1\t10000000\t10000030\tchr1_10000000_10000030" | bedtools getfasta -s -fi my_GRCm38.p5.genome.chr_scaff.fa -bed -  -name | 
	| sed '$!N;s/\n/ /' \
	| perl -lane '$i=0;while($i<=(length($F[1])-28)){print $F[0] , "$i\n", substr $F[1],$i,28 ; $i = $i +1}' \
	| perl -lan -F'_|\)|\(' -e 'if( /^>/){print $F[0],"_",$F[1]+$F[4],"_",$F[2]+$F[4]}else{print @F}' \
	| perl -lanpe 's/>/@/ ; s/^([^@]+)/\1\n+/; if($1){$a="I" x length($1); s/\+/+\n$a/}' #now in fastq format, with perfect sequencing quality