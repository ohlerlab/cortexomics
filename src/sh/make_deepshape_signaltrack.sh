set -x
echo "tpms to reads"
join -1 6 -2 1 <(sort -k 6 $1) <(cat $2 | cut -f 1,2 | sort -k 1) | sort -k 2 > test.bam.withtpm
echo "sum tpms"
datamash -t' ' -g 2 sum 8 <  test.bam.withtpm > test.bam.readtpms
echo "distribute read signal"
join -1 2 -2 1 test.bam.withtpm test.bam.readtpms | awk '{print $0,$8/$9}' | datamash -t' ' -g 2,6 sum 10 > tr_pos_ribosig