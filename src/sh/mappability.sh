#!/bin/bash - 
#===============================================================================
#
#          FILE: mappa_cacl.sh
# 
#         USAGE: ./mappa_cacl.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Vang Le (vql AT rn.dk) 
#  ORGANIZATION: 
#       CREATED: 26/02/16 01:49
#      REVISION:  ---
#      Reference: http://wiki.bits.vib.be/index.php/Create_a_mappability_track
#===============================================================================

set -o nounset                              # Treat unset variables as an error

pref="rn6.softmask.all" # should be changed
reference="rn6.softmask.all.fa" # should be changed
idxpref="${pref}_index"
thr=20; # use 20 cores, change it to the number you want to use
outmappa="rn6.mappa.tsv" # should be changed
#make index for the genome
gem-indexer -T ${thr} -c dna -i ${reference} -o ${idxpref}

# The following calculates index and creates mappability tracks with various kmer lengths.
# it may take long time to finish.
# choose the right kmer length for you.
for kmer in 20 27 31; do

  # compute mappability data
   gem-mappability -T ${thr} -I ${idxpref}.gem -l ${kmer} -o ${pref}_${kmer}
  mpc=$(gawk '{c+=gsub(s,s)}END{print c}' s='!' ${pref}_${kmer}.mappability)
  echo ${pref}_${kmer}"\t"$mpc >> $outmappa
  # convert results to wig and bigwig
  gem-2-wig -I ${idxpref}.gem -i ${pref}_${kmer}.mappability -o ${pref}_${kmer}
  wigToBigWig ${pref}_${kmer}.wig ${pref}.sizes ${pref}_${kmer}.bw

done