

import numpy as np 
import pandas as pd
import pyranges as pr
from threading import Thread

from ipdb import set_trace

import matplotlib.pyplot as plt

from itertools import * 
from pathlib import Path

import multiprocessing

  # def __init__(self, b=None):
  #   Base.__init__(self)
  #   if b:
  #     self.__dict__.update(b.__dict__) # copy instance variables
  #   self.y = "some more variables in derived class"

class ExtPyRanges(pr.PyRanges):
	def __init__(self,tocast):
		pr.PyRanges.__init__(self)
		self.__dict__.update(tocast.__dict__)
	# def __init__(self,*argv,**kwargs):
	# 	coerce = kwargs.get('coerce',False)
	# 	print(coerce)
	# 	if not coerce:
	# 		pr.PyRanges(*argv,**kwargs)
	# 	else:
	# 		print('coercing')
	# 		self.coerce=None
	# 		pr.PyRanges.__init__(self)
	# 		self.__dict__.update(coerce.__dict__)	
	def move(self,dist):
		self.Start += np.where(self.Strand!='-',dist,-dist)
		self.End += np.where(self.Strand!='-',dist,-dist)
		return self
		#
	def cast(cls,ranges: pr.PyRanges) -> 'ExtPyRanges':
		print(type(cls))
		print(cls)
		outranges = cls.__init__(ranges)
		outranges.__dict__.update(ranges.__dict__)
		return outranges
#
#
from Bio import SeqIO
import re 

if __name__ == '__main__':


	bam = pr.read_bam('pipeline/star/data/E13_ribo_1/E13_ribo_1.bam')
	bam = bam[bam.Chromosome]
	#gtf = pr.read_gtf('pipeline/my_gencode.vM12.annotation.gtf')
	pr.readers.read_gtf_restricted('pipeline/my_gencode.vM12.annotation.gtf',annotation='start_codon')

	mybamheadext = ExtPyRanges(mybamhead)

	fasta = Path('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')

	
	bam.fp = 'NNNN'
	bam.tp = 'NNNN'

	processes = [Thread(target = lambda x: results.append(x) ,args=(i,)) for i in range(3) ]
	
	[p.start() for p in processes]
	results

	[len(record) for record in SeqIO.parse(fasta, "fasta")]

	for record in SeqIO.parse(fasta, "fasta"):
		continue



	print(record.id)
	reads = bam[bam.Chromosome==record.id] 
	starts = reads.Start.astype(int)
	ends = reads.End.astype(int)
	isneg = reads.Strand == '-'
	bam[bam.Chromosome==record.id].fp = (record.seq[e-2:e+2].reverse_complement() if st else record.seq[s-2:s+2] for s,e,st in zip(starts,ends,isneg))
	bam[bam.Chromosome==record.id].tp = (record.seq[s-2:s+2].reverse_complement() if st else record.seq[e-2:e+2] for s,e,st in zip(starts,ends,isneg))

	with ProcessPoolExecutor(num_processes) as p_exec:
  		results=p_exec.map(add1, [1,2,3])




	############



	record = next(SeqIO.parse(fasta,"fasta"))
	record[reads[1]]

	reads.iloc[1]

	f = plt.figure()
	plt.plot(range(10), range(10), "o")
	plt.show()
	f.savefig("foo.svg", bbox_inches='tight')
	print(pathlib.Path('foo.svg').absolute())




#functions that I use a lot - resize, shift, countOverlaps (with different modes)
#
#functions that I use a lot - 'iscompatible with splicing'