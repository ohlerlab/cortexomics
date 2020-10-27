import sys
import os
import string
import sklearn
import numpy as np
import numpy
import pandas as pd
import_folder = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/python/'
sys.path = [import_folder]+list(set(sys.path)-set(import_folder)) # this tells python to look in `import_folder` for imports
from scipy.sparse import find
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.preprocessing import normalize
from sklearn.preprocessing import LabelEncoder
import arraychoose
def sparse_choosecols(spmat):
	sptmatnorm = normalize(spmat,axis=1,norm='l1')
	spcols, sprows, spvals = find(sptmatnorm.transpose())
	return arraychoose.choose_colssp(sprows, spcols, spvals)


def readcount2TPM(transcript_readcount, transcript_CDS_len):
	transcript_TPMs = {}
	total_reads = 0
	for key, value in transcript_readcount.items():
		transcript_TPMs[key] = transcript_readcount[key]/transcript_CDS_len[key]*1000
		total_reads += transcript_TPMs[key]
	for key, value in transcript_TPMs.items():
		transcript_TPMs[key] = transcript_TPMs[key]/total_reads*1000000
	return transcript_TPMs

# M step
def M_step(transcript_readcount, transcript_CDS_len, transcript_TPMs):
	transcript_TPMs_new = readcount2TPM(transcript_readcount, transcript_CDS_len)
	TPM_diff = {}
	for key, value in transcript_TPMs.items():
		if value != 0:
			TPM_diff[key] = (transcript_TPMs_new[key] - value)
		else:
			TPM_diff[key] = 0
	return [transcript_TPMs_new, TPM_diff]

# Output TPM_diff and TPM_diff_relative of each loop
def write_runlog(f_out, loopnum, transcript_TPMs, TPM_diff):
	TPM_list = [value for key, value in transcript_TPMs.items()]
	f_out.write('Loop ' + str(loopnum) + ' TPM:\t' + '\t'.join(map(str, TPM_list)) + '\n')
	f_out.write('Loop ' + str(loopnum) + ' TPM_diff:\t' + '\t'.join(map(str, TPM_diff)) + '\n')
	TPM_diff_square = sum([numpy.power(TPM_diff[i],2) for i in range(0, len(TPM_diff))])
	TPM_diff_mean_square = TPM_diff_square/len(TPM_diff)
	f_out.write('Loop ' + str(loopnum) + ' quadratic sum & mean square:\t' + str(TPM_diff_square) + '\t' + str(TPM_diff_mean_square) + '\n')

def write_runlog_1(f_out, f_summary_out, loopnum, transcript_TPMs, TPM_diff,transcript_readcount):
	for key, value in transcript_TPMs.items():
		f_out.write(key + '\t' + str(value) + '\t' + str(TPM_diff[key]) + '\t'+ str(transcript_readcount[key]) + '\n')
	TPM_diff_square = sum([numpy.power(value,2) for key, value in TPM_diff.items()])
	TPM_diff_mean_square = TPM_diff_square/len(TPM_diff)
	f_out.write('Loop ' + str(loopnum) + ' quadratic sum & mean square:\t' + str(TPM_diff_square) + '\t' + str(TPM_diff_mean_square) + '\n')




# Main

# path_out_dir = sys.argv[3]


# numloops = int(sys.argv[4])

# numloops=10
# path_ref='/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa'
# print('load read data')

# df = pd.read_table(
# 	'pipeline/deepshapebamdata/E13_ribo_1.bam.reformat.absolute',
# 	names=['read_name','chrom','strand','pos','cdspos','tr_id','g_id']
# 	)


# df = df[['read_name','tr_id']]

# #assign iniitial tpms
# tpms = df[['tr_id']].drop_duplicates().assign(TPM=1000000/len(df.tr_id.unique()))

# #put these on the reads
# df = df.merge(tpms)

# #get read total tpms
# rptmsums = df.groupby('read_name')['TPM'].apply(lambda x: x/sum(x))

# rptmsums = rptmsums.reset_index(drop=False)

# df.merge(rptmsums)

# tpms['TPM'] = 1000000/len(df.tr_id.unique())
# df['TPM'] = 1000000/len(df.tr_id.unique())


def E_step_sp(readgraph,trenc):
		#normalize rows (mapping TPMS), then sum over columns(transcript readcounts)
		mat_readcount = normalize(readgraph,axis=1,norm='l1').sum(axis=0)
		mat_readcount = pd.Series(np.squeeze(np.array(mat_readcount)))		
		mat_readcount.index=pd.Series(trenc.inverse_transform(range(0,readgraph.shape[1])))
		return mat_readcount

def ribo_EM(readdf, transcript_CDS_len , numloops, verbose=True):
	# try:
		# path_out_dir = 'emondiskcheck'
	# except FileExistsError:
		# pass
	# f_summary_out = open(path_out_dir+"/summary.txt", 'w')
	readdf.sort_values(['read_name','tr_id'],inplace=True)

	trenc = LabelEncoder()
	trenc.fit(readdf.tr_id)
	renc = LabelEncoder()
	renc.fit(readdf.read_name)
	rlabs = renc.transform(readdf.read_name)
	tlabs = trenc.transform(readdf.tr_id)
	readdf = readdf[['read_name','tr_id']]

	readgraph = csr_matrix(
		([1]*readdf.shape[0],(rlabs, tlabs)),
		shape=(len(renc.classes_), len(trenc.classes_))
	)

	transcript_TPMs = readcount2TPM(readdf.tr_id.value_counts(),transcript_CDS_len)

	for loopnum in range(0, numloops):
		print("loop: " + str(loopnum))

		# E step: allot reads to transcripts according to transcript TPM distribution.
		# transcript_readcount = E_step(readdf[['read_name','tr_id']], transcript_TPMs)
		transcript_readcount = E_step_sp(readgraph, trenc)
		# M step: update transcript TPM distribution according to reads number.
		[transcript_TPMs, TPM_diff] = M_step(transcript_readcount, transcript_CDS_len, transcript_TPMs)		#now reasign the tpms to the data of the sparse matrix
		readgraph.data = pd.Series(transcript_TPMs)[pd.Series(trenc.inverse_transform(readgraph.indices))]

		#write_runlog(f_out, loopnum, transcript_TPMs, TPM_diff)
		# f_out = open(path_out_dir + "/runlog_"+str(loopnum)+".txt", 'w')
		#if(loopnum(numloops-1)): 
		# write_runlog_1(f_out, f_summary_out, loopnum, transcript_TPMs, TPM_diff, transcript_readcount)
		diffdf = pd.concat([pd.Series(TPM_diff).abs(),pd.Series(transcript_TPMs)],axis=1)
		diffdf.columns=['tpm_diff','TPM']
		diffdf['fracchange'] = diffdf.tpm_diff/diffdf.TPM
		diffdf = diffdf.query('TPM>1')
		if((diffdf.fracchange.max()<0.01)&(loopnum!=0)&(verbose)): 
			print('minimum change less than 1% - convergence')
			print(diffdf.sort_values(['fracchange'],ascending=False).head(5))
			break
		else:
			print(diffdf.sort_values(['fracchange'],ascending=False).head(5))

	print('assign reads using the TPMs')
	trinds = sparse_choosecols(readgraph)

	sampreads=pd.concat([
		pd.Series(renc.inverse_transform(range(0,trinds.shape[0]))),
		pd.Series(trenc.inverse_transform(trinds))
		],axis=1)
	sampreads.columns = ['read_name','tr_id']

	return sampreads, transcript_TPMs, TPM_diff, transcript_readcount

if False:
	import sklearn
	from scipy.sparse import find
	from scipy.sparse import csr_matrix
	from scipy.sparse.csgraph import connected_components

	testdf = bamdf

	sptmat = csr_matrix(np.matrix([[0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0],
        [1, 1, 1, 0, 0],
        [0, 0, 0, 0, 1]]))
	testmat = sptmat.todense()
	tvect = np.array(range(0,5)).reshape([1,5])

	trcongraph = testmat.transpose().dot(testmat)
	n_comps, labels = connected_components(trcongraph,return_labels=True)
	#normalize(testmat,axis=1,norm='l1') makes rows sums add to 1
	
	testdf = bamdf
	trenc = LabelEncoder()
	trenc.fit(testdf.tr_id)
	renc = LabelEncoder()
	renc.fit(testdf.read_name)

	readgraph = csr_matrix(
		(
			[1]*testdf.TPM,
			(renc.transform(testdf.read_name),trenc.transform(testdf.tr_id))),
		shape=(len(renc.classes_),len(trenc.classes_)))


	trinds = sparse_choosecols(readgraph)

	choicedf=pd.concat([
		pd.Series(renc.inverse_transform(range(0,trinds.shape[0]))),
		pd.Series(trenc.inverse_transform(trinds))
		],axis=1)
	choicedf.columns = ['read_name','tr_id']

	bamdf.merge(choicedf,how='inner')

	counts = pd.Series(np.squeeze(np.array(trinds))).value_counts()
	counts.index = trenc.inverse_transform(counts.index)

	readgraph = csr_matrix(
		(
			[1]*testdf.shape[0],
			(renc.transform(testdf.read_name),trenc.transform(testdf.tr_id))),
		shape=(len(renc.classes_),len(trenc.classes_)))
	trcongraph = readgraph.transpose().dot(readgraph)
	n_comps, labels = connected_components(trcongraph,return_labels=True)
	print(n_comps)

	labelcounts = pd.Series(labels).value_counts()
	multlabels = list(labelcounts[labelcounts>1].index)
	nreadgraph = normalize(readgraph,axis=1,norm='l1')



	testdf = bamdf
	trenc = LabelEncoder()
	trenc.fit(testdf.tr_id)
	renc = LabelEncoder()
	renc.fit(testdf.read_name)

	readgraph = csr_matrix(
		(
			[1]*testdf.TPM,
			(renc.transform(testdf.read_name),trenc.transform(testdf.tr_id))),
		shape=(len(renc.classes_),len(trenc.classes_)))


	trinds = sparse_choosecols(readgraph)
	trinds


	bamdf = bamdf.merge(pd.Series(transcript_TPMs,name='TPM'),left_on='tr_id',right_index=True)

def weightedsamplenp(bamdf):
	choices = []
	for i, g in bamdf.groupby('read_name')['TPM']:
		g = np.cumsum(np.array(g))
		g = g/g[-1]
		choice=np.searchsorted(g,np.random.uniform(0,1))
		choices.append(i[choice])
	return choices

# %prun weightedsamplenp(bamdf.head(100000))



def weightedsamplecython(bamdf):
	testdf = bamdf
	trenc = LabelEncoder()
	trenc.fit(testdf.tr_id)
	renc = LabelEncoder()
	renc.fit(testdf.read_name)

	readgraph = csr_matrix(
		(
			[1]*testdf.TPM,
			(renc.transform(testdf.read_name),trenc.transform(testdf.tr_id))),
		shape=(len(renc.classes_),len(trenc.classes_)))


	trinds = sparse_choosecols(readgraph)
	return trinds

# %prun -l 10 weightedsamplenp(bamdf.head(1000000))
# %prun -l 10 weightedsamplecython(bamdf.head(1000000))

# testdf = bamdf
# trenc = LabelEncoder()
# trenc.fit(testdf.tr_id)
# renc = LabelEncoder()
# renc.fit(testdf.read_name)

# readgraph = csr_matrix(
# 	(
# 		[1]*testdf.TPM,
# 		(renc.transform(testdf.read_name),trenc.transform(testdf.tr_id))),
# 	shape=(len(renc.classes_),len(trenc.classes_)))


# trinds = sparse_choosecols(readgraph)
# trinds




