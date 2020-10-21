import pysam
import sys
import textwrap
import multiprocessing
import numpy as np
import pandas as pd
import pybedtools as pbt
from textwrap import wrap
from scipy import sparse
# from scikit_ribo import asite_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel, RFECV
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score

# from sklearn import preprocessing, svm, tree
# from sklearn.metrics import roc_curve, auc
# from sklearn.preprocessing import label_binarize
# from sklearn.multiclass import OneVsRestClassifier
# from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.feature_selection import SelectFromModel, RFECV
# import glmnet_py.dataprocess as dataprocess

from src.Ribotransformer.ribo_EM import ribo_EM
from src.Ribotransformer.asite_predict import PredictAsite
import argparse
import pandas as pd

# Initialize abundance of all transcript to 1/#transcripts


# # read in the sequence
# cdsfasta = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa'
# # transcript bam
# bam = 'pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.sort.bam'


#Note - this function needs a fasta file with teh following format - (coordinates 1 indexed)
#>ENSMUST00000070533.4|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|Xkr4-001|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634|
def transcript_initialize(cdsfasta):

    f_ref = open(cdsfasta)
    # transcript_TPMs = {}
    transcript_CDS_len = {}
    transcript_cdsstart = {}
    transcript_cdsend = {}
    transcript_gene = {}
    exonseq = {}
    for line in f_ref:
        if line.startswith('>'):
            line = line.lstrip('>').rstrip().rstrip('|')
            linearray = line.split('|')
            transID = linearray[0]
            transcript_gene[transID] = linearray[1]
            exonseq[transID] = ''
            # transcript_TPMs[transID] = 1
            for a in linearray:
                if a.startswith("CDS:"):
                    start, stop = a.lstrip("CDS:").split('-')
                    transcript_CDS_len[transID] = int(stop) - int(start) + 1
                    transcript_cdsstart[transID] = int(start)-1
                    transcript_cdsend[transID] = int(stop)-1
        else:
            exonseq[transID] += line.strip()
    # transNum = len(transcript_TPMs)
    # for key, value in transcript_TPMs.items():
        # transcript_TPMs[key] = 1000000.0/transNum
    f_ref.close()

    # return [transcript_TPMs, transcript_CDS_len, transcript_cdsstart, transcript_cdsend, exonseq]
    return [transcript_CDS_len, transcript_cdsstart, transcript_cdsend, exonseq, transcript_gene]


def make_bamDF(bam, mapq=-1, minRL=20, maxRL=35):
    # create a template/header from the input bam file
    inBam = pysam.AlignmentFile(bam, "rb")
    # read a bam file and extract info
    # cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
    cigar_to_exclude = set([1, 2, 3, 4, 5])  # set(['I','D','S','H'])
    i = 0
    readdf = []
    for read in inBam.fetch():
        i += 1
        if(i==10000):break
        cigars = set([c[0] for c in read.cigartuples])
        if read.mapping_quality > mapq and \
                minRL <= read.query_length <= maxRL and \
                not cigars.intersection(cigar_to_exclude) and \
                not read.is_reverse:
                # The start is 0 indexed, the end is 1 indexed! Make them
                # both 0 indexed, so if I index a string with the trseq using the cols
                # I get the first and last bp of the read
            readdf.append((read.query_name, read.reference_name.split(
                '|')[0], read.reference_start, read.reference_end-1, read.query_length))

    readdf = pd.DataFrame.from_records(
        readdf, columns=['read_name', 'tr_id', 'start', 'end', 'read_length'])

    return readdf


def add_seqinfo(bamdf, transcript_cdsstart, exonseq, trlengths):
    # add sequence
    uniquereadpos = bamdf[['tr_id', 'start', 'end']].drop_duplicates()
    bamtrlens = trlengths[uniquereadpos.tr_id].reset_index(drop=True)
    uniquereadpos = uniquereadpos.reset_index(
        drop=True)[uniquereadpos.end.add(1).reset_index(drop=True) != bamtrlens]
    uniquereadpos = uniquereadpos[uniquereadpos.start != 0]
    # get the edge sequences of our unique read pos
    edgeseq = [(exonseq[tr][st-1], exonseq[tr][st], exonseq[tr][e], exonseq[tr][e+1]) for
               i, tr, st, e in uniquereadpos.itertuples()]
    # concat these to the unique read df
    seqcols = ['nt_-1', 'nt_0', 'nt_n-1', 'nt_n']
    edgeseq = pd.DataFrame.from_records(edgeseq, columns=seqcols)
    uniquereadpos = pd.concat(
        [uniquereadpos.reset_index(drop=True), edgeseq], axis=1)
    # add offsets

    cdsstarts = pd.Series(transcript_cdsstart)[
        uniquereadpos.tr_id].reset_index(drop=True)
    uniquereadpos['cdspos'] = (uniquereadpos.start - cdsstarts).values
    # uniquereadpos.assign(cdspos =  (uniquereadpos.start - cdsstarts).values)
    # uniquereadpos = uniquereadpos.assign(cdspos=lambda x:x.start - pd.Series(transcript_cdsstart)[x.tr_id].reset_index(drop=True))

    uniquereadpos['cdsendpos'] = uniquereadpos.start - \
        (pd.Series(transcript_cdsend) -
         2)[uniquereadpos.tr_id].reset_index(drop=True)

    uniquereadpos['5_offset'] = uniquereadpos['cdspos'] % 3
    uniquereadpos['3_offset'] = (
        uniquereadpos['cdspos']+(uniquereadpos['end']-uniquereadpos['start']+1)-1) % 3

    # merge in to the nonredundant read df
    bamdf = pd.merge(
        bamdf, uniquereadpos,
        left_on=['tr_id', 'start', 'end'],
        right_on=['tr_id', 'start', 'end']
    )
    return bamdf


def makePredTraining(bamdf):

    training = bamdf[
        (bamdf['cdspos'] <= 0) &
        (bamdf['cdspos'] >= -bamdf['read_length'])
    ]

    training['asite'] = 3 - training['cdspos']

    training = (training.query('asite >= 9').
                query('asite <= 18').
                query('asite >= (read_length / 2 - 1)')
                )
    training = training[['read_length', 'nt_-1', 'nt_0',
                         'nt_n-1', 'nt_n', '5_offset', '3_offset', 'asite']]
    return training



def get_coddf(transcript_cdsstart, transcript_cdsend, exonseq, trlength ):
    coddfs = []
    # def tcode():
    for tr in trstouse:
        futrcodpos = list(reversed(range(transcript_cdsstart[tr]-3, -1, -3)))
        cdscodpos = list(range(transcript_cdsstart[tr], transcript_cdsend[tr], 3))
        tputrcodpos = list(range(transcript_cdsend[tr]+1, trlength[tr]-2, 3))
        codstarts = futrcodpos + cdscodpos + tputrcodpos
        codidx = list(reversed(range(-1, -len(futrcodpos)-1, -1)))
        codidx += list(range(0, len(cdscodpos)))
        codidx += list(range(len(cdscodpos), len(cdscodpos)+len(tputrcodpos)))
        cods = wrap(exonseq[tr][codstarts[0]:codstarts[-1]+3], 3)
        coddf = pd.DataFrame(zip(codstarts, codidx, cods))
        coddf.columns = ['start', 'codon_idx', 'codon']
        coddf['end'] = coddf['start']+3
        coddf['chrom'] = tr
        coddfs.append(coddf)

    # %prun -l 10 tcode()
    coddf = pd.concat(coddfs,axis=0)
    coddf['gene_strand']='+'
    coddf['gene'] = coddf['chrom']
    return coddfs



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input bam file",default = "pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.sort.bam")
    parser.add_argument("-f", help="Fasta file", default="ext_data/gencode.vM12.pc_transcripts.fa")
    parser.add_argument("-o", help="ribotrans_test", default="ribotranstest")
    args = parser.parse_args()
    
    bam = args.i
    cdsfasta = args.f

    transcript_CDS_len, transcript_cdsstart, transcript_cdsend, exonseq, transript_gene = transcript_initialize(
        cdsfasta)
    trlength = pd.Series(dict([(s, len(exonseq[s])) for s in exonseq]))

    bamdf = make_bamDF(bam)

    #count only reads overlapping the cds
    bamdf = (bamdf.
        merge(pd.Series(transcript_cdsstart,name='cdsstart'),left_on='tr_id',right_index=True).
        query('start > cdsstart').
        drop(['cdsstart'],axis=1).
        merge(pd.Series(transcript_cdsend,name='cdsend'),left_on='tr_id',right_index=True).
        query('end < cdsend').
        drop(['cdsend'],axis=1)
    )


    print('bam loaded')

    sampreads, transcript_TPMs, TPM_diff, transcript_readcount = ribo_EM(
        bamdf[['read_name', 'tr_id']],
        transcript_CDS_len,
        numloops=100
    )

    sampreads = sampreads[['read_name', 'tr_id']]
    sampreadsfull = sampreads.merge(bamdf, how='inner')

    sampreadsfull = add_seqinfo(sampreadsfull, transcript_cdsstart, exonseq, trlength)

    gtpms = pd.concat([pd.Series(transript_gene,name='gene'),pd.Series(dict(transcript_TPMs),name='TPM')],axis=1)
    trstouse = gtpms.sort_values('TPM',ascending=False).groupby('gene').head(1)['TPM'].index.to_series()

    coddf = get_coddf(transcript_cdsstart, transcript_cdsend, exonseq, trlength)
    coddf = coddf.merge(pd.Series(transcript_TPMs,name='TPM'),left_on='chrom',right_index=True).assign(pair_prob=0)

    training = makePredTraining(sampreadsfull[sampreadsfull.tr_id.isin(trstouse)])

    model = PredictAsite(training, sampreadsfull, 'rf', False,cdsIdxDf = coddf)
    print('fit A sites')
    model.rfFit()
    model.rfPredict()
    print("[execute]\tlocalize the a-site codon and create coverage df", file=sys.stderr)
    model.cds['gene_strand'] = '+'
    model.cds['strand'] = '+'
    model.cds = model.cds.rename(columns={'tr_id': 'chrom'})
    dataFrame = model.recoverAsite()

    dataFrame.to_csv(args.o+'.csv')


# #looks like a lot of reads are extra-ORF, but this isn't due to frame shifting (chceck)


# tpmLB = 1
# unmap=None
# out=None
# # start model fitting
# print("[execute]\tStart the modelling of translation efficiency (TE)", file=sys.stderr)
# mod = ModelTE(dataFrame, unmap, out, tpmLB)
# print("[execute]\tLoading data", file=sys.stderr)
# mod.loadDat()
# print("[execute]\tFiltering the df", file=sys.stderr)
# mod.filterDf()
# print("[execute]\tScaling the variables", file=sys.stderr)
# mod.varScaling()
# print("[execute]\tFitting the GLM", file=sys.stderr)
# X, y, offsets, numCodons, numGenes, varsNames = mod.glmnetArr()
# mod.glmnetFit(X, y, offsets, numCodons, numGenes, varsNames, lambda_min = 0.13)
# mod.codonBetas.to_csv('pipeline/trtestcodondts.csv')

# fplen = 0
# ncods = 4
# tplen = 3
# transcript_cdsstart['foo'] = fplen
# transcript_cdsend['foo'] = transcript_cdsstart['foo']+(3*(ncods-1))+2
# trlength['foo'] = transcript_cdsend['foo']+1+tplen
# exonseq['foo'] = ('A'*fplen)+('C'*(3*ncods))+('G'*tplen)
# tr = 'foo'
# exonseq.pop('foo', None)

# transcript_cdsstart['foo']
# exonseq[tr][transcript_cdsend['foo']+1]



# ################################################################################
# # now do asite model
# ################################################################################


# # model = PredictAsite(training, bamdf, 'rf', False)

# model.rfFit()

# model.rfPredict()

# readdf = model.cds

# readdf.asite.value_counts()

# readdf.codon = readdf.cdspos + readdf.asite

# readdf = readdf.query('cdsendpos<=0').query('read_length > ( - cdspos)')

# readdf = readdf.assign(codidx=(np.floor((readdf.cdspos + readdf.asite) / 3)))

# countdf = readdf.groupby(['tr_id', 'codidx'], observed=True).size(
# ).reset_index().rename(columns={0: 'ribo_count'})

# # for each tr, we want to go back as many CODONS as we can.
# cdsstarts = pd.Series(transcript_cdsstart)
# cdsends = pd.Series(transcript_cdsend)
# cdslens = pd.Series(transcript_CDS_len)
# fpcodonbps = (np.floor((cdsstarts-1)/3)*3)
# tpcodonbps = (3*np.floor((trlens - cdsends) / 3))

# testtr = trlens.index[0]

# codondfs = []
# for tr in readdf['tr_id'].unique()[0:4]:
#     codons = textwrap.wrap(
#         exonseq[tr][int((cdsstarts[tr] - fpcodonbps[tr]) - 1)                    :int(cdsends[tr]+tpcodonbps[tr])],
#         3
#     )
#     lcodon = -int(fpcodonbps[tr]/3)
#     rcodon = int((cdslens[tr]+tpcodonbps[tr])/3)
#     indices = list(range(lcodon, rcodon))
#     #
#     codondf = pd.concat([pd.Series(codons), pd.Series(indices)], axis=1)
#     codondf = codondf.assign(tr_id=tr)
#     codondf.columns = ['codon', 'codidx', 'tr_id']
#     codondfs.append(codondf)

# codondf = pd.concat(codondfs, axis=0)

# codondf = codondf.merge(countdf, how='left')

# codondf['ribo_count'] = codondf['ribo_count'].fillna(0)

# numloops = 20
# path_out_dir = 'emtest'


# # now run EM to distribute the reads appropriately

# bamdf

#      query('asite >= 9').
#             query('asite <= 18').
#             query('asite >= (read_length / 2 - 1)')
