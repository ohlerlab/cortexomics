import os
import re
import numpy as np
import pandas as pd
from inspect import getmembers, isfunction
import glob
import matplotlib.pyplot as plt
matplotlib.use('Agg')

FOLDCHANGETHRESH=0.1
PVALTHRESH=0.05

[o for o in getmembers(pd) if isfunction(o[1])]
[o for o in getmembers(os.path) if isfunction(o[1])]

from sklearn.mixture import GaussianMixture


root = os.path.expanduser('~/projects/cortexomics/')

riboseqresglob = os.path.join(root,'exploration','tables','riboseqres*')
riboseqresfiles =  glob.glob(riboseqresglob)

riboresfile=riboseqresfiles[0]

def get_te_tbl(riboresfile):
	time = re.findall('_(.*?).txt',os.path.basename(riboresfile))[0]
	riboseqtbl=pd.read_table(riboresfile,usecols=[0,3,6])
	riboseqrename=dict(zip(list(riboseqtbl.columns),['feature_id','adj_p_value','log2fc']))
	riboseqtbl=riboseqtbl.rename(columns=riboseqrename)
	riboseqtbl['time']=time
	riboseqtbl['assay']='TE'
	return riboseqtbl

all_te_data = pd.concat(map(get_te_tbl,riboseqresfiles))

all_te_data = all_te_data.groupby('feature_id')
#filter the data for significant fold changes
def rowfilt(x,FOLDCHANGETHRESH=FOLDCHANGETHRESH,PVALTHRESH=PVALTHRESH): 
	return ((x['log2fc']>FOLDCHANGETHRESH) & (x['adj_p_value']<PVALTHRESH)).any()
all_te_data = all_te_data.filter(rowfilt)
#now spread our data
all_te_data_spread = all_te_data.pivot(index='feature_id',columns='time',values='log2fc')


gmm_list = [ GaussianMixture(n_components=i,covariance_type='diag').fit(all_te_data_spread)  for i in range(1,12) ]
plt=pd.Series([gmm.bic(all_te_data_spread) for gmm in gmm_list]).plot()
plt.savefig('myfile.pdf')

gmm.predict_proba(all_te_data_spread)
np.log(gmm.predict_proba(all_te_data_spread)).sum()

dir(gmm)


all_te_data_spread

all_te_data_spread


from sklearn.naive_bayes import GaussianNB

clf = GaussianNB()
clf.fit(all_te_data)

dir(GaussianNB)

GaussianNB(priors=None)
print(clf.predict([[-0.8, -1]]))
[1]
clf_pf = GaussianNB()
clf_pf.partial_fit(X, Y, np.unique(Y))
GaussianNB(priors=None)
print(clf_pf.predict([[-0.8, -1]]))