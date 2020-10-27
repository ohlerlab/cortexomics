import pyximport
pyximport.install(reload_support=True)
#
import numpy as np
import sys
from sklearn.preprocessing import normalize 
from scipy.sparse import csr_matrix 
from scipy.sparse import find
import_folder = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/python/'
sys.path = [import_folder]+list(set(sys.path)-set(import_folder)) # this tells python to look in `import_folder` for imports
import arraychoose
#
sptmat = csr_matrix(np.matrix([[1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 1, 1, 1, 1],
        [0, .2, 0, 0, 5]]))
#
def sparse_choosecols(spmat):
	sptmatnorm = normalize(spmat,axis=1,norm='l1')
	spcols,sprows,spvals = find(sptmatnorm.transpose())
	return arraychoose.choose_colssp(sprows,spcols,spvals)

sparse_choosecols(sptmat).shape


# ndarraychoose.__file__
# import subprocess
# from importlib import reload

# del sys.modules['ndarraychoose']
# ndarraychoose = reload(ndarraychoose)
# ndarraychoose.choose_on
