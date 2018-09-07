import tensorflow as tf
import edward as edw
import edward as ed
import numpy as np
import pandas as pd
from edward.models import OneHotCategorical,Categorical, Dirichlet, InverseGamma, \
    MultivariateNormalDiag, Normal, ParamMixture, \
    Gamma, Mixture, DirichletProcess
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2

n_genes = 3
n_time = 5

# def define_model(n_genes,n_time):

# ribomsbase = tf.constant(todimn(100.,n_genes))




n_genes = 3 
n_time = 5

# vals = define_model(n_genes,n_time)


deg = tf.placeholder(tf.float32,shape=n_genes)

protein0 = tf.placeholder(tf.float32,shape=n_genes)

ribo = tf.placeholder(tf.float32,shape=(n_time,n_genes))

#define our protein
protein= [ (ribo[0,:] + (protein0 * deg))]
for  i in range(1,n_time):
	protein= protein + [ribo[len(protein)-1,:] + (protein[-1] * deg)]
protein = tf.stack(protein)

#define our 
riboreads = Normal(
	tf.stack([ribo,ribo],axis=2),
	tf.sqrt(tf.stack([ribo,ribo],axis=2))
)

#ms is similiar to riboreads for now - but with higher variance
ms = Normal(
	tf.stack([protein]*3,axis=2),
	tf.sqrt(tf.stack([protein]*3,axis=2))
)

vals = ( protein, riboreads, ms )
	# return( protein, riboreads, ms )

vars = tf.trainable_variables()
	
with tf.Session() as sess:

	ms_sample=sess.run(
		fetches = vals,
		feed_dict = {
			protein0:[10000]*n_genes,
			deg:[0.5]*n_genes,
			ribo: np.stack([[1]*n_time]*n_genes,axis=-1)
		}
	)


#This function will stack up the variables our model outputs for eas plotting
def stackarray(d2stack,varname,indexnamelist):
	if(len(d2stack.shape) != 3) : d2stack = np.expand_dims(d2stack,axis = -1)
	indexnamelist = indexnamelist[0:len(d2stack.shape)]
	# indexnamelist = indexnamelist[0:len(d2stack.shape)]

	index = pd.MultiIndex.from_product([range(s)for s in d2stack.shape], names=indexnamelist)
	# df = pd.DataFrame({'d2stack': d2stack.flatten()}, index=index)['d2stack']
	df = pd.DataFrame({'value': d2stack.flatten()}, index=index)
	df['var'] = varname
	return df

varnames = ['riboreads','protein','ms']
dfs = [stackarray(d,nm,varnames) for d,nm in zip(ms_sample,varnames)]

for i in range(0,len(dfs)):
	dfs[i]['var'] = varnames[i]

alldata = pd.concat(dfs,sort=False)

alldata.to_csv('tmp.tsv',sep='\t')

n_riboreps = 2
n_msreps = 3








# #from the linear mixed effects bit in the edward tutorials
# with tf.variable_scope('b'):
# 	deg = tf.get_variable('deg',[n_genes])
# 	protein0 = tf.sigmoid(tf.get_variable('protein0',[n_genes]))
# 	ribo = tf.get_variable('ribo',[n_time,n_genes])

# #define our protein
# protein= [ (ribo[0,:] + (protein0 * deg))]
# for  i in range(1,n_time):
# 	protein= protein + [ribo[len(protein)-1,:] + (protein[-1] * deg)]
# protein = tf.stack(protein)

# #define our 
# riboreads = Normal(
# 	tf.stack([ribo,ribo],axis=2),
# 	tf.sqrt(tf.stack([ribo,ribo],axis=2))
# )

# #ms is similiar to riboreads for now - but with higher variance
# ms = Normal(
# 	tf.stack([protein]*3,axis=2),
# 	tf.sqrt(tf.stack([protein]*3,axis=2))
# )

# vals = ( protein, riboreads, ms )
# 	# return( protein, riboreads, ms )

# vars = tf.trainable_variables()
	
# with tf.Session() as sess:

# 	ms_sample=sess.run(
# 		fetches = vals,
# 		feed_dict = {
# 			protein0:[10000]*n_genes,
# 			deg:[0.5]*n_genes,
# 			ribo: np.stack([[1]*n_time]*n_genes,axis=-1)
# 		}
# 	)




deg = tf.Variable(tf.float32)

protein0 = tf.Variable(tf.float32,shape=n_genes)

ribo = tf.Variable(tf.float32,shape=(n_time,n_genes))

#
riboreadvals = ms_sample[1]
msvals = ms_sample[2]
#define our approximations for inference

q_ribo = Normal(
    loc=tf.get_variable("q_ribo/loc", [n_time,n_genes]),
    scale=tf.nn.softplus(tf.get_variable("q_ribo/scale", [n_time,n_genes])))

q_protein = Normal(
    loc=tf.get_variable("q_protein/loc", [n_time,n_genes]),
    scale=tf.nn.softplus(tf.get_variable("q_protein/scale", [n_time,n_genes])))

q_protein0 = Normal(
    loc=tf.get_variable("q_protein0/loc", [n_genes]),
    scale=tf.nn.softplus(tf.get_variable("q_protein0/scale", [n_genes])))

q_deg = Normal(
    loc=tf.get_variable("q_deg/loc", [n_genes]),
    scale=tf.nn.softplus(tf.get_variable("q_deg/scale", [n_genes])))


inference = ed.KLqp(
		{deg: q_deg},
	 	{riboreads: riboreadvals, ms: msvals }
	)

inference.initialize(n_print=20, n_iter=100)





# # plt.plot()
# # plt.savefig('tmp.png')

# r = robjects.r
# df = r("read.table('tmp.tsv')")

# r(r"""
# 	pdf("tmp.pdf"%T>%{normalizePath(.)%>%message});
# 	print(read_tsv("tmp.tsv")%>%
# 		ggplot(aes(x=time,color=as.factor(gene),y=value))+
# 		facet_grid(var~.)+
# 		theme_bw()+
# 		geom_point()+
# 		ggtitle('expression_patterns'))
# 		dev.off()""")


