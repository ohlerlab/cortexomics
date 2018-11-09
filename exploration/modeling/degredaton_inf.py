import tensorflow as tf
import edward as ed
import numpy as np
import pandas as pd
from edward.models import OneHotCategorical,Categorical, Dirichlet, InverseGamma, \
    MultivariateNormalDiag, Normal, ParamMixture, \
    Gamma, Mixture, DirichletProcess
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2



###Basic testing of how edward works - estimating a normal variable

with tf.variable_scope(tf.get_variable_scope(),reuse=tf.AUTO_REUSE):
	N = 100
	mu = np.random.normal()
	x = np.random.normal(loc=mu, scale=1, size=(N, 1)).astype(np.float32)
	px = ed.models.Normal(loc=tf.Variable(tf.zeros([1])) * tf.ones([N, 1]), scale=tf.ones([1]))
	inf = ed.KLqp(data={px: x})
	inf.run(n_iter=250)
	print(mu, x.mean(), ed.get_session().run(px.mean()[0, 0]))
	print(ed.evaluate('mean_squared_error', data={px: x, : y_test}))

tf.reset_default_graph()
tf.get_variable_scope().reuse_variables()
sess=tf.InteractiveSession()


with tf.variable_scope(tf.get_variable_scope(),reuse=tf.AUTO_REUSE):
	N = 100
	np.random.seed(0)

	mu = tf.placeholder(tf.float32,[1])
	x = ed.models.Normal(loc=mu, scale=tf.ones([1]))

	mu = np.random.normal()
	x = np.random.normal(loc=mu, scale=1, size=(N, 1)).astype(np.float32)
	px = ed.models.Normal(loc=tf.Variable(tf.zeros([1])) * tf.ones([N, 1]), scale=tf.ones([1]))
	inf = ed.KLqp(data={px: x})
	inf.run(n_iter=250)
	print(mu, x.mean(), ed.get_session().run(px.mean()[0, 0]))
	ed.evaluate()



with tf.variable_scope(tf.get_variable_scope(),reuse=tf.AUTO_REUSE):
	
	N = 100#number of reps
	np.random.seed(0)#seed

	mumu = tf.Variable(tf.zeros([1]))
	musig = tf.Variable(tf.ones([1.]))
	mu = ed.models.Normal(mumu,tf.nn.softplus(musig))
	mu = ed.models.Normal(mumu* tf.ones([N, 1]),musig)
	# mu = tf.Variable(tf.ones([1.]))

	#now draw samples
	mureal = np.random.normal()#get some mean 
	x_train = np.random.normal(loc=mureal, scale=1, size=(N, 1)).astype(np.float32)

	#define an approximate distribution over mu
	qmu = ed.models.Normal(
		loc=tf.get_variable('qmu/loc',[1])* tf.ones([N, 1]), 
		scale=tf.get_variable('qmu/scale',[1]))

	dx = ed.models.Normal(loc=mu , scale=tf.ones([1]))

	#run inference
	inf = ed.KLqp({mu:qmu},data={dx: x_train})
	inf.run()

	#fetch our parameters
	vars2get = (qmu.mean()[0],qmu.stddev()[0])
	print(mureal, x.mean(),ed.get_session().run(vars2get) )


#define a simple riboseq/protein system with degredation, initial protein levels and riboseq
#as pointestimates
n_genes = 3 
n_time = 5
n_riboreps = 2
n_msreps = 4

# vals = define_model(n_genes,n_time)
deg = tf.placeholder(tf.float32,shape=n_genes)
protein0 = tf.placeholder(tf.float32,shape=n_genes)
ribo = tf.placeholder(tf.float32,shape=(n_time,n_genes))

#define our protein
protein= [ (ribo[0,:] + (protein0))]
for  i in range(1,n_time):
	protein= protein + [ribo[len(protein)-1,:] + (protein[-1] * deg)]
protein = tf.stack(protein)

#define our 
riboreads = Normal(
	tf.stack([ribo]*n_riboreps,axis=2),
	tf.sqrt(tf.stack([ribo]*n_riboreps,axis=2))
)

#ms is similiar to riboreads for now - but with higher variance
ms = Normal(
	tf.stack([protein]*n_msreps,axis=2),
	tf.sqrt(tf.stack([protein]*n_msreps,axis=2))
)

vals = ( protein, riboreads, ms )
	# return( protein, riboreads, ms )

vars = tf.trainable_variables()


#take samples from this system - feeding in rriboseq, protein levels and degredation levels
#with riboseq set to almost zero
with tf.Session() as sess:

	ms_sample=sess.run(
		fetches = vals,
		feed_dict = {
			# protein0:[10000]*n_genes,
			protein0:[1000,3000,5000],
			deg:[0.5]*n_genes,
			ribo: np.stack([[0]*n_time]*n_genes,axis=-1)
		}
	)

ms_sample[0][1,0]
ms_sample[1][2,:,:]
ms_sample[1][0,2,:]
ms_sample[2].shape

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

#format this data
varnames = ['time','gene','replicate']
dfs = [stackarray(d,nm,varnames) for d,nm in zip(ms_sample,varnames)]
#name vairable
for i in range(0,len(dfs)):
	dfs[i]['var'] = varnames[i]
#stack
alldata = pd.concat(dfs,sort=False)
#and save
alldata.to_csv('tmp.tsv',sep='\t')
alldata.head()




#from the linear mixed effects bit in the edward tutorials
with tf.variable_scope('b'):
	deg = tf.get_variable('deg',[n_genes])
	protein0 = tf.sigmoid(tf.get_variable('protein0',[n_genes]))
	ribo = tf.get_variable('ribo',[n_time,n_genes])

#define our protein
protein= [ (ribo[0,:] + (protein0 * deg))]
for  i in range(1,n_time):
	protein= protein + [ribo[len(protein)-1,:] + (protein[-1] * deg)]
protein = tf.stack(protein)

#define our 
riboreads = Normal(
	tf.stack([protein]*n_riboreps,axis=2),
	tf.sqrt(tf.stack([protein]*n_riboreps,axis=2))
)

#ms is similiar to riboreads for now - but with higher variance
ms = Normal(
	tf.stack([protein]*n_msreps,axis=2),
	tf.sqrt(tf.stack([protein]*n_msreps,axis=2))
)


vals = ( protein, riboreads, ms )
	# return( protein, riboreads, ms )

vars = tf.trainable_variables()
	
	#so 
with tf.Session() as sess:

	ms_sample=sess.run(
		fetches = vals,
		feed_dict = {
			protein0:[10000]*n_genes,
			deg:[0.5]*n_genes,
			ribo: np.stack([[1]*n_time]*n_genes,axis=-1)
		}
	)




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




#Now use R to plot
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


