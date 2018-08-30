

#cluster selection vector
#cluster number k
# the number of genes
import tensorflow as tf
import edward as edw
import edward as ed
import numpy as np
import pandas as pd
from edward.models import OneHotCategorical,Categorical, Dirichlet, InverseGamma, \
    MultivariateNormalDiag, Normal, ParamMixture, \
    Gamma, Mixture, DirichletProcess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt





#lets do a very simple thing, and first just loo




############################
############################
####################################

# ####DOn't RUN me oing math 
# M0 = 12000
# D = 2



# M(0) = M0 * exp(-D_int(t))
# 	 = M0 * exp(-D_int(0))
# 	 = M0 * 1
# 	 = M0

# M() = M0 * exp(-D_int(t)) + P()


# dMS/dt = Prod - Degredation
# 	   = C*Ribo - D*MS




# logfc = log((dMS/dt + MS) / MS) 
# 	  = log(((Prod - Degredation) + MS) / MS)
# 	  = log( C*Ribo - D (MS) + MS ) -log(MS)   



# # M0 = 1000
# # D = 2

# def Dtild(t):
# 	return D*t

# def Intupto(t,myfun):
# 	vlist = [myfun(ti) for ti in range(0,t)]
# 	return sum(vlist)

# def P(x):
# 	return 500 

# def M(t,M0 = 100):
# 	return np.exp(-Dtild(t)) * (
# 		M0 + (
# 			Intupto(t,
# 				lambda x:
# 					P(x) * np.exp(Dtild(x))
# 				)
# 			)		
# 	)


# tmin,tmax = 0,20
# mvals = [M(t) for t in range(tmin,tmax)]


# pllimits = [
# 	np.floor(np.min(tmin)),
# 	np.ceil(np.max(tmax)),
# 	0,
# 	np.ceil(np.max(mvals))
# ]
# plt.axis(pllimits)
# plt.scatter(list(range(tmin,tmax)),mvals)
# plt.title("MS values")
# plt.savefig('tmp.png')
# plt.close()




# logfc = log (ms2 / ms1)
# #Degredation is complete
# logfc = log ((P2) / (P1))
# #PRoduction is zero
# logfc = log (M0*exp(Dt) / (M0*exp(D(n-1))))
# #
# logfc = log (M0*exp(Dt) / (M0*exp(D(n-1))))
# 	  = log (M0*exp(Dt)) - log((M0*exp(D(n-1))))
# 	  = log (M0*exp(Dt)) - log((M0*exp(D(n-1))))

# = log (M0*exp(Dt)) - log((M0*exp(D(n-1))))



#so some amount M0, some time previously
#MS5 = Ribo5*c + MS4*deg
#MS4 = Ribo4*c + MS3+ 
#Final MS consists of newly produced which is 


#okay so this is a working probabalistic system, albeit very stripped down...

print('here')

def todimn(cont,n,step = 100):
	return [cont * i * step for i in range(1,n+1)]

print('here')

#now all of this can probably be compressed...
#I'll want to run a performance evaluation of sequentially
#defined stuff before I use it tho

# ribo1=tf.constant(todimn(1.,n_genes))
# ribo2=tf.constant(todimn(1.1,n_genes))
# ribo3=tf.constant(todimn(3.1,n_genes))
# ribo4=tf.constant(todimn(3.2,n_genes))
# ribo5=tf.constant(todimn(4.,n_genes))

# ribo1=tf.constant(todimn(1.,n_genes))
# ribo2=tf.constant(todimn(1.1,n_genes))
# ribo3=tf.constant(todimn(3.1,n_genes))
# ribo4=tf.constant(todimn(3.2,n_genes))
# ribo5=tf.constant(todimn(4.,n_genes))


# protein1 = (ribo1 + (protein0 * deg))
# protein2 = (ribo2 + (protein1 * deg))
# protein3 = (ribo3 + (protein2 * deg))
# protein4 = (ribo4 + (protein3 * deg))
# protein5 = (ribo5 + (protein4 * deg))

print('here')

# riboreads1=Normal(ribo1,.1)
# riboreads2=Normal(ribo2,.1)
# riboreads3=Normal(ribo3,.1)
# riboreads4=Normal(ribo4,.1)
# riboreads5=Normal(ribo5,.1)
# 

n_genes = 3
n_time = 5

# def define_model(n_genes,n_time):

# ribomsbase = tf.constant(todimn(100.,n_genes))




n_genes = 3 
n_time = 5

# vals = define_model(n_genes,n_time)


with tf.Session() as sess:

	deg = tf.placeholder(tf.float32,shape=n_genes)

	protein0 = tf.placeholder(tf.float32,shape=n_genes)

	ribo = tf.placeholder(tf.float32,shape=(5,n_genes))

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
	
	test=sess.run(
		fetches = vals,
		feed_dict = {
			protein0:[10000]*n_genes,
			deg:[0.5]*n_genes,
			ribo: np.stack([[1000]*n_time]*n_genes,axis=-1)
		}
	)






# #run our model!
# ms_sample = sess.run(
# 	fetches = [
# 		[riboreads1.value(), riboreads2.value(), riboreads3.value(), riboreads4.value(),riboreads5.value()],
# 		[protein1, protein2, protein3, protein4,protein5],
# 		[ms1.value(), ms2.value(), ms3.value(), ms4.value(),ms5.value()]
# 	],
# 	feed_dict = {
# 		protein0:[100]*n_genes,
# 		deg:[0.5]*n_genes,
# 		ribo: np.stack([[100]*n_time]*n_genes,axis=-1)
# 	}
# )

# def stackarray(d2stack,varname,indexnamelist):
# 	if(len(d2stack.shape) != 3) : d2stack = np.expand_dims(d2stack,axis = -1)
# 	indexnamelist = indexnamelist[0:len(d2stack.shape)]
# 	# indexnamelist = indexnamelist[0:len(d2stack.shape)]

# 	index = pd.MultiIndex.from_product([range(s)for s in d2stack.shape], names=indexnamelist)
# 	# df = pd.DataFrame({'d2stack': d2stack.flatten()}, index=index)['d2stack']
# 	df = pd.DataFrame({'value': d2stack.flatten()}, index=index)
# 	df['var'] = varname
# 	return df

# varnames = ['riboreads','protein','ms']

# dfs = [stackarray(np.stack(d),nm,['time','gene','rep']) for d,nm in zip(ms_sample,varnames)]

# for i in range(0,len(dfs)):
# 	dfs[i]['var'] = varnames[i]

# alldata = pd.concat(dfs,sort=False)

# alldata.to_csv('tmp.tsv',sep='\t')

# plt.plot()
# plt.savefig('tmp.png')


# import rpy2.robjects as robjects
# import rpy2.robjects.lib.ggplot2 as ggplot2

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

