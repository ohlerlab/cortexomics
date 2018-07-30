

#cluster selection vector
#cluster number k
# the number of genes
import tensorflow as tf
import edward as edw
import edward as ed
import numpy as np
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

M0 = 12000
D = 2



M(0) = M0 * exp(-D_int(t))
	 = M0 * exp(-D_int(0))
	 = M0 * 1
	 = M0

M() = M0 * exp(-D_int(t)) + P()


dMS/dt = Prod - Degredation
	   = C*Ribo - D*MS




logfc = log(dMS/dt + MS / MS) 
	  = log( ((Prod - Degredation) + MS) / (MS))
	  = log( (C*Ribo - D (MS) ) / (MS))
	  = log((C*Ribo - D(MS)))  - log(MS)



M0 = 1000
D = 2

def Dtild(t):
	return D*t

def Intupto(t,myfun):
	vlist = [myfun(ti) for ti in range(0,t)]
	return sum(vlist)

def P(x):
	return 500 

def M(t,M0 = 100):
	return np.exp(-Dtild(t)) * (
		M0 + (
			Intupto(t,
				lambda x:
					P(x) * np.exp(Dtild(x))
				)
			)		
	)


tmin,tmax = 0,20
mvals = [M(t) for t in range(tmin,tmax)]


pllimits = [
	np.floor(np.min(tmin)),
	np.ceil(np.max(tmax)),
	0,
	np.ceil(np.max(mvals))
]
plt.axis(pllimits)
plt.scatter(list(range(tmin,tmax)),mvals)
plt.title("MS values")
plt.savefig('tmp.png')
plt.close()



logfc = log (ms2 / ms1)
#Degredation is complete
logfc = log ((P2) / (P1))
#PRoduction is zero
logfc = log (M0*exp(Dt) / (M0*exp(D(n-1))))
#
logfc = log (M0*exp(Dt) / (M0*exp(D(n-1))))
	  = log (M0*exp(Dt)) - log((M0*exp(D(n-1))))
	  = log (M0*exp(Dt)) - log((M0*exp(D(n-1))))

= log (M0*exp(Dt)) - log((M0*exp(D(n-1))))



#so some amount M0, some time previously
#MS5 = Ribo5*c + MS4*deg
#MS4 = Ribo4*c + MS3+ 
#Final MS consists of newly produced which is 


#okay so this is a working probabalistic system, albeit very stripped down...

def todimn(cont,n,step = 100):
	return [cont * i*step for i in range(1,n+1)]

ribomsbase = tf.constant(todimn(100.,2))
deg = tf.constant([0.5,0.5])

ms0 = tf.placeholder(tf.float32,shape=2)

ribo1=tf.constant(todimn(1.,2))
ribo2=tf.constant(todimn(1.1,2))
ribo3=tf.constant(todimn(3.1,2))
ribo4=tf.constant(todimn(3.2,2))
ribo5=tf.constant(todimn(4.,2))

ms1	= Normal(ribo1*ribomsbase + ms0 * deg,.1)
ms2 = Normal(ribo2*ribomsbase + ms1 * deg,.1)
ms3 = Normal(ribo3*ribomsbase + ms2 * deg,.1)
ms4 = Normal(ribo4*ribomsbase + ms3 * deg,.1)
ms5 = Normal(ribo5*ribomsbase + ms4 * deg,.1)


sess = tf.Session()

ms_sample = sess.run(
	fetches = [ms1.value(), ms2.value(), ms3.value(), ms4.value(),ms5.value()],
	feed_dict = {ms0:[100,100]}
)






