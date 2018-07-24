#I'm still not sure what the fuck a dirichlet process is.
#chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://cran.r-project.org/web/packages/dirichletprocess/vignettes/dirichletprocess.pdf
#okay I get that a dirchlet process is a probability distribution over proabilty distributions,
#what I don't get is how the inputs to the DP() - the G0 distribution, and the concentration
#parameter, interact.
#Okay I think I"M getting it - it.. the dirchlet process just gives gives you catagories
#it's a dirchlet process mixture model that one generally wants, in which the 
#hmm but the wiki page suggests that it's more like... picking 
#https://en.wikipedia.org/wiki/Dirichlet_process


#Need to thinka bout all this more, let's just go with a constant number of clusters for now

#let's play with edward!

#a constant number of clusters + multivariate gaussian on the fold changes should be doable
#another option would be a gaussian process - this would constrain the covariance matrix in a
#more reasonable way...
#could I use cell content or something rather than time for the GP?
#In the mean time I"m pretty sure I could just specify something identical to the 



#https://www.statlect.com/probability-distributions/normal-distribution-linear-combinations
#formulas for covariances between sums of linear combiantions
#so the variances add, 
#what about the covariances?
#How do I make all this happen in edward?
#I can imagine just doing it timepoint by timepoint without explicitly specifying the 
#covariance matrix for a whole cluster...
#can I in fact represent conditional independance using a simple multivariate gaussian and
#covariance matrix? Actually maybe not.... but that could be fine


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





#lets do a very simple thing, and first just look at doing ribodiff with edward.





k=2
concentration_vector = tf.ones(k)/k
#so these give us k dimensional vectors
clustvects = OneHotCategorical(concentration_vector)
#note that clustdist.sample(3).eval(session=sess).shape gives us (4,8) - rows are samples.
clustvects = Categorical(shape=(n,k))#not sure if this call is right...

#draw all of these (this these are the ones we want to infer)


#prior distribution of means 
fc_prior_mean_variance=1000.
N = 1000
k = 2
D= 2
clustmeans = [[i*100]*D for i in list(range(0,k))]

#make data from a normal mixture
x_sim=np.random.normal(0,1,size=100) + np.array(([100]*N/k)+([20]*N/k))

meanvect = Normal(tf.zeros(k),tf.ones(k)*fc_prior_mean_variance)

clustprobs = Dirichlet(tf.ones(k)/k)
clustvects = OneHotCategorical(clustprobs)
clustvects = tf.cast(clustvects,tf.float32)
sess = tf.Session()

N = 1000
K = 3
D= 2
clustmeans = np.array([[i*100]*D for i in list(range(0,K))])

#betas have this prior distriubtion
beta = Normal(loc=tf.zeros([K, D])+clustmeans, scale=tf.ones([K, D]) )
#the clusters have this latent variable
z = Categorical(logits=tf.zeros([N, K]))
#the observed variable has this distribution
x = MultivariateNormalDiag(loc=tf.gather(beta, z), scale_diag=tf.ones(D))


w2 = MultivariateNormalDiag(loc=tf.zeros([D,D]),scale_diag = tf.ones(D)*10)


allw2 = MultivariateNormalDiag(loc=tf.zeros([D,D]),scale_diag = tf.ones(D)*10,batch_shape=3)
w2 = tf.gather()

b2 = Normal(loc=tf.zeros([N,D]),scale = tf.ones([N,D])*10)
x2 = MultivariateNormalDiag(loc=tf.matmul(x,w2) + b2, scale_diag=tf.ones(D))

w3 = Normal(loc=tf.zeros([D,D]),scale = tf.ones([D,D])*10)

w3 = MultivariateNormalDiag(loc=tf.zeros([D,D]),scale_diag = tf.ones(D)*10)
b3 = Normal(loc=tf.zeros([N,D]),scale = tf.ones([N,D])*10)
x3 = MultivariateNormalDiag(loc=tf.matmul(x2,w3) + b3, scale_diag=tf.ones(D))

w4 = MultivariateNormalDiag(loc=tf.zeros([D,D]),scale_diag = tf.ones(D)*10)
b4 = Normal(loc=tf.zeros([N,D]),scale = tf.ones([N,D])*10)
x4 = MultivariateNormalDiag(loc=tf.matmul(x3,w4) + b4, scale_diag=tf.ones(D))

X = tf.stack([x,x2])
x_test = ed.copy([x,x2,x3,x4])
x_test = ed.copy(X)

x_test.value()

x_samp = x.sample(1).eval(session=sess)

X = tf.ones([100,2])+np.array(list(range(0,200))).reshape([100,2])


X = tf.placeholder(tf.float32, [N, D])
w = Normal(loc=tf.zeros(D), scale=tf.ones(D))
b = Normal(loc=tf.zeros(1), scale=tf.ones(1))
y = Normal(loc=ed.dot(X, w) + b, scale=tf.ones(N))


X_sample = sess.run([x.value(), x2.value(), x3.value(), x4.value(), w2,w3,w4])



X_sample[-3]


# fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
# ax.plot([0,1,2], [10,20,3])
# fig.savefig('testimage.png')   # save the figure to file
# plt.close(fig)   


#this code plots the data as sampled from my 
pllimits = [
	np.ceil(np.max(x_samp)/10)*10,
	np.floor(np.min(x_samp)/10)*10,
	np.ceil(np.max(x_samp)/10)*10,
	np.floor(np.min(x_samp)/10)*10
]
plt.scatter(x_samp[0,:, 0], x_samp[0,:, 1])
plt.axis(pllimits)
plt.title("Predicted cluster assignments")
plt.savefig('tmp.png')
plt.close()






# x = Normal(loc=tf.gather(beta, z), scale=tf.ones([N, D]))


clustmeans = edw.dot(clustvects,meanvect)



x = Normal(clustmeans,Normal(varprior))

#prior distribution of variances
varpriordist Inversegamma(variance_prior_paramns)
#weight prior dist
fccorweights ~ Normal(0,fc_prior_weight_variance)#something vauge

#matrix of means for fold changes
#note this is a k,ntrack,ntps tensor
fcmeans ~ meanpriordist(shape=(k,ntracks,ntps))
fcvars ~ varpriordist(shape=(k,ntracks,ntps))
fccorweights ~ varpriordist(shape=(k,ntracks,ntps-1))

#first fold changes for first tracks ay and first tp

#this should get the means for each sample, using it's cluster membershiop
genemeans ~ tf.matmult(clustervects,fcmeans[:,:,:])
genevars ~tf.matmult(clustervects,fcvars[:,:,:])#and the vars
#we also need to learn the weights

#distribution of the first timespiont
tp1 ~ Normal(genemeans[:,:,1],genevars[:,:,1],observed = X[:,:,1])
tp1 ~ Normal(genemeans[:,:,2],genevars[:,:,2],observed = X[:,:,1])

#now for the second
X[:,:,2] ~ Normal(
	tf.matmult(genemeans[:,:,1],t(genemeans[:,:,2])),
	tf.matmult(genevars[:,:,1],t(genevars[:,:,2]))
)

#now for the second
#note that the structure of the weight matrix here will determine whether we allow
#covariances
#the weights between each timestage would be a K*K matrix, diagonal if we're not
#allowing inter-track depenencies.
#note that I could also be doing timepoint indexign with matrix multiplication here...
X[:,:,2] ~ Normal(
	tf.matmult(genemeans[:,:,1],t(genemeans[:,:,2])),
	tf.matmult(genevars[:,:,1],t(genevars[:,:,2]))
)
#....and so on...

#Now I THINK this should work, or something like it. 
#not too sure about how i"ve connected the distirbutions to the data. I need to work on that...

#then I need to figure out fi all this can be optimized.

#The above could also be done with use of tf.where

#there are several componenets to the final thing I want, and several approximations I want to hit on the 
#way. I would like to do clusteirng - i.e. fit a mixture distribution. 



#Things that I think
#I THINK that the timeless model is kiiiiind of similiar to a gaussian process. 
#There will stil be covariance between distant timestages, the conditionality structure
#simply makes the constraint that it can only decrease over time, and allows the timeless model
#to learn what the distances (covariances) between timepoints ought to be, instead of assuming 
#that the clock time has any kind of meaning, or common meanign between clusters.


