from __future__ import print_function
import math
import os
import torch
import torch.distributions.constraints as constraints
import pyro
import pandas as pd
from pyro.optim import Adam
from pyro.infer import SVI, Trace_ELBO
import pyro.distributions as dist


# this is for running the notebook in our testing framework
smoke_test = ('CI' in os.environ)
n_steps = 2 if smoke_test else 2000

# enable validation (e.g. validate parameters of distributions)
# assert pyro.__version__.startswith('0.3.1')
pyro.enable_validation(True)

# clear the param store in case we're in a REPL
pyro.clear_param_store()

# create some data with 6 observed heads and 4 observed tails
data = []
for _ in range(6):
    data.append(torch.tensor(1.0))
for _ in range(4):
    data.append(torch.tensor(0.0))

def model(data):
    # define the hyperparameters that control the beta prior
    alpha0 = torch.tensor(10.0)
    beta0 = torch.tensor(10.0)
    # sample f from the beta prior
    f = pyro.sample("latent_fairness", dist.Beta(alpha0, beta0))
    # loop over the observed data
    for i in range(len(data)):
        # observe datapoint i using the bernoulli likelihood
        pyro.sample("obs_{}".format(i), dist.Bernoulli(f), obs=data[i])

def guide(data):
    # register the two variational parameters with Pyro
    # - both parameters will have initial value 15.0.
    # - because we invoke constraints.positive, the optimizer
    # will take gradients on the unconstrained parameters
    # (which are related to the constrained parameters by a log)
    alpha_q = pyro.param("alpha_q", torch.tensor(15.0),
                         constraint=constraints.positive)
    beta_q = pyro.param("beta_q", torch.tensor(15.0),
                        constraint=constraints.positive)
    # sample latent_fairness from the distribution Beta(alpha_q, beta_q)
    pyro.sample("latent_fairness", dist.Beta(alpha_q, beta_q))

# setup the optimizer
adam_params = {"lr": 0.0005, "betas": (0.90, 0.999)}
optimizer = Adam(adam_params)

# setup the inference algorithm
svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

# do gradient steps
for step in range(n_steps):
    svi.step(data)
    if step % 100 == 0:
        print('.', end='')

# grab the learned variational parameters
alpha_q = pyro.param("alpha_q").item()
beta_q = pyro.param("beta_q").item()

# here we use some facts about the beta distribution
# compute the inferred mean of the coin's fairness
inferred_mean = alpha_q / (alpha_q + beta_q)
# compute inferred standard deviation
factor = beta_q / (alpha_q * (1.0 + alpha_q + beta_q))
inferred_std = inferred_mean * math.sqrt(factor)

print("\nbased on the data and our prior belief, the fairness " +
      "of the coin is %.3f +- %.3f" % (inferred_mean, inferred_std))


##Attempt to implement this for my data



testgeneset = ("Acadvl", "Ak1", "Asna1", "Cald1", "Cst3", "Ctnnd2", "Cul2",
"Dclk1", "Dpysl3", "Epb41l1", "Flna", "Hspa12a", "Igsf21", "Nos1",
"Orc3", "Phactr4", "Rap1b", "Rasa3", "Rps6ka5", "Satb2", "Srcin1",
"Stx1b", "Tbc1d24", "Trmt1l", "Zc2hc1a", "Zmym3")

# clear the param store in case we're in a REPL


exprdata = pd.read_csv('exprdata/transformed_data.txt',sep='\t')
exprdata = exprdata.loc[exprdata.gene_name.isin(testgeneset)]
ribodata = exprdata.loc[:,exprdata.columns.str.contains('ribo')]
msdata = exprdata.loc[:,exprdata.columns.str.contains('MS')].fillna(method='ffill')

ribotensor = 2**torch.tensor(ribodata.values,dtype=torch.float64).reshape(len(testgeneset),5,2)
mstensor = 2**torch.tensor(msdata.values,dtype=torch.float64).reshape(len(testgeneset),5,3)

msmean=mstensor.mean()
mssd=mstensor.std()

ms_mu_p = msmean.log()
ms_sg_p = mssd.log()*2.
mssig_mu_p = mssd.log()
mssig_sg_p = mssd.log()*10


#first simple model, 
data=mstensor

pyro.clear_param_store()


def model(data):
    # define the hyperparameters that control the beta prior
    
    i=0;t=0;r=0
    # sample f from the beta prior    # loop over the observed data
    for i in range(data.shape[0]):
        sig = pyro.sample("log_sig_{}".format(i),dist.Normal(mssig_mu_p,mssig_sg_p)).exp()
        for t in range(data.shape[1]):
            prot = pyro.sample("prot_{}_{}".format(i,t),dist.Normal(ms_mu_p,mssig_mu_p))
            # observe datapoint i using the bernoulli likelihood
            for r in range(data.shape[2]):
                pyro.sample("ms.{}.{}.{}".format(i,t,r),dist.Normal(prot,sig),obs=data[i,t,r])


model(data=mstensor[0:2])


def guide(data):    
    # sample f from the beta prior    # loop over the observed data
    for i in range(data.shape[0]):
        for t in range(data.shape[1]):
            qloc = pyro.param("qloc_{}_{}".format(i,t),dist.Normal(ms_mu_p,mssig_mu_p))
            qscale = pyro.param("qscale_{}_{}".format(i,t),dist.Normal(ms_mu_p,mssig_mu_p))
            prot = pyro.sample("prot_{}_{}".format(i,t),dist.Normal(qloc,qscale))
            # observe datapoint i using the bernoulli likelihood

guide(mstensor[0:2])
# guide(mstensor)

# setup the optimizer
adam_params = {"lr": 0.0005, "betas": (0.90, 0.999)}
optimizer = Adam(adam_params)

# setup the inference algorithm
svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

svi.step(data)

# do gradient steps
for step in range(n_steps):
    svi.step(data)
    if step % 100 == 0:
        print('.', end='')
# grab the learned variational parameters
alpha_q = pyro.param("alpha_q").item()
beta_q = pyro.param("beta_q").item()


################################################################################
########Now without the loops
################################################################################
    


def model(ribotensor,mstensor):
    # define the hyperparameters that control the beta prior
    ribotensor 
    # sample f from the beta prior    # loop over the observed data


def guide(data):    
    # sample f from the beta prior    # loop over the observed data
    for i in range(data.shape[0]):
        for t in range(data.shape[1]):
            qloc = pyro.param("qloc_{}_{}".format(i,t),dist.Normal(ms_mu_p,mssig_mu_p))
            qscale = pyro.param("qscale_{}_{}".format(i,t),dist.Normal(ms_mu_p,mssig_mu_p))
            prot = pyro.sample("prot_{}_{}".format(i,t),dist.Normal(qloc,qscale))
            # observe datapoint i using the bernoulli likelihood









