import pandas as pd
from observations import insteval
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
tf.reset_default_graph()

data, metadata = insteval("~/data")
data = pd.DataFrame(data, columns=metadata['columns'])

# s - students - 1:2972
# d - instructors - codes that need to be remapped
# dept also needs to be remapped
data['s'] = data['s'] - 1
data['dcodes'] = data['d'].astype('category').cat.codes
data['deptcodes'] = data['dept'].astype('category').cat.codes
data['y'] = data['y'].values.astype(float)

train = data.sample(frac=0.8)
test = data.drop(train.index)


s_train = train['s'].values
d_train = train['dcodes'].values
dept_train = train['deptcodes'].values
y_train = train['y'].values
service_train = train['service'].values
n_obs_train = train.shape[0]

s_test = test['s'].values
d_test = test['dcodes'].values
dept_test = test['deptcodes'].values
y_test = test['y'].values
service_test = test['service'].values
n_obs_test = test.shape[0]

n_s = max(s_train) + 1  # number of students
n_d = max(d_train) + 1  # number of instructors
n_dept = max(dept_train) + 1  # number of departments
n_obs = train.shape[0]  # number of observations

print("Number of students: {}".format(n_s))
print("Number of instructors: {}".format(n_d))
print("Number of departments: {}".format(n_dept))
print("Number of observations: {}".format(n_obs))


# Set up placeholders for the data inputs.
s_ph = tf.placeholder(tf.int32, [None])
d_ph = tf.placeholder(tf.int32, [None])
dept_ph = tf.placeholder(tf.int32, [None])
service_ph = tf.placeholder(tf.float32, [None])


# Set up fixed effects.
mu = tf.get_variable("mu", [])
service = tf.get_variable("service", [])
musig = tf.get_variable("musig", [])
mumu = tf.get_variable("mumu", [])
mu = Normal(mumu,musig)

sigma_s = tf.sqrt(tf.exp(tf.get_variable("sigma_s", [])))
sigma_d = tf.sqrt(tf.exp(tf.get_variable("sigma_d", [])))
sigma_dept = tf.sqrt(tf.exp(tf.get_variable("sigma_dept", [])))

# Set up random effects.
eta_s = Normal(loc=tf.zeros(n_s), scale=sigma_s * tf.ones(n_s))
eta_d = Normal(loc=tf.zeros(n_d), scale=sigma_d * tf.ones(n_d))
eta_dept = Normal(loc=tf.zeros(n_dept), scale=sigma_dept * tf.ones(n_dept))

yhat = (tf.gather(eta_s, s_ph) +
        tf.gather(eta_d, d_ph) +
        tf.gather(eta_dept, dept_ph) +
        mu + service * service_ph)

y = Normal(loc=yhat, scale=tf.ones(n_obs))

q_eta_s = Normal(
    loc=tf.get_variable("q_eta_s/loc", [n_s]),
    scale=tf.nn.softplus(tf.get_variable("q_eta_s/scale", [n_s])))
q_eta_d = Normal(
    loc=tf.get_variable("q_eta_d/loc", [n_d]),
    scale=tf.nn.softplus(tf.get_variable("q_eta_d/scale", [n_d])))
q_eta_dept = Normal(
    loc=tf.get_variable("q_eta_dept/loc", [n_dept]),
    scale=tf.nn.softplus(tf.get_variable("q_eta_dept/scale", [n_dept])))
q_mu = Normal(
    loc=tf.get_variable("q_mu/loc", [1]),
    scale=tf.nn.softplus(tf.get_variable("q_mu/scale", [1])))

latent_vars = {
    eta_s: q_eta_s,
    eta_d: q_eta_d,
    eta_dept: q_eta_dept,   q_mu}

data = {
    y: y_train,
    s_ph: s_train,
    d_ph: d_train,
    dept_ph: dept_train,
    service_ph: service_train}

inference = ed.KLqp(latent_vars, data)

yhat_test = ed.copy(yhat, {
    eta_s: q_eta_s.mean(),
    eta_d: q_eta_d.mean(),
    eta_dept: q_eta_dept.mean()})

inference.initialize(n_print=20, n_iter=100)
tf.global_variables_initializer().run()

for _ in range(inference.n_iter):
  # Update and print progress of algorithm.
  info_dict = inference.update()
  inference.print_progress(info_dict)

  t = info_dict['t']
  if t == 1 or t % inference.n_print == 0:
    # Make predictions on test data.
    yhat_vals = yhat_test.eval(feed_dict={
        s_ph: s_test,
        d_ph: d_test,
        dept_ph: dept_test,
        service_ph: service_test})

mumu.eval()
musig.eval()
mu.mean().eval()
# with tf.Session() as sess:
# 	##Toy examples
# 	# Set up probability model.
# 	mu = Normal(loc=0.0, scale=1.0)
# 	x = Normal(loc=mu, scale=1.0, sample_shape=50)

# 	# Set up posterior approximation.
# 	qmu_loc = tf.Variable(tf.random_normal([]))
# 	qmu_scale = tf.nn.softplus(tf.Variable(tf.random_normal([])))
# 	qmu = Normal(loc=qmu_loc, scale=qmu_scale)

# 	inference = ed.Inference({mu: qmu}, data={x: tf.zeros(50)})

