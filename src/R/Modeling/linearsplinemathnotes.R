

###############################################################################
#######Math  - this time linear space splines
###############################################################################
So I think this is all fine if we just use linear space splines
IF we can express a time constant function as a simple combination of splines.
Now, our linear thing is a 


time = 0:23 + 10
sig = 0:23 + 10

dbs%>%args

fit=lm(sig ~ 0+ bs(time,intercept=TRUE))
fit
fit_x2=lm(sig*2 ~ bs(time,intercept=TRUE))
fit_x2


pdf('tmp.pdf')
bs(time,intercept=TRUE)%>%matplot
dev.off()
normalizePath('tmp.pdf')

dbs(time,intercept=TRUE) %*% (fit$coef)
bs(time,intercept=TRUE) %*% (fit$coef)
ibs(time,intercept=TRUE) %*% (fit$coef)

#so if we want our RNA to be 1+4 dimensions for an arbitrary fit
#then we can have 
#and we want our prot to be 1+4 dimensions for arbitrary fit
inter_dbasis = cbind(1,dbs(time,df=4))
integ_inter_dbasis = cbind(1,1:(1*length(time)),bs(time,df=4))





#So weith default df = 3
#so note that with an intercept, bs has only 3
princomp(bs(time))%>%summary
princomp(bs(time,intercept=TRUE))%>%summary
#And dbs only has two
princomp(dbs(time,intercept=TRUE))%>%summary
#so to fit e.g. 5 timepoint in our rate of change, we'd need 6 
mydf = 6
#we would like there to be a meaningful intercept at both levels actually
princomp(dbs(time,df=mydf,intercept=TRUE))%>%summary

let bs be the linear spline basis
px = ([m0,a1,a2,a3])
so that bs %*% pv, bs is a txD matrix
dbs is now also txD
let zv be teh coefficients of a linear fit in our spline basis

dP(t) = Ks P(t) - Kd (P(t))



zv = lm(time ~ 0 + inter_dbasis)$coef
zv = lm(1:24 ~ 0 + integ_inter_dbasis[,,-5])$coef


note that px1 is t+1 for a perfect  


dbasis = inter_dbasis
mybs = integ_inter_dbasis



###actually hadamard producct

#Preeeety sure this works

t,D    D,1    =   t,D     D,1 -  t,D     D,1  %.%  t,D1   D1,1
dbasis %*% px = dbasis %*% yv  -  dbasis %*% zv %.%  mybs %*% px1

dbasis %*% px = dbasis %*% yv -  dbasis %*% zv %.%  mybs %*% px1
#add
dbasis %*% px + dbasis %*% zv %.%  mybs %*% px = dbasis %*% yv
#flilp
dbasis %*% yv = dbasis %*% px + dbasis %*% zv %.%  mybs %*% px1 

#Okay now can I factorize this???
#note the hadamard (elementwise) product there
#t,D      D,1 = t,D    %*% D,1 +     t,D   . D,1        t,D1 %*% D1,1
dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

#Okay so no, this isn't a linear operation, but it's lcose
#t,D      D,1 =  t,D  %*% D,1 +  t,D   .    D,1  %*%    
dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 
#Now let's assume that there's a left inverse to dbasis - i.e. it has linearly independent columns
dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 

#now we end up with this.......
yv = px + dbinv %*% diag(mybs %*% px1) %*% dbasis %*% t(zv) 


#We COULD also just express our zv in terms of the protein level basis
dbasis %*% yv = dbasis %*% px +  (mybs %*% zv) %.% (mybs %*% px1)

dbasis %*% yv = dbasis %*% px +  mybs %*%  (zv %.% px1)

#again assume linear independence in the columns of dbasis
dbasis %*% yv = dbasis %*% px +  mybs %*% zv  %.% mybs %*% zv



#what about ks...



dbasis %*% yv = dbasis %*% yv +  mybs %*%  (zv %.% px1)
D,1 =D,1 +       D,t   %*% t,D1  %*%   D1,1

#this is WAY nicer

                    t,D1    D1,1
yv = px + invdbasis %*% mybs %*%  (zv %.% px1)

#waht if we assume that zv[1] = 0?

yv = px + invdbasis %*% mybs %*%  (zv %.% px)


hadarmard product can be expressed as linear product with a diagnonal matrix

yv = px + invdbasis %*% mybs %*%  diag_z %*%  px

yv = px + invdbasis %*% mybs %*%  diag_z %*%  px

yv = px + invdbasis %*% mybs %*%  diag_z %*%  px

yv =  (1+invdbasis %*% mybs %*%  diag_z) %*% px

yv =  (1+invdbasis %*% mybs %*%  diag_z) %*% px








now in orthogonolaizing our basis... we need to think about how the coefficients are effected.

hang on, if one of our basis terms is the intercept in dbs space, then every term but one in


So my problem is just optimizing such that the RNA is always a positive value.
That means I need M splines
I need to devise a basis, then, show i can get parameter values for this basis that, when I apply the equations I have,
lead to reasonablee values for thee other measure.

okay well the log spline is working okay now, though the prot trajectory got weird when I increased the lower limit boundary


I should maybe think about a linear spline
which means parametrizing the whole thing in terms of the log fold change
And then the degredation.
