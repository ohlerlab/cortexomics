

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
dbasis %*% yv = dbasis %*% px +  mybs %*% px  %.% mybs %*% zv

0_dbasis %*% ms0_yv = dbasis %*% px +  mybs %*% px  %.% mybs %*% zv


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


#Okay so trying for piecewise linear
problem, when I construct px, how do I do it? The protein implies a given amount of synthesis and degredation.
But the synthesis depends on the average amount between two points, not the amount. So 4 degrees of freedom, 5 ribo points
ribo isnt specified.

edbasis = cbind(0,mydbs)

edbasis = cbind(0,mydbs)


#This doesn't work cos edbasis not invertible...
edbasis %*% ms0_yv = edbasis %*% px +  mybs %*% px  %.% mybs %*% zv
edbasis %*% ms0_yv = edbasis %*% px +  mybs %*% px  %*%  diag_z
edbasis %*% ms0_yv = edbasis %*% px +  diag_z %*% mybs %*% px
edbasis %*% ms0_yv = (edbasis +  diag_z %*% mybs )  %*% px  


#This somehow doesn't wok....
#But the one on the right might be no?
5,6			6			5,6			5,5.    5,6.        6

edbasis %*% ms0_yv = (edbasis +  diag_z %*% mybs )  %*% px
ginv((edbasis +  diag_z %*% mybs )) %*% edbasis %*% ms0_yv = px


myTE=20
ms0=4000
ms0_yv <- c(4000,ribo*20)
ribo
ntps=length(time)
diag_z = diag(rep(mydeg,ntps))
#woooooot this is invertible
(mydbs) %>% {. %*% ginv(.)}
(edbasis +  diag_z %*% mybs) %>% {. %*% ginv(.)}
txtplot(ms0 + ginv((edbasis +  diag_z %*% mybs )) %*% edbasis %*% ms0_yv)

###sooo othat doesn't work
#Our differnce in dimensionality is not just a nuisanc eheree
#We can work out our yv given px but this is a problem
#Because we want a test set of px
#This is why they haad boundaries mebbe
#

5,5.      5.     5,6.      5,5.      5,6.      6
mydbs %*% yv = (edbasis +  diag_z %*% mybs )  %*% px
yv = ginv(edbasis) %*% ((edbasis +  diag_z %*% mybs )) %*% px

#Maybe I need to fix the bounadary asa well?
#Like the rate of change at th first boundaary is just 
#ribo0*myTE - ms0*mydeg
#synth0 - ms0*mydeg

#aha yes. 
splfit = lm(ribo ~ 0+mydbs)
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=FALSE)
prot <- num_ode_res$P
startchange <- num_ode_res$ribo[1]*myrTE - prot[1]*mydeg
endchange <- num_ode_res$ribo[ntps]*myrTE - prot[5]*mydeg
px_lm = lm(c(startchange,prot) ~ 0+ rbind(edbasis[1,],mybs))
px = px_lm$coef
txtplot(predict(px_lm)[-1])
txtplot(predict(px_lm)[-1],prot)


#but applying this gives us bullshit:
diag_z = diag(rep(mydeg,nrow(mybs)))
yv = ginv(mydbs) %*% ((edbasis +  diag_z %*% mybs )) %*% px
mydbs %*% yv
txtplot(yv)
txtplot(yv,ribo)

#what if we instead require non negativity?
#Nah how. 
px_lm = lm(c(startchange,prot) ~ 0+ edbasis[1,],mybs)
px = px_lm$coef



#One night later I'm an idiot - mybs needs to be one parameter higher than the tp, then it determines
dP and thus RNA

#Okay sooooo, I wondder if theere's. problem with how I 

################################################################################
########NOW with loog
################################################################################
	
splfit = lm(log(ribo) ~ 0+mydbs)
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=TRUE)
prot <- num_ode_res$P
startchange <- (num_ode_res$ribo[1]*myrTE - prot[1]*mydeg) / prot[1]
txtplot(num_ode_res$ribo)
txtplot(num_ode_res$P)
px_lm = lm(c(startchange,log(prot)) ~ 0  + rbind(edbasis[1,],mybs))
px = px_lm$coef

txtplot(predict(px_lm)[-1])
txtplot(predict(px_lm)[-1],log(prot))
#we face the basic prooblem that our method is not propeerly coontraainted.

(mybs %*% px) + log(edbasis%*%(px + mydeg))






thalf = log(2)/deg
deg = log(2)/thalf

#what if we parametrize instead with
steady = Ks/deg




dP = R Ks - deg P
dP/deg = R Ks / deg - P
dP/deg = R steady - P
dP/deg = R steady - P
dP/deg = R steady - P
dP/deg + P = R steady
R steady = dP/deg + P
log(R steady) = log( dP/deg + P )
log(R) + log(steady) = log(dP/deg + P)
log(R)  = log(dP/deg + P) - log(steady)
log(R)  = log(P * dlP/deg + P) - log(steady)
log(R)  = log(P * (dlP/deg + 1)) - log(steady)
log(R)  = log(P * (dlP/deg + 1)) - log(steady)

#IF deg is high then
log(R)  = log(P * (dlP/deg + 1)) - log(steady)
		~= log(P * (0 + 1)) - log(steady)
		~= log(P) - log(steady)

#But if deg is low then
log(R)  = log(P) + log(dlP/deg + 1) - log(steady)
log(R)  = log(P * (dlP/deg + 1)) - log(steady)


#or with half life
log(R)  = log(P * (dlP*thalf/log(2)) + 1) - log(steady)

log(R)  = log(P) + log(dlP*(thalf/log(2)) + 1) - log(steady)

#
log(R)  = log(P) + log(dlP*(thalf/log(2)) + 1) - log(steady)

#
log(R)  = log(P) + log(dlP*(thalf/log(2)) + 1) - log(steady)


