expr<-fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/exprdata/transformed_data.txt')
ribo <- expr%>%filter(gene_name=='Satb2')%>%gather(dataset,signal)%>%tail(-1)%>%separate(dataset,c('time','assay','rep'))%>%filter(assay=='ribo')%>%group_by(time)%>%summarise(signal=median(as.numeric(signal)))%>%.$signal%>%{2^.}

deg=0.1
#ribo=c(100,120,500,600,650)
  # identity
  #sample(5)
# prot=c((ribo[1]*rTE/deg)*exp(runif(1,-2,2)),0,0,0,0)
prot=c((ribo[1]*rTE/deg),0,0,0,0)
rTE=100

degs = seq(0,1,by=0.1)
#
for(t in 2:5){
  prot[t]<-prot[t-1] + (ribo[t]*rTE) - (prot[t-1]*deg)
}
#
#
# prot = (ribo-(ribo[1]*0.3))*rTE

simdat<-data_frame(ribo,prot,time=seq_along(ribo))%>%
  mutate(ribofc = ribo/lag(ribo))%>%
  mutate(protfc = prot/lag(prot))%>%
  mutate(sprotfc = protfc-ribofc)

simdat%>%ggplot(data=.,aes(x=ribofc,y=sprotfc))+geom_point()


simdat%>%select(time,prot,ribo)%>%gather(assay,signal,-time)%>%group_by(assay)%>%mutate(signal=signal/median(na.rm = T,(signal)))%>%
    ggplot(aes(x=time,y=log10(signal),color=assay))+geom_line()+geom_point()


ribo <- ts(simdat$ribo)
prot <- ts(simdat$prot)

exprts <- cbind(ribo,prot)

plot(exprts)

VARselect(exprts, lag.max = 1, type = "none")


mod1 <- VAR(exprts, p = 1, type = "none",lag.max=1)

summary(mod1)

################################################################################
########dMod?
################################################################################
#install('dMod')
library(dMod)
library(ggplot2)

# Reactions
f <- NULL
f <- addReaction(f, 
                 from = "Enz + Sub", 
                 to = "Compl", 
                 rate = "k1*Enz*Sub",
                 description = "production of complex")
f <- addReaction(f, 
                 from = "Compl", 
                 to = "Enz + Sub", 
                 rate = "k2*Compl",
                 description = "decay of complex")
f <- addReaction(f, 
                 from = "Compl", 
                 to = "Enz + Prod", 
                 rate = "k3*Compl",
                 description = "production of product")
f <- addReaction(f, 
                 from = "Enz", 
                 to = ""     , 
                 rate = "k4*Enz",
                 description = "enzyme degradation")

# ODE model
model <- odemodel(f, modelname = "enzymeKinetics")

# Prediction function
x <- Xs(model)


observables <- eqnvec(
  product = "Prod", 
  substrate = "(Sub + Compl)", 
  enzyme = "(Enz + Compl)"
)

# Generate observation functions
g <- Y(observables, x, compile = TRUE, modelname = "obsfn", attach.input = FALSE)


# Get all parameters
innerpars <- getParameters(g*x)
# Identity transformation
trafo <- repar("x~x", x = innerpars)
# Set some initial value parameters
trafo <- repar("x~0", x = c("Compl", "Prod"), trafo)
# Explicit log-transform of all parameters
trafo <- repar("x~exp(x)", x = innerpars, trafo)

## Split transformation into two
trafo1 <- trafo2<- trafo

# Set the degradation rate in the first condition to 0
trafo1["k4"] <- "0"

# Generate parameter transformation functions
p <- NULL
p <- p + P(trafo1, condition = "noDegradation")
p <- p + P(trafo2, condition = "withDegradation")


# Initialize with randomly chosen parameters
set.seed(1)
outerpars <- getParameters(p)
pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)
times <- 0:100

modelfit <- (g*x*p)(times, pouter)

plot((g*x*p)(times, pouter))


data <- datalist(
  noDegradation = data.frame(
    name = c("product", "product", "product", "substrate", "substrate", "substrate"),
    time = c(0, 25, 100, 0, 25, 100),
    value = c(0.0025, 0.2012, 0.3080, 0.3372, 0.1662, 0.0166),
    sigma = 0.02),
  withDegradation = data.frame(
    name = c("product", "product", "product", "substrate", "substrate", "substrate", "enzyme", "enzyme", "enzyme"),
    time = c(0, 25, 100, 0, 25, 100, 0, 25, 100),
    value = c(-0.0301,  0.1512, 0.2403, 0.3013, 0.1635, 0.0411, 0.4701, 0.2001, 0.0383),
    sigma = 0.02)
)

timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))

# Compare data to prediction
plot(data) + geom_line()


plot((g*x*p)(times, pouter), data)



# Define prior values for parameters
prior <- structure(rep(0, length(pouter)), names = names(pouter))

# Set up objective function
obj <- normL2(data, g*x*p) + constraintL2(mu = prior, sigma = 10)

# Optimize the objective function
myfit <- trust(obj, pouter, rinit = 1, rmax = 10)

plot((g*x*p)(times, myfit$argument), data)