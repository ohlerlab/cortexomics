library(tidyverse)
library(data.table)
library(magrittr)
library(stats4)
library(parallel)
library(ggExtra)


#modify a function so that some of it's arguments are taken from the elements of a 'par' vector, for optim
with_vectargs <- function (FUN, ... ) {
   vect <- map_chr(rlang::ensyms(...),as.character)
   stopifnot(all(vect %in% names(formals(FUN))))
   nonvectargs <- setdiff(vect,names(formals(FUN)))

   #modify the function as having a vector as it's first arg
   .FUN <- FUN
   formals(.FUN) <- as.pairlist(c(list(par=quote(NULL)),as.list(formals(.FUN))))

   #now assign the default for each of our vector arguments, as that element of the vector
   for(i in seq_along(vect)){
       vect_argname <- vect[i]
       formals(.FUN)[[vect_argname]] <- parse(text=paste0('par[',i,']',collapse=''))[[1]]
   }
   
   .FUN
}k

mlemod <- function(minuslogl, start = formals(minuslogl), method = "BFGS",
    fixed = list(), nobs, ...){
    call <- match.call()
    n <- names(fixed)
    fullcoef <- formals(minuslogl)
    if (any(!n %in% names(fullcoef)))
        stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
    fullcoef[n] <- fixed
    if (!missing(start) && (!is.list(start) || is.null(names(start))))
        stop("'start' must be a named list")
    start[n] <- NULL
    start <- unlist(lapply(start, eval.parent))
    nm <- names(start)
    oo <- match(nm, names(fullcoef))
    if (anyNA(oo))
        stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
    start <- start[order(oo)]
    nm <- names(start)
    f <- function(p) {
        l <- as.list(p)
        names(l) <- nm
        l[n] <- fixed
        do.call("minuslogl", l)
    }
    oout <- if (length(start))
        optim(start, f, method = method, hessian = TRUE, ...)
    else list(par = numeric(), value = f(start))
    coef <- oout$par
    vcov <- if (length(coef))
        solve(oout$hessian)
    else matrix(numeric(), 0L, 0L)
    min <- oout$value
    fullcoef[nm] <- coef
    new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef),
        vcov = vcov, min = min, details = oout, minuslogl = minuslogl,
        nobs = if (missing(nobs))
            NA_integer_
        else nobs, method = method)
}

nLL_model <- function(ldeg,prot0,ribo,MS,rTE,ms_params){
   #our degredation constant gets optimized in log space, but used in linear space
   deg = exp(ldeg)   
   #Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
   prot = MS[,,1,drop=FALSE]
   dim(prot)=dim(prot)[1:2]
   prot[] <- 0
   prot[,1] <- prot0
   #build up our protein array tp by tp
   i=2
   for (i in 2:ncol(ribo)){
      prot[,i] <- prot[,i-1,drop=FALSE] + (rTE*ribo[,i,drop=FALSE]) - (prot[,i-1,drop=FALSE]*deg)
   }
   #transform our MS to log scale, using our parameters
   l2ms_params<-get_l2_ms_params(prot[1,],ms_params$cvslope,ms_params$cvint)

   #finally get our log likelihood
   out <- 
   -sum(
      dnorm(log2(MS+1),mean=l2ms_params$u,sd=l2ms_params$s,log=TRUE)
   ,na.rm=TRUE)
   # if(!is.finite(out)) browser()
   # browser()
   out
}


#the model but with variance as a free parameter to fit
nLL_model_svar <- function(ldeg,prot0,ribo,MS,rTE,msvar){
   #our degredation constant gets optimized in log space, but used in linear space
   deg = exp(ldeg)   
   #Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
   prot = MS[,,1,drop=FALSE]
   dim(prot)=dim(prot)[1:2]
   prot[] <- 0
   prot[,1] <- prot0
   #build up our protein array tp by tp
   i=2
   for (i in 2:ncol(ribo)){
      prot[,i] <- prot[,i-1,drop=FALSE] + (rTE*ribo[,i,drop=FALSE]) - (prot[,i-1,drop=FALSE]*deg)
   }

   #finally get our log likelihood
   out <- 
   -sum(
      dnorm(log2(MS+1),mean=log2(prot),sd=msvar,log=TRUE)
   ,na.rm=TRUE)
   # if(!is.finite(out)) browser()
   # browser()
   out
}

#Model's plotting 
nLL_model_plot <- function(ldeg,prot0,ribo,MS,rTE,ms_params){
   #our degredation constant gets optimized in log space, but used in linear space
   deg = exp(ldeg)   
   #Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
   prot = MS[,,1,drop=FALSE]
   dim(prot)=dim(prot)[1:2]
   prot[] <- 0
   prot[,1] <- prot0
   #build up our protein array tp by tp
   i=2
   for (i in 2:ncol(ribo)){
      prot[,i] <- prot[,i-1,drop=FALSE] + (rTE*ribo[,i,drop=FALSE]) - (prot[,i-1,drop=FALSE]*deg)
   }
   #transform our MS to log scale
   # prot <- log2(prot+1)
   #the sd will be a linear func of the mean, as per plots
   # ms_sd <- sweep(prot,2,ms_sd,FUN='*')
   #finally get our log likelihood
   #transform our MS to log scale, using our parameters
   l2ms_params<-get_l2_ms_params(prot[1,],ms_params$cvslope,ms_params$cvint)

  
   melt(MS)%>%set_colnames(c('gene_id','time','rep','signal'))%>%
   mutate(var='MS',pred=FALSE,signal=log2(signal+1))%>%
   {bind_rows(.,data.frame(gene_id=1,time=1:length(l2ms_params$u),rep=1,var='MS',pred=TRUE,signal=l2ms_params$u) )}%>%
   {bind_rows(.,data.frame(gene_id=1,time=1:length(ribo),pred=FALSE,rep=1,var='ribo',signal=log2(ribo[1,]) ))}%>%
   mutate(ymin=ifelse(pred==TRUE,signal-2*l2ms_params$s,NA),ymax=ifelse(pred==TRUE,signal+2*l2ms_params$s,NA))

   # ggsave(file='../plots/modelling/example_ms_fit.pdf'%TRUE>%{normalizePath(.)%>%message})
}


#the model but with variance as a free parameter to fit
nLL_model_svar_plot <- function(ldeg,prot0,ribo,MS,rTE,msvar){
   #our degredation constant gets optimized in log space, but used in linear space
   deg = exp(ldeg)   
   #Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
   prot = MS[,,1,drop=FALSE]
   dim(prot)=dim(prot)[1:2]
   prot[] <- 0
   prot[,1] <- prot0
   #build up our protein array tp by tp
   i=2
   for (i in 2:ncol(ribo)){
      prot[,i] <- prot[,i-1,drop=FALSE] + (rTE*ribo[,i,drop=FALSE]) - (prot[,i-1,drop=FALSE]*deg)
   }

   melt(MS)%>%set_colnames(c('gene_id','time','rep','signal'))%>%
   mutate(var='MS',pred=FALSE,signal=log2(signal+1))%>%
   {bind_rows(.,data.frame(gene_id=1,time=1:length(log2(prot[1,])),rep=1,var='MS',pred=TRUE,signal=log2(prot[1,])) )}%>%
   {bind_rows(.,data.frame(gene_id=1,time=1:length(ribo),pred=FALSE,rep=1,var='ribo',signal=log2(ribo[1,]) ))}%>%
   mutate(ymin=ifelse(pred==TRUE,signal-2*msvar,NA),ymax=ifelse(pred==TRUE,signal+2*msvar,NA))
}



#########Load up data
exprdata_all<-fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/exprdata/transformed_data.txt')
mscols <- exprdata_all%>%colnames%>%str_subset('MS_')
ribocols <- exprdata_all%>%colnames%>%str_subset('ribo_')


ms_tall<-exprdata_all%>%
   select(one_of(mscols),gene_name)%>%
   select(1:3,gene_name)%>%
   gather(dataset,value,-gene_name)%>%
   separate(dataset,c('time','assay','rep'))%>%
   group_by(time,gene_name)

ms_meansd <- ms_tall%>% summarise(mean=mean(na.omit(value)),sd=sd(na.omit(value)))

ms_meansd%>%filter(sd>0.5)%>%ungroup%>%arrange(gene_name)%>%as.data.frame%>%left_join(ms_tall)
#there is some higher var at the low end of the protein scores as well....
ms_meansd%>%ggplot(aes(x=sd))+geom_histogram()+coord_cartesian(xlim=c(0,0.5))
ms_meansd%>%ggplot(aes(x=mean,y=sd))+geom_point()+geom_smooth()

exprdatareshape<-exprdata_all%>%select(gene_name,one_of(c(ribocols,mscols)))%>%
   gather(dataset,value,-gene_name)%>%
   mutate(value = 2^value)%>% ####NOTE I'm de-logging it here... probably not good
   separate(dataset,c('time','assay','rep'))%>%
   identity%>%
   group_by(gene_name,time)%>%mutate(datname=paste0(assay,rep))%>%
   select(-assay,-rep)%>%
   spread(datname,value)

'Satb2' %in% exprdata_all$gene_name
'Satb2' %in% exprdatareshape$gene_name

#Function to get the ribo vector we use to infer the protein levels, can mod later with splines
get_trans_vect <- function(ribodata){
   ribodata%>%as.matrix%>%apply(1,mean,na.rm=TRUE)
}
exprdata <- exprdatareshape
exprdata$ribo <- get_trans_vect(exprdata%>%ungroup%>%select(matches('ribo')))
exprdata$gene_id <-exprdata$gene_name

mscols <- colnames(exprdatareshape)%>% str_subset('MS')

#expr_array <-
expr_array=exprdata%>%ungroup%>%select(gene_id,one_of(mscols))%>%
   # head(10)%>%
   group_by(gene_id)%>%
   nest%>%mutate(data=map(data,as.matrix))%>%
   # {setNames(.$data,.$gene_id)}%>%
   pluck('data')%>%
   {array(unlist(.),dim=c(nrow(.[[1]]),ncol(.[[1]]),length(.)))}%>%
   aperm(c(3,1,2))

ribo_matrix <- exprdata%>%{split(.$ribo,.$gene_id)}%>%do.call(what=rbind)

satb2ind <- exprdata$gene_name%>%unique%>%`==`('Satb2')%>%which







#Calculates for MS variance

#Load our LFQ values, ge the log2 mean and log2 sd for the complete cases
sdplotdf<-'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/ms_tables/ms_LFQ_total_ms_tall.tsv'%>%
   read_tsv%>%
   filter(fraction=='total')%>%
   group_by(time,fraction,Protein_IDs)%>%
   filter(length(signal)==3)%>%
   filter(!any(is.na(signal)))%>%
   summarise(mean_sig=log2(mean(signal)),sd_signal = log2(sd(signal)))

sdplotdftrans<-exprdata%>%ungroup%>%
   select(gene_id,time,one_of(mscols))%>%
   gather(rep,signal,-time,-gene_id)%>%
   group_by(time,gene_id)%>%
   filter(all(!is.na(signal)))%>%
   summarise(mean_sig=mean(signal),sd_signal = sd(signal))

#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{
{
   sdplotdf%>% 
   ggplot(aes(x=log2(mean_sig),y=log2(sd_signal)))+
   # geom_bkde2d()+
   facet_wrap(~time)+
   # geom_point(alpha=I(0.1))+
   geom_smooth(method = "lm", se = FALSE)+
   stat_density_2d(aes(fill = ..level..), geom = "polygon")+
   geom_abline(slope=1)+
   theme_bw()
}
})

pdfexpr('../plots/modelling/mass_spec_sdplot_trans.pdf',{
{
   sdplotdftrans%>% 
   ggplot(aes(x=mean_sig,y=sd_signal))+
   # geom_bkde2d()+
   facet_wrap(~time)+
   # geom_point(alpha=I(0.1))+
   geom_smooth(method = "lm", se = FALSE)+
   stat_density_2d(aes(fill = ..level..), geom = "polygon")+
   geom_abline(slope=1)+
   theme_bw()
}
})



#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{
{
   sdplotdf%>% 
   ggplot(aes(x=(mean_sig),y=(sd_signal)))+
   # geom_bkde2d()+
   facet_wrap(~time)+
   # geom_point(alpha=I(0.1))+
   geom_smooth(method = "lm", se = FALSE,data=sdplotdf%>%filter(between(mean_sig,27,32) ))+
   stat_density_2d(aes(fill = ..level..), geom = "polygon")+
   geom_abline(slope=1)+
   theme_bw()
}
})



#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot_all.pdf',{
{
   sdplotdf%>% 
   ggplot(aes(x=(mean_sig),y=(sd_signal)))+
   # geom_bkde2d()+
   geom_point(alpha=I(0.1))+
   geom_smooth(method = "lm", se = TRUE)+
   # stat_density_2d(aes(fill = ..level..), geom = "polygon")+
   geom_abline(slope=1)+
   theme_bw()
}
})

#seems to be a predictable enough linear function, let's get those fits

sd_mean_lmfits <- sdplotdf%>%group_by(time)%>%nest%>%mutate(fit = map(data, ~ lm (data=. , sd_signal ~ mean_sig)) )
meansiglopes <- sd_mean_lmfits%>%mutate(l2mean_sig_slp =  fit%>%map_dbl(~ .$coef['mean_sig']))%>%{setNames(.$l2mean_sig_slp,.$time)}
meansigint <- sd_mean_lmfits%>%mutate(l2mean_sig_slp =  fit%>%map_dbl(~ .$coef['Intercept']))%>%{setNames(.$l2mean_sig_slp,.$time)}

get_l2_ms_params <- function(m,cvslope,cvint){
   l2m = log2(m)
   sds = 2^((l2m*cvslope)+cvint)
   s = sqrt(log2(1 + (sds^2) / (m^2) ))
   u = log2(m/(sqrt(1+(sds)/(m^2))))
   list(u=u,s=s)
}

ms_params <- 
 list(
   cvslope=sd_mean_lmfits$fit%>%map_dbl(~.$coef['mean_sig'])%>%setNames(sd_mean_lmfits$time),
   cvint=sd_mean_lmfits$fit%>%map_dbl(~.$coef['(Intercept)'])%>%setNames(sd_mean_lmfits$time)
)



# cvslope = sd_mean_lmfits$fit[[1]]$coef['mean_sig']
# cvint = sd_mean_lmfits$fit[[1]]$coef['(Intercept)']
# l2means = (sdplotdf%>%filter(time=='E13')%>%.$mean_sig)

# means = 2^(sdplotdf%>%filter(time=='E13')%>%.$mean_sig)
# sds = 2^((l2means*cvslope)+cvint)
# s = sqrt(log2(1 + (sds^2) / (means^2) ))
# u = log2(means/(sqrt(1+(sds)/(means^2))))
# simdata=replicate(3,rnorm(n=length(u),mean=u,sd=s))
# # simdatamean = simdata%>%apply(1,mean)
# simdatasd = 2^simdata%>%apply(1,sd)
#get_l2_ms_params(means[1:5],ms_params$cvslope,ms_params$cvint)








#Now test out model calls. 
igenevect=c(1:2)
exprdata%>%filter(gene_name=='Satb2')
expr_array[3021,,]

nLL_model(
   ribo=ribo_matrix[igenevect,,drop=FALSE],
   MS=expr_array[igenevect,,,drop=FALSE],
   ldeg=rep(log(0.5),length(igenevect)),
   rTE=100,
   prot0=expr_array[igenevect,1,1],
   ms_params=ms_params
)

nLL_model(
   ribo=ribo_matrix[10,,drop=FALSE],
   MS=expr_array[10,,,drop=FALSE],
   ldeg=c(log(0.5),log(0.1))[1],
   rTE=100,
   prot0=rep(30986,2)[1],
   ms_params=ms_params
)

optim%>%args

#optimize with free rTE
fTEfit <- optim(
   fn=with_vectargs(nLL_model,ldeg,prot0,rTE),
   par=c(log(0.1),expr_array[satb2ind,,][[1]],100),
   ribo=ribo_matrix[satb2ind,,drop=FALSE],
   MS=expr_array[satb2ind,,,drop=FALSE],
   ms_params=ms_params,
   method='L-BFGS-B',
   lower = c(-Inf,1,1e-12),
   upper= c(0,+Inf,+Inf),
   control=list(trace=1)
)

#optimize constrained rTE
constr_fit<-optim(
   fn=with_vectargs(nLL_model,ldeg,prot0),
   par=c(log(0.1),expr_array[satb2ind,,][[1]]),
   ribo=ribo_matrix[satb2ind,,drop=FALSE],
   MS=expr_array[satb2ind,,,drop=FALSE],
   rTE=100,
   ms_params=ms_params,
   method='L-BFGS-B',
   lower = c(-Inf,1),
   upper= c(0,+Inf)
)

plot_LL<-function(.){
   ggplot(.,aes(y=signal,x=time))+
   geom_point(data=filter(.,!pred))+
   #facet_grid(scale='free',var ~ .)
   geom_ribbon(data=filter(.,pred),aes(y=signal,ymax=ymax,ymin=ymin,alpha=0.1))+
   geom_line(data=filter(.,pred),alpha=0.1)+
   facet_grid(scale='free',var ~. )
}

#plot the free TE fit
with_vectargs(nLL_model_plot,ldeg,prot0,rTE)(
   par=fTEfit$par,
   ribo=ribo_matrix[satb2ind,,drop=FALSE],
   MS=expr_array[satb2ind,,,drop=FALSE],
   ms_params=ms_params,
)%>%plot_LL+ggtitle('free fit',sub=paste0(c('ldeg','prot0','rTE'),' = ',round(fTEfit$par,3),collapse=';'))


with_vectargs(nLL_model_plot,ldeg,prot0,rTE)(
   par=c(constr_fit[1:2,]),
   ribo=ribo_matrix[satb2ind,,drop=FALSE],
   MS=expr_array[satb2ind,,,drop=FALSE],
   ms_params=ms_params,
)

with_vectargs(nLL_model_plot,ldeg,prot0)(
   par=constr_fit$par,
   ribo=ribo_matrix[satb2ind,,drop=FALSE],
   MS=expr_array[satb2ind,,,drop=FALSE],
   rTE=100,
   ms_params=ms_params,
)%>%plot_LL+ggtitle('constrained fit rTE=100',sub=paste0(c('ldeg','prot0'),' = ',round(fTEfit$par,3),collapse=';'))




############Now optimize with free TE for all genes
#optimize with free rTE
ginds = 1:nrow(expr_array)
ginds = sample(ginds,10)
gnames <- exprdata$gene_id%>%unique%>%.[ginds]


###With optim, this works
fTEfits <- mclapply( ginds ,function(gind) {safely(optim)(
   fn=with_vectargs(nLL_model,ldeg,prot0,rTE),
   par=c(log(0.1),expr_array[gind,,][[1]],100),
   ribo=ribo_matrix[gind,,drop=FALSE],
   MS=expr_array[gind,,,drop=FALSE],
   ms_params=ms_params,
   method='L-BFGS-B',
   lower = c(-20,1,1e-12),
   upper= c(0,+Inf,+Inf),
   hessian=TRUE,
   control=list(trace=1)
   )
})

###With optim, and free variance
fTEfits_svar <- mclapply( ginds ,function(gind) {safely(optim)(
   fn=with_vectargs(nLL_model_svar,ldeg,prot0,rTE,msvar),
   par=c(log(0.1),expr_array[gind,,][[1]],1,1),
   ribo=ribo_matrix[gind,,drop=FALSE],
   MS=expr_array[gind,,,drop=FALSE],
   method='L-BFGS-B',
   lower = c(-20,1,1e-12,1e-12),
   upper= c(0,+Inf,+Inf,10),
   hessian=TRUE,
   control=list(trace=1)
   )
})

#mle works, but then I can't get confints
fTEfitsmle <- mclapply( ginds ,function(gind) {
   safely(mle)(
   minuslogl=nLL_model,
   start=list(ldeg=log(0.1),prot0=expr_array[gind,,][[1]],rTE=100),
   fixed=list(
      ribo=ribo_matrix[gind,,drop=FALSE],
      MS=expr_array[gind,,,drop=FALSE],
      ms_params=ms_params
   ),
   method='L-BFGS-B',
   lower = c(-20,1,1e-12),
   upper= c(0,+Inf,+Inf),
   )
})
fTEfitsmle%<>%setNames(gnames)

#but how to get confints???
confint(fTEfitsmle$Satb2$result)
fitconfints <- map(fTEfits,'result')%>%keep(Negate(is.null))%>%map(safely(confint))
fitconfints%>%map('result')%>%keep(Negate(is.null))

solve(fTEfits$Satb2$result$hessian)

failedfits <- fTEfits%>%map('error')%>%map_lgl(Negate(is.null))%>%which

fTEfits <- mclapply( ginds[1] ,function(gind) {
   (bbmle::mle2)(
   minuslogl=nLL_model,
   start=list(ldeg=log(0.1),prot0=expr_array[gind,,][[1]],rTE=100),
   fixed=list(
      ribo=ribo_matrix[gind,,drop=FALSE],
      MS=expr_array[gind,,,drop=FALSE],
      ms_params=ms_params
   ),
   method='L-BFGS-B',
   lower = c(-20,1,1e-12),
   upper= c(0,+Inf,+Inf),
   )
})
options(error=recover)
fTEfits%<>%setNames(gnames)


#####analyzing results..

#get estimates of parameter values
freeest_df<-fTEfits%>%
   map('result')%>%
   keep(Negate(is.null))%>%
   map('par')%>%
   map(setNames,c('ldeg','prot0','rTE'))%>%
   map(enframe,name='var',value='est')%>%
   bind_rows(.id='gene_id')%>%
   spread(var,est)%>%
   mutate(Gene_Name=fct_other(gene_id,'Satb2'))%TRUE>%print
 

#get estimates of parameter values
freeest_df<-fTEfits_svar%>%
   map('result')%>%
   keep(Negate(is.null))%>%
   map('par')%>%
   map(setNames,c('ldeg','prot0','rTE','msvar'))%>%
   map(enframe,name='var',value='est')%>%
   bind_rows(.id='gene_id')%>%
   spread(var,est)%>%
   mutate(Gene_Name=fct_other(gene_id,'Satb2'))%TRUE>%print
 
freeest_df$rTE%>%log10%>%keep(~ . > -5) %>% hist(50)

freeest_df %>%  {
      ggMarginal(ggplot(.,aes(y=exp(ldeg),x=log(rTE),color=Gene_Name))+geom_point(alpha=I(0.1)),type='histogram')
}


freeest_df %>%  {
      ggMarginal(ggplot(.,aes(y=log10(prot0),x=log((msvar)^2),color=Gene_Name))+geom_point(alpha=I(0.1)),type='histogram')
}

igeneind <- freeest_df%>%filter(ldeg==0)%>%.$gene_id%>%sample(1)%>%match(gnames)

#plot the free TE fit
igeneind=7

with_vectargs(nLL_model_plot,ldeg,prot0,rTE)(
   par=fTEfits[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
   ms_params=ms_params,
)%>%plot_LL + ggtitle('free fit',sub=paste0(c('ldeg','prot0','rTE'),' = ',round(fTEfits[[igeneind]]$result$par,3),collapse=';'))

#find the one with the weird ribo at timepoint 2
igeneind <- exprdata%>%group_by(gene_name,time)%>%summarise(time_sd=sd(c(ribo1,ribo2)))%>%
   summarise(oratio =  time_sd[time=='E145'] / 4*(mean(time_sd)))%>%ungroup%>%slice(which.max(oratio))%>%
   .$gene_name%>%unique%TRUE>%print%>%
   `==`(gnames,.)%>%which
0
options(error=NULL)
testgene <- gnames[igeneind]
exprdata%>%filter(gene_id==testgene)
expr_array[igeneind,,]
ribo_matrix[igeneind,]
fTEfits[[igeneind]]
fTEfits[igeneind]
freeest_df%>%filter(gene_id==testgene)

0


igeneind = which(gnames==ribodevgenes[3])

with_vectargs(nLL_model_plot,ldeg,prot0,rTE)(
   par=fTEfits[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
   ms_params=ms_params,
)%>%plot_LL + ggtitle('free fit',sub=paste0(c('ldeg','prot0','rTE'),' = ',round(fTEfits[[igeneind]]$result$par,3),collapse=';'))



igeneind = 


with_vectargs(nLL_model_svar_plot,ldeg,prot0,rTE,msvar)(
   par=fTEfits_svar[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
)%>%plot_LL + ggtitle('free fit',sub=paste0(c('ldeg','prot0','rTE'),' = ',round(fTEfits_svar[[igeneind]]$result$par,3),collapse=';'))


freeest_df%>%group_by(gene_id)%>%nest%>%slice(igeneind)

#messing with creating mle objets - but tough due to the call structure
mle
oout <- fTEfits$Satb2$result
call = quote(NULL)
coef = oout$par



ribodevgenes<- exprdata%>%summarize(tmp=abs(ribo1-ribo2)/sum(abs(ribo1)))%>%filter(time=='E145')%>%ungroup%>%arrange(desc(tmp))%>%
   .$gene_name



exprdata%>%select(gene_id,time,ribo1,ribo2)




###Now let's try with a fixed RTE
fTEfits_rtefix <- mclapply( ginds ,function(gind) {safely(optim)(
   fn=with_vectargs(nLL_model_svar,ldeg,prot0,msvar),
   par=c(log(0.1),expr_array[gind,,][[1]],1),
   rTE= -1,
   ribo=ribo_matrix[gind,,drop=FALSE],
   MS=expr_array[gind,,,drop=FALSE],
   method='L-BFGS-B',
   lower = c(-20,1,0.00001),
   upper= c(0,+Inf,10),
   hessian=TRUE,
   control=list(trace=1)
   )
})


###And with the degredation rate fixed at one
fTEfits_degfix <- mclapply( ginds ,function(gind) {safely(optim)(
   fn=with_vectargs(nLL_model_svar,prot0,msvar,rTE),
   par=c(expr_array[gind,,][[1]],1,0.1),
   ldeg= 0,
   ribo=ribo_matrix[gind,,drop=FALSE],
   MS=expr_array[gind,,,drop=FALSE],
   method='L-BFGS-B',
   lower = c(1,0.00001,1,1e-12),
   upper= c(+Inf,10,+Inf),
   hessian=TRUE,
   control=list(trace=1)
   )
})

gind=igeneind

###Quick chisq tests comparing likelihood under fixed deg and 
modelcomppvals <- map(ginds,safely(function(gind){
   ms_degfix <- 2  * fTEfits_degfix[[gind]]$result$value
   ms_freedeg <- 2  * fTEfits_svar[[gind]]$result$value
   Deviance <-  ms_degfix - ms_freedeg
   1 - pchisq(Deviance,1)
}))%>%map('result')%>%setNames(gnames)%>%keep(Negate(is.null))%>%unlist

igene <- modelcomppvals%>%sort%>%tail(3)%>%head(1)%>%names
igeneind <- which(gnames==igene)

#plot the fixed TE fit
with_vectargs(nLL_model_svar_plot,ldeg,prot0,msvar)(
   par=fTEfits_rtefix[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
   rTE=-1
)%>% plot_LL + ggtitle('free fit',sub=paste0(c('ldeg','prot0','msvar'),' = ',round(fTEfits_rtefix[[igeneind]]$result$par,3),collapse=';'))

#plot the free fit
with_vectargs(nLL_model_svar_plot,ldeg,prot0,rTE,msvar)(
   par=fTEfits_svar[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
)%>%plot_LL + ggtitle('free fit',sub=paste0(c('ldeg','prot0','rTE','msvar'),' = ',round(fTEfits_svar[[igeneind]]$result$par,3),collapse=';'))
c

#and the fixed degredation
with_vectargs(nLL_model_svar_plot,prot0,msvar,rTE)(
   par=fTEfits_degfix[[igeneind]]$result$par,
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
   ldeg=0,
)%>%plot_LL + ggtitle('free fit',sub=paste0(c('prot0','msvar','rTE'),' = ',round(fTEfits_degfix[[igeneind]]$result$par,3),collapse=';'))

fTEfits_svar[[igeneind]]$value


#LL would be product, much less than 1 greater than zero
#log LL is real number, very negative



igeneind=3915
#plot the free fit
with_vectargs(nLL_model_svar,ldeg,prot0,rTE,msvar)(
   par=fTEfits_svar[[igeneind]]$result$par%>%{.},
   ribo=ribo_matrix[igeneind,,drop=FALSE],
   MS=expr_array[igeneind,,,drop=FALSE],
)
# dnormcall <- as.list(body(nLL_model))%>%{.[[length(.)-1]][[3]][[2]][[2]]}
# dnormcall<-c( quote(list),dnormcall$mean,dnormcall$sd)

# as.list(body(nLL_model))[[10]][[3]][[2]][[2]]$mean











# D <- matrix(10^c(0:4))
# Di <- matrix(10^c(-(0:4)))
# L <- lower.tri(diag(5),d=TRUE)
# o <- matrix(rep(1,5))





# t(Di) %*% L %*% D

# o  %*%  t(D) %*% L

#First choose a gene
#show it's expression data
#optimize over it's parameters
#show the fit given linear
#also show the fit given curv.
#Quantify the difference

#####So profiling is going to be an issue while I have complex fixed arguments, I think....d

#####Okay, I need to figure out 1) Sane way of parametrizing variance amongs the 

###My estimate of varianc in the ms look good but there's still a relationship between the mean and variance there.....

####Maybe I should even do teh steps on the log scale, so that the degredation factor is a constant, not sure how to do this

####Okay, maybe we need to go back and make totally sure that linear fitting can't work.... Let's find hte guys who are most
#####Strongly loaded on the "MS only" pca, based on the linear fitting..
#####Also, we should 

#I need to simulate some data and show that my tests are actually finding the ones with degredation kinetics...

#IT really looks like fixing degredation at 0 helps! So why then, doesn't optim just do that???
#Why are the values negative???
