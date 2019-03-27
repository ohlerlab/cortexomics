
#Plot our model with select paramters - we get divergences
g2plot=1
allpars <- stanfit%>%as.data.frame%>%colnames
parstoplot <- allpars%>%
  grep(inv=F,val=T,patt=str_interp('\\[${g2plot}\\]|\\,${g2plot}\\]'))%>%
  grep(inv=T,val=T,patt='prot')
parstoplot%<>%append(allpars%>%{.[!str_detect(.,'\\[')]})
parstoplot%<>%setdiff(c('deg','ms0logratio[1]'))
print(pairs(stanfit,pars=parstoplot))

#rte
rTEints <- summary(stanfit)%>%as.data.frame%>%rownames_to_column('parameter')%>%filter(parameter%>%str_detect('rTE\\[')) %>%
  select(summary.mean,summary.2.5.,summary.97.5.) %>%
  mutate(actual = simdata$rTE) %>%
  rowwise%>%
  mutate(accurate=between(actual,summary.2.5.,summary.97.5.))%>%
  ungroup%>%
  mutate(type=names(ribopatterns)[ribopatinds])%>%
  mutate(dit = summary.mean-actual)
#the rTE values are consistently under estimated
print(rTEints)

#Consistently 
ldegints <- summary(stanfit)%>%as.data.frame%>%rownames_to_column('parameter')%>%filter(parameter%>%str_detect('deg'))%>%
  select(summary.mean,summary.2.5.,summary.97.5.)%>%
  mutate(actual = exp(simdata$ldeg)) %>%
  rowwise%>%
  mutate(accurate=between(actual,summary.2.5.,summary.97.5.))%>%
  ungroup%>%
  mutate(type=names(ribopatterns)[ribopatinds])%>%
  mutate(dit = summary.mean-actual)
#Not very accurate either
print(ldegints)

#See what happens when the model is at the actual values I used
stanfit%>%as.data.frame%>%rownames_to_column('parameter')%>%
  filter(between(`ldeg[1]`,simdata$ldeg[1]*1.1,0))%>%
  filter(between(`ldeg[1]`,simdata$ldeg[1]*1.1,simdata$ldeg[1]*0.9))%>%

    filter(between(`MS0[1]`,simdata$prot0[1]*0.9,simdata$prot0[1]*1.1))%>%
  identity
##....they just aren't there. 


