lrTECONST <- 19
  ################################################################################
  ########Draw trajetory plot for saved csv files 
  ################################################################################
  source('src/R/modeling/stan_predict_impute.R')
 
standata<-allstandata
igene='Satb2'

  plotdata <- function(standata,igene){
    allstandata$MS[,,igene]%>%as.data.frame%>%set_colnames(tps)%>%mutate(rep=1:3)%>%gather(time,signal,-rep)%>%select(time,signal)
  }

modelfiles <- Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/stansamples/allhierarch_[0-9].csv')
modelfileslin <- Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/stansamples/alllinear_[0-9].csv')
modelfiles_rna <- Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/stansamples/allhierarch_rna_[0-9].csv')
modelfileslin_rna <- Sys.glob('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/stansamples/alllinear_rna_[0-9].csv')

igenename = "Satb2"
pars <- 'prot'

stanfiles<-modelfiles
geneind=3016

get_genepars <- function(stanfiles,pars,geneind){

  #get the parameter names from the file
  filepars <- fread(str_interp('head -n 100 ${stanfiles[1]} | grep -e "lp__,"   '))%>%colnames
  #pars their indexes
  pardf<-filepars%>%vparse_stan_pars
  pardf$n <- 1:nrow(pardf)

  #select the ones for our gene
  pardf <- pardf%>%filter(replace_na(geneind%in%geneind,TRUE))
  #select the ones for our parameter
  pardf <- pardf%>%filter(parameter %in% pars)

  stopifnot(nrow(pardf)>0)
  

  #how much to skip
  stanheadlength <- system(intern=T,str_interp("sed -n '0,/# Diagonal/p' ${stanfiles[1]} | wc -l"))%>%as.numeric%>%add(1)
  # stanhead <- readLines(stanfile,stanheadlength)%>%str_extract('.{0,20}')
  # stanfile<-stanfiles[1]

  out <- map(stanfiles,~fread(cmd=paste('grep -v "#"',.),select=pardf$n,header=T))%>%
    bind_rows%>%
    mutate(sample=1:nrow(.))%>%
    gather(par,signal,-sample)%>%
    inner_join(pardf,by='par')%>%
    select(signal,time=timeind,parameter,sample)
  
  out%>%head

}

igene='Orc3'
plot_hiearachtraj<-function(igene){

  igeneind = which(allstandata$MS%>%dimnames%>%.[[3]]%>%`==`(igene))
  stopifnot(length(igeneind)==1)
  assay2model <- c('ribo'='Riboseq Data','total'='RNAseq Data')
  model2color <- c('Kinetic'='red','Linear'='blue','MS Data'='black','rTE: 0'='purple','Riboseq Data'='dark green','RNAseq Data'='purple')
    
  #Mass spec actual data
  msdata<-plotdata(standata,igene)
  #estimates form the models
  message('fetching protein estimates from model csvs')
  protsamples <- get_genepars( modelfiles,'prot',igeneind)
  message('fetching protein estimates from linear model csvs')
  protsampleslin <- get_genepars(modelfileslin,'prot',igeneind)

  #and the ocunt data
  seqdataused <- preddf%>%
    filter(gene_name==igene,assay%in%c('total','ribo'),!is.na(predicted_signal_spline_3))%>%
    mutate(signal=2^predicted_signal_spline_3,model=assay2model[assay])
  seqdatapoints <- exprdata%>%
    filter(gene_name==igene,assay%in%c('total','ribo'),!is.na(signal))%>%
    mutate(signal=2^signal,model=assay2model[assay])

  #plot showing trajectories and fits with MCMC samples 
  msdata%<>%mutate(model='MS Data')

  # browser()

  trajectoryplot<-ggplot(
        data=msdata,
      aes(color=model,y=log2(signal),x=as.numeric(as_factor(time)))
    )+
    stat_summary(geom='line',fun.y=mean)+
    geom_point()+
    geom_line(size=I(1),data=protsamples%>%group_by(time)%>%summarise(signal=mean(signal))%>%mutate(model='Kinetic'),linetype=1)+
    geom_line(size=I(1),data=protsampleslin%>%group_by(time)%>%summarise(signal=mean(signal))%>%mutate(model='Linear'),linetype=1)+
    geom_line(alpha=I(0.005),data=protsampleslin%>%mutate(model='Linear'),aes(group=sample))+
    geom_line(alpha=I(0.005),data=protsamples%>%mutate(model='Kinetic'),aes(group=sample))+
    scale_x_continuous(name='Stage',labels=tps)+
    scale_y_continuous(name='Log2 LFQ / Log2 Normalized Counts')+
    scale_color_manual(values = model2color)+
    theme_bw()+
    geom_line(data=seqdataused%>%mutate(time=as.numeric(as_factor(time))))+
    facet_grid(scale='free',ifelse(model %in% assay2model,'Seq Data','MS') ~ . )+
    geom_point(data=seqdatapoints%>%mutate(time=as.numeric(as_factor(time))))+
    ggtitle(label = str_interp('Linear vs. Kinetic Model - ${igene}'),
      sub="Faded Lines represent samples from Posterior\nSolid Line is Median Value")
    #=
  # ggsave(plot=trajectoryplot,file=plotfile,w=7*1.5,h=7)
    #

  message('DEBUG')
  for(v in ls()) assign(v,get(v),.GlobalEnv)


  # trajectoryplot
  gglayout = c(1,1,2,3)%>%matrix(ncol=2)
  

  arrangedplot<-gridExtra::grid.arrange(
    trajectoryplot,
    get_genepars( modelfiles,'lrTE',igeneind)%>%.$signal%>%add(lrTECONST)%>%qplot(bins=50,main=str_interp('Posterior Distribution log rTE\n${igene}'))+theme_bw(),
    get_genepars( modelfiles,'ldeg',igeneind)%>%.$signal%>%qplot(bins=50,main=str_interp('Posterior Distribution - log(degredation / 1.5 days) ${igene}'))+theme_bw(),
  layout_matrix=gglayout)

  plotfile <- here(str_interp('plots/modelling/allhierarch/${igene}.pdf'))
  plotfile%>%dirname%>%dir.create(rec=T)
  ggsave(plot=arrangedplot,file=plotfile,w=7*1.5,h=7)
  message(plotfile)

}

out<-get_genepars(modelfiles,'prot',3016)
get_genepars(modelfiles,'prot',3016)


# plot_hiearachtraj('Satb2')
# plot_hiearachtraj('Rasa3')
# plot_hiearachtraj('Ewsr1')
# plot_hiearachtraj('Cul2')







#' sab2 increasing, linear fit, kind of bimodal posterior now.
#' Orc3 weird zig zag in the riboseq that the model fits by just completely disregarding it.... This may be an instance of the riboseq signal itself being weird, rnaseq is slightly better...
#' Orc3 and now looks fucked up??? Why is linear mod3el so off??
#' being fucked up
#' Ewsr1 - another fucking weird, zig-zaggy riboseq signal, that the model just chucks away. Hierarch looks pretty good.
#' 'Ap1' nice clear signal for rTE, increasing, didn't go this time
#' 'Rasa3' - looks like shit. The linear fit is one step off - as if the MS just takes a single stage to pick up...
#' Rap1b looks fucked' - looks like shit. In this case the degredation model works great but only by fitting zero degredation to what is clearly just unrelated Ribo and RNA - ribo goes way down at P0, protein doesn't.
#' Zc2hc1a  looks kind of like an okay fit, although rTE is a little bimodal, still
#' #' Rap1b looks fucked' - looks like shit. In this case the degredation model works great but only by fitting zero degredation to what is clearly just unrelated Ribo and RNA - ribo goes way down at P0, protein doesn't.
#' Acadvl - rTE sucks, cos again, the Riboseq doesn't really bear much relationship to the mS...
#' Cul2 - once again, bananas!! wtf is happening here. RNAseq maaaayyyyyybe looks like it mtches better but... not that much
#' Zmym3 - again something weird...
#' Trmt1l - wtf???? So little correlation between the MS and the other two... I don't know what the fuck is happening...
#' 'Hspa12a' - is increasing, and looks correct.
#' Rps6ka5 - also increasing, also looks correct
#' Stx1b increasing - a case in which it still kind of looks like the RNAseq fits better, just like, compressed a bit - making me wonder if tp specific factors + linear would work best
#' Srcin1 Same phenomenon - I wonder if these guys tend to have like, 
#' Ak1 Same phenomenon - again
#' Tbc1d24 Same phenomenon - again
#' Cald1 decreases - again it's like, linear fit pretty good, but just underestimates the effect that the Riboseq has on the MS
#' Phactr4 too

