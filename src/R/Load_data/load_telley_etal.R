
#+ setup, include=FALSE, echo=FALSE, eval=T
library(rmarkdown)
library(knitr)
library(here)
knitr::opts_chunk$set(root.dir = here::here(),eval=FALSE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=FALSE)
isknitr<-isTRUE(getOption('knitr.in.progress'))
#if(!isknitr) rmarkdown::render(knitr::spin(here('src/R/load_telley_etal.R'),knit=F),output_dir=here('Reports'),knit_root_dir=here())

source('src/R/Rprofile.R')
library('DeconRNASeq')

telleyfiles<-c(
	twavetbl="/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/telley/aav2522_Data-S2.xlsx",
	tdynamicstbl="/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/telley/aav2522_Data-S4.xlsx",
	tcoretable="/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/telley/aav2522_Data-S6.xlsx",
	tclusters="/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/telley/aav2522_Data-S7.xlsx"
)
library(readxl)

modules::import(attach=TRUE,'dplyr')
conflicted::conflict_prefer('filter','dplyr')
conflicted::conflict_prefer('slice','dplyr')
conflicted::conflict_prefer('matches','dplyr')





t_wavetbl<-telleyfiles[[1]]%>%read_excel%>%rename(gene_name=`Gene symbol`)
t_dynamicstbl<-telleyfiles[[2]]%>%read_excel%>%rename(gene_name=`Gene symbol`)

telleyfiles[[3]]%>%excel_sheets

t_timecoretbl<- telleyfiles[[3]]%>%read_excel(sheet=1)%>%rename(gene_name=`Gene symbol`)%>%rename(Time_core_Specificity=Specificity)
t_diffcoretbl<- telleyfiles[[3]]%>%read_excel(sheet=2)%>%rename(gene_name=`Gene symbol`)%>%rename(Diff_core_Specificity=Specificity)

chronotypic_clusters<-telleyfiles['tclusters']%>%read_excel(sheet='Clusters')%>%rename(gene_name=`Gene symbol`)%>%rename(cluster=`tSNE cluster`)

#get info on TE
gnamevect <- fread(here('pipeline/ids.txt'))%>%{setNames(.$gene_name,.$gene_id)}


xtailfiles <- Sys.glob('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/xtail//*.txt')
names(xtailfiles)<-xtailfiles%>%str_extract('(?<=/xtail_).*?(?=.txt)')
xtaildata<-xtailfiles%>%map(fread)
xtaildata%<>%bind_rows(.id='time')
xtaildata%<>%rename('gene_id'='feature_id')
xtaildata%<>%filter(!str_detect(gene_id,'uORF'))
xtaildata%<>%safe_left_join(fread(here('pipeline/ids.txt')))
0
tegenes <-xtaildata%>%filter(log2fc>0.32,adj_p_value<0.05)%>%.$gene_name

ntegenes <- xtaildata%>%
  filter(adj_p_value>0.05)%>%
  filter(between(log2fc,-0.025,0.025))%>%
  .$gene_name


source('src/R/modeling/stan_predict_impute.R')

assert_that(exists('exprdata'))
assert_that(all(c("gene_name", "dataset", "signal", "time", "assay", "rep", "ribo",
"MS")%in% colnames(exprdata)))



sample_using_empdens <- function(sampv,size=NULL,inpvals,empvals){
		sample(sampv,size,replace=TRUE,prob =
			replace_na(approxfun(density(empvals))(inpvals),0)
			)
}

library(org.Mm.eg.db)
x <- org.Mm.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

symid <- unlist(as.list(org.Mm.egSYMBOL2EG))


getmatchingset<-function(testset,data){
	data$testvar <- data[[1]] %in% testset
	stopifnot(sum(data[['testvar']]) > 10)
	matchformula <- as.formula(paste0('testvar ~ ',colnames(data)%>%str_subset('total')%>%paste0(collapse='+')) )
	matchobject <- MatchIt::matchit(matchformula,data=data)
	data[[1]][matchobject[[1]]%>%as.numeric]
}


matchsettbl <- exprdata%>%filter(assay=='total')%>%select(gene_name,dataset,signal)%>%spread(dataset,signal)
matchsettbl[[1]] <- unlist(symid[matchsettbl[[1]]])


#First try with MatchIT to make conotrol set

matchingnoncoreset <- getmatchingset(symid[t_timecoretbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[t_timecoretbl$gene_name])
fisher.test(table(testentrez %in% symid[t_timecoretbl$gene_name],testentrez %in% symid[tegenes]))

#Also just matching mean
noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[t_timecoretbl$gene_name])
cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[t_timecoretbl$gene_name])
noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
	size=nrow(cormatchsettbl)*3,
	inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
	empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
)
testentrez <- c(symid[t_timecoretbl$gene_name],noncorentrez)
fisher.test(table(testentrez %in% symid[t_timecoretbl$gene_name],testentrez %in% symid[tegenes]))


#Now same thing, but with the diff cor set
##First try with MatchIT to make conotrol set
matchingnoncoreset <- getmatchingset(symid[t_diffcoretbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[t_diffcoretbl$gene_name])
fisher.test(table(testentrez %in% symid[t_diffcoretbl$gene_name],testentrez %in% symid[tegenes]))

##Also just matching mean
noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[t_diffcoretbl$gene_name])
cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[t_diffcoretbl$gene_name])
noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
	size=nrow(cormatchsettbl)*3,
	inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
	empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
)
testentrez <- c(symid[t_diffcoretbl$gene_name],noncorentrez)
fisher.test(table(testentrez %in% symid[t_diffcoretbl$gene_name],testentrez %in% symid[tegenes]))

#Now same thing but with both core sets

##Also just matching mean
coretbl<-t_diffcoretbl%>%bind_rows(t_timecoretbl)

matchingnoncoreset <- getmatchingset(symid[coretbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[coretbl$gene_name])
fisher.test(table(testentrez %in% symid[coretbl$gene_name],testentrez %in% symid[tegenes]))

noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[coretbl$gene_name])
cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[coretbl$gene_name])
noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
	size=nrow(cormatchsettbl)*3,
	inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
	empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
)
testentrez <- c(symid[coretbl$gene_name],noncorentrez)
fisher.test(table(testentrez %in% symid[coretbl$gene_name],testentrez %in% symid[tegenes]))

#Now everything in the waves - actually this kind of gives a result!
matchingnoncoreset <- getmatchingset(symid[t_wavetbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[t_wavetbl$gene_name])
fisher.test(table(testentrez %in% symid[t_wavetbl$gene_name],testentrez %in% symid[tegenes]))

noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[t_wavetbl$gene_name])
cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[t_wavetbl$gene_name])
noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
	size=nrow(cormatchsettbl),
	inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
	empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
)
testentrez <- c(symid[t_wavetbl$gene_name],noncorentrez)
fisher.test(table(testentrez %in% symid[t_wavetbl$gene_name],testentrez %in% symid[tegenes]))

#and everything in the dynamics table
matchingnoncoreset <- getmatchingset(symid[t_dynamicstbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[t_dynamicstbl$gene_name])
fisher.test(table(testentrez %in% symid[t_dynamicstbl$gene_name],testentrez %in% symid[tegenes]))

noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[t_dynamicstbl$gene_name])
cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[t_dynamicstbl$gene_name])
noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
	size=nrow(cormatchsettbl),
	inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
	empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
)
testentrez <- c(symid[t_dynamicstbl$gene_name],noncorentrez)
fisher.test(table(testentrez %in% symid[t_dynamicstbl$gene_name],testentrez %in% symid[tegenes]))


#and finally, both
t_alltbl <- t_dynamicstbl%>%bind_rows(t_wavetbl)
matchingnoncoreset <- getmatchingset(symid[t_alltbl$gene_name],matchsettbl%>%filter(!is.na(gene_name)))
testentrez <- c(matchingnoncoreset,symid[t_alltbl$gene_name])
fisher.test(table(testentrez %in% symid[t_alltbl$gene_name],testentrez %in% symid[tegenes]))


chronotypic_clusters$entrez<-symid[chronotypic_clusters$gene_name]
#Finally, look at clusters
chronotypic_clusters$is_te<- chronotypic_clusters$entrez %in% symid[tegenes]



chronotypic_clusters$cluster%<>%as.factor
chronotypic_clusters$is_te
chronotypic_clusters

medexprs<-matchsettbl%>%select(-1)%>%apply(1,function(x) median(2^x))%>%setNames(.,matchsettbl$gene_name)
#get data with mediane exprs and clusters
expr_clust_data <- chronotypic_clusters%>%left_join(enframe(medexprs,name='entrez',value='med_expr'))


expr_clust_data%>%group_by(gene_name)%>%filter(n()==48)


expr_clust_data%>%filter(med_expr>quantile(med_expr,0.9,na.rm=T))%>%.$is_te%>%table
expr_clust_data%>%filter(med_expr<quantile(med_expr,0.1,na.rm=T))%>%.$is_te%>%table
expr_clust_data%>%filter(between(med_expr,quantile(na.rm=T,med_expr,0.3),quantile(na.rm=T,med_expr,0.6)))%>%.$is_te%>%table


cluster_is_te_fits<-glm(data=expr_clust_data,is_te~ as.factor(cluster) + med_expr)
cluster_is_te_fits%>%summary

clustpredictions<-predict(expr_clust_data,expr_clust_data,se.fit=TRUE)
with(clustpredictions,data_frame(cluster=clustpredictions$fit%>%names,est=fit,hl=fit+(1.96*se.fit),ll=fit-(1.96*se.fit)))

base_is_te_fits<-glm(data=chronotypic_clusters,is_te~1)
basepredictions<-predict(base_is_te_fits,chronotypic_clusters%>%distinct(cluster),se.fit=TRUE)
with(basepredictions,data_frame(cluster=basepredictions$fit%>%names,est=fit,hl=fit+(1.96*se.fit),ll=fit-(1.96*se.fit)))







test_geneset_ctl_expr<-function(gnames2test,testname,tegenes=tegenes,matchsettbl=matchsettbl,symid=symid,plotexprhist=FALSE){
	noncormatchsettbl <- matchsettbl%>%filter(!gene_name %in% symid[gnames2test])
	cormatchsettbl <- matchsettbl%>%filter(gene_name %in% symid[gnames2test])
	noncorentrez<-sample_using_empdens(sampv=noncormatchsettbl$gene_name,
		size=nrow(cormatchsettbl),
		inpvals=noncormatchsettbl%>%select(-gene_name)%>%apply(1,median),
		empvals=cormatchsettbl%>%select(-gene_name)%>%apply(1,median)
	)
	testentrez <- c(symid[gnames2test],noncorentrez)

	med_exprs <- matchsettbl%>%select(-gene_name)%>%apply(1,median)%>%setNames(matchsettbl$gene_name)

	if(plotexprhist==TRUE){
		exprhist <- enframe(	med_exprs)%>%mutate(set=case_when(name%in%noncormatchsettbl[[1]] ~ 'Nontest',name%in%cormatchsettbl[[1]] ~ 'Test',TRUE ~ 'NA' ))%>%
			qplot(data=.,geom='histogram',fill=set,x=value)+facet_grid(set~.)+ggtitle('Expression Distributions For Test and Matching non-test sets')
	}

	ftab <- table(test=testentrez %in% symid[gnames2test],TE_change=testentrez %in% symid[tegenes])

	ftest <- fisher.test(ftab)

	xlabs <- c('Ctrl, Static TE','Ctrl, Dynamic TE',str_interp('${testname}, Static TE'),str_interp('${testname}, Dynamic TE'))

	ggplot(as.data.frame(ftab), aes(x=test:TE_change, y = Freq, fill=TE_change)) + 
    	geom_bar(stat="identity",position='dodge',width=I(0.8))+
    	scale_x_discrete(labels=xlabs)+
    	theme_bw()+
    	theme(axis.text.x=element_text(angle=45,vjust= 0.5))+
    	ggtitle(
    		paste0('Dynamic TE enrichment: ',testname),
    		sub=paste0('Fishers Exact Test OR = ',round(ftest$conf.int[1],2),' - ',round(ftest$conf.int[2],2),' p = ',round(ftest$p.value,3)
    		)
    	)

}

test_geneset_ctl_expr(c(t_dynamicstbl$gene_name,t_wavetbl$gene_name),testname='All_Wave_Genes',tegenes,matchsettbl=matchsettbl,symid=symid)


fisher.test(table(exprdata$gene_name %in% t_diffcoretbl$gene_name,exprdata$gene_name %in% tegenes))



t_wave_expr<-exprdata%>%
	select(gene_name,time,assay,signal)%>%
	group_by(gene_name,assay)%>%mutate(signal = signal/median(signal))%>%
	left_join(t_wavetbl)%>%
	left_join(t_dynamicstbl)%>%
	left_join(t_timecoretbl)%>%
	left_join(t_diffcoretbl,by='gene_name')

t_wave_expr$`Diff_core_Specificity`%>%table


t_wave_expr%<>%mutate_at(vars(matches('wave')),list(~replace_na(.,0)))

plotwaves <- function(data,wavecol) {
	wavecol<-enquo(wavecol)
	# browser()
	tlabs <- unique(data$time)
	data$time%<>%as.factor%>%as.numeric
	data$assay%<>%factor
	#
	newassaynames<- c('RNASeq','Riboseq','Mass Spec')
	assaylvls<-c('total','ribo','MS')%>%setNames(newassaynames)
	data$assay %<>% fct_recode(!!!assaylvls)%>%fct_relevel(newassaynames)
	# data$assay %<>% factor
	browser()
	datagsum<-data%>%group_by(!!wavecol,assay,time,gene_name)%>%summarise(signal=median(signal))
	datasum<-data%>%group_by(!!wavecol,assay,time,gene_name)%>%summarise(signal=median(signal))%>%group_by(!!wavecol,assay,time)%>%summarise(y=median(signal),ymax=quantile(signal,0.75),ymin=quantile(signal,0.25))
	wavecolname<-str_replace_all(str_replace_all(rlang::quo_text(wavecol),'`',''),'`','')
	data[[wavecolname]]=as.factor(data[[wavecolname]])
	#

	data %>%subset(gene_name %in% unique(gene_name))%>%{
		# ggplot(data=.,aes(color=!!wavecol,y=signal,x=time))+
		ggplot(data=.,aes(y=signal,color=as.factor(!!wavecol),x=time))+
    	geom_smooth(method='glm')+
    	# geom_ribbon(alpha=I(0.5),data=datasum,aes(fill=!!wavecol,y=y,ymax=ymax,ymin=ymin))+
    	# geom_smooth() +
    	geom_point(alpha=0.1)+
    	# geom_line(alpha=I(.1),data=datagsum%>%subset(gene_name==gene_name),aes(group=gene_name))+
   		# scale_x_continuous(name='Time Point',breaks=seq_along(tlabs),labels=tlabs)+
   		scale_y_continuous(name='Read Density / LFQ , Scaled by Gene/Assay median')+
    	facet_grid(~assay)+
    	theme_bw()+
    	scale_color_discrete(name=wavecolname)
	}

}
0
t_wave_expr%>%plotwaves(`Time_core_Specificity`)

t_wave_expr%>%plotwaves(`E12 wave`)

t_wave_expr%>%plotwaves(`E13 wave`)

t_wave_expr%>%plotwaves(`E14 wave`)

t_wave_expr%>%plotwaves(`E15 wave`)

t_wave_expr%>%plotwaves(`Dynamics type`)

t_wave_expr%>%plotwaves(`Diff_core_Specificity`)

#read in the core gene models

##read in the diff model

##read in the age model


#read in the 'wave' data.


#maybe come up with some distinction between gene regulators and non-gene regulators
#GO cats will do.


