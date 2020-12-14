base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
rename<-dplyr::rename
first<-dplyr::first
last<-dplyr::last

  gid2gnm<-load_hashmap('gid2gnm.hmp')
  gnm2gid<-load_hashmap('gnm2gid.hmp')
  gnm2gid<-load_hashmap('gnm2gid.hmp')


sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')
ms_metadf<-readRDS('data/ms_metadf.rds')
allvoom <- readRDS(here('data/allvoom.rds'))
prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)

exprdf = bind_rows( 
  allvoom$E%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate')),
  sel_ms_mat[,]%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))
)
exprdf%>%head
exprdf$gene_name = gid2gnm[[exprdf$gene_id]]
prediction_df$gene_name = gid2gnm[[prediction_df$gene_id]]


#' Thi sis a title with inline R code `r foo`

#' First we load the list of protein IDs, handpicked by Matt, using only the small
#' or large subunits - no mitochondrial riboproteins  
#define ambigous protein groups as those which have elements that appear in more than one protein group
# allpgroups <- mstall$Protein_IDs%>%unique
# multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
# all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)
# library(data.table)


rpexprdf = exprdf%>%inner_join(ms_metadf%>%filter(is_rpl|is_rps))
rpexprdf$cat = case_when(
  rpexprdf$gene_id %in% (ms_metadf%>%filter(is_rpl)%>%.$gene_id) ~ 'RPL',
  rpexprdf$gene_id %in% (ms_metadf%>%filter(is_rps)%>%.$gene_id) ~ 'RPS',
  TRUE ~ 'other'
)
rpexprdf$cat%>%table
rpexprdf%<>%bind_rows(rpexprdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (rpexprdf%>%filter(assay=='total')%>%.$signal)))

ggdf = rpexprdf%>%
  arrange(time,assay=='MS',assay=='TE',assay=='ribo')%>%
  mutate(assay=as_factor(assay))%>%
  group_by(gene_name,assay,time,cat)%>%
  summarise(signal=mean(signal,na.rm=TRUE))%>%
    group_by(gene_name,assay,cat)%>%
    mutate(signal = signal - median(signal[time=='E13']))

medggdf = ggdf%>%group_by(assay,cat,time)%>%summarise(signal = median(signal))

#now plot
plotfile<- here('plots/figures/figure1/RP_traj.pdf')
pdf(plotfile,w=14,h=7)
ggdf%>%
    # filter(signal < -3)
  ggplot(.,aes(y=signal,x=as.numeric(as.factor(time)),group=gene_name))+
  geom_line(alpha=I(.8),color=I("grey"))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('Time'),labels = rpexprdf$time%>%unique)+
  scale_y_continuous(paste0('log2(CPM/iBAQ) relative to Mean'),limits=c(-1,1))+
  ggtitle(paste0('Ribosomal Proteins - Expression Trajectory'))+
  facet_grid(cat~assay,scales='free')+
  geom_line(data=medggdf,aes(group=cat),color=I('black'))+
  theme_bw()
dev.off()
normalizePath(plotfile)

stop()

################################################################################
########
################################################################################
  
prediction_df%>%slice(1)%>%left_join(exprdf)

#this contains the 

#' Need to think of a way to show if proDA is working well.
#' and ADDITIONALLY, if our dropout inclusive model works better...
#' an example of that would be great.
#' Also reports that work would be great....
assaynames = c('total'='RNA-seq','ribo'='Ribo-seq','TE'='TE','MS'='Mass-Spec')
stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
tpnames = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))

{
make_traj_plot <- function(genes2plot,myymin,myymax,exprdf,prediction_df){
  plotlistlist = list()
  #
  assert_that(all(genes2plot %in% exprdf$gene_name))
  #get data for that test gene
  ggdf <- exprdf%>%filter(gene_name%in%genes2plot)
  #add TE to our data frame
  ggdf%<>%bind_rows(ggdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (ggdf%>%filter(assay=='total')%>%.$signal)))

  #scaling function
  scaledata = function(x,assay2plot){
    if(assay2plot=='MS'){#use msrescale2lfq
      out = x%>%mutate_at(vars(one_of(c('signal','estimate','CI.R','CI.L'))),list(function(x)x))
    }else{
      out = x 
    }
    rescale = if('signal' %in% colnames(x)){ median(x$signal[x$time=='E13']) } else {median(x$estimate[x$time=='E13'])}
    out = x %>%mutate_at(vars(one_of(c('signal','estimate','CI.R','CI.L'))),list(function(x)x-rescale))
    out
  }
  # plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot='MS')
  # scaledata = identity
  #
  trajfile = './plots/tmp.pdf'
  for(gname in genes2plot){
    assaytitles <- if(gname==first(genes2plot)) assaynames else NULL
    xname <- if(gname==last(genes2plot)) 'Stage' else NULL
    xtpnames <- if(gname==last(genes2plot)) tpnames else NULL
    uidfilt<-.%>%filter(gene_name==gname)
    # ms_id2protein_id%>%filter(uprotein_id==testuid)
    plotdf<-ggdf%>%uidfilt
    assays2plot <- unique(ggdf$assay)%>%sort%>%rev%>%.[order(.=='TE')]
    trajectoryplots<-lapply(assays2plot,function(assay2plot){
      if(assay2plot==assays2plot[1]){
        yname = str_interp('Fold Change ( ${gname} )')
      }else{
        yname =''
      }
      if(assay2plot!='MS'){
        points = geom_point()
        linerange=NULL
      }else{
        # points = geom_point(data%>%filter%>%scaledata('dont unscale'))
        points = geom_point()
        # linerange=geom_linerange(color=I('blue'),data=ggdf_msconf%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(y=signal,ymin=CI.L,ymax=CI.R))
      }
      browser()
      ggplot(
      data = plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot),
      aes(
        x=as.numeric(as_factor(time)),
        y=signal
      ))+
      points+
      # linerange+
      geom_ribbon(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(x=as.numeric(as_factor(time)),y=estimate,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+ 
      geom_line(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),linetype=2,aes(x=as.numeric(as_factor(time)),y=estimate))+  
      # geom_line(data=prediction_df%>%uidfilt,aes(x=time,y=estimate))+  
      scale_x_continuous(name=xname,labels=xtpnames)+
      theme_bw()+
      scale_y_continuous(name=yname)+
      ggtitle(assaytitles[assay2plot])
      # facet_wrap( ~ assay,scales='free')+
    })
    # trajectoryplots<-commonyspanplots(trajectoryplots,breakint=0.25,minrange=8)
    # trajectoryplots[c(1:3)]%<>%lapply(function(plot) plot+ coord_cartesian(ylim=u[[gene2plot]]))
    # trajectoryplots[c(4)]%<>%lapply(function(plot) plot+ coord_cartesian(ylim=c(-3,3)))
    trajectoryplots%<>%lapply(function(plot) plot+ coord_cartesian(ylim=c(myymin,myymax)))
    #
    plotlistlist = append(plotlistlist,list(trajectoryplots))
    # genenamelist = append(genenamelist,gname)
  }
  #
  assert_that(length(plotlistlist) < 9,msg='Too many genes!')
  trajectoryplot<-ggarrange(plotlist=plotlistlist%>%flatten,ncol=4,nrow=length(plotlistlist))
  trajectoryplot<-annotate_figure(trajectoryplot)
  trajectoryplot
}
}




#What should our plot do???


################################################################################
########The app
################################################################################

genechoices <- c('Satb2','Flna')


plotable_genes = exprdf$gene_name%>%unique
# conflict_prefer("Position",'ggplot2')


make_traj_plot_p <- partial(make_traj_plot,
  exprdf=exprdf,
  prediction_df=prediction_df,
)
# make_traj_plot_p('Flna')


if(!shiny::isRunning()) {
  pdf('tmp.pdf')
  print(make_traj_plot_p(c('Satb2','Flna'),myymin=-5,myymax=5))
  dev.off()
  message(normalizePath('tmp.pdf'))
}


stop()
splitgenes <- .%>%str_split(' ')%>%unlist%>%setdiff(c('',' '))

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Behold Cortexomics!"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

  # selectInput('genes2plot', 'Gene For Plotting', genechoices, selected = 'Satb2', multiple = FALSE,
  #     selectize = TRUE, width = NULL, size = NULL)
    textInput('genes2plot', 'Gene For Plotting', value = "Satb2 Flna", width = NULL,placeholder=NULL),
    sliderInput("ExprYmax", "Ymax",
                       min = -10, max = 10, value = 4),
    sliderInput("ExprYmin", "Ymin",
                       min = -10, max = 10, value = -4)
  
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      # plotOutput(outputId = "trajplot",height = plotHeight)
      # plotOutput(outputId = "trajplot",height = '1000px')
      uiOutput("ui")
    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
 
  output$ui <- renderUI({
 
    output$trajplot <- renderPlot({
      
        genes2plot <- input$genes2plot
        genes2plot %<>% splitgenes

        if(! all(genes2plot %in% plotable_genes)){
              nonfoundgenes <- genes2plot%>%setdiff(plotable_genes)
              nonfoundgenes <- paste(sep=',',nonfoundgenes)
              nonfoundgenecol <- paste0(nonfoundgenes,collapse='')
              othergenes <- sample(unique(plotable_genes))
              alternativelist <- othergenes[head(order((adist(toupper(nonfoundgenes[1]),toupper(othergenes)))),n=5)]
              alternativelist <- paste0(paste0(alternativelist,' ?\n'),collapse='')
              notfoundtext <- str_interp('Genes ${nonfoundgenecol} Not Found. By ${nonfoundgenes[1]} Did you mean...\n${alternativelist}\n (Case sensitive matching)')

              print(qplot(1,1,label=notfoundtext,geom='text',size=I(10))+theme_bw())
              # renderText(notfoundtext)

        }else{
          print(suppressWarnings({make_traj_plot_p(genes2plot,myymin=input$ExprYmin,myymax=input$ExprYmax)}))                   
        }
      # browser()
    })
    message(names(input))
    message('......')
    plotn=length(splitgenes(input$genes2plot))
    message(plotn)
    if(! all(splitgenes(input$genes2plot) %in% plotable_genes)) plotn=2
    message(splitgenes(input$genes2plot))
    tagList(
      # plotOutput("trajplot", width = '100%', height = 300*plotn) ,## nn is my scaling factor,
      plotOutput("trajplot", width = '100%', height = 200*plotn) ## nn is my scaling factor
    )

  })
  message('done')

}

shinyApp(ui, server)

#https://www.shinyapps.io/admin/#/dashboard - see for auth
# appfolder <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/Shiny/'
# shiny::runApp(appfolder)
# rsconnect::deployApp(appfolder,appFiles=c('app.R','trajplotobjects.Rdata'),appTitle='Cortexomics Trajectories',appName='cortex_traj')

