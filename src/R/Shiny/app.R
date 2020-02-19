# install.packages(c('magrittr','stringr','ggpubr','data.table','assertthat','tidyverse','dplyr','here','conflicted','ggplot2'))
#install.packages('tidyverse')

suppressMessages(library(magrittr))
suppressMessages(library(shiny))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(here))
library(conflicted)
library(ggplot2)
library(rsconnect)
# conflict_prefer("filter",'dplyr')
# conflict_prefer("last",'dplyr')
# conflict_prefer("setdiff",'dplyr')

if(file.exists('trajplotobjects.Rdata')){
  datafile <- 'trajplotobjects.Rdata'
}else{
  # datafile <- here('data/trajplotobjects.Rdata')
    datafile <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics/data/trajplotobjects.Rdata'

}
if(!exists('exprdf')){
   load(file=datafile)
}

#todo stick gene name on this

#assaynames - prettier names for the 


#currently msrescale2lfq is a problem

#The scale we want is the median LFQ one I guess. So the median normalized protmat, but before it was centered on zero
#for numerical issues in stan
#That means we need to rescale the output of our plots to that.
#We'll do that by simply re-aligning the matrices, another function should do the same for confidence interval dfs.
#
#this contains all chosen ms_ids , median normalized
ms_mat_mednorm <- 
#this is the same but for count data
countmat_mednorm <- 

#This contains the limma predictions, rescaled

#This contains the predictions from the proDA run.

#This contains arbitrary predictions, from a variable number of models. The
# function will add them in as different panels.

################################################################################
########
################################################################################
  


#this contains the 

#' Need to think of a way to show if proDA is working well.
#' and ADDITIONALLY, if our dropout inclusive model works better...
#' an example of that would be great.
#' Also reports that work would be great....

{
make_traj_plot <- function(genes2plot,myymin,myymax,msdf,postmeanmat,postprecsdmat,exprdf,msrescale2lfq,best_uprotein_ids,ci_df,tpnames,assaynames){
  plotlistlist = list()
  # genenamelist = list()
  #
  u = list('Bcl11b'=c(-2,4),'Flna'=c(-4,1),'Nes'=c(-4,1),'Satb2'=c(-1,5))
  #
  assert_that(all(genes2plot %in% exprdf$gene_name))
  testpids <- exprdf%>%subset(gene_name%in%genes2plot)%>%distinct(protein_id,gene_name)%>%{setNames(.$protein_id,.$gene_name)}
  test_uids<-unique(exprdf$uprotein_id)%>%str_subset(testpids%>%paste0(collapse='|'))
  #get data for that test gene
  ggdf <- exprdf%>%filter(uprotein_id%in%test_uids)
  ggdf%<>%bind_rows(ggdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (ggdf%>%filter(assay=='total')%>%.$signal)))
  postmeanmat_scaled <- postmeanmat - msrescale2lfq
  ggdf_msconf <- 
    left_join(
      ((postmeanmat_scaled[test_uids,,drop=F])-(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.L,-uprotein_id)%>%separate(dataset,c('time','assay')),
      ((postmeanmat_scaled[test_uids,,drop=F])+(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.R,-uprotein_id)%>%separate(dataset,c('time','assay'))
    )%>%
    left_join(
      ((postmeanmat_scaled[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%separate(dataset,c('time','assay')),
      )
  # ggdf_msconf$time%<>%factor%>%time
  #get ms-protein id pairs
  test_uids<-ggdf%>%distinct(uprotein_id,gene_name)%>%{setNames(.$uprotein_id,.$gene_name)}
  #scaling function
  scaledata = function(x,assay2plot){
    if(assay2plot=='MS'){#use msrescale2lfq
      out = x%>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x+msrescale2lfq))
    }else{
      out = x 
    }
    rescale = if('signal' %in% colnames(x)){ median(x$signal[x$time=='E13']) } else {median(x$logFC[x$time=='E13'])}
    out = x %>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x-rescale))
    out
  }
  # plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot='MS')
  # scaledata = identity
  #
  trajfile = './plots/tmp.pdf'
  for(testuid in test_uids[genes2plot]){
    assaytitles <- if(testuid==test_uids[genes2plot][1]) assaynames else NULL
    xname <- if(testuid==last(test_uids[genes2plot])) 'Stage' else NULL
    xtpnames <- if(testuid==last(test_uids[genes2plot])) tpnames else NULL
    gname = names(test_uids)[test_uids==testuid]
    uidfilt<-.%>%filter(uprotein_id==testuid)
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
        points = geom_point(data=msdf%>%filter(uprotein_id%in%testuid)%>%scaledata('dont unscale'))
        linerange=geom_linerange(color=I('blue'),data=ggdf_msconf%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(y=signal,ymin=CI.L,ymax=CI.R))
      }
      
      browser()
      ggplot(
      data = plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot),
      aes(
        x=as.numeric(as_factor(time)),
        y=signal
      ))+
      points+linerange+
      geom_ribbon(data=ci_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(x=as.numeric(as_factor(time)),y=logFC,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+ 
      geom_line(data=ci_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),linetype=2,aes(x=as.numeric(as_factor(time)),y=logFC))+  
      # geom_line(data=prediction_df%>%uidfilt,aes(x=time,y=logFC))+  
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
  msdf=msdf,
  postmeanmat=postmeanmat,
  postprecsdmat=postprecsdmat,
  exprdf=exprdf,
  tpnames=mytpnames,
  msrescale2lfq=msrescale2lfq,
  best_uprotein_ids=best_uprotein_ids,
  ci_df=prediction_df,
  assaynames=myassaynames,
)
make_traj_plot_p('Flna')

stop()

if(!shiny::isRunning()) {
  pdf('tmp.pdf')
  print(make_traj_plot_p(c('Satb2','Flna'),myymin=-3,myymax=3))
  dev.off()
}


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

