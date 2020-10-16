base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
rename<-dplyr::rename
first<-dplyr::first
last<-dplyr::last

sel_prodalfcs<-readRDS('data/sel_prodalfcs.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')

prediction_df = bind_rows(
  sel_prodalfcs%>%
    mutate(assay='MS')%>%
    rename('logFC'=diff)%>%
    select(gene_id,time,assay,logFC,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,logFC,CI.L,CI.R)
)

exprdf = bind_rows( 
  allvoom$E%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate')),
  sel_ms_mat[!is.na(rownames(sel_ms_mat)),]%>%as.data.frame%>%
    rownames_to_column('gene_id')%>%
    filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))
)

exprdf$gene_name = gid2gnm[[exprdf$gene_id]]
prediction_df$gene_name = gid2gnm[[prediction_df$gene_id]]

################################################################################
########
################################################################################
  


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
      out = x%>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x))
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
      geom_ribbon(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(x=as.numeric(as_factor(time)),y=logFC,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+ 
      geom_line(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),linetype=2,aes(x=as.numeric(as_factor(time)),y=logFC))+  
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

