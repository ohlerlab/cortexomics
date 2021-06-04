# install.packages(c('magrittr','stringr','ggpubr','data.table','assertthat','tidyverse','dplyr','here','conflicted','ggplot2'))
#install.packages('tidyverse')

suppressMessages(library(BiocManager))
options(repos = BiocManager::repositories())
suppressMessages(library(shiny))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(here))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(Gviz))
suppressMessages(library(magrittr))
library(ggplot2)
library(rsconnect)

#Note this function contains the expression data, and predictions needed
#to plot the trajectory plots - so it's very much a closure, rather than
#a function
onserver = (getwd()%>%str_detect('/srv/connect/apps'))
datafile='data/make_trajplots_arglist.rds'


make_trajplots_arglist <- readRDS(datafile)


tpcols <- c(E13='#214098',E145='#2AA9DF',E16='#F17E22',E175='#D14E28',P0='#ED3124')

make_trajplots<-function(gnames,make_trajplot,...){
  ngenes <- length(gnames)
  ggarrange(plotlist=map(gnames,make_trajplot,...),nrow=ngenes)
}
make_trajplots <- do.call(what=partial,args=c(list(make_trajplots),make_trajplots_arglist))

get_gviztracks <- function(igene,shinytrackfolder='data/Shiny_track_data'){
  tracks<-readRDS(file.path(shinytrackfolder,paste0(igene,'.rds')))
  tracks
}

if(getwd()%>%str_detect('AG_Ohler')) { 
  #now plot
  plotfile<- here(paste0('plots/','tmp','.pdf'))
  pdf(plotfile)
  print(make_trajplots(c('Satb2','Flna'),myymin=-1,myymax=4))
  dev.off()
  plotfile<- here(paste0('plots/','tmplocus','.pdf'))
  pdf(plotfile)
  gviztracks <- get_gviztracks('Satb2')
  print(Gviz::plotTracks(gviztracks,col = tpcols))
  dev.off()
  message(normalizePath(plotfile))
  stop()
  if(TRUE){
  #https://www.shinyapps.io/admin/#/dashboard - see for auth
  appfolder <- here('src/R/Shiny/')
  # shiny::runApp(appfolder)
  rsconnect::deployApp(appfolder,
      appFiles=c('App.R',basename(datafile),'data'),
      appTitle='Cortexomics Trajectories',
      appName='cortex_traj')
  }
}

################################################################################
########The app
################################################################################

genechoices <- c('Satb2','Flna')
# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Mouse Neocortex Developmental Translatome"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

  # selectInput('genes2plot', 'Gene For Plotting', genechoices, selected = 'Satb2', multiple = FALSE,
  #     selectize = TRUE, width = NULL, size = NULL)
    textInput('genes2plot', 'Gene For Plotting', value = "Flna", width = NULL,placeholder=NULL),
    sliderInput("ExprYmax", "Ymax",
                       min = -10, max = 10, value = 1),
    sliderInput("ExprYmin", "Ymin",
                       min = -10, max = 10, value = -4),
    imageOutput("sideimage")
    # img('~/Desktop/timage1.svg')
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      #
      # uiOutput("panelimages"),
  
      # Output: Histogram ----
      # plotOutput(outputId = "trajplot",height = plotHeight)
      # plotOutput(outputId = "trajplot",height = '1000px')
      uiOutput("ui")
    )
  )
)



splitgenes <- .%>%str_split(' ')%>%unlist%>%setdiff(c('',' '))
plotable_genes <- make_trajplots_arglist$exprdf$gene_name%>%unique
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
              # alternativelist <- get_alternativelist()
              notfoundtext <- str_interp('Genes ${nonfoundgenecol} Not Found. By ${nonfoundgenes[1]} Did you mean...\n${alternativelist}\n (Case sensitive matching, see table S1 for complete list)')

              print(qplot(1,1,label=notfoundtext,geom='text',size=I(10))+theme_bw())
              # renderText(notfoundtext)

        }else{
          print(suppressWarnings({make_trajplots(genes2plot,myymin=input$ExprYmin,myymax=input$ExprYmax)}))                   
        }
    })
    output$locusplot <- renderPlot({
        genes2plot <- input$genes2plot
        genes2plot %<>% splitgenes
        if(length(genes2plot)>1){
          print(qplot(1,1,label='enter 1 gene name to display compressed locus plot',geom='text',size=I(10))+theme_bw())
        }else{
            gviztracks <- possibly(get_gviztracks,NULL)(genes2plot[1])
            if(is.null(gviztracks)){
                print(qplot(1,1,label=str_interp('No locusplot data saved for ${genes2plot[1]}'),geom='text',size=I(10))+theme_bw())
            }else{
              gviztracks <- head(gviztracks,5)
              gviztracks[[3]]@name <- 'log2(iBAQ)'
              print(Gviz::plotTracks(gviztracks,col = tpcols,col.title='black',main=str_interp("${genes2plot[1]} plotted 5'->3' ")))
            }
        }
    })
    # output$panelimages <- renderUI(
    #     renderImage(list(src = "4_web-02.svg", width = '25%', height = 200)),
    #     renderImage(list(src = "4_web-03.svg", width = '25%', height = 200)),
    #     renderImage(list(src = "4_web-04.svg", width = '25%', height = 200)),
    #     renderImage(list(src = "4_web-05.svg", width = '25%', height = 200))
    # )
    message(names(input))

    message('......')
    plotn=length(splitgenes(input$genes2plot))
    message(plotn)
    if(! all(splitgenes(input$genes2plot) %in% plotable_genes)) plotn=2
    message(splitgenes(input$genes2plot))
    tagList(
      {if(length(splitgenes(input$genes2plot))>1){
        "Enter 1 gene name to show locus plot"
      }else{
        plotOutput("locusplot", width = '100%', height = 200*3)
      }}
      ,
      fluidRow(plotOutput("trajplot", width = '100%', height = 200*plotn)), ## nn is my scaling factor
      fluidRow(
        column(3,imageOutput("panelimages")),
        column(3,imageOutput("panelimages3")),
        column(3,imageOutput("panelimages4")),
        column(3,imageOutput("panelimages2"))
      ),
    )

  })
  output$panelimages <- renderImage({ list(src ='data/4_web2-02.svg',contentType = 'image/svg+xml', width='100%',height=160, alt = paste(""))},deleteFile=FALSE)
  output$panelimages2 <- renderImage({ list(src ='data/4_web2-03.svg',contentType = 'image/svg+xml',  width='100%',height=160,alt = paste(""))},deleteFile=FALSE)
  output$panelimages3 <- renderImage({ list(src ='data/4_web2-04.svg',contentType = 'image/svg+xml',  width='100%',height=160,alt = paste(""))},deleteFile=FALSE)
  output$panelimages4 <- renderImage({ list(src ='data/4_web2-05.svg',contentType = 'image/svg+xml',  width='100%',height=160,alt = paste(""))},deleteFile=FALSE)

  output$sideimage <- renderImage(list(src='data/4_web2-01.svg',contentType = 'image/svg+xml',width='100%',height='100%'),deleteFile=FALSE)
  message('done')

}


shinyApp(ui, server)
