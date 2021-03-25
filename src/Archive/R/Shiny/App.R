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
library(ggplot2)
library(rsconnect)

#Note this function contains the expression data, and predictions needed
#to plot the trajectory plots - so it's very much a closure, rather than
#a function
onserver = (getwd()%>%str_detect('/srv/connect/apps'))
datafile <- 'make_trajplots_arglist.rds'

if(!onserver) {
  datafile = here('data',datafile)
}


if(onserver){
  message(datafile)
  message(getwd())
  message(paste0(list.files(full=T,rec=T),sep='\n'))
}

make_trajplots_arglist <- readRDS(datafile)


tpcols <- c(E13='#214098',E145='#2AA9DF',E16='#F17E22',E175='#D14E28',P0='#ED3124')

make_trajplots<-function(gnames,make_trajplot,...){
  ngenes <- length(gnames)
  ggarrange(plotlist=map(gnames,make_trajplot,...),nrow=ngenes)
}
make_trajplots <- do.call(what=partial,args=c(list(make_trajplots),make_trajplots_arglist))

get_gviztracks <- function(igene,shinytrackfolder='data/Shiny_track_data'){
  readRDS(file.path(shinytrackfolder,paste0(igene,'.rds')))
}

if(!onserver) { 
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
  titlePanel("Mouse Neocortical Translation Trajectories\nE12.5 - P0"),

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
              alternativelist <- get_alternativelist()
              notfoundtext <- str_interp('Genes ${nonfoundgenecol} Not Found. By ${nonfoundgenes[1]} Did you mean...\n${alternativelist}\n (Case sensitive matching)')

              print(qplot(1,1,label=notfoundtext,geom='text',size=I(10))+theme_bw())
              # renderText(notfoundtext)

        }else{
          print(suppressWarnings({make_trajplots(genes2plot,myymin=input$ExprYmin,myymax=input$ExprYmax)}))                   
        }
    })
    output$locusplot <- renderPlot({
      
        genes2plot <- input$genes2plot
        genes2plot %<>% splitgenes

        if(! all(genes2plot %in% plotable_genes)){
          print(qplot(1))
        }else{
          gviztracks <- get_gviztracks(genes2plot[1])
          print(Gviz::plotTracks(gviztracks,col = tpcols))
        }
    })
    message(names(input))
    message('......')
    plotn=length(splitgenes(input$genes2plot))
    message(plotn)
    if(! all(splitgenes(input$genes2plot) %in% plotable_genes)) plotn=2
    message(splitgenes(input$genes2plot))
    tagList(
      plotOutput("trajplot", width = '100%', height = 200*plotn), ## nn is my scaling factor
      plotOutput("locusplot", width = '100%', height = 200*1)
    )

  })
  message('done')

}

shinyApp(ui, server)
