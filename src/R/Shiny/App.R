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
library(ggplot2)
library(rsconnect)

#Note this function contains the expression data, and predictions needed
#to plot the trajectory plots - so it's very much a closure, rather than
#a function
datafile <- 'make_trajplots_arglist.rds'

if(!shiny::isRunning()) { datafile = here('data',datafile)}

make_trajplots_arglist <- readRDS(datafile)

make_trajplots<-function(gnames,make_trajplot,...){
  ngenes <- length(gnames)
  ggarrange(plotlist=map(gnames,make_trajplot,...),nrow=ngenes)
}
make_trajplots <- do.call(what=partial,args=c(list(make_trajplots),make_trajplots_arglist))

if(!shiny::isRunning()) { 
  #now plot
  plotfile<- here(paste0('plots/','tmp','.pdf'))
  pdf(plotfile)
  print(make_trajplots(c('Satb2','Flna'),myymin=-1,myymax=4))
  dev.off()
  message(normalizePath(plotfile))
  if(TRUE){
  #https://www.shinyapps.io/admin/#/dashboard - see for auth
  appfolder <- here('src/R/Shiny/')
  file.copy(datafile,appfolder)
  shiny::runApp(appfolder)
  rsconnect::deployApp(appfolder,
      appFiles=c('App.R',datafile),
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



splitgenes <- .%>%str_split(' ')%>%unlist%>%setdiff(c('',' '))
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
