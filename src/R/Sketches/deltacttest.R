### R code from vignette source 'rtPCR-usage.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=200)


###################################################
### code chunk number 2: exec ddCt inside R (eval = FALSE)
###################################################
## scFile <- system.file("scripts/ddCt.R", package="ddCt")
## inputFile <- system.file("extdata/Experiment1.txt", package="ddCt")
## source(scFile)
## ddCtExec(list(inputFile=inputFile, ...))


###################################################
### code chunk number 3: load ddct
###################################################
source(system.file("scripts", "ddCt.R", package="ddCt"))


###################################################
### code chunk number 4: specifyFiles
###################################################
params <- list(loadPath = "ext_data/tRNA_data//Raw data/",
               savePath = getwd())  
###################################################
### code chunk number 5: The Data
###################################################
params$confFile <- system.file("scripts", "ddCt.conf", package="ddCt")
# params$inputFile <- c("Experiment1.txt")
params$inputFile <- list.files(rec=T,params$loadPath,pattern='txt',full=F)%>%.[[1]]
params$sampleAnnotationFile = NULL
params$columnForGrouping    = NULL
params$onlyFromAnnotation   = FALSE

testrawdata <- list.files(rec=T,params$loadPath,pattern='txt',full=T)%>%.[[1]]
ReadqPCR::read.qPCR(testrawdata)

testrawdata <- list.files(rec=T,params$loadPath,pattern='txt',full=T)%>%.[[1]]
ReadqPCR::read.qPCR(testrawdata)
library(HTqPCR)








raw <-ctobject
pdf('tmp.pdf')
g <- featureNames(ctobject)[1:10]
plotCtOverview(ctobject, genes = g, xlim = c(0, 50), groups = mygroups,
conf.int = TRUE, ylim = c(0, 55))
plotCtCard(raw, col.range = c(10, 35), well.size = 2.6)
dev.off()
normalizePath('tmp.pdf')

plotCtReps(raw, card = 2, percent = 20)

plotCtOverview(ctobject, genes = g, xlim = c(0, 50), groups = mygroups,
calibrator = "Control")
dev.off()

readLines(file.path(path,files$File[1]))
c('flag', 'feature', 'position', 'type' ,'Ct')

###################################################
### code chunk number 6: transform
###################################################
#params$geneAlias   = c("Gene1"="Gene A",
#                       "Gene4"="Gene B",
#			"Gene5"="Gene C",
#			"Gene2"="HK1",
#			"Gene3"="HK2")
#params$sampleAlias = NULL


###################################################
### code chunk number 7: housekeeping
###################################################
params$referenceGene <- c("Gene1", "Gene2")
params$referenceSample <- c("Sample1")


###################################################
### code chunk number 8: threshold
###################################################
params$threshold <- 40


###################################################
### code chunk number 9: mode
###################################################
params$mode = "median"


###################################################
### code chunk number 10: effizienzen synta (eval = FALSE)
###################################################
## params$efficiencies       = c("Gene1"=1.9, "Gene2"=2, ...)


###################################################
### code chunk number 11: effizienzen
###################################################
#params$efficiencies       = c("Gene A"=1.9,"Gene B"=1.8,"HK1"=2,"Gene C"=2,"HK2"=2)
#params$efficienciesError  = c("Gene A"=0.01,"Gene B"=0.1,"HK1"=0.05,"Gene C"=0.01,"HK2"=0.2)


###################################################
### code chunk number 12: PlotMode
###################################################
params$plotMode = c("level","Ct")


###################################################
### code chunk number 13: remaining
###################################################
#REMAIN
params$genesRemainInOutput   = NULL
params$samplesRemainInOutput = NULL
#NOT
params$genesNotInOutput   = NULL
params$samplesNotInOutput = NULL


###################################################
### code chunk number 14: grouping
###################################################
params$groupingBySamples = FALSE


###################################################
### code chunk number 15: plot0
###################################################
params$plotPerObject = TRUE


###################################################
### code chunk number 16: plot1
###################################################
#params$groupingForPlot = NULL
#params$groupingColor   = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")


###################################################
### code chunk number 17: plot2
###################################################
#params$cutoff = NULL


###################################################
### code chunk number 18: plot3
###################################################
#params$brewerColor = c("Set3","Set1","Accent","Dark2","Spectral","PuOr","BrBG")


###################################################
### code chunk number 19: plot4
###################################################
params$legende  = TRUE


###################################################
### code chunk number 20: ttest1
###################################################
params$ttestPerform = FALSE


###################################################
### code chunk number 21: ttest
###################################################
#params$ttestCol = NULL


###################################################
### code chunk number 22: ttest3
###################################################
#params$pairsCol = NULL


###################################################
### code chunk number 23: ttest4
###################################################
#params$samplesRemainInTTest	 = NULL
#params$samplesNotInTTest	 = NULL


###################################################
### code chunk number 24: correlation
###################################################
#params$samplesRemainInCor = NULL
#params$samplesNotInCor    = NULL


###################################################
### code chunk number 25: execute ddct
###################################################
ddCtExec(params)


###################################################
### code chunk number 26: Sweave (eval = FALSE)
###################################################
## Sweave("rtPCR-usage.Rnw")


###################################################
### code chunk number 27: End of session
###################################################
toLatex(sessionInfo())

