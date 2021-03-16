print('foo')
library(tidyverse)
library(here)

#now plot
plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)
mtcars%>%
	ggplot(.,aes(y=mpg,x=gear))+
	geom_point()+
	scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0('xname'))+
	scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('title'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

