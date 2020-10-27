library(splines)
seq = seq(0,1,1/1e4)


nsbasis <- ns(seq,4)
svdnsbasis <- svd(nsbasis)$u


matplot(nsbasis,type='l')

matplot(svdnsbasis,type='l')
library(tidyverse)

svdnsbasis %*% runif(4,0,1) %>% matplot(type='l')


data <- c(16,16.2,17,18,18.2)

fit <- lm( data ~ svd(ns(1:5,3))$u)

plot(data)
points(predict(fit),type='l')



