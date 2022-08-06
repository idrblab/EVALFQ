## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library('EVALFQ')
set.seed(135)

## ----kable1-------------------------------------------------------------------
load("../data/my_spiked.rda")
dim(my_spiked)
head(my_spiked[1:4])
load("../data/spiked_data.rda")
dim(spiked_data)
head(spiked_data[1:5])

## ----kable--------------------------------------------------------------------
load("../data/allranks.rda")
head(allranks)

