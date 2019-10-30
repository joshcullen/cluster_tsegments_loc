rm(list=ls(all=TRUE))
set.seed(1)

library('Rcpp')
library('MCMCpack')
library(dplyr)
library(tidyr) #for gather function
library(ggnewscale) #for multiple fill scales in ggplot2
library(pals) # for more color palettes


sourceCpp('aux1.cpp')
source('gibbs functions.R') #for clustering
source('gibbs sampler2.R')
source('helper functions.R')



# Load data
dat=read.csv('Snail Kite Gridded Data_AC.csv', header =T, sep = ",")
obs=get.summary.stats_obs(dat)
dat=dat[which(apply(dat,1,sum)>10),]