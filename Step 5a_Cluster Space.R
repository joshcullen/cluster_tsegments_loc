set.seed(1)

library('Rcpp')
library('MCMCpack')
library(dplyr)
library(purrr)
library(tidyr) #for gather function
library(ggplot2)
library(ggnewscale) #for multiple fill scales in ggplot2
library(pals) # for more color palettes
library(progress) #for progress bar


sourceCpp('aux1.cpp')
source('gibbs functions.R') #for clustering
source('gibbs sampler2.R')
source('helper functions.R')



# Load data
dat=read.csv('Snail Kite Gridded Data_AC.csv', header =T, sep = ",")
dat.list=df.to.list(dat)
obs=get.summary.stats_ac(dat)
obs.list=df.to.list(obs) %>% lapply(function(x) {as.matrix(x[,!is.na(colSums(x))])})  #remove ACs not used and convert to matrix from DF
# obs.list2<- map(obs.list, function(x) {x[which(apply(x,1,sum)>10),]})  #remove time segs w low obs



#################################
#### Run Gibbs Sampler by ID ####
#################################

#priors
psi=0.01
gamma1=0.1

#basic settings
ngibbs = 1000
nburn=ngibbs/2




## ID 1 ##

#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)


dat1.res<- gibbs.clust.space(dat = obs.list$`1`[,-1], ngibbs = ngibbs, nburn = nburn,
                             nclustmax = 6)

plot(dat1.res$loglikel,type='l')
plot(dat1.res$phi[ngibbs,],type='h')
plot(dat1.res$z[ngibbs,],type='h')

MAP1<- which(dat1.res$loglikel==max(dat1.res$loglikel))  # iteration 620 of MAP
# MAP1<- dat1.res$loglikel %>% order(decreasing = T) %>% subset(. > 500) %>% first() #iteration 541
tbsp.clust1<- dat1.res$z[MAP1,]
table(tbsp.clust1)  # 6 clusters

time.seg<- 1:nrow(obs.list$`1`)
tbsp.clust1<- cbind(tbsp.clust1,time.seg) %>% data.frame()
dat.list$`1`<- left_join(dat.list$`1`, tbsp.clust1, by="time.seg")


## Plot heatmap of clusters

#format data
colnames(obs.list$`1`)[-1]=1:ncol(obs.list$`1`[,-1]) %>% as.character()
nobs=nrow(obs.list$`1`)
nloc=ncol(obs.list$`1`[,-1])
obs1.long<- obs.list$`1`[,-1] %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs1.long$key<- as.factor(obs1.long$key)
levels(obs1.long$key)<- 1:nloc
obs1.long$key<- as.numeric(obs1.long$key)

tbsp.clust1[,1]<- tbsp.clust1[,1] %>% as.numeric()
tbsp.clust1[,2]<- tbsp.clust1[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust1$tbsp.clust1)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust1")
rect.lims.new<- data.frame(rect.lims.new)

#plot
ggplot() +
  geom_tile(data=obs1.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs1.long$key) + 0.5,
                                    ymax = max(obs1.long$key) + 0.75, fill = tbsp.clust1),
                                    color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## ID 12 ##

#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)


dat12.res<- gibbs.clust.space(dat = obs.list$`12`[,-1], ngibbs = ngibbs, nburn = nburn,
                             nclustmax = 6)

plot(dat12.res$loglikel,type='l')
plot(dat12.res$phi[ngibbs,],type='h')
plot(dat12.res$z[ngibbs,],type='h')

MAP12<- which(dat12.res$loglikel==max(dat12.res$loglikel))  # iteration 966 of MAP
# MAP12<- dat1.res$loglikel %>% order(decreasing = T) %>% subset(. > 500) %>% first() 
tbsp.clust12<- dat12.res$z[MAP1,]
table(tbsp.clust12)  # 6 clusters

time.seg<- 1:nrow(obs.list$`12`)
tbsp.clust12<- cbind(tbsp.clust12,time.seg) %>% data.frame()
dat.list$`12`<- left_join(dat.list$`12`, tbsp.clust12, by="time.seg")


## Plot heatmap of clusters

#format data
colnames(obs.list$`12`)[-1]=1:ncol(obs.list$`12`[,-1]) %>% as.character()
nobs=nrow(obs.list$`12`)
nloc=ncol(obs.list$`12`[,-1])
obs12.long<- obs.list$`12`[,-1] %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs12.long$key<- as.factor(obs12.long$key)
levels(obs12.long$key)<- 1:nloc
obs12.long$key<- as.numeric(obs12.long$key)

tbsp.clust12[,1]<- tbsp.clust12[,1] %>% as.numeric()
tbsp.clust12[,2]<- tbsp.clust12[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust12$tbsp.clust12)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust12")
rect.lims.new<- data.frame(rect.lims.new)

#plot
ggplot() +
  geom_tile(data=obs12.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs12.long$key) + 0.5,
                                    ymax = max(obs12.long$key) + 0.75, fill = tbsp.clust12),
            color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## ID 19 ##

#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)


dat19.res<- gibbs.clust.space(dat = obs.list$`19`[,-1], ngibbs = ngibbs, nburn = nburn,
                              nclustmax = 6)

plot(dat19.res$loglikel,type='l')
plot(dat19.res$phi[ngibbs,],type='h')
plot(dat19.res$z[ngibbs,],type='h')

MAP19<- which(dat19.res$loglikel==max(dat19.res$loglikel))  # iteration 511 of MAP
# MAP19<- dat1.res$loglikel %>% order(decreasing = T) %>% subset(. > 500) %>% first() 
tbsp.clust19<- dat19.res$z[MAP1,]
table(tbsp.clust19)  # 6 clusters

time.seg<- 1:nrow(obs.list$`19`)
tbsp.clust19<- cbind(tbsp.clust19,time.seg) %>% data.frame()
dat.list$`19`<- left_join(dat.list$`19`, tbsp.clust19, by="time.seg")


## Plot heatmap of clusters

#format data
colnames(obs.list$`19`)[-1]=1:ncol(obs.list$`19`[,-1]) %>% as.character()
nobs=nrow(obs.list$`19`)
nloc=ncol(obs.list$`19`[,-1])
obs19.long<- obs.list$`19`[,-1] %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs19.long$key<- as.factor(obs19.long$key)
levels(obs19.long$key)<- 1:nloc
obs19.long$key<- as.numeric(obs19.long$key)

tbsp.clust19[,1]<- tbsp.clust19[,1] %>% as.numeric()
tbsp.clust19[,2]<- tbsp.clust19[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust19$tbsp.clust19)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust19")
rect.lims.new<- data.frame(rect.lims.new)

#plot
ggplot() +
  geom_tile(data=obs19.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs19.long$key) + 0.5,
                                    ymax = max(obs19.long$key) + 0.75, fill = tbsp.clust19),
            color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## ID 27 ##

#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)


dat27.res<- gibbs.clust.space(dat = obs.list$`27`[,-1], ngibbs = ngibbs, nburn = nburn,
                              nclustmax = 3)

plot(dat27.res$loglikel,type='l')
plot(dat27.res$phi[ngibbs,],type='h')
plot(dat27.res$z[ngibbs,],type='h')

# MAP27<- which(dat27.res$loglikel==max(dat27.res$loglikel))  # iteration 966 of MAP
MAP27<- dat1.res$loglikel %>% order(decreasing = T) %>% subset(. > 500) %>% first() #iteration 541
tbsp.clust27<- dat27.res$z[MAP1,]
table(tbsp.clust27)  # 6 clusters

time.seg<- 1:nrow(obs.list$`27`)
tbsp.clust27<- cbind(tbsp.clust27,time.seg) %>% data.frame()
dat.list$`27`<- left_join(dat.list$`27`, tbsp.clust27, by="time.seg")


## Plot heatmap of clusters

#format data
colnames(obs.list$`27`)[-1]=1:ncol(obs.list$`27`[,-1]) %>% as.character()
nobs=nrow(obs.list$`27`)
nloc=ncol(obs.list$`27`[,-1])
obs27.long<- obs.list$`27`[,-1] %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs27.long$key<- as.factor(obs27.long$key)
levels(obs27.long$key)<- 1:nloc
obs27.long$key<- as.numeric(obs27.long$key)

tbsp.clust27[,1]<- tbsp.clust27[,1] %>% as.numeric()
tbsp.clust27[,2]<- tbsp.clust27[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust27$tbsp.clust27)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust27")
rect.lims.new<- data.frame(rect.lims.new)

#plot
ggplot() +
  geom_tile(data=obs27.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs27.long$key) + 0.5,
                                    ymax = max(obs27.long$key) + 0.75, fill = tbsp.clust27),
            color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))
