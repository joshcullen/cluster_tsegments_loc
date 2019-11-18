

library('Rcpp')
library('MCMCpack')
library(dplyr)
library(purrr)
library(tidyr) #for gather function
library(ggplot2)
library(ggnewscale) #for multiple fill scales in ggplot2
library(pals) # for more color palettes
library(progress) #for progress bar
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(rgdal)



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

#basic settings
ngibbs = 1000
nburn=ngibbs/2

#priors
psi=0.01
gamma1=0.1



## ID 1 ##

dat1.res<- gibbs.clust.space(dat = obs.list$`1`[,-1], ngibbs = ngibbs, nburn = nburn,
                             nclustmax = ncol(obs.list$`1`)-1)

plot(dat1.res$loglikel,type='l')
plot(dat1.res$phi[ngibbs,],type='h')
plot(dat1.res$z[ngibbs,],type='h')

MAP1<- find.MAP(dat = dat1.res, nburn = nburn) #iteration 541
tbsp.clust1<- dat1.res$z[MAP1,]
table(tbsp.clust1)

time.seg<- 1:nrow(obs.list$`1`)
tbsp.clust1<- cbind(tbsp.clust1,time.seg) %>% data.frame()
dat.list$`1`<- left_join(dat.list$`1`, tbsp.clust1, by="time.seg")
nclust<- tbsp.clust1$tbsp.clust1 %>% unique() %>% length()



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
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = viridis(n=9)[7], size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs1.long$key) + 0.5,
                                    ymax = max(obs1.long$key) + 1.0, fill = tbsp.clust1),
                                    color = NA, size = 1.5) +
  scale_fill_gradientn("Cluster", colours = ocean.amp(nclust)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


#Plot geographic map of clusters
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat.list$`1`, aes(x = utmlong, y = utmlat), size=0.35) +
  geom_point(data = dat.list$`1`, aes(utmlong, utmlat, fill=as.numeric(tbsp.clust1)), size=2,
             pch=21, stroke=0.25) +
  scale_fill_gradientn("Cluster", colours = ocean.amp(nclust)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  facet_wrap(~tbsp.clust1)





## ID 12 ##

dat12.res<- gibbs.clust.space(dat = obs.list$`12`[,-1], ngibbs = ngibbs, nburn = nburn,
                             nclustmax = ncol(obs.list$`12`)-1)

plot(dat12.res$loglikel,type='l')
plot(dat12.res$phi[ngibbs,],type='h')
plot(dat12.res$z[ngibbs,],type='h')

MAP12<- find.MAP(dat = dat12.res, nburn = nburn)  #iteration 666
tbsp.clust12<- dat12.res$z[MAP12,]
table(tbsp.clust12)

time.seg<- 1:nrow(obs.list$`12`)
tbsp.clust12<- cbind(tbsp.clust12,time.seg) %>% data.frame()
dat.list$`12`<- left_join(dat.list$`12`, tbsp.clust12, by="time.seg")





## ID 19 ##

dat19.res<- gibbs.clust.space(dat = obs.list$`19`[,-1], ngibbs = ngibbs, nburn = nburn,
                              nclustmax = ncol(obs.list$`19`)-1)

plot(dat19.res$loglikel,type='l')
plot(dat19.res$phi[ngibbs,],type='h')
plot(dat19.res$z[ngibbs,],type='h')

MAP19<- find.MAP(dat = dat19.res, nburn = nburn)  #iteration 945
tbsp.clust19<- dat19.res$z[MAP19,]
table(tbsp.clust19)

time.seg<- 1:nrow(obs.list$`19`)
tbsp.clust19<- cbind(tbsp.clust19,time.seg) %>% data.frame()
dat.list$`19`<- left_join(dat.list$`19`, tbsp.clust19, by="time.seg")




## ID 27 ##

dat27.res<- gibbs.clust.space(dat = obs.list$`27`[,-1], ngibbs = ngibbs, nburn = nburn,
                              nclustmax = ncol(obs.list$`27`)-1)

plot(dat27.res$loglikel,type='l')
plot(dat27.res$phi[ngibbs,],type='h')
plot(dat27.res$z[ngibbs,],type='h')

MAP27<- find.MAP(dat = dat27.res, nburn = nburn)  # iteration 633
tbsp.clust27<- dat27.res$z[MAP27,]
table(tbsp.clust27)

time.seg<- 1:nrow(obs.list$`27`)
tbsp.clust27<- cbind(tbsp.clust27,time.seg) %>% data.frame()
dat.list$`27`<- left_join(dat.list$`27`, tbsp.clust27, by="time.seg")


