rm(list=ls(all=TRUE))
set.seed(1)

library('Rcpp')
library('MCMCpack')
library(tidyr) #for gather function
library(ggnewscale) #for multiple fill scales in ggplot2
library(pals) # for more color palettes


sourceCpp('aux1.cpp')
source('gibbs functions.R') #for clustering


## ID 1 ##

dat=read.csv('ID1 Seg x Clust.csv',header =T, sep = ",")
dat=data.matrix(dat)
dat=dat[which(apply(dat,1,sum)>10),]
n=rowSums(dat)
nobs=nrow(dat)
nloc=ncol(dat)
lo=0.000000000000001
  
#priors
psi=0.01
gamma1=0.1

#starting values
nclustmax=10
z=sample(1:nclustmax,size=nobs,replace=T)
theta=matrix(1/nloc,nclustmax,nloc)
phi=rep(1/nclustmax,nclustmax)

#store results
ngibbs=1000
store.phi=matrix(NA,ngibbs,nclustmax)
store.z=matrix(NA,ngibbs,nobs)
store.theta=matrix(NA,ngibbs,nclustmax*nloc)
store.loglikel=matrix(NA,ngibbs,1)

#gibbs sampler
nburn=ngibbs/2
for (i in 1:ngibbs){
  print(i)

  #occasionally re-order this
  if (i<nburn & i%%50==0){
    ind=order(phi,decreasing=T)
    theta=theta[ind,]
    phi=phi[ind]
    znew=z
    for (j in 1:nclustmax){
      znew[z==ind[j]]=j
    }
    z=znew
  }
  
  #draw samples from FCD's
  z=sample.z(dat=dat,theta=theta,phi=phi,
             nobs=nobs,nclustmax=nclustmax,nloc=nloc,z=z,n=n)
  # z=z.true
  
  v=sample.v(z=z,nclustmax=nclustmax,gamma1=gamma1)
  phi=GetPhi(vec=c(v,1),nclustmax=nclustmax)

  theta=sample.theta(dat=dat,nclustmax=nclustmax,nloc=nloc,z=z,psi=psi)
  #to avoid numerical issues
  theta[theta<lo]=lo
  # theta=theta.true
  
  #get logl
  tmp=sum(dat*log(theta)[z,])+sum(dbeta(v,1,gamma1,log=T))+sum((psi-1)*log(theta))
  
  #store results
  store.loglikel[i]=tmp
  store.theta[i,]=theta
  store.phi[i,]=phi
  store.z[i,]=z
}



plot(store.loglikel,type='l')

plot(store.phi[ngibbs,],type='h')
plot(store.theta[ngibbs,],type='h')
plot(store.z[ngibbs,],type='h')

#by col or row?
new.theta=matrix(store.theta[ngibbs,], nclustmax, nloc)
image(new.theta)
new.theta=matrix(store.theta[ngibbs,], nclustmax, nloc, byrow = T)
image(new.theta)

MAP1<- which(store.loglikel==max(store.loglikel))  #107 iteration is MAP





tbsp.clust<- store.z[MAP1[1],]
time.seg<- 1:nobs
tbsp.clust<- cbind(tbsp.clust,time.seg) %>% data.frame()
tbsp.clust$tbsp.clust<- as.factor(tbsp.clust$tbsp.clust)
levels(tbsp.clust$tbsp.clust)<- 1:length(levels(tbsp.clust$tbsp.clust))
dat1$time.seg<- as.factor(dat1$time.seg)
tbsp.clust$time.seg<- as.factor(tbsp.clust$time.seg)

dat1<- left_join(dat1, tbsp.clust, by="time.seg")



colnames(dat)=1:ncol(dat)
obs1.breakpts<- data.frame(breaks=obs1.breakpts)
obs1.long<- dat %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs1.long$key<- as.factor(obs1.long$key)
levels(obs1.long$key)<- 1:nloc
obs1.long$key<- as.numeric(obs1.long$key)

tbsp.clust[,1]<- tbsp.clust[,1] %>% as.numeric()
tbsp.clust[,2]<- tbsp.clust[,2] %>% as.numeric()



rect.lims<- rle(tbsp.clust$tbsp.clust)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust")
rect.lims.new<- data.frame(rect.lims.new)

ggplot() +
  geom_tile(data=obs1.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(breaks = 1:6, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = 6.5,
                                    ymax = 6.75, fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Spatial Cluster") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## ID 12 ##

dat=read.csv('ID12 Seg x Loc.csv',header =T, sep = ",")
dat=data.matrix(dat)
n=rowSums(dat)
nobs=nrow(dat)
nloc=ncol(dat)
lo=0.000000000000001

#priors
psi=0.01
gamma1=0.1

#starting values
nclustmax=nobs
z=sample(1:nclustmax,size=nobs,replace=T)
theta=matrix(1/nloc,nclustmax,nloc)
phi=rep(1/nclustmax,nclustmax)

#store results
ngibbs=1000
store.phi=matrix(NA,ngibbs,nclustmax)
store.z=matrix(NA,ngibbs,nobs)
store.theta=matrix(NA,ngibbs,nclustmax*nloc)
store.loglikel=matrix(NA,ngibbs,1)

#gibbs sampler
nburn=ngibbs/2
for (i in 1:ngibbs){
  print(i)
  
  #occasionally re-order this
  if (i<nburn & i%%50==0){
    ind=order(phi,decreasing=T)
    theta=theta[ind,]
    phi=phi[ind]
    znew=z
    for (j in 1:nclustmax){
      znew[z==ind[j]]=j
    }
    z=znew
  }
  
  #draw samples from FCD's
  z=sample.z(dat=dat,theta=theta,phi=phi,
             nobs=nobs,nclustmax=nclustmax,nloc=nloc,z=z,n=n)
  # z=z.true
  
  v=sample.v(z=z,nclustmax=nclustmax,gamma1=gamma1)
  phi=GetPhi(vec=c(v,1),nclustmax=nclustmax)
  
  theta=sample.theta(dat=dat,nclustmax=nclustmax,nloc=nloc,z=z,psi=psi)
  #to avoid numerical issues
  theta[theta<lo]=lo
  # theta=theta.true
  
  #get logl
  tmp=sum(dat*log(theta)[z,])+sum(dbeta(v,1,gamma1,log=T))+sum((psi-1)*log(theta))
  
  #store results
  store.loglikel[i]=tmp
  store.theta[i,]=theta
  store.phi[i,]=phi
  store.z[i,]=z
}



plot(store.loglikel,type='l')

plot(store.phi[ngibbs,],type='h')
plot(store.theta[ngibbs,],type='h')
plot(store.z[ngibbs,],type='h')










## ID 19 ##

dat=read.csv('ID19 Seg x Loc.csv',header =T, sep = ",")
dat=data.matrix(dat)
n=rowSums(dat)
nobs=nrow(dat)
nloc=ncol(dat)
lo=0.000000000000001

#priors
psi=0.01
gamma1=0.1

#starting values
nclustmax=nobs
z=sample(1:nclustmax,size=nobs,replace=T)
theta=matrix(1/nloc,nclustmax,nloc)
phi=rep(1/nclustmax,nclustmax)

#store results
ngibbs=1000
store.phi=matrix(NA,ngibbs,nclustmax)
store.z=matrix(NA,ngibbs,nobs)
store.theta=matrix(NA,ngibbs,nclustmax*nloc)
store.loglikel=matrix(NA,ngibbs,1)

#gibbs sampler
nburn=ngibbs/2
for (i in 1:ngibbs){
  print(i)
  
  #occasionally re-order this
  if (i<nburn & i%%50==0){
    ind=order(phi,decreasing=T)
    theta=theta[ind,]
    phi=phi[ind]
    znew=z
    for (j in 1:nclustmax){
      znew[z==ind[j]]=j
    }
    z=znew
  }
  
  #draw samples from FCD's
  z=sample.z(dat=dat,theta=theta,phi=phi,
             nobs=nobs,nclustmax=nclustmax,nloc=nloc,z=z,n=n)
  # z=z.true
  
  v=sample.v(z=z,nclustmax=nclustmax,gamma1=gamma1)
  phi=GetPhi(vec=c(v,1),nclustmax=nclustmax)
  
  theta=sample.theta(dat=dat,nclustmax=nclustmax,nloc=nloc,z=z,psi=psi)
  #to avoid numerical issues
  theta[theta<lo]=lo
  # theta=theta.true
  
  #get logl
  tmp=sum(dat*log(theta)[z,])+sum(dbeta(v,1,gamma1,log=T))+sum((psi-1)*log(theta))
  
  #store results
  store.loglikel[i]=tmp
  store.theta[i,]=theta
  store.phi[i,]=phi
  store.z[i,]=z
}



plot(store.loglikel,type='l')

plot(store.phi[ngibbs,],type='h')
plot(store.theta[ngibbs,],type='h')
plot(store.z[ngibbs,],type='h')










## ID 27 ##

dat=read.csv('ID27 Seg x Loc.csv',header =T, sep = ",")
dat=data.matrix(dat)
n=rowSums(dat)
nobs=nrow(dat)
nloc=ncol(dat)
lo=0.000000000000001

#priors
psi=0.01
gamma1=0.1

#starting values
nclustmax=nobs
z=sample(1:nclustmax,size=nobs,replace=T)
theta=matrix(1/nloc,nclustmax,nloc)
phi=rep(1/nclustmax,nclustmax)

#store results
ngibbs=1000
store.phi=matrix(NA,ngibbs,nclustmax)
store.z=matrix(NA,ngibbs,nobs)
store.theta=matrix(NA,ngibbs,nclustmax*nloc)
store.loglikel=matrix(NA,ngibbs,1)

#gibbs sampler
nburn=ngibbs/2
for (i in 1:ngibbs){
  print(i)
  
  #occasionally re-order this
  if (i<nburn & i%%50==0){
    ind=order(phi,decreasing=T)
    theta=theta[ind,]
    phi=phi[ind]
    znew=z
    for (j in 1:nclustmax){
      znew[z==ind[j]]=j
    }
    z=znew
  }
  
  #draw samples from FCD's
  z=sample.z(dat=dat,theta=theta,phi=phi,
             nobs=nobs,nclustmax=nclustmax,nloc=nloc,z=z,n=n)
  # z=z.true
  
  v=sample.v(z=z,nclustmax=nclustmax,gamma1=gamma1)
  phi=GetPhi(vec=c(v,1),nclustmax=nclustmax)
  
  theta=sample.theta(dat=dat,nclustmax=nclustmax,nloc=nloc,z=z,psi=psi)
  #to avoid numerical issues
  theta[theta<lo]=lo
  # theta=theta.true
  
  #get logl
  tmp=sum(dat*log(theta)[z,])+sum(dbeta(v,1,gamma1,log=T))+sum((psi-1)*log(theta))
  
  #store results
  store.loglikel[i]=tmp
  store.theta[i,]=theta
  store.phi[i,]=phi
  store.z[i,]=z
}



plot(store.loglikel,type='l')

plot(store.phi[ngibbs,],type='h')
plot(store.theta[ngibbs,],type='h')
plot(store.z[ngibbs,],type='h')
