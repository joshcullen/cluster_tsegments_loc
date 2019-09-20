rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(3)

nobs=1000
nloc=150

nclust=10
z.true=z=sample(1:nclust,size=nobs,replace=T)
theta.true=theta=rdirichlet(nclust,rep(0.01,nloc))
image(theta.true)
n=round(runif(nobs,min=50,max=250))

obs=matrix(NA,nobs,nloc)
for (i in 1:nobs){
  tmp=theta[z[i],]
  obs[i,]=rmultinom(1,size=n[i],prob=tmp)
}
image(obs)
rowSums(obs)[1:5]
n[1:5]

setwd('U:\\GIT_models\\cluster_tsegments_loc')
write.csv(obs,'fake data.csv',row.names=F)
