rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\GIT_models\\cluster_tsegments_loc')
library('Rcpp')
sourceCpp('aux1.cpp')
source('gibbs functions.R')

dat=read.csv('fake data.csv',as.is=T)
dat=data.matrix(dat)
n=rowSums(dat)
nobs=nrow(dat)
nloc=ncol(dat)
lo=0.000000000000001
  
#priors
psi=0.01
gamma1=0.1

#starting values
nclustmax=20
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
