sample.theta=function(dat,nclustmax,nloc,z,psi){
  theta=matrix(NA,nclustmax,nloc)
  for (i in 1:nclustmax){
    cond=z==i
    soma=sum(cond)
    if (soma==0) tmp=rep(0,nloc)
    if (soma==1) tmp=dat[cond,]
    if (soma >1) tmp=colSums(dat[cond,])
    theta[i,]=rdirichlet(1,tmp+psi)
  }
  theta
}
#----------------------------------------------
sample.v=function(z,nclustmax,gamma1){
  #get n
  n=rep(0,nclustmax)
  tmp=table(z)
  n[as.numeric(names(tmp))]=tmp

  #get ngreater
  seq1=nclustmax:1
  tmp=cumsum(n[seq1])[seq1]
  ngreater=tmp[-1]
  
  #get v's
  rbeta(nclustmax-1,n[-nclustmax]+1,ngreater+gamma1)
}
#----------------------------------------------
sample.z=function(dat,theta,phi,nobs,nclustmax,nloc,z,n){
  #pre-calculate some useful quantities
  ltheta=log(theta)
  lphi=log(phi)
  
  #determine the number of locations in each group
  tab=rep(0,nclustmax)
  tmp=table(z)
  tab[as.numeric(names(tmp))]=tmp
  
  #calculate log-probability
  tmp=matrix(NA,nobs,nclustmax)
  for (i in 1:nclustmax){
    ltheta1=matrix(ltheta[i,],nobs,nloc,byrow=T)
    tmp[,i]=rowSums(dat*ltheta1)+lphi[i]
  }

  #sample z
  for (i in 1:nobs){
    tab[z[i]]=tab[z[i]]-1
    prob=rep(NA,nclustmax)
    cond=tab==0
    prob[ cond]=lphi[z[i]]+lgamma(nloc*psi)-nloc*lgamma(psi)+sum(lgamma(dat[i,]+psi))-lgamma(n[i]+nloc*psi) #log probability for a new group 
    prob[!cond]=tmp[i,!cond]
   
    #get normalized probs
    tmp1=prob-max(prob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1
    
    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    z[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  z
}
