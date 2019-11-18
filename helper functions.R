get.summary.stats_ac=function(dat){  #dat must have time.seg assigned; for all IDs
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  
  #calculate # of obs in each AC by time.seg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$time.seg)
    nloc=length(unique(dat.ind$ac))
    res=matrix(0, ntseg, nloc)
    colnames(res)=1:nloc
    
    #Re-label ACs to be consecutive numbers
    dat.ind$ac<- as.factor(dat.ind$ac)
    levels(dat.ind$ac)<- 1:nloc
    dat.ind$ac<- as.numeric(dat.ind$ac)
    
    for (j in 1:ntseg){
      ind=dat.ind %>% filter(time.seg==j) %>% group_by(ac) %>% count()
      res[j,ind$ac]=ind$n #takes count of each cluster within given time segment
    }
    
    id<- rep(unique(dat.ind$id), ntseg)
    res=cbind(id, res) %>% data.frame()
    obs.list[[i]]=res
  }
  #obs<- do.call(rbind.data.frame, obs.list)
  obs<- map_dfr(obs.list, `[`)
  obs
}
#---------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#---------------------------------------------
find.MAP=function(dat, nburn) {  #select MAP value that is beyond burn-in phase
  if (length(max(dat$loglikel)) > 1) {
    stop("> 1 MAP value; inspect likelihood vector")
  } else if (which(dat$loglikel==max(dat$loglikel)) < nburn) {
  MAP<- dat$loglikel %>% order(decreasing = T) %>% subset(. > nburn) %>% first()
  } else {
  MAP<- which(dat$loglikel==max(dat$loglikel))
  }
  return(MAP)
}
