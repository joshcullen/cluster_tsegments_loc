get.summary.stats_obs=function(dat){  #dat must have time.seg assigned; for all IDs
  
  #change values of grid cells for easy manipulation
  dat$grid.cell<- as.factor(dat$grid.cell)
  levels(dat$grid.cell)<- 1:length(levels(dat$grid.cell))
  dat$grid.cell<- as.numeric(dat$grid.cell)
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  
  #calculate # of obs in each grid.cell by time.seg
  for (i in 1:length(dat.list)) {
    ntseg=max(dat.list[[i]]$time.seg)
    nloc=length(unique(dat$grid.cell))
    res=matrix(0, ntseg, nloc)
    colnames(res)=1:nloc
    
    for (j in 1:ntseg){
      ind=dat.list[[i]] %>% filter(time.seg==j) %>% group_by(grid.cell) %>% count()
      res[j,ind$grid.cell]=ind$n #takes count of each cluster within given time segment
    }
    id<- rep(unique(dat.list[[i]]$id), ntseg)
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