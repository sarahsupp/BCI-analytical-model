# Functions for UBC Working Group on assessing change in biodiversity manuscript
# modified by Sarah Supp, from Fangliang He, based on He and Legendre 20002, and He 2012
# May 4-7, 2015

divmetrics.main = function(data0,data03,data05,data08) {
  # calculate diversity change using: 
  # log SR, log Shannon, log Simpson, slope SR, slope abundance, slope Shannon, slope Simpson, 
  # Jaccard, Bray-Curtis
  
  effect.size03=numeric()
  effect.size05=numeric()
  effect.size08=numeric()
  
  Hshannon00=numeric()
  Hshannon03=numeric() 
  Hshannon05=numeric()
  Hshannon08=numeric()
  
  Hsimpson00=numeric()
  Hsimpson03=numeric()
  Hsimpson05=numeric()
  Hsimpson08=numeric()
  
  bc03=numeric()
  bc05=numeric()
  bc08=numeric()
  
  data0=data0[,-1]
  data03=data03[,-1]
  data05=data05[,-1]
  data08=data08[,-1]
  
  ncell=dim(data0)[1]
  
  print(ncell)
  
  for (i in 1:ncell) {
    sp00=data0[i,]
    sp03=data03[i,]
    sp05=data05[i,]
    sp08=data08[i,]
    
    # richness effect size
    nxsp00=length(sp00[sp00>0])
    nxsp03=length(sp03[sp03>0])
    nxsp05=length(sp05[sp05>0])
    nxsp08=length(sp08[sp08>0])
    
    effect.size03[i]=log(nxsp03/nxsp00)
    effect.size05[i]=log(nxsp05/nxsp00)
    effect.size08[i]=log(nxsp08/nxsp00)
    
    # Hill's Shannon
    p00=sp00/sum(sp00)
    p03=sp03/sum(sp03)
    p05=sp05/sum(sp05)
    p08=sp08/sum(sp08)
    
    Hshannon00[i]=exp(-sum(p00[p00>0]*log(p00[p00>0])))
    Hshannon03[i]=exp(-sum(p03[p03>0]*log(p03[p03>0]))) 
    Hshannon05[i]=exp(-sum(p05[p05>0]*log(p05[p05>0])))
    Hshannon08[i]=exp(-sum(p08[p08>0]*log(p08[p08>0])))
    
    # Hill's Simpson
    Hsimpson00[i]=1/sum(p00^2)
    Hsimpson03[i]=1/sum(p03^2) 
    Hsimpson05[i]=1/sum(p05^2)
    Hsimpson08[i]=1/sum(p08^2)
    
    # Bray-Curtis
    bc03[i]=sum(abs(sp00-sp03))/sum(sp00+sp03)
    bc05[i]=sum(abs(sp00-sp05))/sum(sp00+sp05)
    bc08[i]=sum(abs(sp00-sp08))/sum(sp00+sp08)
    
    print(i)
  }
  
  stress1=c(rep(0.3,ncell),rep(0.5,ncell),rep(0.8,ncell))
  stress2=c(rep(0,ncell),rep(0.3,ncell),rep(0.5,ncell),rep(0.8,ncell))
  
  effect.size.out=data.frame(stress1,effect.size=c(effect.size03,effect.size05,effect.size08))
  Hshannon.out=data.frame(stress2,Hshannon=c(Hshannon00,Hshannon03,Hshannon05,Hshannon08))
  Hsimpson.out=data.frame(stress2,Hsimpson=c(Hsimpson00,Hsimpson03,Hsimpson05,Hsimpson08))
  bc.out=data.frame(stress1,bc=c(bc03,bc05,bc08))
  
  par(mfrow=c(1,4))
  
  boxplot(effect.size~stress1, data=effect.size.out,main="Richness effect size")
  boxplot(Hshannon~stress2, data=Hshannon.out,main="Hill Shannon")
  boxplot(Hsimpson~stress2, data=Hsimpson.out, main="Hill Simpson")
  boxplot(bc~stress1, data=bc.out, main="Bray-Curtis")
  
  # return(effect.size.out)
}