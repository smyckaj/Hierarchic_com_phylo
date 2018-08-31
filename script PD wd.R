#PD GAM
#NTND is PD in this script
############################
#prepare data
############################


setwd("~/Desktop/Spatial/ML species trees")
library(ape)
library(parallel)
files=list.files()
trees=list()
for (i in 1:length(files)){
  tr=read.tree(files[i])
  trees[[i]]=tr
}


setwd("~/Desktop/Spatial")

##getting distribution of mntd

#load data
data=read.csv("data.csv")
row.names(data)=data$X
data=data[,-1]
spdata=read.csv("spdata.csv")

#angiosperms only
m=match(c("Lycopodium_clavatum",
          "Diphasiastrum_alpinum",
          "Diphasiastrum_oellgaardii",
          "Huperzia_selago", 
          "Isoetes_echinospora",
          "Selaginella_selaginoides",
          "Pinus_cembra",
          "Pinus_mugo",
          "Larix_decidua",
          "Juniperus_communis",
          "Juniperus_sabina",
          "Botrychium_simplex",
          "Equisetum_arvense",
          "Cryptogramma_crispa",
          "Asplenium_fissum",
          "Asplenium_viride",
          "Cystopteris_alpina",
          "Cystopteris_montana",
          "Woodsia_alpina",
          "Woodsia_pulchella",
          "Athyrium_distentifolium",
          "Polystichum_lonchitis",
          "Dryopteris_villarii"),names(data))
data=data[,-m]

angonly=function(tree){
  tree=drop.tip(tree,c("Lycopodium_clavatum",
                       "Diphasiastrum_alpinum",
                       "Diphasiastrum_oellgaardii",
                       "Huperzia_selago", 
                       "Isoetes_echinospora",
                       "Selaginella_selaginoides",
                       "Pinus_cembra",
                       "Pinus_mugo",
                       "Larix_decidua",
                       "Juniperus_communis",
                       "Juniperus_sabina",
                       "Botrychium_simplex",
                       "Equisetum_arvense",
                       "Cryptogramma_crispa",
                       "Asplenium_fissum",
                       "Asplenium_viride",
                       "Cystopteris_alpina",
                       "Cystopteris_montana",
                       "Woodsia_alpina",
                       "Woodsia_pulchella",
                       "Athyrium_distentifolium",
                       "Polystichum_lonchitis",
                       "Dryopteris_villarii"))
}
mltree=angonly(mltree)
trees=mclapply(trees,angonly,mc.cores = 10)

library(PhyloMeasures)

mntds=function(tree){
  datareshuf=data[ ,match(tree$tip.label,names(data))]
  mntd=pd.query(tree,datareshuf,standardize=T)
  mntd
}

mntds=mclapply(trees,mntds,mc.cores = 10)
mntdsar=simplify2array(mntds)

#distance matrix
dist=as.matrix(dist(cbind(spdata$LON,spdata$LAT),diag=T,upper=T))

library(mgcv)

spdata$meanaltkm=spdata$newmeanalt/1000
spdata$rangekm=spdata$newrange/1000

jm=jagam(sspda~s(LON,LAT)+meanaltkm+rangekm+rock+rocksil+caref+siref,data=spdata,file="jagam.bug")
jmnew=list(n=jm$jags.data$n,
           S1=jm$jags.data$S1,
           X=jm$jags.data$X,
           zero=jm$jags.data$zero,
           mntd=mntdsar)

#!!!manually modifyjagam to jagamnew

#run model
library(rjags)
load.module("glm")

model.mntdgam=jags.model("jagamnew.bug",data=jmnew,n.chains=5)
update(model.mntdgam,20000)
sims.mntdgam=jags.samples(model.mntdgam,c("b","rho","mntdmean","lambda"),n.iter=50000,n.thin=50)
msims.mntdgam=lapply(sims.mntdgam,as.mcmc.list)
summary(msims.mntdgam$b)
plot(msims.mntdgam$b)
gelman.diag(msims.mntdgam$b)

save(sims.mntdgam,file="PDgamwd.RData")

#plot
mntdscores=summary(msims.mntdgam$mntdmean)[[1]][,1]
mntdcoefs=summary(msims.mntdgam$b)[[1]][,1]
mntdsmooth=jmnew$X[,8:36]%*%mntdcoefs[8:36]
#lin res mntdscoresres=mntdscores-mntdcoefs[3]*spdata$rangekm-mntdcoefs[4]*spdata$rock

#smooth res mntdscoresres=mntdscores-mntdcoefs[3]*spdata$rangekm-mntdcoefs[4]*spdata$rock-mntdsmooth
#smooth res by categories 
mntdscoresres=mntdscores-mntdcoefs[3]*spdata$rangekm-mntdcoefs[4]*spdata$rock-mntdcoefs[5]*spdata$rocksil-spdata$siref*mean(mntdsmooth[which(spdata$siref==1)])-spdata$caref*mean(mntdsmooth[which(spdata$caref==1)])-mean(mntdsmooth)
#non logit mntdscoresres=mntdscores/((1-mntdscores)*exp(mntdcoefs[3]*spdata$rangekm+mntdcoefs[4]*spdata$rock)+mntdscores)

par(mfrow=c(1,2))

plot(mntdscoresres~spdata$meanaltkm,xlab="mean elevation[km]",ylab="phylogenetic diversity",cex=0, cex.lab=0.6,cex.axis=0.6)

library(gtools)

points(spdata$meanaltkm[which(spdata$caref==0 & spdata$siref==0)],mntdscoresres[which(spdata$caref==0 & spdata$siref==0)],cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==1 & spdata$siref==0)],mntdscoresres[which(spdata$caref==1 & spdata$siref==0)],col="blue", pch=2,cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==0 & spdata$siref==1)],mntdscoresres[which(spdata$caref==0 & spdata$siref==1)],col="red",pch=6,cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==1 & spdata$siref==1)],mntdscoresres[which(spdata$caref==1 & spdata$siref==1)],col="purple",pch=0,cex=0.5)

curve((mntdcoefs[1]+x*(mntdcoefs[2])),add=T,lwd=3)
curve((mntdcoefs[1]+mntdcoefs[6]+x*(mntdcoefs[2])),add=T,col="blue",lwd=3, lty=6)
curve((mntdcoefs[1]+mntdcoefs[7]+x*(mntdcoefs[2])),add=T,col="red", lwd=3,lty=5)

text(2.4,2.5,labels="Calc. refugia",col="blue", cex=0.6)
points(1.95,2.5, pch=2, col="blue",cex=0.7)

text(2.4,2.2,labels="Silic. refugia",col="red", cex=0.6)
points(1.95,2.2, pch=6, col="red",cex=0.7)

text(2.4,1.9,labels="Both refugia",col="purple",cex=0.6)
points(1.95,1.9, pch=0, col="purple",cex=0.7)

text(2.4,1.6,labels="No refugia",col="black",cex=0.6)
points(1.95,1.6, col="black",cex=0.7)

#spline
jgs=sim2jam(sims.mntdgam,jm$pregam)
plot(jgs,se=F,scheme=2,main="",xlab="LON",ylab="LAT",cex.lab=0.6,cex.axis=0.6,hcolors=rev(heat.colors(80)))

points(5.7245,45.1885,pch=16)
text(5.7245+1,45.1885-0.3,labels="Grenoble",cex=0.6)

points(11.4041,47.2692,pch=16)
text(11.4041+1,47.2692-0.3,labels="Innsbruck",cex=0.6)

points(9.1859,45.4654,pch=16)
text(9.1859+1,45.4654-0.3,labels="Milan",cex=0.6)

points(14.5058,46.0569,pch=16)
text(14.5058+0.5,46.0569-0.3,labels="Ljubljana",cex=0.6)

