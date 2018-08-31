#rarity GAM impdet

############################
#prepare data
############################


setwd("~/Desktop/Spatial")
library(parallel)

#load data
spdata=read.csv("spdata.csv")

data=read.csv("data.csv")
data=data[,-1]

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

#FIdata
FIdata=read.csv("FI.csv")
FIdatafilt=FIdata[match(names(data),FIdata$Taxon_FA),]

#endemics only
dataend=data[,which(FIdatafilt$Endemicity==1)]
dataend=dataend[which(spdata$end>0) ,]

#rarity data
rartab=matrix(0,nrow=248,ncol=358)

for(i in 1:358){
  for(j in 1:248){
    rartab[j,i]=(1/sum(dataend[,j]))*dataend[i,j]
  }
}

rartabsort=matrix(0,nrow=248,ncol=358)
for(k in 1:358){
  rartabsort[,k]=sort(rartab[,k],decreasing=T)
}

lrartabsort=log(rartabsort)

#jagam
library(mgcv)

spdata$meanaltkm=spdata$newmeanalt/1000
spdata$rangekm=spdata$newrange/1000
spdata$lrar=log(spdata$rar)

spdatafilt=spdata[which(spdata$end!=0),]

jm=jagam((lrar)~s(LON,LAT)+meanaltkm+rangekm+rock+rocksil+caref+siref,data=spdatafilt,file="jagamrar.bug")

jmnew=list(n=jm$jags.data$n,
           S1=jm$jags.data$S1,
           X=jm$jags.data$X,
           zero=jm$jags.data$zero,
           nend=colSums(sign(rartabsort)),
           lrar=lrartabsort)

#run model
library(rjags)
load.module("glm")


model.rargam=jags.model("jagamrarnewln.bug",data=jmnew,n.chains=5)
update(model.rargam,20000)
sims.rargam=jags.samples(model.rargam,c("b","rho","rarmean","lambda"),n.iter=50000,n.thin=50)
msims.rargam=lapply(sims.rargam,as.mcmc.list)
summary(msims.rargam$b)
plot(msims.rargam$b)
gelman.diag(msims.rargam$b)


save(sims.rargam,file="rargamwd.RData")

#plot
rarscores=summary(msims.rargam$rarmean)[[1]][,1]
rarcoefs=summary(msims.rargam$b)[[1]][,1]
rarsmooth=jmnew$X[,8:36]%*%rarcoefs[8:36]
#linear res rarscoresres=rarscores-rarcoefs[3]*spdatafilt$rangekm-rarcoefs[4]*spdatafilt$rock
#smooth res rarscoresres=rarscores-rarcoefs[3]*spdata$rangekm-rarcoefs[4]*spdata$rock-rarsmooth
#smooth res by categories 
rarscoresres=rarscores-rarcoefs[3]*spdatafilt$rangekm-rarcoefs[4]*spdatafilt$rock-rarcoefs[5]*spdatafilt$rocksil-spdatafilt$siref*mean(rarsmooth[which(spdatafilt$siref==1)])-spdatafilt$caref*mean(rarsmooth[which(spdatafilt$caref==1)])-mean(rarsmooth)
#non logit rarscoresres=rarscores/((1-rarscores)*exp(rarcoefs[3]*spdata$rangekm+rarcoefs[4]*spdata$rock)+rarscores)


par(mfrow=c(1,2))

plot(rarscoresres~spdatafilt$meanaltkm,xlab="mean elevation[km]",ylab="endemics rarity",cex=0, cex.lab=0.6,cex.axis=0.6)

library(gtools)

points(spdatafilt$meanaltkm[which(spdatafilt$caref==0 & spdatafilt$siref==0)],rarscoresres[which(spdatafilt$caref==0 & spdatafilt$siref==0)],cex=0.5)
points(spdatafilt$meanaltkm[which(spdatafilt$caref==1 & spdatafilt$siref==0)],rarscoresres[which(spdatafilt$caref==1 & spdatafilt$siref==0)],col="blue", pch=2,cex=0.5)
points(spdatafilt$meanaltkm[which(spdatafilt$caref==0 & spdatafilt$siref==1)],rarscoresres[which(spdatafilt$caref==0 & spdatafilt$siref==1)],col="red",pch=6,cex=0.5)
points(spdatafilt$meanaltkm[which(spdatafilt$caref==1 & spdatafilt$siref==1)],rarscoresres[which(spdatafilt$caref==1 & spdatafilt$siref==1)],col="purple",pch=0,cex=0.5)

curve((rarcoefs[1]+x*(rarcoefs[2])),add=T)
curve((rarcoefs[1]+rarcoefs[6]+x*(rarcoefs[2])),add=T,col="blue",lwd=3, lty=6)
curve((rarcoefs[1]+rarcoefs[7]+x*(rarcoefs[2])),add=T,col="red",lty=5)

text(2.4,-4.9,labels="Calc. refugia",col="blue", cex=0.6)
points(1.95,-4.9, pch=2, col="blue",cex=0.7)

text(2.4,-4.97,labels="Silic. refugia",col="red", cex=0.6)
points(1.95,-4.97, pch=6, col="red",cex=0.7)

text(2.4,-5.05,labels="Both refugia",col="purple",cex=0.6)
points(1.95,-5.05, pch=0, col="purple",cex=0.7)

text(2.4,-5.12,labels="No refugia",col="black",cex=0.6)
points(1.95,-5.12, col="black",cex=0.7)

#spline
jgs=sim2jam(sims.rargam,jm$pregam)
plot(jgs,se=F,scheme=2,main="",xlab="LON",ylab="LAT",cex.lab=0.6,cex.axis=0.6,hcolors=rev(heat.colors(80)))

points(5.7245,45.1885,pch=16)
text(5.7245+1,45.1885-0.3,labels="Grenoble",cex=0.6)

points(11.4041,47.2692,pch=16)
text(11.4041+1,47.2692-0.3,labels="Innsbruck",cex=0.6)

points(9.1859,45.4654,pch=16)
text(9.1859+1,45.4654-0.3,labels="Milan",cex=0.6)

points(14.5058,46.0569,pch=16)
text(14.5058+0.5,46.0569-0.3,labels="Ljubljana",cex=0.6)

