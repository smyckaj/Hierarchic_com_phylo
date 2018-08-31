#endemism GAM

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

#endemism
FIdata=read.csv("FI.csv")
FIdatafilt=FIdata[match(names(data),FIdata$Taxon_FA),]

#endemics only
dataendnf=data[,which(FIdatafilt$Endemicity==1)]

spdata$endnew=rowSums(dataendnf)

#diversity
spdata$div=rowSums(data)

#jagam
library(mgcv)
library(gtools)

spdata$meanaltkm=spdata$newmeanalt/1000
spdata$rangekm=spdata$newrange/1000

jm=jagam(inv.logit(relend+0.0001)~s(LON,LAT)+meanaltkm+rangekm+rock+rocksil+caref+siref,data=spdata,file="jagamend.bug")
jmnew=list(n=jm$jags.data$n,
           S1=jm$jags.data$S1,
           X=jm$jags.data$X,
           zero=jm$jags.data$zero,
           end=spdata$endnew,
           div=spdata$div)

#!!!manually modifyjagam to jagamnew

#run model
library(rjags)
load.module("glm")

model.endgam=jags.model("jagamendnewbin.bug",data=jmnew,n.chains=5)
update(model.endgam,20000)
sims.endgam=jags.samples(model.endgam,c("b","rho","mu","lambda"),n.iter=50000,n.thin=50)
msims.endgam=lapply(sims.endgam,as.mcmc.list)
summary(msims.endgam$b)
plot(msims.endgam$b)
gelman.diag(msims.endgam$b)

save(sims.endgam,file="endgamwd.RData")

#plot
library(gtools)
relendscores=logit(summary(msims.endgam$mu)[[1]][,1])
relendcoefs=summary(msims.endgam$b)[[1]][,1]
relendsmooth=jmnew$X[,8:36]%*%relendcoefs[8:36]
#smooth with linear terms relendscoresres=relendscores-relendcoefs[3]*spdata$rangekm-relendcoefs[4]*spdata$rock
#smooth res relendscoresres=relendscores-relendcoefs[3]*spdata$rangekm-relendcoefs[4]*spdata$rock-relendsmooth
#smooth res by categories 
relendscoresres=relendscores-relendcoefs[3]*spdata$rangekm-relendcoefs[4]*spdata$rock-relendcoefs[5]*spdata$rocksil-spdata$siref*mean(relendsmooth[which(spdata$siref==1)])-spdata$caref*mean(relendsmooth[which(spdata$caref==1)])-mean(relendsmooth)
#non logit relendscoresres=relendscores/((1-relendscores)*exp(relendcoefs[3]*spdata$rangekm+relendcoefs[4]*spdata$rock)+relendscores)

par(mfrow=c(2,2))

plot(relendscoresres~spdata$meanaltkm,xlab="mean elevation[km]",ylab="logit(endemics/species)",cex=0, cex.lab=0.6,cex.axis=0.6)

points(spdata$meanaltkm[which(spdata$caref==0 & spdata$siref==0)],relendscoresres[which(spdata$caref==0 & spdata$siref==0)],cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==1 & spdata$siref==0)],relendscoresres[which(spdata$caref==1 & spdata$siref==0)],col="blue", pch=2,cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==0 & spdata$siref==1)],relendscoresres[which(spdata$caref==0 & spdata$siref==1)],col="red",pch=6,cex=0.5)
points(spdata$meanaltkm[which(spdata$caref==1 & spdata$siref==1)],relendscoresres[which(spdata$caref==1 & spdata$siref==1)],col="purple",pch=0,cex=0.5)

#make sure if this should indeed be inverse logit
curve((relendcoefs[1]+x*(relendcoefs[2])),add=T,lwd=3)
curve((relendcoefs[1]+relendcoefs[6]+x*(relendcoefs[2])),add=T,col="blue",lwd=3, lty=6)
curve((relendcoefs[1]+relendcoefs[7]+x*(relendcoefs[2])),add=T,col="red",lwd=3,lty=5)

text(2.4,-3.2,labels="Calc. refugia",col="blue", cex=0.6)
points(1.95,-3.2, pch=2, col="blue",cex=0.7)

text(2.4,-3.35,labels="Silic. refugia",col="red", cex=0.6)
points(1.95,-3.35, pch=6, col="red",cex=0.7)

text(2.4,-3.5,labels="Both refugia",col="purple",cex=0.6)
points(1.95,-3.5, pch=0, col="purple",cex=0.7)

text(2.4,-3.65,labels="No refugia",col="black",cex=0.6)
points(1.95,-3.65, col="black",cex=0.7)

#spline
jgs=sim2jam(sims.endgam,jm$pregam)
plot(jgs,se=F,scheme=2,main="",xlab="LON",ylab="LAT",cex.lab=0.6,cex.axis=0.6, hcolors=rev(heat.colors(80)))

points(5.7245,45.1885,pch=16)
text(5.7245+1,45.1885-0.3,labels="Grenoble",cex=0.6)

points(11.4041,47.2692,pch=16)
text(11.4041+1,47.2692-0.3,labels="Innsbruck",cex=0.6)

points(9.1859,45.4654,pch=16)
text(9.1859+1,45.4654-0.3,labels="Milan",cex=0.6)

points(14.5058,46.0569,pch=16)
text(14.5058+0.5,46.0569-0.3,labels="Ljubljana",cex=0.6)





