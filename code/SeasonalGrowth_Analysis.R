################################################################################
################################################################################
##
## Analysis script for ...
##
################################################################################
################################################################################

################################################################################
## SETUP
################################################################################
## Load required packages
library(FSAdata)   # for Bonito and Mosquitofish data
library(FSA)       # for Somers function function
library(nlstools)  # for nls model bootstrapping

## Create a Somers function
vbSO <- vbFuns("Somers")

## Create the Pauly et al. function
vbSCGF <- vbFuns("Pauly")

## Set the random number seed so that the bootstraps stay reproducible
set.seed(730987)

################################################################################
## Create the demo figure of the Somers model
################################################################################
pdf("results/Figure_1.PDF",width=4,height=4)
par(xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1)
ts <- 0.05; Linf <- 30; K <- 0.3; t0 <- -0.1
C <- c(0,0.5,1,2)
clr <- col2rgbt(rep("black",4),1/(1:4))
t <- seq(0,3,length.out=299)
plot(vbSO(t,Linf=c(Linf,K,t0,C=C[1],ts))~t,type="l",lwd=2,col=clr[1],
     ylim=c(0,20),ylab="Mean Length",xlab="Age (years)",xaxt="n")
axis(1,0:3)
lines(vbSO(t,Linf=c(Linf,K,t0,C=C[2],ts))~t,lwd=2,col=clr[2])
lines(vbSO(t,Linf=c(Linf,K,t0,C=C[3],ts))~t,lwd=2,col=clr[3])
lines(vbSO(t,Linf=c(Linf,K,t0,C=C[4],ts))~t,lwd=2,col=clr[4])
legend("topleft",paste("C",C,sep="="),lwd=2,col=clr,bty="n",inset=-0.02)
dev.off()


################################################################################
## Create the demo figure of the Pauly et al. (1992) function
################################################################################
pdf("results/Figure_2.PDF",width=4,height=4)
par(xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1)
Kpr <- 0.35; NGT <- 0.3
Pcf <- c(Linf,Kpr,t0,ts,NGT)
t <- seq(-0.1,3.3,length.out=499)
PL <- vbSCGF(t,Linf=Pcf)
plot(c(-1,-1),xlim=c(-0.1,3.5),ylim=c(0,22),xaxt="n",ylab="Length",xlab="Age (years)")
WPs <- 0:2+ts+0.5
LatWPs <- vbSCGF(WPs,Linf=Pcf)
SNGTs <- WPs-NGT/2
ENGTs <- WPs+NGT/2
for (i in 1:length(SNGTs))
  polygon(c(SNGTs[i],SNGTs[i],ENGTs[i],ENGTs[i]),c(0,LatWPs[i],LatWPs[i],0),
          col=col2rgbt("black",1/20),border=NA)
arrows(WPs,LatWPs-2,WPs,LatWPs-0.2,lwd=2,length=0.1)
arrows(SNGTs,LatWPs-1.2,ENGTs,LatWPs-1.2,lwd=2,length=0.025,angle=90,code=3)
lines(PL~t,lwd=2)
tss <- 0:3+ts
Lattss <- vbSCGF(tss,Linf=Pcf)
points(tss,Lattss,pch=16,col="gray50",cex=1.1)
axis(1,0:3)
axis(1,SNGTs,tcl=-0.2)
axis(1,ENGTs,tcl=-0.2)
# makes inside ticks
axis(1,at=c(0:3,SNGTs,ENGTs),labels=NA,tcl=0.2)
axis(1,at=1:3,labels=FSA:::iCalc_tpr(1:3,ts,NGT),tcl=0.2,line=-1.7,lwd=0)
axis(1,at=SNGTs,labels=FSA:::iCalc_tpr(SNGTs,ts,NGT),tcl=0.2,line=-1.7,lwd=0)
axis(1,at=ENGTs,labels=FSA:::iCalc_tpr(ENGTs,ts,NGT),tcl=0.2,line=-1.7,lwd=0)
axis(1,at=3.3,labels="<= t",tick=FALSE)
axis(1,at=3.3,labels="<= t'",tick=FALSE,line=-1.7,lwd=0)
dev.off()


################################################################################
## Example analysis with the Bonito data
################################################################################
data(Bonito)

## 1. Fit Somers function with C<=1 constraint (as in Stewart et al. (2013))
##    Linf=71.9, K=0.27, t0=-1.92, C=1, and ts=0.09 ...
##    they all match (within rounding)
Slwrbnd <- c(Linf=0,K=0,t0=-Inf,C=0,ts=0)
Suprbnd <- c(Linf=Inf,K=Inf,t0=Inf,C=1,ts=1)
SsvBon <- vbStarts(fl~age,data=Bonito,param="Somers",fixed=list(C=0.6,ts=0.2))
SfitBon <- nls(fl~vbSO(age,Linf,K,t0,C,ts),data=Bonito,
               start=SsvBon,lower=Slwrbnd,upper=Suprbnd,
               algorithm="port",control=list(maxiter=100))
SbootBon <- nlsBoot(SfitBon)
ScfBon <- cbind(Est=coef(SfitBon),confint(SbootBon))

## 2. Fit new Pauly et al. (1992) function
Plwrbnd <- c(Linf=0,Kpr=0,t0=-Inf,ts=0,NGT=0)
Puprbnd <- c(Linf=Inf,Kpr=Inf,t0=Inf,ts=1,NGT=1)
PsvBon <- vbStarts(fl~age,data=Bonito,param="Pauly",fixed=list(ts=0.25,NGT=0.2))
PfitBon <- nls(fl~vbSCGF(age,Linf,Kpr,t0,ts,NGT),data=Bonito,
             start=PsvBon,lower=Plwrbnd,upper=Puprbnd,
             algorithm="port",control=list(maxiter=100))
PbootBon <- nlsBoot(PfitBon)
PcfBon <- cbind(Est=coef(PfitBon),confint(PbootBon))

## 3. Summary results
ScfBon <- rbind(ScfBon,c(sum(summary(SfitBon)$residuals^2),NA,NA),
                c(AIC(SfitBon),NA,NA))
PcfBon <- rbind(PcfBon,c(sum(summary(PfitBon)$residuals^2),NA,NA),
                c(AIC(PfitBon),NA,NA))
rownames(ScfBon)[6:7] <- rownames(PcfBon)[6:7] <- c("RSS","AIC")
print(round(ScfBon,2),na.print="-")
print(round(PcfBon,2),na.print="-")


################################################################################
## Example analysis with the Bonito data -- Site 2
################################################################################
data(Mosquitofish)
mqf2 <- subset(Mosquitofish,sitenum==2)

## 1. Fit Somers function with C unconstrained (as in Carmona-Catot et al. (2014))
##    Linf=35.85, K=2.012, t0=-0.02, C=1.95, and ts=-0.118 ...
##    they all match (within rounding) except ts but it is off by +1
Suprbnd <- c(Linf=Inf,K=Inf,t0=Inf,C=Inf,ts=1)
Ssvmqf2 <- vbStarts(sl~age2,data=mqf2,param="Somers",fixed=list(C=1.5,ts=0.9))
Sfitmqf2 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf2,
                start=Ssvmqf2,lower=Slwrbnd,upper=Suprbnd,
                algorithm="port",control=list(maxiter=100))
Sbootmqf2 <- nlsBoot(Sfitmqf2)
Scfmqf2 <- cbind(Est=coef(Sfitmqf2),confint(Sbootmqf2))

## 2. Fit new Pauly et al. (1992) function
Psvmqf2 <- vbStarts(sl~age2,data=mqf2,param="Pauly",fixed=list(ts=0.9,NGT=0.5))
Pfitmqf2 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf2,
                start=Psvmqf2,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf2 <- nlsBoot(Pfitmqf2)
Pcfmqf2 <- cbind(Est=coef(Pfitmqf2),confint(Pbootmqf2))


## 3. Summary results
Scfmqf2 <- rbind(Scfmqf2,c(sum(summary(Sfitmqf2)$residuals^2),NA,NA),
                c(AIC(Sfitmqf2),NA,NA))
Pcfmqf2 <- rbind(Pcfmqf2,c(sum(summary(Pfitmqf2)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf2),NA,NA))
rownames(Scfmqf2)[6:7] <- rownames(Pcfmqf2)[6:7] <- c("RSS","AIC")
print(round(Scfmqf2,2),na.print="-")
print(round(Pcfmqf2,2),na.print="-")


################################################################################
## Example analysis with the Bonito data -- Site 4
################################################################################
mqf4 <- subset(Mosquitofish,sitenum==4)

Ssvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Somers",fixed=list(K=0.6,C=1.6,ts=0.9),plot=TRUE)
Sfitmqf4 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf4,
                start=Ssvmqf4,lower=Slwrbnd,upper=Suprbnd,
                algorithm="port",control=list(maxiter=1000))
Sbootmqf4 <- nlsBoot(Sfitmqf4)
Scfmqf4 <- cbind(Est=coef(Sfitmqf4),confint(Sbootmqf4))

## 2. Fit new Pauly et al. (1992) function
Psvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Pauly",fixed=list(ts=0.9,NGT=0.5))
Pfitmqf4 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf4,
                start=Psvmqf4,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf4 <- nlsBoot(Pfitmqf4)
Pcfmqf4 <- cbind(Est=coef(Pfitmqf4),confint(Pbootmqf4))

## 3. Summary results
Scfmqf4 <- rbind(Scfmqf4,c(sum(summary(Sfitmqf4)$residuals^2),NA,NA),
                 c(AIC(Sfitmqf4),NA,NA))
Pcfmqf4 <- rbind(Pcfmqf4,c(sum(summary(Pfitmqf4)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf4),NA,NA))
rownames(Scfmqf4)[6:7] <- rownames(Pcfmqf4)[6:7] <- c("RSS","AIC")
print(round(Scfmqf4,2),na.print="-")
print(round(Pcfmqf4,2),na.print="-")


################################################################################
## Example analysis with the Bonito data -- Site 9
################################################################################
mqf9 <- subset(Mosquitofish,sitenum==9)

Ssvmqf9 <- vbStarts(sl~age2,data=mqf9,param="Somers",fixed=list(t0=-2.5,C=0.6,ts=0.85))
Sfitmqf9 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf9,
                start=Ssvmqf9,lower=Slwrbnd,upper=Suprbnd,
                algorithm="port",control=list(maxiter=100))
Sbootmqf9 <- nlsBoot(Sfitmqf9)
Scfmqf9 <- cbind(Est=coef(Sfitmqf9),confint(Sbootmqf9))

## 2. Fit new Pauly et al. (1992) function
Psvmqf9 <- vbStarts(sl~age2,data=mqf9,param="Pauly",fixed=list(ts=0.85,NGT=0.05))
Pfitmqf9 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf9,
                start=Psvmqf9,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf9 <- nlsBoot(Pfitmqf9)
Pcfmqf9 <- cbind(Est=coef(Pfitmqf9),confint(Pbootmqf9))

## 3. Summary results
Scfmqf9 <- rbind(Scfmqf9,c(sum(summary(Sfitmqf9)$residuals^2),NA,NA),
                 c(AIC(Sfitmqf9),NA,NA))
Pcfmqf9 <- rbind(Pcfmqf9,c(sum(summary(Pfitmqf9)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf9),NA,NA))
rownames(Scfmqf9)[6:7] <- rownames(Pcfmqf9)[6:7] <- c("RSS","AIC")
print(round(Scfmqf9,2),na.print="-")
print(round(Pcfmqf9,2),na.print="-")


################################################################################
## Summary graphic
################################################################################
pdf("results/Figure_3.PDF",width=8,height=8)
par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1)
plot(fl~age,data=Bonito,pch=19,col=col2rgbt("black",1/4),
     xlab="Age (years)",ylab="Fork Length (mm)",xlim=c(0,4),ylim=c(0,70))
curve(vbSCGF(x,coef(PfitBon)),from=0,to=4,lwd=4,add=TRUE)
curve(vbSO(x,coef(SfitBon)),from=0,to=4,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"A",cex=1.5)

plot(sl~age2,data=mqf2,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbSCGF(x,coef(Pfitmqf2)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf2)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"B",cex=1.5)

plot(sl~age2,data=mqf4,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbSCGF(x,coef(Pfitmqf4)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf4)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"C",cex=1.5)

plot(sl~age2,data=mqf9,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbSCGF(x,coef(Pfitmqf9)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf9)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"D",cex=1.5)
dev.off()
