################################################################################
################################################################################
##
## Analysis script for ... Ogle, DH.  201X.  An Algorithm for the von Bertalanffy
##   Seasonal Cessation in Growth Function of Pauly et al. (1992).  Submitted to
##   Fisheries Research.
##
## Need to be patient with bootstrapping functions.  May also need to create a
##   directory called "results" in your current working directory to hold the
##   figures produced by the pdf() function (or not run the pdf() functions to
##   simple produce the figures on a local device.)
##
################################################################################
################################################################################

################################################################################
## SETUP
## 
## Requires "latest" versions of FSA (>=0.8.8) and FSAdata (>=0.3.3) from CRAN.
##   These can be installed as follows:
##
##   install.packages("FSA")
##   install.packages("FSAdata")
##
################################################################################
## Load required packages
library(FSAdata)   # for Bonito and Mosquitofish data
library(FSA)       # for Somers and Pauly function functions
library(nlstools)  # for nls model bootstrapping

## Create a Typical function
vbTyp <- vbFuns("Typical")

## Create a Somers function
vbSO <- vbFuns("Somers")

## Create the Pauly et al. function
( vbSCGF <- vbFuns("Pauly") )

## Note that vbSCGF uses an internal function for computing the t-prime values.
##  The next line displays this internal function (Step comments correspond to
##  the steps outlined in the manuscript).
FSA:::iCalc_tpr

## Set the random number seed so that the bootstraps stay reproducible
set.seed(730987)

################################################################################
## Create the demo figure of the Somers model
################################################################################
pdf("results/Figure_1.PDF",width=4,height=4)
par(xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1,cex=0.9)
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
par(xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1,cex=0.9)
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

## 1. Fit typical model
Tlwrbnd <- c(Linf=0,K=0,t0=-Inf)
Tuprbnd <- c(Linf=Inf,K=Inf,t0=Inf)
TsvBon <- vbStarts(fl~age,data=Bonito,param="Typical")
TfitBon <- nls(fl~vbTyp(age,Linf,K,t0),data=Bonito,
               start=TsvBon,lower=Tlwrbnd,upper=Tuprbnd,
               algorithm="port",control=list(maxiter=100))
TbootBon <- nlsBoot(TfitBon)
TcfBon <- cbind(Est=coef(TfitBon),confint(TbootBon))

## 2. Fit Somers function with C<=1 constraint (as in Stewart et al. (2013))
##    Linf=71.9, K=0.27, t0=-1.92, C=1, and ts=0.09 ...
##    they all match (within rounding)
Slwrbnd <- c(Linf=0,K=0,t0=-Inf,C=0,ts=0)
SuprbndB <- c(Linf=Inf,K=Inf,t0=Inf,C=1,ts=1)
SsvBon <- vbStarts(fl~age,data=Bonito,param="Somers",fixed=list(C=0.6,ts=0.2))
SfitBon <- nls(fl~vbSO(age,Linf,K,t0,C,ts),data=Bonito,
               start=SsvBon,lower=Slwrbnd,upper=SuprbndB,
               algorithm="port",control=list(maxiter=100))
SbootBon <- nlsBoot(SfitBon)
ScfBon <- cbind(Est=coef(SfitBon),confint(SbootBon))

## 3. Fit new Pauly et al. (1992) function
Plwrbnd <- c(Linf=0,Kpr=0,t0=-Inf,ts=0,NGT=0)
Puprbnd <- c(Linf=Inf,Kpr=Inf,t0=Inf,ts=1,NGT=1)
PsvBon <- vbStarts(fl~age,data=Bonito,param="Pauly",fixed=list(ts=0.25,NGT=0.2))
PfitBon <- nls(fl~vbSCGF(age,Linf,Kpr,t0,ts,NGT),data=Bonito,
             start=PsvBon,lower=Plwrbnd,upper=Puprbnd,
             algorithm="port",control=list(maxiter=100))
PbootBon <- nlsBoot(PfitBon)
PcfBon <- cbind(Est=coef(PfitBon),confint(PbootBon))

## 4. Summary results
TcfBon <- rbind(TcfBon,c(sum(summary(TfitBon)$residuals^2),NA,NA),
                c(AIC(TfitBon),NA,NA))
ScfBon <- rbind(ScfBon,c(sum(summary(SfitBon)$residuals^2),NA,NA),
                c(AIC(SfitBon),NA,NA))
PcfBon <- rbind(PcfBon,c(sum(summary(PfitBon)$residuals^2),NA,NA),
                c(AIC(PfitBon),NA,NA))
rownames(TcfBon)[4:5] <- rownames(ScfBon)[6:7] <- rownames(PcfBon)[6:7] <- c("RSS","AIC")
print(round(TcfBon,2),na.print="-")
print(round(ScfBon,2),na.print="-")
print(round(PcfBon,2),na.print="-")


################################################################################
## Example analysis with the Mosquitofish data -- Site 2
################################################################################
data(Mosquitofish)
mqf2 <- subset(Mosquitofish,sitenum==2)

## 1. Fit typical model
Tsvmqf2 <- vbStarts(sl~age2,data=mqf2,param="Typical")
Tfitmqf2 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf2,
               start=Tsvmqf2,lower=Tlwrbnd,upper=Tuprbnd,
               algorithm="port",control=list(maxiter=100))
Tbootmqf2 <- nlsBoot(Tfitmqf2)
Tcfmqf2 <- cbind(Est=coef(Tfitmqf2),confint(Tbootmqf2))

## 2. Fit Somers function with C unconstrained (as in Carmona-Catot et al. (2014))
##    Linf=35.85, K=2.012, t0=-0.02, C=1.95, and ts=-0.118 ...
##    they all match (within rounding) except ts but it is off by +1
SuprbndM <- c(Linf=Inf,K=Inf,t0=Inf,C=Inf,ts=1)
Ssvmqf2 <- vbStarts(sl~age2,data=mqf2,param="Somers",fixed=list(C=1.5,ts=0.9))
Sfitmqf2 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf2,
                start=Ssvmqf2,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=100))
Sbootmqf2 <- nlsBoot(Sfitmqf2)
Scfmqf2 <- cbind(Est=coef(Sfitmqf2),confint(Sbootmqf2))

## 3. Fit new Pauly et al. (1992) function
Psvmqf2 <- vbStarts(sl~age2,data=mqf2,param="Pauly",fixed=list(ts=0.9,NGT=0.5))
Pfitmqf2 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf2,
                start=Psvmqf2,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf2 <- nlsBoot(Pfitmqf2)
Pcfmqf2 <- cbind(Est=coef(Pfitmqf2),confint(Pbootmqf2))

## 4. Summary results
Tcfmqf2 <- rbind(Tcfmqf2,c(sum(summary(Tfitmqf2)$residuals^2),NA,NA),
                c(AIC(Tfitmqf2),NA,NA))
Scfmqf2 <- rbind(Scfmqf2,c(sum(summary(Sfitmqf2)$residuals^2),NA,NA),
                c(AIC(Sfitmqf2),NA,NA))
Pcfmqf2 <- rbind(Pcfmqf2,c(sum(summary(Pfitmqf2)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf2),NA,NA))
rownames(Tcfmqf2)[4:5] <- rownames(Scfmqf2)[6:7] <- rownames(Pcfmqf2)[6:7] <- c("RSS","AIC")
print(round(Tcfmqf2,2),na.print="-")
print(round(Scfmqf2,2),na.print="-")
print(round(Pcfmqf2,2),na.print="-")


################################################################################
## Example analysis with the Mosquitofish data -- Site 4
################################################################################
mqf4 <- subset(Mosquitofish,sitenum==4)

## 1. Fit typical model
# Tsvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Typical",plot=TRUE) # did not work
Tsvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Typical",plot=TRUE,
                    fixed=list(Linf=45,K=1.2,t0=-0.2))
Tfitmqf4 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf4,
                start=Tsvmqf4,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
Tbootmqf4 <- nlsBoot(Tfitmqf4)
Tcfmqf4 <- cbind(Est=coef(Tfitmqf4),confint(Tbootmqf4))

## 2. Fit Somers function with C unconstrained (as in Carmona-Catot et al. (2014))
Ssvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Somers",fixed=list(K=0.6,C=1.6,ts=0.9),plot=TRUE)
Sfitmqf4 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf4,
                start=Ssvmqf4,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=1000))
Sbootmqf4 <- nlsBoot(Sfitmqf4)
Scfmqf4 <- cbind(Est=coef(Sfitmqf4),confint(Sbootmqf4))

## 3. Fit new Pauly et al. (1992) function
Psvmqf4 <- vbStarts(sl~age2,data=mqf4,param="Pauly",fixed=list(ts=0.9,NGT=0.5))
Pfitmqf4 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf4,
                start=Psvmqf4,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf4 <- nlsBoot(Pfitmqf4)
Pcfmqf4 <- cbind(Est=coef(Pfitmqf4),confint(Pbootmqf4))

## 4. Summary results
Tcfmqf4 <- rbind(Tcfmqf4,c(sum(summary(Tfitmqf4)$residuals^2),NA,NA),
                 c(AIC(Tfitmqf4),NA,NA))
Scfmqf4 <- rbind(Scfmqf4,c(sum(summary(Sfitmqf4)$residuals^2),NA,NA),
                 c(AIC(Sfitmqf4),NA,NA))
Pcfmqf4 <- rbind(Pcfmqf4,c(sum(summary(Pfitmqf4)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf4),NA,NA))
rownames(Tcfmqf4)[4:5] <- rownames(Scfmqf4)[6:7] <- rownames(Pcfmqf4)[6:7] <- c("RSS","AIC")
print(round(Tcfmqf4,2),na.print="-")
print(round(Scfmqf4,2),na.print="-")
print(round(Pcfmqf4,2),na.print="-")


################################################################################
## Example analysis with the Mosquitofish data -- Site 9
################################################################################
mqf9 <- subset(Mosquitofish,sitenum==9)

## 1. Fit typical model
Tsvmqf9 <- vbStarts(sl~age2,data=mqf9,param="Typical")
Tfitmqf9 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf9,
                start=Tsvmqf9,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
Tbootmqf9 <- nlsBoot(Tfitmqf9)
Tcfmqf9 <- cbind(Est=coef(Tfitmqf9),confint(Tbootmqf9))

## 2. Fit Somers function with C unconstrained (as in Carmona-Catot et al. (2014))
Ssvmqf9 <- vbStarts(sl~age2,data=mqf9,param="Somers",fixed=list(t0=-2.5,C=0.6,ts=0.85))
Sfitmqf9 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf9,
                start=Ssvmqf9,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=100))
Sbootmqf9 <- nlsBoot(Sfitmqf9)
Scfmqf9 <- cbind(Est=coef(Sfitmqf9),confint(Sbootmqf9))

## 3. Fit new Pauly et al. (1992) function
Psvmqf9 <- vbStarts(sl~age2,data=mqf9,param="Pauly",fixed=list(ts=0.85,NGT=0.05))
Pfitmqf9 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf9,
                start=Psvmqf9,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Pbootmqf9 <- nlsBoot(Pfitmqf9)
Pcfmqf9 <- cbind(Est=coef(Pfitmqf9),confint(Pbootmqf9))

## 4. Summary results
Tcfmqf9 <- rbind(Tcfmqf9,c(sum(summary(Tfitmqf9)$residuals^2),NA,NA),
                 c(AIC(Tfitmqf9),NA,NA))
Scfmqf9 <- rbind(Scfmqf9,c(sum(summary(Sfitmqf9)$residuals^2),NA,NA),
                 c(AIC(Sfitmqf9),NA,NA))
Pcfmqf9 <- rbind(Pcfmqf9,c(sum(summary(Pfitmqf9)$residuals^2),NA,NA),
                 c(AIC(Pfitmqf9),NA,NA))
rownames(Tcfmqf4)[4:5] <- rownames(Scfmqf9)[6:7] <- rownames(Pcfmqf9)[6:7] <- c("RSS","AIC")
print(round(Tcfmqf9,2),na.print="-")
print(round(Scfmqf9,2),na.print="-")
print(round(Pcfmqf9,2),na.print="-")


################################################################################
## Summary graphic
################################################################################
pdf("results/Figure_3.PDF",width=8,height=8)
par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(3,3,0.6,0.6),mgp=c(1.7,.4,0),tcl=-0.2,las=1)
plot(fl~age,data=Bonito,pch=19,col=col2rgbt("black",1/4),
     xlab="Age (years)",ylab="Fork Length (mm)",xlim=c(0,4),ylim=c(0,70))
curve(vbTyp(x,coef(TfitBon)),from=0,to=4,lwd=4,add=TRUE,col="gray25",lty=2)
curve(vbSCGF(x,coef(PfitBon)),from=0,to=4,lwd=4,add=TRUE)
curve(vbSO(x,coef(SfitBon)),from=0,to=4,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"A",cex=1.5)

plot(sl~age2,data=mqf2,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbTyp(x,coef(Tfitmqf2)),from=0,to=3,lwd=4,add=TRUE,col="gray25",lty=2)
curve(vbSCGF(x,coef(Pfitmqf2)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf2)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"B",cex=1.5)

plot(sl~age2,data=mqf4,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbTyp(x,coef(Tfitmqf4)),from=0,to=3,lwd=4,add=TRUE,col="gray25",lty=2)
curve(vbSCGF(x,coef(Pfitmqf4)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf4)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"C",cex=1.5)

plot(sl~age2,data=mqf9,pch=19,col=col2rgbt("black",1/6),xaxt="n",
     xlab="Age (years)",ylab="Standard Length (mm)",xlim=c(0,2.4),ylim=c(0,50))
axis(1,0:2)
curve(vbTyp(x,coef(Tfitmqf9)),from=0,to=3,lwd=4,add=TRUE,col="gray25",lty=2)
curve(vbSCGF(x,coef(Pfitmqf9)),from=0,to=3,lwd=4,add=TRUE)
curve(vbSO(x,coef(Sfitmqf9)),from=0,to=3,lwd=2,add=TRUE,col="gray50")
text(grconvertX(0.08,"npc","user"),grconvertY(0.92,"npc","user"),"D",cex=1.5)
dev.off()


################################################################################
## Testing different starting values for model fits ... per reviewer request
##   Just checking for convergence and relationship to parameter estimates
##   from the fits above.
################################################################################

# Bonito
TsvBon1 <- list(Linf=60,K=0.3,t0=0)
TfitBon1 <- nls(fl~vbTyp(age,Linf,K,t0),data=Bonito,
                start=TsvBon1,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
TsvBon2 <- list(Linf=90,K=0.6,t0=0)
TfitBon2 <- nls(fl~vbTyp(age,Linf,K,t0),data=Bonito,
                start=TsvBon2,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
TsvBon3 <- list(Linf=80,K=0.1,t0=-2)
TfitBon3 <- nls(fl~vbTyp(age,Linf,K,t0),data=Bonito,
                start=TsvBon3,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
cbind(TsvBon,TsvBon1,TsvBon2,TsvBon3)  # The starting values
round(cbind(coef(TfitBon),coef(TfitBon1),coef(TfitBon2),
            coef(TfitBon2)),4)   # Parameter estimates look good

SsvBon1 <- list(Linf=60,K=0.3,t0=-2,C=0.3,ts=0.1)
SfitBon1 <- nls(fl~vbSO(age,Linf,K,t0,C,ts),data=Bonito,
                start=SsvBon1,lower=Slwrbnd,upper=SuprbndB,
                algorithm="port",control=list(maxiter=100))
SsvBon2 <- list(Linf=40,K=0.1,t0=0,C=0.7,ts=0.3)
SfitBon2 <- nls(fl~vbSO(age,Linf,K,t0,C,ts),data=Bonito,
                start=SsvBon2,lower=Slwrbnd,upper=SuprbndB,
                algorithm="port",control=list(maxiter=100))
SsvBon3 <- list(Linf=60,K=0.5,t0=0,C=0.5,ts=0.5)
SfitBon3 <- nls(fl~vbSO(age,Linf,K,t0,C,ts),data=Bonito,
                start=SsvBon3,lower=Slwrbnd,upper=SuprbndB,
                algorithm="port",control=list(maxiter=100))
cbind(SsvBon,SsvBon1,SsvBon2,SsvBon3)  # The starting values
round(cbind(coef(SfitBon),coef(SfitBon1),coef(SfitBon2),
            coef(SfitBon3)),4)         # Parameter estimates look good

PsvBon1 <- list(Linf=60,Kpr=0.2,t0=0,ts=0.1,NGT=0.1)
PfitBon1 <- nls(fl~vbSCGF(age,Linf,Kpr,t0,ts,NGT),data=Bonito,
                start=PsvBon1,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
PsvBon2 <- list(Linf=40,Kpr=0.7,t0=-1,ts=0.3,NGT=0.5)
PfitBon2 <- nls(fl~vbSCGF(age,Linf,Kpr,t0,ts,NGT),data=Bonito,
                start=PsvBon2,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
PsvBon3 <- list(Linf=70,Kpr=0.5,t0=-1,ts=0.5,NGT=0.5)
PfitBon3 <- nls(fl~vbSCGF(age,Linf,Kpr,t0,ts,NGT),data=Bonito,
                start=PsvBon3,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
cbind(PsvBon,PsvBon1,PsvBon2,PsvBon3)   # Starting values
round(cbind(coef(PfitBon),coef(PfitBon1),coef(PfitBon2),
            coef(PfitBon3)),4)          # Parameter estimates look good

# Mosquitofish Site 2
Tsvmqf21 <- list(Linf=60,K=0.3,t0=0)
Tfitmqf21 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf2,
                start=Tsvmqf21,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
Tsvmqf22 <- list(Linf=90,K=0.6,t0=0)
Tfitmqf22 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf2,
                start=Tsvmqf22,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
Tsvmqf23 <- list(Linf=80,K=0.1,t0=-2)
Tfitmqf23 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf2,
                start=Tsvmqf23,lower=Tlwrbnd,upper=Tuprbnd,
                algorithm="port",control=list(maxiter=100))
cbind(Tsvmqf2,Tsvmqf21,Tsvmqf22,Tsvmqf23)  # The starting values
round(cbind(coef(Tfitmqf2),coef(Tfitmqf21),coef(Tfitmqf22),
            coef(Tfitmqf22)),4)   # Parameter estimates look good

Ssvmqf21 <- list(Linf=60,K=0.3,t0=-2,C=1,ts=0.9)
Sfitmqf21 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf2,
                start=Ssvmqf21,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=100))
Ssvmqf22 <- list(Linf=40,K=0.1,t0=0,C=1.9,ts=0.7)
Sfitmqf22 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf2,
                start=Ssvmqf22,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=100))
Ssvmqf23 <- list(Linf=40,K=0.5,t0=0,C=1.4,ts=0.7)
Sfitmqf23 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf2,
                start=Ssvmqf23,lower=Slwrbnd,upper=SuprbndM,
                algorithm="port",control=list(maxiter=100))
cbind(Ssvmqf2,Ssvmqf21,Ssvmqf22,Ssvmqf23)  # The starting values
round(cbind(coef(Sfitmqf2),coef(Sfitmqf21),coef(Sfitmqf22),
            coef(Sfitmqf23)),4)         # Parameter estimates look good
# Some starting values did produce different estimates, but the fit was worse
# e.g.,  Ssvmqf23 <- list(Linf=60,K=0.5,t0=-2,C=1.3,ts=0.5)
#        deviance(Sfitmqf23)

Psvmqf21 <- list(Linf=60,Kpr=0.2,t0=0,ts=0.7,NGT=0.1)
Pfitmqf21 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf2,
                start=Psvmqf21,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Psvmqf22 <- list(Linf=40,Kpr=0.4,t0=-1,ts=0.7,NGT=0.5)
Pfitmqf22 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf2,
                start=Psvmqf22,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
Psvmqf23 <- list(Linf=60,Kpr=0.5,t0=0,ts=0.5,NGT=0.3)
Pfitmqf23 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf2,
                start=Psvmqf23,lower=Plwrbnd,upper=Puprbnd,
                algorithm="port",control=list(maxiter=100))
cbind(Psvmqf2,Psvmqf21,Psvmqf22,Psvmqf23)   # Starting values
round(cbind(coef(Pfitmqf2),coef(Pfitmqf21),coef(Pfitmqf22),
            coef(Pfitmqf23)),4)          # Parameter estimates look good

# Mosquitofish Site 4
Tsvmqf41 <- list(Linf=60,K=0.3,t0=0)
Tfitmqf41 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf4,
                 start=Tsvmqf4,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
Tsvmqf42 <- list(Linf=90,K=0.6,t0=0)
Tfitmqf42 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf4,
                 start=Tsvmqf42,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
Tsvmqf43 <- list(Linf=80,K=0.1,t0=-2)
Tfitmqf43 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf4,
                 start=Tsvmqf43,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
cbind(Tsvmqf4,Tsvmqf41,Tsvmqf42,Tsvmqf43)  # The starting values
round(cbind(coef(Tfitmqf4),coef(Tfitmqf41),coef(Tfitmqf42),
            coef(Tfitmqf42)),4)   # Parameter estimates look good (i.e., comparable)

Ssvmqf41 <- list(Linf=60,K=0.3,t0=0,C=1,ts=0.9)
Sfitmqf41 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf4,
                 start=Ssvmqf41,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
Ssvmqf42 <- list(Linf=60,K=0.9,t0=0,C=2,ts=0.8)
Sfitmqf42 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf4,
                 start=Ssvmqf42,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
Ssvmqf43 <- list(Linf=40,K=0.9,t0=0,C=1,ts=0.8)
Sfitmqf43 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf4,
                 start=Ssvmqf43,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
cbind(Ssvmqf4,Ssvmqf41,Ssvmqf42,Ssvmqf43)  # The starting values
round(cbind(coef(Sfitmqf4),coef(Sfitmqf41),coef(Sfitmqf42),
            coef(Sfitmqf43)),4)         # Parameter estimates look good
# Seems very sensitive to choice of ts (in terms of convergence)

Psvmqf41 <- list(Linf=60,Kpr=0.4,t0=-2,ts=0.7,NGT=0.5)
Pfitmqf41 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf4,
                 start=Psvmqf41,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
Psvmqf42 <- list(Linf=40,Kpr=0.3,t0=0,ts=0.9,NGT=0.7)
Pfitmqf42 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf4,
                 start=Psvmqf42,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
Psvmqf43 <- list(Linf=60,Kpr=0.2,t0=-2,ts=0.7,NGT=0.7)
Pfitmqf43 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf4,
                 start=Psvmqf43,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
cbind(Psvmqf4,Psvmqf41,Psvmqf42,Psvmqf43)   # Starting values
round(cbind(coef(Pfitmqf4),coef(Pfitmqf41),coef(Pfitmqf42),
            coef(Pfitmqf43)),4)          # Parameter estimates look good

# Mosquitofish Site 9
Tsvmqf91 <- list(Linf=60,K=0.3,t0=0)
Tfitmqf91 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf9,
                 start=Tsvmqf9,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
Tsvmqf92 <- list(Linf=90,K=0.6,t0=0)
Tfitmqf92 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf9,
                 start=Tsvmqf92,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
Tsvmqf93 <- list(Linf=80,K=0.1,t0=-2)
Tfitmqf93 <- nls(sl~vbTyp(age2,Linf,K,t0),data=mqf9,
                 start=Tsvmqf93,lower=Tlwrbnd,upper=Tuprbnd,
                 algorithm="port",control=list(maxiter=100))
cbind(Tsvmqf9,Tsvmqf91,Tsvmqf92,Tsvmqf93)  # The starting values
round(cbind(coef(Tfitmqf9),coef(Tfitmqf91),coef(Tfitmqf92),
            coef(Tfitmqf92)),4)   # Parameter estimates look good (i.e., comparable)

Ssvmqf91 <- list(Linf=60,K=0.3,t0=-1,C=0.4,ts=0.8)
Sfitmqf91 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf9,
                 start=Ssvmqf91,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
Ssvmqf92 <- list(Linf=30,K=0.1,t0=-1,C=0.4,ts=0.9)
Sfitmqf92 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf9,
                 start=Ssvmqf92,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
Ssvmqf93 <- list(Linf=60,K=0.2,t0=-1,C=0.2,ts=0.6)
Sfitmqf93 <- nls(sl~vbSO(age2,Linf,K,t0,C,ts),data=mqf9,
                 start=Ssvmqf93,lower=Slwrbnd,upper=SuprbndM,
                 algorithm="port",control=list(maxiter=100))
cbind(Ssvmqf9,Ssvmqf91,Ssvmqf92,Ssvmqf93)  # The starting values
round(cbind(coef(Sfitmqf9),coef(Sfitmqf91),coef(Sfitmqf92),
            coef(Sfitmqf93)),4)         # Parameter estimates look good
# Difficult time finding starting values that led to convergence

Psvmqf91 <- list(Linf=40,Kpr=0.1,t0=0,ts=0.65,NGT=0.25)
Pfitmqf91 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf9,
                 start=Psvmqf91,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
Psvmqf92 <- list(Linf=60,Kpr=0.3,t0=-1,ts=0.95,NGT=0.5)
Pfitmqf92 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf9,
                 start=Psvmqf92,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
Psvmqf93 <- list(Linf=40,Kpr=0.2,t0=-1,ts=0.5,NGT=0.5)
Pfitmqf93 <- nls(sl~vbSCGF(age2,Linf,Kpr,t0,ts,NGT),data=mqf9,
                 start=Psvmqf93,lower=Plwrbnd,upper=Puprbnd,
                 algorithm="port",control=list(maxiter=100))
cbind(Psvmqf9,Psvmqf91,Psvmqf92,Psvmqf93)   # Starting values
round(cbind(coef(Pfitmqf9),coef(Pfitmqf91),coef(Pfitmqf92),
            coef(Pfitmqf93)),4)          # Parameter estimates look good
# Some starting values did produce different estimates, but the fit was slightly worse
# e.g.,  Psvmqf94 <- list(Linf=60,Kpr=0.3,t0=0,ts=0.95,NGT=0.5)
#        deviance(Pfitmqf94)
