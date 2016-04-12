################################################################################
## This is an attempt at a function to create the Pauly et al. (1992)
## seasonal cessation Von B growth model
##
##   Linf, K, t0 as usual
##   WP = "Winter Period" (point where growth=0) = ts+0.5 (ts from Pauly et al.)
##   NGT = "No Growth Time" = "fraction of a year where no growth occurs
##
##   tpr = "t-prime" = actual age (t) minus cumulative NGT prior to t
##   Q is as defined in Pauly et al.
##   qt and qtr are intermediate values to make Pauly SC look similar to the
##     Somers and Somers2 models from vbFuns() in FSA
##
##   The final line is basically Equation 4 from Pauly et al.
################################################################################
paulySC <- function(t,Linf,K,t0,WP,NGT) {
  ## Shift time so that first NGT is at 1
  # NGT starts at WP
  shift <- 1-WP
  tmp.t <- t+shift
  ## Find fraction of year on this new time scale
  tmp.t2 <- tmp.t-floor(tmp.t)
  ## Adjust this fraction for no growth
  for (i in 1:length(tmp.t2)) {
    if (tmp.t2[i]<=NGT) tmp.t2[i] <- 0
    else tmp.t2[i] <- tmp.t2[i]-NGT
  }
  tmp.t <- floor(tmp.t)*(1-NGT)+tmp.t2
  ## Shift back to get tprime
  tpr <- tmp.t-shift+NGT/2  # why is NGT/2 needed here
  ## return model values
  Q <- (2*pi)/(1-NGT)
  qt <- (K/Q)*sin(Q*(tpr-WP+0.5))
  qt0 <- (K/Q)*sin(Q*(t0-WP+0.5))
  Linf*(1-exp(-K*(tpr-t0)-qt+qt0))
}


################################################################################
## Some Checks?
################################################################################
## Does this basically reproduce Figure 1.4 in Magnifico ... IT DOES NOT
par(xaxs="i",yaxs="i")
t <- seq(-1,5,0.01)
WP <- 0.2; NGT=0.5
L <- paulySC(t,Linf=100,K=1,t0=-0.2,WP=WP,NGT=NGT)
plot(L~t,type="l",lwd=2,ylim=c(0,100))
abline(v=(0:5)+WP,lty=3,col="red")
abline(v=(0:5)+WP+NGT,lty=3,col="blue")
grid()


## Does this plot make sense ... should cease growth in late fall for half year
##   then commences a half year later
t <- seq(0,5,0.01)
WP <- 0.8; NGT=0.5
L <- paulySC(t,Linf=100,K=1,t0=-0.2,WP=WP,NGT=NGT)
plot(L~t,type="l",lwd=2,ylim=c(0,100))
abline(v=(0:5)+WP,lty=3,col="red")
abline(v=(0:5)+WP+NGT,lty=3,col="blue")
grid()


## Do the model fits match those in Pauly's Figure 2 ... NO, but they seem reasonable
##    could be that the data are not actual .. i.e., I guessed at the data??
setwd("C:/aaaWork/Consulting/R_MiscOther/Newton_Erica_PaulyVBGF")
f2 <- read.csv("Fig2.csv")
plot(len~age,data=f2,pch=19,ylim=c(0,16),xlim=c(0,2))
ages2 <- seq(0,2,0.01)
predL2 <- paulySC(ages2,Linf=18,K=1.07,t0=-0.03,WP=0.71,NGT=0.05)
lines(predL2~ages2,col="gray70")

fit2 <- nls(len~paulySC(age,Linf,K,t0,WP,NGT),data=f2,
            start=list(Linf=18,K=1.07,t0=-0.03,WP=0.71,NGT=0.05))
predL2f <- predict(fit2,data.frame(age=ages2))
lines(predL2f~ages2,col="red")
grid()
coef(fit2)
#confint(fit2)


## Do the model fits match those in Pauly's Figure 3 ... NO, but they seem reasonable
##    could be that the data are not actual??
f3 <- read.csv("Fig3.csv")
plot(len~age,data=f3,pch=19,ylim=c(0,120),xlim=c(0,2.5))
ages3 <- seq(0,2.5,0.01)
predL3 <- paulySC(ages3,Linf=156,K=0.76,t0=0.10,WP=0.66,NGT=0.27)
lines(predL3~ages3,col="gray70")

fit3 <- nls(len~paulySC(age,Linf,K,t0,WP,NGT),data=f3,
            start=list(Linf=156,K=0.76,t0=0.10,WP=0.66,NGT=0.27))
predL3f <- predict(fit3,data.frame(age=ages3))
lines(predL3f~ages3,col="red")
grid()
cbind(coef(fit3),confint(fit3))



################################################################################
# Fit some actual data
################################################################################
# Make Somers2 model
library(FSA)
( vbS <- vbFuns("Somers2") )

library(FSAdata)

### Araucanian Herring from Chile
data(AHerringChile)
## Fit Somers model
# find starting values
plot(len~age,data=AHerringChile,pch=19,col=col2rgbt("black",1/4))
#curve(vbS(x,16.3,0.5,-1,0.9,0.9),from=0.2,to=3.1,add=TRUE,col="red")
H.st <- list(Linf=16,K=0.5,t0=-1,C=0.9,WP=0.9)
fitH.s <- nls(len~vbS(age,Linf,K,t0,C,WP),data=AHerringChile,start=H.st) 
curve(vbS(x,coef(fitH.s)),from=0.2,to=3.1,add=TRUE,col="blue",lwd=2)
coef(fitH.s)
# Fit Pauly SC model
fitH.p <- nls(len~paulySC(age,Linf,K,t0,WP,NGT),data=AHerringChile,
              start=list(Linf=16,K=0.5,t0=-1,WP=0.9,NGT=0.01))
curve(paulySC(x,coef(fitH.p)[1],coef(fitH.p)[2],coef(fitH.p)[3],coef(fitH.p)[4],
              coef(fitH.p)[5]),from=0.2,to=3.1,add=TRUE,col="green",lwd=2)
coef(fitH.p)
##!!! Notice slight parameter differences but essentially the same model fit!!


### Anchoveta from Chile
data(AnchovetaChile)
AnchovetaChile$age <- AnchovetaChile$age.mon/12
## Fit Somers model
# find starting values
plot(tl.cm~age,data=AnchovetaChile,pch=19,col=col2rgbt("black",1/4))
#curve(vbS(x,20,0.6,0.2,0.9,0.9),from=0.2,to=4.5,add=TRUE,col="red")
A.st <- list(Linf=20,K=0.6,t0=0.2,C=0.9,WP=0.9)
fitA.s <- nls(tl.cm~vbS(age,Linf,K,t0,C,WP),data=AnchovetaChile,start=A.st) 
curve(vbS(x,coef(fitA.s)),from=0.2,to=4.5,add=TRUE,col="blue",lwd=2)
coef(fitA.s)
# Fit Pauly SC model
fitA.p <- nls(tl.cm~paulySC(age,Linf,K,t0,WP,NGT),data=AnchovetaChile,
              start=list(Linf=16,K=0.5,t0=-1,WP=0.9,NGT=0.01))
curve(paulySC(x,coef(fitA.p)[1],coef(fitA.p)[2],coef(fitA.p)[3],coef(fitA.p)[4],
              coef(fitA.p)[5]),from=0.2,to=4.5,add=TRUE,col="green",lwd=2)
coef(fitA.p)
##!!! Notice slight parameter differences but essentially the same model fit!!


################################################################################
# TO DO
#
# 1. Hammer out the mathematics for t-prime.
# 2. Provide a better derivation of the model then what is in Pauly et al.
# 3. Why do Somers and PaulySC give basically the same fit.
# 4. Simulate some data to see about parameter interpretations.  Do Somers and
#    PaulySC fit the same if NGT is longer.
# 5. Clean up the PaulySC function (put into FSA if it seems correct).
# 6. Find some real data that is a more interesting example.
# 7. Write a note????