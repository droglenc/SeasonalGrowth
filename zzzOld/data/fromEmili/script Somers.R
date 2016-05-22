library(FSA);library(NCStats);library(nlstools);library(nls2)
library(car);library(qpcR)

setwd("D:/Emili2/R3 projectes/projectes MEC General/MEC 2006/Gambúsia latitud 2007/creix estacional/resultats diferenciant cohorts")
file9<-read.csv("fig3.csv",sep = ";")

file9$Females<-factor(file9$Females)
file9$Sitenum<-factor(file9$Sitenum)
sapply(file9, data.class)
names(file9)
summary(file9)
head(file9)

scatterplot(SL~ Age2  | Site, data=file9,smooth=T,reg.line=F,boxplots=F)

####################    VBGF   #########################


table(file9$Site,file9$Sitenum)

# Altea Bourdigou Ebro Fluvià Millars Orb Segura Ter Vistre Xuquer
# Select only 1; they are ordered by number/latitude
SELECTION="Vistre"
SELECTION="Orb"
SELECTION="Bourdigou"
SELECTION="Fluvià"
SELECTION="Ter"
SELECTION="Ebro"
SELECTION="Millars"
SELECTION="Xúquer"
SELECTION="Altea"
SELECTION="Segura"

file8=file9[file9$Site==SELECTION, ]

file8$len=file8$SL;file8$age=file8$Age2
summary(file8$len);summary(file8$age)
vbso <- vbFuns("Somers")

scatterplot(SL~ Age2, data=file8,smooth=T,reg.line=F,boxplots=F)

start <- data.frame(Linf = c(20, 50), K  = c(0, 1), t0 = c(-1, +1), C = c(0, 2), ts = c(0, +1))
vbso1 <- nls2(len ~ vbso(age, Linf, K, t0, C, ts), data = file8, start = start, trace=T)
summary(vbso1, correlation = TRUE)

# confint(vbso1)
Rsq(vbso1)
ages <- seq(0,3,0.01) # choose age range
plen <- predict(vbso1,data.frame(age=ages)) # predict length @ age
plot(len~age,data=file8,pch=19,xlab="Age (years)",ylab="Standard Length (mm)",main=SELECTION)
lines(ages,plen,lwd=2,col="red")

remove(vbso1)

