#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day three materials.
#
########################################################
# ONeal/Russett data & examples

library(RCurl)
library(foreign)
library(survival)

ORURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2015-git/master/Data/OR.csv"
temp<-getURL(ORURL)
OR<-read.csv(textConnection(temp))

summary(OR)

# Surv object...

OR.S<-Surv(OR$start,OR$stop,OR$dispute,type=c('counting'))
OR.KM<-survfit(OR.S~1)

pdf("ORKM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(OR.KM,mark.time=FALSE,lwd=c(2,1,1),
     xlab="Time (in years)",ylab="Survival Probability")
dev.off()

# Cox model w/OR data:

ORCox.br<-coxph(OR.S~allies+contig+capratio+growth+democracy+trade,
                data=OR,na.action=na.exclude, method="breslow")

# Scaling covariates:

OR$growthPct<-OR$growth*100
summary(coxph(OR.S~allies+contig+capratio+growthPct+democracy+trade,
                data=OR,na.action=na.exclude, method="breslow"))

# Baseline (cumulative) hazard:

OR.BH<-basehaz(ORCox.br,centered=FALSE)

pdf("ORCoxBaseH.pdf",10,5)
par(mar=c(4,4,2,2))
plot(OR.BH$time,OR.BH$hazard,t="l",lwd=4,col="red",
     xlab="Time (in years)",ylab="Baseline Integrated Hazard")
lines(abline(lm(OR.BH$hazard~0+OR.BH$time),lty=2,lwd=2))
legend("bottomright",inset=0.02,bty="n",
       c("Baseline Hazard","Linear Fit"),lty=c(1,2),
       lwd=c(4,2),col=c("red","black"))
dev.off()

# Comparing survival curves:

FakeContig<-as.data.frame(t(c(mean(OR$allies),1,mean(OR$capratio),mean(OR$growth),
              mean(OR$democracy),mean(OR$trade))))
FakeApart<-as.data.frame(t(c(mean(OR$allies),0,mean(OR$capratio),mean(OR$growth),
              mean(OR$democracy),mean(OR$trade))))
colnames(FakeContig)<-c("allies","contig","capratio","growth",
                        "democracy","trade")
colnames(FakeApart)<-c("allies","contig","capratio","growth",
                        "democracy","trade")

FCHat<-survfit(ORCox.br,FakeContig)
FAHat<-survfit(ORCox.br,FakeApart)

pdf("SurvCompare.pdf",6,5)
par(mar=c(4,4,2,2))
plot(FAHat,col=c(rep("black",times=3)),lwd=c(3,1,1),lty=c(1,2,2),
     xlab="Time (in years)",ylab="Survival Probability",
     mark.time=FALSE)
par(new=TRUE)
plot(FCHat,col=c(rep("red",times=3)),lwd=c(3,1,1),lty=c(2,2,2),
     mark.time=FALSE)
legend("bottomleft",inset=0.02,bty="n",
       c("Non-Contiguous","Contiguous"),lty=c(1,2),
       lwd=c(3,3),col=c("black","red"))
dev.off()

# Ties:

set.seed(7222009)
Data<-as.data.frame(cbind(c(rep(1,times=400)),
            c(rep(c(0,1),times=200))))
colnames(Data)<-c("C","X")
Data$T<-rexp(400,exp(0+1*Data$X)) # B = 1.0
Data.S<-Surv(Data$T,Data$C)
plot(survfit(Data.S~Data$X),col=c("black","red"),
     xlab="Time",ylab="Survival") # plot

D.br<-coxph(Data.S~X,data=Data,method="breslow")
D.ef<-coxph(Data.S~X,data=Data,method="efron")
D.ex<-coxph(Data.S~X,data=Data,method="exact")

D.Bs<-c(D.br$coefficients,D.ef$coefficients,D.ex$coefficients)
Dlab<-c("Breslow","Efron","Exact")

# Plot:

pdf("SimNoTies.pdf",12,5)
par(mar=c(4,4,2,2))
dotchart(D.Bs,labels=Dlab,pch=19,cex=1.8,xlab="Estimated Beta")
dev.off()

# Now add ties via rounding (nearest multiple of 5):

Data$Tied<-round(Data$T,0)

DataT.S<-Surv(Data$Tied,Data$C)

DT.br<-coxph(DataT.S~X,data=Data,method="breslow")
DT.ef<-coxph(DataT.S~X,data=Data,method="efron")
DT.ex<-coxph(DataT.S~X,data=Data,method="exact")
DT.Bs<-c(DT.br$coefficients,DT.ef$coefficients,DT.ex$coefficients)

#Plot:

pdf("SimTies.pdf",12,5)
par(mar=c(4,4,2,2))
dotchart(DT.Bs,labels=Dlab,pch=19,xlab="Estimated Beta",
         cex=1.8)
abline(v=D.ex$coefficients,lty=2,lwd=2)
dev.off()
