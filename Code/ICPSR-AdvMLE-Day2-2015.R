#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day two materials. 
#
# Note: Please see the comments and recommendations on
# the last slide from class.
########################################################
# KABL data & examples

library(RCurl)
library(foreign)
library(survival)
library(eha)
library(rms)
library(flexsurv)

KABLURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2015-git/master/Data/KABL.csv"
temp<-getURL(KABLURL)
KABL<-read.csv(textConnection(temp))

KABL.S<-Surv(KABL$durat,KABL$ciep12)

# KM plot:

KABL.KapM<-survfit(KABL.S~1)

pdf("KABLKapM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.KapM,lwd=c(3,1,1),mark.time=FALSE,
     ylab="Survival Probability",xlab="Time (in months)")
dev.off()

# Covariates:

xvars<-c("fract","polar","format","invest","numst2","eltime2","caretk2")
MODEL<-as.formula(paste(paste("KABL.S ~ ", paste(xvars,collapse="+"))))

# Exponential:

KABL.exp.AFT<-survreg(MODEL,data=KABL,dist="exponential")

KABL.exp.PH<-(-KABL.exp.AFT$coefficients)

KABL.exp.HRs<-exp(-KABL.exp.AFT$coefficients)

# Comparing predicted survival curves:
#
# Refit model using flexsurvreg...

KABL.exp<-flexsurvreg(MODEL,data=KABL,dist="exp")
FakeInvest<-t(c(mean(KABL$fract),mean(KABL$polar),mean(KABL$format),1,
            mean(KABL$numst2),mean(KABL$eltime2),mean(KABL$caretk2)))
FakeNoInvest<-t(c(mean(KABL$fract),mean(KABL$polar),mean(KABL$format),0,
                mean(KABL$numst2),mean(KABL$eltime2),mean(KABL$caretk2)))
colnames(FakeInvest)<-xvars
colnames(FakeNoInvest)<-xvars

pdf("ExpSurvCompare.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.exp,FakeInvest,mark.time=FALSE,col.obs="black",
     lty.obs=c(0,0,0),xlab="Time (in months)",ylab="Survival Probability")
lines(KABL.exp,FakeNoInvest,mark.time=FALSE,col.obs="black",
      lty.obs=c(0,0,0),col=c(rep("green",times=3)))
legend("topright",inset=0.05,bty="n",
       c("Investiture Requirement","No Investiture Requirement"),
       lty=c(1,1),lwd=c(2,2),col=c("red","green"))
dev.off()


# Weibulls...
#
# Plot of various hazard shapes:

t<-cbind(1:60,1:60,1:60)
P<-c(0.5,1,2)
WeibullHs<-t(apply(t,1,function(t) 0.02*P*((0.02*t)^(P-1))))
WeibullSs<-t(apply(t,1,function(t) (exp(-0.02*t))^P))

pdf("WeibullHSims.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],WeibullHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Hazard",ylim=c(0,0.08))
lines(t[,2],WeibullHs[,2],t="l",lwd=3,lty=2,col="black")
lines(t[,3],WeibullHs[,3],t="l",lwd=3,lty=3,col="red")
legend("topright",inset=.02,
       c("p = 0.5","p = 1.0","p = 2.0"),
       lty=c(1,2,3),lwd=c(3,3,3),col=c("green","black","red"),
       cex=1.2,bty="n")
dev.off()

pdf("WeibullSSims.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],WeibullSs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Survival Probability",ylim=c(0,1))
lines(t[,2],WeibullSs[,2],t="l",lwd=3,lty=2,col="black")
lines(t[,3],WeibullSs[,3],t="l",lwd=3,lty=3,col="red")
legend("bottomleft",inset=.02,
       c("p = 0.5","p = 1.0","p = 2.0"),
       lty=c(1,2,3),lwd=c(3,3,3),col=c("green","black","red"),
       cex=1.2,bty="n")
dev.off()

# Weibull KABL:

KABL.weib.AFT<-survreg(MODEL,data=KABL,dist="weibull")

KABL.weib.PH<-(-KABL.weib.AFT$coefficients)/(KABL.weib.AFT$scale) 

KABL.weib.HRs<-exp(KABL.weib.PH)

# Comparing Weibull survival curves

KABL.weib.Ihat<-predict(KABL.weib.AFT,newdata=as.data.frame(FakeInvest),
                type="quantile",se.fit=TRUE,p=seq(.01,.99,by=.01))

KABL.weib.NoIhat<-predict(KABL.weib.AFT,newdata=as.data.frame(FakeNoInvest),
                type="quantile",se.fit=TRUE,p=seq(.01,.99,by=.01))

pdf("WeibSurvCompare.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.weib.NoIhat$fit,seq(.99,.01,by=-.01),t="l",lwd=3,col="green",
     xlab="Time (in months)",ylab="Survival Probability")
lines(KABL.weib.Ihat$fit,seq(.99,.01,by=-.01),lwd=3,col="red")
lines(KABL.weib.NoIhat$fit+1.96*(KABL.weib.NoIhat$se),
      seq(.99,.01,by=-.01),lty=2,lwd=1,col="green")
lines(KABL.weib.NoIhat$fit-1.96*(KABL.weib.NoIhat$se),
      seq(.99,.01,by=-.01),lty=2,lwd=1,col="green")
lines(KABL.weib.Ihat$fit+1.96*(KABL.weib.NoIhat$se),
      seq(.99,.01,by=-.01),lty=2,lwd=1,col="red")
lines(KABL.weib.Ihat$fit-1.96*(KABL.weib.NoIhat$se),
      seq(.99,.01,by=-.01),lty=2,lwd=1,col="red")
legend("topright",inset=0.05,bty="n",
       c("Investiture Requirement","No Investiture Requirement"),
       lty=c(1,1),lwd=c(2,2),col=c("red","green"))
dev.off()

# Gompertz sims:

t<-cbind(1:1000,1:1000,1:1000)
t<-t/100
G<-c(-0.5,0,0.1)
GompertzHs<-t(apply(t,1,function(t) exp(0.2)*exp(G*t)))

pdf("GompertzHSims.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],GompertzHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Hazard",ylim=c(0,4))
lines(t[,2],GompertzHs[,2],lwd=3,lty=2,col="black")
lines(t[,3],GompertzHs[,3],lwd=3,lty=3,col="red")
legend("topleft",inset=.02,
       c("gamma = -0.5","gamma = 0","gamma = 0.2"),
       lty=c(1,2,3),lwd=c(3,3,3),col=c("green","black","red"),
       cex=1.2,bty="n")
dev.off()

# KABL Gompertz:

library(flexsurv)

KABL.Gomp<-flexsurvreg(MODEL,data=KABL,dist="gompertz")

# Log-logistic sims:

t<-cbind(1:100,1:100,1:100)
t<-t/2
p<-c(0.5,1.0,3.0)
LogLogHs<-t(apply(t,1,function(t) (0.05*p*((0.05*t)^(p-1)))/(1+(0.05*t)^p)))

pdf("LogLogHSims.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],LogLogHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Hazard",ylim=c(0,0.15))
lines(t[,2],LogLogHs[,2],lwd=3,lty=2,col="black")
lines(t[,3],LogLogHs[,3],lwd=3,lty=3,col="red")
legend("topright",inset=.02,
       c("p = 0.5","p = 1.0","p = 3.0"),
       lty=c(1,2,3),lwd=c(3,3,3),col=c("green","black","red"),
       cex=1.2,bty="n")
dev.off()

# KABL Log-Logistic:

KABL.loglog<-survreg(MODEL,data=KABL,dist="loglogistic")

# KABL Log-Normal:

KABL.logN<-survreg(MODEL,data=KABL, dist="lognormal")
