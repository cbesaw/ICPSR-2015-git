#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day One materials.
#
########################################################
# King et al. (1990) data examples

library(RCurl)
library(foreign)
library(survival)

KABLURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2015-git/master/Data/KABL.csv"
temp<-getURL(KABLURL)
KABL<-read.csv(textConnection(temp))

KABL.S<-Surv(KABL$durat,KABL$ciep12)

# Survival curve estimates:

KABL.fit<-survfit(KABL.S~1)

pdf("KABLKM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.fit,xlab="Time (in months)",ylab="Survival Estimate",
     lwd=c(3,1,1))
dev.off()

# Estimates of H(t):

KABL.FH<-survfit(KABL.S~1,type="fleming")

pdf("KABLNA.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.FH$time,-log(KABL.FH$surv),t="l",lwd=3,lty=1,
     xlab="Time (in months)",ylab="Cumulative Hazard",
     ylim=c(0,2.5))
lines(KABL.FH$time,-log(KABL.FH$upper),t="l",lwd=2,lty=2)
lines(KABL.FH$time,-log(KABL.FH$lower),t="l",lwd=2,lty=2)
points(KABL.FH$time[KABL.FH$n.censor>0],
       -log(KABL.FH$surv[KABL.FH$n.censor>0]),pch=19)
dev.off()

# Alternatively, w/no CIs...:
#
# NAalen<-cumsum(KABL.fit$n.event / KABL.fit$n.risk) 
# Ht<- -log(KABL.fit$surv) # Alternative estimate
# 
# pdf("KABLNA.pdf",6,5)
# par(mar=c(4,4,2,2))
# plot(NAalen~KABL.fit$time,t="l",lwd=3,xlab="Time (in months)",
#      ylab="Cumulative Hazard Estimate")
# points(KABL.fit$time[KABL.fit$n.censor>0],
#        NAalen[KABL.fit$n.censor>0],pch=4)
# dev.off()

# Log-rank test:

survdiff(KABL.S~invest,data=KABL,rho=0)

# Survival curve comparison:

KABL.fit2<-survfit(KABL.S~invest,data=KABL)

pdf("KABL-By-Invest.pdf",6,5)
par(mar=c(4,4,2,2))
plot(KABL.fit2,conf.int=TRUE,lty=c(1,2),lwd=c(3,1,1,3,1,1),
     col=c("black","red"),ylab="Survival Estimate",
     xlab="Time (in months)", mark.time=FALSE)
legend("topright",inset=.02,
       c("No Investiture Requirement","Investiture Requirement"),
       lty=c(1,2),lwd=c(3,3),col=c("black","red"),bty="n")
dev.off()
