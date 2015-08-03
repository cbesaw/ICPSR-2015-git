#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day five materials.
#
########################################################

library(RCurl)
library(foreign)
library(gtools)
library(survival)
library(flexsurv)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

# Illustrating proportional hazards

T<-1:30
HE<-rep(0.02,times=30)
HR<-seq(from=0.01,to=0.03,length.out=30)
HD<-seq(from=0.025,to=0.005,length.out=30)

pdf("PropHazs.pdf",7,5)
par(mfrow=c(1,3))
plot(T,HE,t="l",lwd=3,lty=1,xlab="Time",ylab="Hazard",main="Constant Hazards",
     ylim=c(0,0.1))
lines(T,2.718*HE,lwd=3,lty=2)
legend("topleft",bty="n",c("h(t|X=0)","h(t|X=1)"),lwd=c(3,3),lty=c(1,2))
plot(T,HR,t="l",lwd=3,lty=1,xlab="Time",ylab="Hazard",main="Rising Hazards",
     ylim=c(0,0.1))
lines(T,2.718*HR,lwd=3,lty=2)
legend("topleft",bty="n",c("h(t|X=0)","h(t|X=1)"),lwd=c(3,3),lty=c(1,2))
plot(T,HD,t="l",lwd=3,lty=1,xlab="Time",ylab="Hazard",main="Declining Hazards",
     ylim=c(0,0.1))
lines(T,2.718*HD,lwd=3,lty=2)
legend("topleft",bty="n",c("h(t|X=0)","h(t|X=1)"),lwd=c(3,3),lty=c(1,2))
dev.off()

# Crossing hazards:

pdf("CrossHaz.pdf",6,5)
par(mar=c(4,4,2,2))
plot(T,2.718*HR,t="l",lwd=3,lty=1,xlab="Time",ylab="Hazard",
         ylim=c(0,0.1))
lines(T,2.718*HD,lwd=3,lty=2)
legend("topright",bty="n",c("h(t|X=0)","h(t|X=1)"),lwd=c(3,3),lty=c(1,2))
dev.off()



# Real example: SCOTUS

scotusURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2015-git/master/Data/scotus.csv"
temp<-getURL(scotusURL)
scotus<-read.csv(textConnection(temp))
rm(temp)

scotus.S<-Surv((scotus$service-1),scotus$service,scotus$retire)

pdf("scotusKM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(scotus.S~1),mark.time=F,lwd=c(3,1,1),
     xlab="Time (in years)",ylab="Survival")
dev.off()

scotus.Cox<-coxph(scotus.S~age+pension+pagree,data=scotus,ties="efron")

# log-log-S vs. log(T) plots:

pdf("loglogS.pdf",10,5)
par(mfrow=c(1,3))
plot(survfit(scotus.S~quantcut(age),data=scotus),lwd=c(2,2,2,2),
     fun="cloglog",mark.time=F,col=c("black","red","green","blue"),
     main="Age",xlab="log-Time",ylab="log(-log(S))")
legend("topleft",bty="n",c("Bottom Quartile","Second Quartile",
                           "Third Quartile","Top Quartile"),
       lwd=c(2,2,2,2),col=c("black","red","green","blue"))
plot(survfit(scotus.S~pension,data=scotus),
     fun="cloglog",mark.time=F,lwd=c(2,2),lty=c(1,2),col=c("blue","red"),
     main="Pension Eligibility",xlab="log-Time",ylab="log(-log(S))")
legend("topleft",bty="n",c("Not Pension Eligible","Pension Eligible"),
       lwd=c(2,2),lty=c(1,2),col=c("blue","red"))
plot(survfit(scotus.S~pagree,data=scotus),
     fun="cloglog",mark.time=F,lwd=c(2,2),lty=c(1,2),col=c("blue","red"),
     main="Party Agreement",xlab="log-Time",ylab="log(-log(S))")
legend("topleft",bty="n",c("No Party Agreement","Party Agreement"),
       lwd=c(2,2),lty=c(1,2),col=c("blue","red"))
dev.off()

# Residuals...
#
# Note: martingale residuals are also already in scotus.Cox$residuals...

scotus$mgres<-residuals(scotus.Cox,type="martingale")

print(scotus[scotus$justice==69,])
print(scotus[scotus$justice==49,])

scotus.schres<-residuals(scotus.Cox,type="schoenfeld")
scotus.scares<-residuals(scotus.Cox,type="scaledsch")

# Proportional hazards...

PHtest<-cox.zph(scotus.Cox)

pdf("scotusPHplots.pdf",10,5)
par(mar=c(4,4,2,2))
par(mfrow=c(1,3))
plot(PHtest,df=2)
dev.off()

# log-T Interactions

scotus$lnT<-log(scotus$service)
scotus$ageLnT<-scotus$age*(scotus$lnT)
scotus.NPH<-coxph(scotus.S~age+pension+pagree+ageLnT,
            data=scotus,ties="efron")
summary(scotus.NPH)

PHtest2<-cox.zph(scotus.NPH)

ageBs<-seq(from=scotus.NPH$coefficient[1]+(scotus.NPH$coefficient[4]*log(min(scotus$service))),
             to=scotus.NPH$coefficient[1]+(scotus.NPH$coefficient[4]*log(max(scotus$service))),
             length.out=(max(scotus$service)-min(scotus$service))+1)

pdf("scotusAgeNPHs.pdf",6,5)
par(mar=c(4,4,2,2))
scotusAgeNPHs<-plot(seq(min(scotus$service),max(scotus$service)),exp(ageBs),
                    t="l",lwd=3,lty=1,xlab="Court Tenure",
                    ylab="Hazard Ration for Age")
abline(h=1,lty=2)
dev.off()

# Unobserved heterogeneity plot:

pdf("UnobsHet.pdf",6,5)
par(mar=c(4,4,2,2))
plot(seq(0:20),rep(0.05,times=21),pch=NA,ylim=c(0.015,0.055),
     xlim=c(0,20),xlab="Time",ylab="Hazard")
abline(h=0.02,lwd=3,col="red")
abline(h=0.05,lwd=3,col="black")
abline(0.045,-0.001,lwd=3,lty=2,col="blue")
text(18,0.017,label=c("h(t) | Z=1"),col="red")
text(18,0.053,label=c("h(t) | Z=0"),col="black")
text(15.5,0.035,label=c("Estimated Hazard"),col="blue")
dev.off()

# Heterogeneity sim:

set.seed(7222009)
W<-rnorm(500)
X<-rnorm(500)
Z<-rnorm(500)
T<-rexp(500,rate=(exp(0+0.5*W+0.5*X-0.6*Z)))
C<-rep(1,times=500)
S<-Surv(T,C)

summary(survreg(S~W,dist="weibull"))
summary(survreg(S~W+X,dist="weibull"))
summary(survreg(S~W+X+Z,dist="weibull"))

M1<-survreg(S~W,dist="weibull")
M2<-survreg(S~W+X,dist="weibull")
M3<-survreg(S~W+X+Z,dist="weibull")

# Plot:

t<-cbind(1:60,1:60,1:60)
P<-c(1/(M1$scale),1/(M2$scale),1/(M3$scale))
DurDepHs<-t(apply(t,1,function(t) 0.02*P*((0.02*t)^(P-1))))

pdf("DurDepHs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],DurDepHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Weibull Hazard",ylim=c(0.01,0.04))
lines(t[,2],DurDepHs[,2],t="l",lwd=3,lty=2,col="blue")
lines(t[,3],DurDepHs[,3],t="l",lwd=3,lty=3,col="red")
abline(h=0.02,lty=4,lwd=2)
legend("topright",inset=.02,
       c("One Covariate","Two Covariates","Correct Specification","True Hazard"),
       lty=c(1,2,3,4),lwd=c(3,3,3,2),col=c("green","blue","red","black"),
       cex=1.2,bty="n")
dev.off()

# Parameterized duration dependence:

ct.weib<-flexsurvreg(scotus.S~age+pension+pagree,
                       data=scotus,dist="weibull")
ct.weib.DD<-flexsurvreg(scotus.S~age+pension+pagree+shape(age),
                         data=scotus,dist="weibull")

# Plots (this code is a hotter mess than usual, and could use
# a solid smack with a few Hadleyverse tools):

t<-1:max(scotus$service)

Age<-seq(min(scotus$age),max(scotus$age),by=1)
P.vary<-exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*Age))

pdf("PbyAge.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Age,P.vary,t="l",lwd=3,lty=1,ylab="Estimate of p")
abline(h=1,lty=2,lwd=2)
dev.off()            

age55<-55
age75<-75

P<-c(exp(scotus.weib$coefficients[1]),
  exp(ct.weib$coefficients[1]),
  exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age55)),
  exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age75)))
XB55<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age55)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB75<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age75)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB55DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age55)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB75DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age75)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB<-c(XB55,XB75,XB55DD,XB75DD)

h55<-dweibull(t,P[1],scale=XB[1]) / (1 - pweibull(t,P[1],scale=XB[1]))
h75<-dweibull(t,P[2],scale=XB[2]) / (1 - pweibull(t,P[2],scale=XB[2]))
h55DD<-dweibull(t,P[3],scale=XB[3]) / (1 - pweibull(t,P[3],scale=XB[3]))
h75DD<-dweibull(t,P[4],scale=XB[4]) / (1 - pweibull(t,P[4],scale=XB[4]))

pdf("DDHazards.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t,h75DD,t="l",lwd=3,lty=1,col="red",ylim=c(0.2,0.8),
     xlab="Time (in years)",ylab="Hazard")
lines(t,h55DD,lwd=3,lty=1,col="blue")
lines(t,h75,lwd=3,lty=2,col="red")
lines(t,h55,lwd=3,lty=2,col="blue")
legend("topleft",inset=0.02,
       c("Age 55 (p varying)","Age 75 (p varying)","Age 55 (p fixed)",
         "Age 75 (p fixed)"),lty=c(1,1,2,2),lwd=c(3,3,3,3),
       col=c("blue","red","blue","red"),bty="n")
dev.off()
