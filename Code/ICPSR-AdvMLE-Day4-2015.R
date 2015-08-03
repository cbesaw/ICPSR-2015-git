#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day four materials.
#
########################################################
# SCOTUS data & examples

library(RCurl)
library(foreign)
library(survival)

##############################
# ONeal-Russett data redux...

ORURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2015-git/master/Data/OR.csv"
temp<-getURL(ORURL)
OR<-read.csv(textConnection(temp))

# Logit models of disputes...

OR.logit<-glm(dispute~allies+contig+capratio+growth+democracy+trade,
              data=OR,na.action=na.exclude,family="binomial")

OR$duration<-OR$stop

OR.trend<-glm(dispute~allies+contig+capratio+growth+democracy+trade
              +duration,data=OR,na.action=na.exclude,family="binomial")

OR$d2<-OR$duration^2*0.1
OR$d3<-OR$duration^3*0.01
OR$d4<-OR$duration^4*0.001

OR.P4<-glm(dispute~allies+contig+capratio+growth+democracy+trade
              +duration+d2+d3+d4,data=OR,na.action=na.exclude,
              family="binomial")

P4test<-anova(OR.logit,OR.P4,test="Chisq")
P4test

OR.dummy<-glm(dispute~allies+contig+capratio+growth+democracy+trade
           +as.factor(duration),data=OR,na.action=na.exclude,
           family="binomial")

Test.Dummies<-anova(OR.logit,OR.dummy,test="Chisq")
Test.Dummies


# Predicted probabilities:

Xhats<-as.data.frame(t(c(mean(OR$allies),mean(OR$contig),mean(OR$capratio),
          mean(OR$growth),mean(OR$democracy),mean(OR$trade))))
Xhats<-Xhats[rep(1:nrow(Xhats),each=max(OR$duration)),]
Xhats$duration<-1:max(OR$duration)        
Xhats$d2<-Xhats$duration^2*0.1
Xhats$d3<-Xhats$duration^3*0.01
Xhats$d4<-Xhats$duration^4*0.001
colnames(Xhats)<-c("allies","contig","capratio","growth","democracy",
                 "trade","duration","d2","d3","d4")

Hat.logit<-predict(OR.logit,Xhats,type="response")
Hat.trend<-predict(OR.trend,Xhats,type="response")
Hat.P4<-predict(OR.P4,Xhats,type="response")
Hat.dummy<-predict(OR.dummy,Xhats,type="response")

pdf("DiscreteHats.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Xhats$duration,Hat.logit,ylim=c(0,0.04),t="l",lwd=3,col="black",
     xlab="Time (in years)",ylab="Predicted Probability")
lines(Xhats$duration,Hat.trend,lwd=3,lty=2,col="blue")
lines(Xhats$duration,Hat.P4,lwd=3,lty=3,col="green")
lines(Xhats$duration,Hat.dummy,lwd=3,lty=4,col="red")
legend("topright",inset=0.05,bty="n",
       c("No Duration Dependence","Linear Dependence",
         "Fourth-Degree Polynomial","Duration Dummies"),lty=c(1,2,3,4),
       lwd=c(3,3,3,3),col=c("black","blue","green","red"))
dev.off()

# Cox/Poisson equivalence:

OR.Cox<-coxph(Surv(OR$start,OR$stop,OR$dispute)~allies+contig+capratio+
                growth+democracy+trade,data=OR,method="breslow")

OR.Poisson<-glm(dispute~allies+contig+capratio+growth+democracy+trade
                +as.factor(duration),data=OR,na.action=na.exclude,
                family="poisson")

pdf("CoxPoisson.pdf",10,5)
plot(OR.Cox$coefficients[1:6],OR.Poisson$coefficients[2:7],pch=19,
     xlab="Cox Estimates",ylab="Poisson Estimates",
     xlim=c(-4,1.5),ylim=c(-4,1.5))
abline(0,1)
text(OR.Cox$coefficients[1:6],OR.Poisson$coefficients[2:7],
     labels=colnames(Xhats[1:6]),pos=4,cex=0.8)
legend("bottomright",bty="n",inset=0.02,c("Line is 45-degree line"))
dev.off()

