library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)
library(xtable)
library(dplyr)
library(ggplot2)
library(gplots)
library(rcompanion)
#### Data load ####
#LG
load("C:/Users/com/Desktop/대학원/1. 박사 논문/2. 논문준비 (Quantile residual life)/4. Real data/spec.RData")
#Mac
#load("~/Desktop/Q/Paper/Real data/spec.RData")
# Select amalgam restoration
subDat = subset(spec, restorationType==1)
# Select molar
subDat = subset(subDat, toothType=="molar1"|toothType=="molar23")

#### Subset setting ####
set.seed(2017313045)
unique.id = as.numeric(as.character(unique(subDat$ID)))
data1002 = subDat[FALSE,]
i=1
for (i in 1:length(unique.id)){
  temp = subDat[subDat$ID %in% unique.id[i], ]
  temp_event = temp[temp$event==1,]
  if (nrow(temp_event)>=1){
    data1002 = rbind(data1002, temp_event[which.min(temp_event$survTime),])  
  }
  else{
    data1002 = rbind(data1002, temp[sample(nrow(temp), 1), ])
  }
}
data1020 = subset(data1002[,c(2,3,4,10,11,17,19)])

# Make dummy variable for factors
dummy_Male = as.numeric(data1020$gender == "Male")
dummy_cohort2 = as.numeric(data1020$cohort == 2)
dummy_cohort3 = as.numeric(data1020$cohort == 3)
dummy_cohort4 = as.numeric(data1020$cohort == 4)
dummy_cohort5 = as.numeric(data1020$cohort == 5)
dummy_cohort6 = as.numeric(data1020$cohort == 6)
dummy_Grad = as.numeric(data1020$provider == "Grad")
dummy_Predoc = as.numeric(data1020$provider == "Predoc")
dummy_Private = as.numeric(data1020$Ins == "Private")
dummy_XIX = as.numeric(data1020$Ins == "XIX")
dummy=cbind(dummy_Male,data1020$age,dummy_cohort2,dummy_cohort3,dummy_cohort4,dummy_cohort5,dummy_cohort6,dummy_Grad,dummy_Predoc,dummy_Private,dummy_XIX,data1020[,6:7])
dummy=as.data.frame(dummy)
colnames(dummy) = c("Male","Age","cohort2","cohort3","cohort4","cohort5","cohort6","Grad","Predoc","Private","XIX","event","survTime")

#### Descriptive statistics ####
# covariate 1 : Gender
table(data1002$gender)
gender_mean = groupwiseMean(survTime ~ gender, data=data1002, conf=0.95, digits=3)
gender_median = groupwiseMedian(survTime ~ gender, data=data1002, conf=0.95, digits=3)

# covariate 2 : ae 
table(cut(data1002$age, breaks=seq(0,100,by=10)))

# covariate 3 : cohort
table(data1002$cohort)
cohort_mean = groupwiseMean(survTime ~ cohort, data=data1002, conf=0.95, digits=3)
cohort_median = groupwiseMedian(survTime ~ cohort, data=data1002, conf=0.95, digits=3)

# covariate 4 : provider
table(data1002$provider)
provider_mean = groupwiseMean(survTime ~ provider, data=data1002, conf=0.95, digits=3)
provider_median = groupwiseMedian(survTime ~ provider, data=data1002, conf=0.95, digits=3)

# covariate 5 : Ins(Payment method)
table(data1002$Ins)
Ins_mean = groupwiseMean(survTime ~ Ins, data=data1002, conf=0.95, digits=3)
Ins_median = groupwiseMedian(survTime ~ Ins, data=data1002, conf=0.95, digits=3)

#### Estimation functions ####
# 1. Crq
crq.est = function(Y, nc, covariate, D, t_0, Q){
  n = length(Y)
  data = matrix(NA, n, nc+8)
  # Y
  data[,1] = Y-t_0
  data[,2] = log(data[,1])
  data[,3:(nc+2)] = covariate
  data[,(nc+3)] = D
  
  # Delete data that smaller than t_0
  data = data[complete.cases(data[,1:(nc+3)]),]
  data = data[order(data[,1]),]
  data = as.data.frame(data)
  
  colnames(data)[1:2]=c("Y-t_0", "log(y-t_0)")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[3:(nc+2)]=covar
  colnames(data)[(nc+3)] = c("delta")
  
  cov = as.matrix(data[,3:(nc+2)])
  crq.fit = crq(Surv(data[,2],data[,(nc+3)]) ~ cov,method='Portnoy')
  a = summary(crq.fit,taus = c(0.01, Q))
  beta.sd = matrix(NA, nc+1, 2)
  for (j in 1:(nc+1)){
    beta.sd[j,] = a[2][[1]]$coefficient[j,c(1,4)]
  }
  print(beta.sd)
}

# 2. ISMB
ismb.est = function(Y, nc, covariate, D, t_0, Q, ne){
  n = length(Y)
  data = matrix(NA, n, nc+8)
  # Y
  data[,1] = Y-t_0
  data[,2] = log(data[,1])
  data[,3:(nc+2)] = covariate
  data[,(nc+3)] = D
  # Delete data that smaller than t_0
  data = data[complete.cases(data[,1:(nc+3)]),]
  data = data[order(data[,1]),]
  data = as.data.frame(data)
  
  # Weight
  m = nrow(data)
  for (i in 1:m){
    data[i,(nc+4)] = sum(data[,1] >= data[i,1])
    data[i,(nc+5)] = sum(data[,1] == data[i,1])
    if (data[i,(nc+3)] == 0){
      data[i,(nc+6)] = (data[i,(nc+4)]-data[i,(nc+5)])/data[i,(nc+4)]
    } else {
      data[i,(nc+6)] = 1
    }
    if (i == 1){
      data[i,(nc+7)] = data[i,(nc+6)]
    } else {
      data[i,(nc+7)] = data[i-1,(nc+7)]*data[i,(nc+6)]
    }
    # weight denominator = 0
    if (data[i,(nc+7)] == 0){
      if (data[i,(nc+3)] == 0){
        data[i,(nc+8)] = 0
      } else {
        data[i,(nc+8)] = data[i-1,(nc+8)]}
    } else {
      data[i,(nc+8)] = data[i,(nc+3)]/data[i,(nc+7)]
    }
  }
  colnames(data)[1:2]=c("Y-t_0", "log(y-t_0)")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[3:(nc+2)]=covar
  colnames(data)[(nc+3):(nc+8)] = c("delta","# at risk","# event","s/d","G_KM","Weight")
  
  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(X) %*% (W*(pnorm((logT-X%*%beta)/sqrt(diag(X %*% G %*% t(X))))-1)+Q)
  }
  
  #### revised object equation ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*X) %*% (W*(pnorm((logT-X%*%beta)/sqrt(diag(X %*% G %*% t(X))))-1)+Q)
  }
  
  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,m)),data[,3:(nc+2)]))
  W = data[,(nc+8)]
  logT = data[,2]
  G = diag(1/n, nc+1, nc+1)
  
  betastart = c(-2,rep(0,nc))
  is.fit = nleqslv(betastart, objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(m,1)
      result = t(eta*X)%*%(W*(pnorm((logT-X%*%solbeta)/sqrt(diag(X %*% G %*% t(X))))-1)+Q)
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% G %*% t(X)))))*X)%*%(X/sqrt(diag(X %*% G %*% t(X))))
    inva.beta = qr.solve(a.beta)
    sigma = inva.beta %*% v %*% t(inva.beta)
    sd = sqrt(diag(sigma))
    beta.sd = cbind(solbeta, sd) 
    print(beta.sd)
  } else {
    solbeta = c(NA,NA)
    sd = c(NA,NA)
    beta.sd = cbind(solbeta, sd) 
    print(beta.sd)
  }
}

#### Test ####
# Mixed model
covariate=as.matrix(dummy[,1:11])
result.crq = list()
result.ismb = list()
for (t in 1:5){
  result.crq[[t]] = crq.est(dummy$survTime, 11, covariate, dummy$event, 0, t*0.05)
  result.ismb[[t]] = ismb.est(dummy$survTime, 11, covariate, dummy$event, 0, t*0.05 ,500)
}

#### Draw 95% Confidence Interval for quantile ####
# Covariate 1 : Gender
center.gender = c(result.ismb[[1]][2,1],result.ismb[[2]][2,1],result.ismb[[3]][2,1],result.ismb[[4]][2,1],result.ismb[[5]][2,1])
width.gender = 1.96*c(result.ismb[[1]][2,2],result.ismb[[2]][2,2],result.ismb[[3]][2,2],result.ismb[[4]][2,2],result.ismb[[5]][2,2])
plotCI(center.gender, y=NULL, uiw = width.gender, liw = width.gender, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "Gender",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.gender = center.gender + width.gender
lower.gender = center.gender - width.gender
points(upper.gender, type="b", col="blue", lwd=2)
points(lower.gender, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 2 : age
center.age = c(result.ismb[[1]][3,1],result.ismb[[2]][3,1],result.ismb[[3]][3,1],result.ismb[[4]][3,1],result.ismb[[5]][3,1])
width.age = 1.96*c(result.ismb[[1]][3,2],result.ismb[[2]][3,2],result.ismb[[3]][3,2],result.ismb[[4]][3,2],result.ismb[[5]][3,2])
plotCI(center.age, y=NULL, uiw = width.age, liw = width.age, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "age",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.age = center.age + width.age
lower.age = center.age - width.age
points(upper.age, type="b", col="blue", lwd=2)
points(lower.age, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 3 : cohort2
center.cohort2 = c(result.ismb[[1]][4,1],result.ismb[[2]][4,1],result.ismb[[3]][4,1],result.ismb[[4]][4,1],result.ismb[[5]][4,1])
width.cohort2 = 1.96*c(result.ismb[[1]][4,2],result.ismb[[2]][4,2],result.ismb[[3]][4,2],result.ismb[[4]][4,2],result.ismb[[5]][4,2])
plotCI(center.cohort2, y=NULL, uiw = width.cohort2, liw = width.cohort2, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort2",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort2 = center.cohort2 + width.cohort2
lower.cohort2 = center.cohort2 - width.cohort2
points(upper.cohort2, type="b", col="blue", lwd=2)
points(lower.cohort2, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 4 : cohort3
center.cohort3 = c(result.ismb[[1]][5,1],result.ismb[[2]][5,1],result.ismb[[3]][5,1],result.ismb[[4]][5,1],result.ismb[[5]][5,1])
width.cohort3 = 1.96*c(result.ismb[[1]][5,2],result.ismb[[2]][5,2],result.ismb[[3]][5,2],result.ismb[[4]][5,2],result.ismb[[5]][5,2])
plotCI(center.cohort3, y=NULL, uiw = width.cohort3, liw = width.cohort3, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort3",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort3 = center.cohort3 + width.cohort3
lower.cohort3 = center.cohort3 - width.cohort3
points(upper.cohort3, type="b", col="blue", lwd=2)
points(lower.cohort3, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 5 : cohort4
center.cohort4 = c(result.ismb[[1]][6,1],result.ismb[[2]][6,1],result.ismb[[3]][6,1],result.ismb[[4]][6,1],result.ismb[[5]][6,1])
width.cohort4 = 1.96*c(result.ismb[[1]][6,2],result.ismb[[2]][6,2],result.ismb[[3]][6,2],result.ismb[[4]][6,2],result.ismb[[5]][6,2])
plotCI(center.cohort4, y=NULL, uiw = width.cohort4, liw = width.cohort4, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort4",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort4 = center.cohort4 + width.cohort4
lower.cohort4 = center.cohort4 - width.cohort4
points(upper.cohort4, type="b", col="blue", lwd=2)
points(lower.cohort4, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 6 : cohort5
center.cohort5 = c(result.ismb[[1]][7,1],result.ismb[[2]][7,1],result.ismb[[3]][7,1],result.ismb[[4]][7,1],result.ismb[[5]][7,1])
width.cohort5 = 1.96*c(result.ismb[[1]][7,2],result.ismb[[2]][7,2],result.ismb[[3]][7,2],result.ismb[[4]][7,2],result.ismb[[5]][7,2])
plotCI(center.cohort5, y=NULL, uiw = width.cohort5, liw = width.cohort5, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort5",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort5 = center.cohort5 + width.cohort5
lower.cohort5 = center.cohort5 - width.cohort5
points(upper.cohort5, type="b", col="blue", lwd=2)
points(lower.cohort5, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 7 : cohort6
center.cohort6 = c(result.ismb[[1]][8,1],result.ismb[[2]][8,1],result.ismb[[3]][8,1],result.ismb[[4]][8,1],result.ismb[[5]][8,1])
width.cohort6 = 1.96*c(result.ismb[[1]][8,2],result.ismb[[2]][8,2],result.ismb[[3]][8,2],result.ismb[[4]][8,2],result.ismb[[5]][8,2])
plotCI(center.cohort6, y=NULL, uiw = width.cohort6, liw = width.cohort6, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort6",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort6 = center.cohort6 + width.cohort6
lower.cohort6 = center.cohort6 - width.cohort6
points(upper.cohort6, type="b", col="blue", lwd=2)
points(lower.cohort6, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 8 : Grad
center.Grad = c(result.ismb[[1]][9,1],result.ismb[[2]][9,1],result.ismb[[3]][9,1],result.ismb[[4]][9,1],result.ismb[[5]][9,1])
width.Grad = 1.96*c(result.ismb[[1]][9,2],result.ismb[[2]][9,2],result.ismb[[3]][9,2],result.ismb[[4]][9,2],result.ismb[[5]][9,2])
plotCI(center.Grad, y=NULL, uiw = width.Grad, liw = width.Grad, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "Grad",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.Grad = center.Grad + width.Grad
lower.Grad = center.Grad - width.Grad
points(upper.Grad, type="b", col="blue", lwd=2)
points(lower.Grad, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 9 : Predoc
center.Predoc = c(result.ismb[[1]][10,1],result.ismb[[2]][10,1],result.ismb[[3]][10,1],result.ismb[[4]][10,1],result.ismb[[5]][10,1])
width.Predoc = 1.96*c(result.ismb[[1]][10,2],result.ismb[[2]][10,2],result.ismb[[3]][10,2],result.ismb[[4]][10,2],result.ismb[[5]][10,2])
plotCI(center.Predoc, y=NULL, uiw = width.Predoc, liw = width.Predoc, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "Predoc",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.Predoc = center.Predoc + width.Predoc
lower.Predoc = center.Predoc - width.Predoc
points(upper.Predoc, type="b", col="blue", lwd=2)
points(lower.Predoc, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 10 : private
center.Private = c(result.ismb[[1]][11,1],result.ismb[[2]][11,1],result.ismb[[3]][11,1],result.ismb[[4]][11,1],result.ismb[[5]][11,1])
width.Private = 1.96*c(result.ismb[[1]][11,2],result.ismb[[2]][11,2],result.ismb[[3]][11,2],result.ismb[[4]][11,2],result.ismb[[5]][11,2])
plotCI(center.Private, y=NULL, uiw = width.Private, liw = width.Private, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "Private",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.Private = center.Private + width.Private
lower.Private = center.Private - width.Private
points(upper.Private, type="b", col="blue", lwd=2)
points(lower.Private, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 11 : XIX
center.XIX = c(result.ismb[[1]][12,1],result.ismb[[2]][12,1],result.ismb[[3]][12,1],result.ismb[[4]][12,1],result.ismb[[5]][12,1])
width.XIX = 1.96*c(result.ismb[[1]][12,2],result.ismb[[2]][12,2],result.ismb[[3]][12,2],result.ismb[[4]][12,2],result.ismb[[5]][12,2])
plotCI(center.XIX, y=NULL, uiw = width.XIX, liw = width.XIX, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "XIX",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.XIX = center.XIX + width.XIX
lower.XIX = center.XIX - width.XIX
points(upper.XIX, type="b", col="blue", lwd=2)
points(lower.XIX, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 12 : Base (Female, age=0, Cohort1, Faculty, Self Pay)
center.base = c(result.ismb[[1]][1,1],result.ismb[[2]][1,1],result.ismb[[3]][1,1],result.ismb[[4]][1,1],result.ismb[[5]][1,1])
width.base = 1.96*c(result.ismb[[1]][1,2],result.ismb[[2]][1,2],result.ismb[[3]][1,2],result.ismb[[4]][1,2],result.ismb[[5]][1,2])
plotCI(center.base, y=NULL, uiw = width.base, liw = width.base, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "base",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.base = center.base + width.base
lower.base = center.base - width.base
points(upper.base, type="b", col="blue", lwd=2)
points(lower.base, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

mean.age = mean(dummy$Age)

#### Kaplan Meier fitting for median survival time####
female.fit = survfit(Surv(survTime, event) ~ gender, data=data1020, subset = gender == "Female", conf.type = "none")
male.fit = survfit(Surv(survTime, event) ~ gender, data=data1020, subset = gender == "Male", conf.type = "none")
cohort1.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 1, conf.type = "none")
cohort2.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 2, conf.type = "none")
cohort3.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 3, conf.type = "none")
cohort4.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 4, conf.type = "none")
cohort5.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 5, conf.type = "none")
cohort6.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = cohort == 6, conf.type = "none")
faculty.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = provider == "Faculty", conf.type = "none")
grad.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = provider == "Grad", conf.type = "none")
predoc.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = provider == "Predoc", conf.type = "none")
selfpay.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = Ins == "Self Pay", conf.type = "none")
private.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = Ins == "Private", conf.type = "none")
XIX.fit = survfit(Surv(survTime, event) ~ 1, data=data1020, subset = Ins == "XIX", conf.type = "none")

# Cox modeling
cox.fit = coxph(Surv(survTime, event) ~ ., data=dummy_data1020)
median.crq = crq.est(dummy$survTime, 11, covariate, dummy$event, 0, 0.5)
median.ismb = ismb.est(dummy$survTime, 11, covariate, dummy$event, 0, 0.5 ,500)




#### change of t_0 ####
# t_0 = 0 #
result.ismb = list()
for (t in 1:5){
  result.ismb[[t]] = ismb.est(dummy$survTime, 11, covariate, dummy$event, 0, t*0.05 ,500)
}

# t_0 = 1 #
result.ismb_1 = list()
for (t in 1:5){
  result.ismb_1[[t]] = ismb.est(dummy$survTime, 11, covariate, dummy$event, 1, t*0.05 ,500)
}

# t_0 = 2 #
#dummy_2 = dummy
#dummy_2$survTime = dummy$survTime-2
#covariate_2=as.matrix(dummy_2[,1:11])
result.ismb_2 = list()
for (t in 1:5){
  result.ismb_2[[t]] = ismb.est(dummy$survTime, 11, covariate, dummy$event, 2, t*0.05 ,500)
}

# t_0 = 3 #
#dummy_3 = dummy
#dummy_3$survTime = dummy$survTime-3
#covariate_3=as.matrix(dummy_3[,1:11])
result.ismb_3 = list()
for (t in 1:5){
  result.ismb_3[[t]] = ismb.est(dummy$survTime, 11, covariate, dummy$event, 3, t*0.05 ,500)
}

#### Draw 95% Confidence Interval for t_0 w/ 10% quantile####
# Covariate 1 : Gender
center.gender = c(result.ismb[[2]][2,1],result.ismb_1[[2]][2,1],result.ismb_2[[2]][2,1],result.ismb_3[[2]][2,1])
width.gender = 1.96*c(result.ismb[[2]][2,2],result.ismb_1[[2]][2,2],result.ismb_2[[2]][2,2],result.ismb_3[[2]][2,2])
plotCI(center.gender, y=NULL, uiw = width.gender, liw = width.gender, axes=FALSE, xlab = "t_0", ylab = "beta", main = "Gender",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.gender = center.gender + width.gender
lower.gender = center.gender - width.gender
points(upper.gender, type="b", col="blue", lwd=2)
points(lower.gender, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 2 : age
center.age = c(result.ismb[[2]][3,1],result.ismb_1[[2]][3,1],result.ismb_2[[2]][3,1],result.ismb_3[[2]][3,1])
width.age = 1.96*c(result.ismb[[2]][3,2],result.ismb_1[[2]][3,2],result.ismb_2[[2]][3,2],result.ismb_3[[2]][3,2])
plotCI(center.age, y=NULL, uiw = width.age, liw = width.age, axes=FALSE, xlab = "t_0", ylab = "beta", main = "age",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.age = center.age + width.age
lower.age = center.age - width.age
points(upper.age, type="b", col="blue", lwd=2)
points(lower.age, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 3 : cohort2
center.cohort2 = c(result.ismb[[2]][4,1],result.ismb_1[[2]][4,1],result.ismb_2[[2]][4,1],result.ismb_3[[2]][4,1])
width.cohort2 = 1.96*c(result.ismb[[2]][4,2],result.ismb_1[[2]][4,2],result.ismb_2[[2]][4,2],result.ismb_3[[2]][4,2])
plotCI(center.cohort2, y=NULL, uiw = width.cohort2, liw = width.cohort2, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort2",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.cohort2 = center.cohort2 + width.cohort2
lower.cohort2 = center.cohort2 - width.cohort2
points(upper.cohort2, type="b", col="blue", lwd=2)
points(lower.cohort2, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 4 : cohort3
center.cohort3 = c(result.ismb[[2]][5,1],result.ismb_1[[2]][5,1],result.ismb_2[[2]][5,1],result.ismb_3[[2]][5,1])
width.cohort3 = 1.96*c(result.ismb[[2]][5,2],result.ismb_1[[2]][5,2],result.ismb_2[[2]][5,2],result.ismb_3[[2]][5,2])
plotCI(center.cohort3, y=NULL, uiw = width.cohort3, liw = width.cohort3, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort3",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.cohort3 = center.cohort3 + width.cohort3
lower.cohort3 = center.cohort3 - width.cohort3
points(upper.cohort3, type="b", col="blue", lwd=2)
points(lower.cohort3, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 5 : cohort4
center.cohort4 = c(result.ismb[[2]][6,1],result.ismb_1[[2]][6,1],result.ismb_2[[2]][6,1],result.ismb_3[[2]][6,1])
width.cohort4 = 1.96*c(result.ismb[[2]][6,2],result.ismb_1[[2]][6,2],result.ismb_2[[2]][6,2],result.ismb_3[[2]][6,2])
plotCI(center.cohort4, y=NULL, uiw = width.cohort4, liw = width.cohort4, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort4",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.cohort4 = center.cohort4 + width.cohort4
lower.cohort4 = center.cohort4 - width.cohort4
points(upper.cohort4, type="b", col="blue", lwd=2)
points(lower.cohort4, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 6 : cohort5
center.cohort5 = c(result.ismb[[2]][7,1],result.ismb_1[[2]][7,1],result.ismb_2[[2]][7,1],result.ismb_3[[2]][7,1])
width.cohort5 = 1.96*c(result.ismb[[2]][7,2],result.ismb_1[[2]][7,2],result.ismb_2[[2]][7,2],result.ismb_3[[2]][7,2])
plotCI(center.cohort5, y=NULL, uiw = width.cohort5, liw = width.cohort5, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort5",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.cohort5 = center.cohort5 + width.cohort5
lower.cohort5 = center.cohort5 - width.cohort5
points(upper.cohort5, type="b", col="blue", lwd=2)
points(lower.cohort5, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 7 : cohort6
center.cohort6 = c(result.ismb[[2]][8,1],result.ismb_1[[2]][8,1],result.ismb_2[[2]][8,1],result.ismb_3[[2]][8,1])
width.cohort6 = 1.96*c(result.ismb[[2]][8,2],result.ismb_1[[2]][8,2],result.ismb_2[[2]][8,2],result.ismb_3[[2]][8,2])
plotCI(center.cohort6, y=NULL, uiw = width.cohort6, liw = width.cohort6, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort6",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.cohort6 = center.cohort6 + width.cohort6
lower.cohort6 = center.cohort6 - width.cohort6
points(upper.cohort6, type="b", col="blue", lwd=2)
points(lower.cohort6, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 8 : Grad
center.Grad = c(result.ismb[[2]][9,1],result.ismb_1[[2]][9,1],result.ismb_2[[2]][9,1],result.ismb_3[[2]][9,1])
width.Grad = 1.96*c(result.ismb[[2]][9,2],result.ismb_1[[2]][9,2],result.ismb_2[[2]][9,2],result.ismb_3[[2]][9,2])
plotCI(center.Grad, y=NULL, uiw = width.Grad, liw = width.Grad, axes=FALSE, xlab = "t_0", ylab = "beta", main = "Grad",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.Grad = center.Grad + width.Grad
lower.Grad = center.Grad - width.Grad
points(upper.Grad, type="b", col="blue", lwd=2)
points(lower.Grad, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 9 : Predoc
center.Predoc = c(result.ismb[[2]][10,1],result.ismb_1[[2]][10,1],result.ismb_2[[2]][10,1],result.ismb_3[[2]][10,1])
width.Predoc = 1.96*c(result.ismb[[2]][10,2],result.ismb_1[[2]][10,2],result.ismb_2[[2]][10,2],result.ismb_3[[2]][10,2])
plotCI(center.Predoc, y=NULL, uiw = width.Predoc, liw = width.Predoc, axes=FALSE, xlab = "t_0", ylab = "beta", main = "Predoc",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.Predoc = center.Predoc + width.Predoc
lower.Predoc = center.Predoc - width.Predoc
points(upper.Predoc, type="b", col="blue", lwd=2)
points(lower.Predoc, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 10 : private
center.Private = c(result.ismb[[2]][11,1],result.ismb_1[[2]][11,1],result.ismb_2[[2]][11,1],result.ismb_3[[2]][11,1])
width.Private = 1.96*c(result.ismb[[2]][11,2],result.ismb_1[[2]][11,2],result.ismb_2[[2]][11,2],result.ismb_3[[2]][11,2])
plotCI(center.Private, y=NULL, uiw = width.Private, liw = width.Private, axes=FALSE, xlab = "t_0", ylab = "beta", main = "Private",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.Private = center.Private + width.Private
lower.Private = center.Private - width.Private
points(upper.Private, type="b", col="blue", lwd=2)
points(lower.Private, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 11 : XIX
center.XIX = c(result.ismb[[2]][12,1],result.ismb_1[[2]][12,1],result.ismb_2[[2]][12,1],result.ismb_3[[2]][12,1])
width.XIX = 1.96*c(result.ismb[[2]][12,2],result.ismb_1[[2]][12,2],result.ismb_2[[2]][12,2],result.ismb_3[[2]][12,2])
plotCI(center.XIX, y=NULL, uiw = width.XIX, liw = width.XIX, axes=FALSE, xlab = "t_0", ylab = "beta", main = "XIX",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.XIX = center.XIX + width.XIX
lower.XIX = center.XIX - width.XIX
points(upper.XIX, type="b", col="blue", lwd=2)
points(lower.XIX, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 12 : Base (Female, age=0, Cohort1, Faculty, Self Pay)
center.base = c(result.ismb[[2]][1,1],result.ismb_1[[2]][1,1],result.ismb_2[[2]][1,1],result.ismb_3[[2]][1,1])
width.base = 1.96*c(result.ismb[[2]][1,2],result.ismb_1[[2]][1,2],result.ismb_2[[2]][1,2],result.ismb_3[[2]][1,2])
plotCI(center.base, y=NULL, uiw = width.base, liw = width.base, axes=FALSE, xlab = "t_0", ylab = "beta", main = "base",type = "o")
axis(side=2)
axis(side=1,at=1:4,label=c("0","1","2","3"))
upper.base = center.base + width.base
lower.base = center.base - width.base
points(upper.base, type="b", col="blue", lwd=2)
points(lower.base, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)


# For latex #

hoy_1 = matrix(NA, 12, 8)
rownames(hoy_1) = c("Base", "Male", "Age", "Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Grad","Predoc","Private","XIX")
hoy_2 = matrix(NA, 12, 8)
rownames(hoy_2) = c("Base", "Male", "Age", "Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Grad","Predoc","Private","XIX")
hoy_3 = matrix(NA, 12, 8)
rownames(hoy_3) = c("Base", "Male", "Age", "Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Grad","Predoc","Private","XIX")
hoy_4 = matrix(NA, 12, 8)
rownames(hoy_4) = c("Base", "Male", "Age", "Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Grad","Predoc","Private","XIX")
hoy_5 = matrix(NA, 12, 8)
rownames(hoy_5) = c("Base", "Male", "Age", "Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Grad","Predoc","Private","XIX")

hoy_1[,1:2] = result.ismb[[1]]
hoy_2[,1:2] = result.ismb[[2]]
hoy_3[,1:2] = result.ismb[[3]]
hoy_4[,1:2] = result.ismb[[4]]
hoy_5[,1:2] = result.ismb[[5]]

hoy_1[,3:4] = result.ismb_1[[1]]
hoy_2[,3:4] = result.ismb_1[[2]]
hoy_3[,3:4] = result.ismb_1[[3]]
hoy_4[,3:4] = result.ismb_1[[4]]
hoy_5[,3:4] = result.ismb_1[[5]]

hoy_1[,5:6] = result.ismb_2[[1]]
hoy_2[,5:6] = result.ismb_2[[2]]
hoy_3[,5:6] = result.ismb_2[[3]]
hoy_4[,5:6] = result.ismb_2[[4]]
hoy_5[,5:6] = result.ismb_2[[5]]

hoy_1[,7:8] = result.ismb_3[[1]]
hoy_2[,7:8] = result.ismb_3[[2]]
hoy_3[,7:8] = result.ismb_3[[3]]
hoy_4[,7:8] = result.ismb_3[[4]]
hoy_5[,7:8] = result.ismb_3[[5]]

xtable(hoy_1, digits=4)
xtable(hoy_2, digits=4)
xtable(hoy_3, digits=4)
xtable(hoy_4, digits=4)
xtable(hoy_5, digits=4)