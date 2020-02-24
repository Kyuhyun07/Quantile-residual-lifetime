# 2020-02-07 by SW
# - Some changes in creating a subset of the data & covariates
# 1. Only Cohorts 4 - 6 and age >= 65 will be used
#   (following Caplan et al. (2019)'s construction of data )
# 2. This results in 1551 unique patients(same as that in Caplan et al.)
# 3. Event: any restoration, extraction or endodontic therapy on any surface
# 4. Survival time: time from initial restoration (any type) to event
# 5. Time to first event will be chosen as before
# 6. If there is no event then randomly choose one observation as before
# 7. Restoration type and tooth type will be considered as covariates
# 8. Some covariates were re-defined (based on Caplan et al.'s paper)

# Not all packages are needed....
library(quantreg)
library(survival)
library(nleqslv)
library(ggplot2)
library(gplots)
# Setting the working directory
#setwd("C:\\Users\\sangw\\Dropbox\\Sangwook\\Research\\Projects\\SQRL\\R")

#### Data load ####
#load("data/spec.RData")
# Mac
load("/Users/kyuhyunkim/Desktop/Q/Paper/4. Real data/spec.RData")
# PC
# load("C:/Users/com/Desktop/???п?/1. ?ڻ? ?���/2. ?���?غ? (Quantile residual life)/4. Real data/spec.RData")
#temp <- subset(spec, !is.na(spec$restTooth))


#### Subsetting data ####
# Added by SW (2020-02-06): cohort 4 - 6 & restoration placed >=65 years
# Note: criterion in the original paper (Caplan et al. 2019)
# 9184 observations with 1551 unique individuals
subDat <- subset(spec, cohort %in% paste(4:6) & age >= 65)

#### Subset setting ####
# - Selecting the only one observation per person;
# - Selecting the very first event time with amalgam;
# - When there is no event, select one molar randomly (censoring)
set.seed(1)
unique.id = as.numeric(as.character(unique(subDat$ID)))
data1002 = subDat[FALSE,]
#i <- 2
for (i in 1:length(unique.id)){ # length(unique.id): # of subjects
  temp = subDat[subDat$ID %in% unique.id[i], ]
  temp_event = temp[temp$event==1,]
  if (nrow(temp_event)>=1){
    data1002 = rbind(data1002, temp_event[which.min(temp_event$survTime),])
  }
  else{
    data1002 = rbind(data1002, temp[sample(nrow(temp), 1), ])
  }
}

data1020 = data1002

# Make dummy variable for factors
## For Gender (Female as reference)
dummy_Male = as.numeric(data1020$gender == "Male")

## For Cohort (6 as reference)
dummy_cohort4 = as.numeric(data1020$cohort == 4)
dummy_cohort5 = as.numeric(data1020$cohort == 5)
# For provider (Grad + Faculty as reference)
dummy_Predoc = as.numeric(data1020$provider == "Predoc")

# For insurance (Self-pay as reference)
dummy_Private = as.numeric(data1020$Ins == "Private")
dummy_XIX = as.numeric(data1020$Ins == "XIX")

# For tooth type (Anterior as reference)
dummy_Molar = as.numeric(data1020$toothType=="molar1"|
                        data1020$toothType=="molar23")
dummy_PreMolar = as.numeric(data1020$toothType=="premolar")

# For restoration type (crown or bridge as reference)
dummy_Amalgam = as.numeric(data1020$restorationType == 1)
dummy_Composite = as.numeric(data1020$restorationType == 2)
dummy_GIC = as.numeric(data1020$restorationType == 3)

# Only "Age"
dummy=cbind(dummy_Male,data1020$age-65,dummy_cohort4,dummy_cohort5,dummy_Predoc,dummy_Private,dummy_XIX,dummy_Molar,dummy_PreMolar,dummy_Amalgam,dummy_Composite,dummy_GIC,data1020[,c(17,19)])
dummy=as.data.frame(dummy)
colnames(dummy) = c("Male","Age","cohort4","cohort5","Predoc","Private","XIX","Molar","Premolar","Amalgam","Composite","GIC","event","survTime")


# Add "Age^2"
#dummy=cbind(dummy_Male,data1020$age-65,(data1020$age-65)^2,dummy_cohort4,dummy_cohort5,dummy_Predoc,dummy_Private,dummy_XIX,dummy_Molar,dummy_PreMolar,dummy_Amalgam,dummy_Composite,dummy_GIC,data1020[,c(17,19)])
#dummy=as.data.frame(dummy)
#colnames(dummy) = c("Male","Age","Age^2","cohort4","cohort5","Predoc","Private","XIX","Molar","Premolar","Amalgam","Composite","GIC","event","survTime")

#### 1. Delete tied data in Z scale ####
# round up survtime
dummy$survTime = round(dummy$survTime,digits = 2)
unique.survTime = as.numeric(as.character(unique(dummy$survTime)))
data_survTime = dummy[FALSE,]
# one data per each survtime (no matter censored or not)
for (i in 1:length(unique.id)){ # length(unique.id): # of subjects
  temp = dummy[dummy$survTime %in% unique.survTime[i], ]
  temp_event = temp[temp$event==1,] 
  if (nrow(temp_event)>=1){
    data_survTime = rbind(data_survTime, temp[temp$event==0,], temp_event[sample(nrow(temp_event), 1), ])
  } else {
    data_survTime = rbind(data_survTime, temp[temp$event==0,])
  }
}
data_survTime = data_survTime[order(data_survTime$survTime),]

# # save data
# save(data_survTime, file = "data_survTime.RData")

# Weight comparison in real data => w2 is best
load("/Users/kyuhyunkim/Desktop/Q/Paper/4. Real data/data_survTime.RData")
covar=as.matrix(data_survTime[,1:12])
w1 = weight1(data_survTime[,14], 12, covar, data_survTime[,13], 0)
w2 = weight2(data_survTime[,14], 12, covar, data_survTime[,13], 0)
w3 = weight3(data_survTime[,14], 12, covar, data_survTime[,13], 0)
result = cbind(w1, w2, w3, data_survTime$event)

# #### Example parameter ####
# Z=data_survTime$survTime
# nc=12
# covariate=as.matrix(data_survTime[,1:12])
# D=data_survTime$event
# t_0=1
# Q=0.05
# ne=100

#### Estimation functions ####
# 1. Crq
crq.est = function(Z, nc, covariate, D, t_0, Q){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  # Y
  data[,1] = Z
  data[,2] = log(data[,1]-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.nan(data[,2]),2]=NA
  data[data[,2]==-Inf,2]=NA
  data[,(nc+4)] = D
  data = as.data.frame(data)
  
  # Weight
  # Kaplan-Meier estimator for censoring
  fit = survfit(Surv(data[,1],(1-data[,(nc+4)])) ~ 1)
  for (i in 1:length(fit$surv))
  {
    data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  }
  data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  data[data[,(nc+5)]==0,(nc+6)]=0
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[T>=C]")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[4:(nc+3)]=covar
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  
  cov = as.matrix(data[,4:(nc+3)])
  crq.fit = crq(Surv(data[,2],data[,(nc+4)]) ~ cov,method='Portnoy')
  a = summary(crq.fit,taus = c(0.01, Q))
  beta.sd = matrix(NA, nc+1, 2)
  for (j in 1:(nc+1)){
    beta.sd[j,] = a[2][[1]]$coefficient[j,c(1,4)]
  }
  print(beta.sd)
}

# 2. rq with weight
rq.est = function(Z, nc, covariate, D, t_0, Q){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  # Y
  data[,1] = Z
  data[,2] = log(data[,1]-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.nan(data[,2]),2]=NA
  data[data[,2]==-Inf,2]=NA
  data[,(nc+4)] = D
  data = as.data.frame(data)
  
  # Weight
  # Kaplan-Meier estimator for censoring
  fit = survfit(Surv(data[,1],(1-data[,(nc+4)])) ~ 1)
  for (i in 1:length(fit$surv))
  {
    data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  }
  data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  data[data[,(nc+5)]==0,(nc+6)]=0
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[T>=C]")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[4:(nc+3)]=covar
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  

  cov = as.matrix(data[,4:(nc+3)])
  rq.fit = rq(data[,2] ~ cov, tau=Q, weight = data[,(nc+6)])
  rst = summary.rq(rq.fit, se="boot")
  print(rst$coefficients[,c(1,2)])
}

# 3. ISMB
ismb.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  # Y
  data[,1] = Z
  data[,2] = log(data[,1]-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.nan(data[,2]),2]=-10
  data[,(nc+4)] = D
  data = as.data.frame(data)
  
  # Weight
  # Kaplan-Meier estimator for censoring
  fit = survfit(Surv(data[,1],(1-data[,(nc+4)])) ~ 1)
  for (i in 1:length(fit$surv))
  {
    data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  }
  data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  data[data[,(nc+5)]==0,(nc+6)]=0
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[T>=C]")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[4:(nc+3)]=covar
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  
  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(I*X*W) %*% {Q - (pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
  }
  
  #### revised object equation ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*X*I*W) %*% {Q - (pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
  }
  
  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+6)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)
  
  betastart = c(-2,rep(0,nc))
  is.fit = nleqslv(betastart,objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I) %*% {Q - W*(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(I*X*W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
    inva.beta = qr.solve(a.beta)
    sigma = t(inva.beta) %*% v %*% inva.beta
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



#### Model fitting (20200224) ####
# 1. different quantile
covariate=as.matrix(data_survTime[,1:12])
result.crq = list()
result.rq = list()
result.ismb = list()
for (t in 1:6){
  result.crq[[t]] = crq.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, t*0.05)
  result.rq[[t]] = rq.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, t*0.05)
  result.ismb[[t]] = ismb.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, t*0.05 ,100)
}

# 2. different t_0
covariate=as.matrix(data_survTime[,1:12])
tresult.crq = list()
tresult.rq = list()
tresult.ismb = list()

for (t in 0:4){
  tresult.crq[[t+1]] = crq.est(data_survTime$survTime, 12, covariate, data_survTime$event, t*0.5, 0.20)
  tresult.rq[[t+1]] = rq.est(data_survTime$survTime, 12, covariate, data_survTime$event, t*0.5, 0.20)
  tresult.ismb[[t+1]] = ismb.est(data_survTime$survTime, 12, covariate, data_survTime$event, t*0.5, 0.20 ,100)
}
#### time comparison ###
tic()
crq.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, 0.40)
toc()

tic()
rq.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, 0.99)
toc()

tic()
ismb.est(data_survTime$survTime, 12, covariate, data_survTime$event, 0, 0.37 ,100)
toc()

#### Draw CI with different quantiles ####
# Covariate 1 : gender
center.gender = c(result.ismb[[1]][2,1],result.ismb[[2]][2,1],result.ismb[[3]][2,1],result.ismb[[4]][2,1],result.ismb[[5]][2,1])
width.gender = 1.96*c(result.ismb[[1]][2,2],result.ismb[[2]][2,2],result.ismb[[3]][2,2],result.ismb[[4]][2,2],result.ismb[[5]][2,2])
plotCI(center.gender, y=NULL, uiw = width.gender, liw = width.gender, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "gender",type = "o")
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

# Covariate 3 : cohort 4
center.cohort4 = c(result.ismb[[1]][4,1],result.ismb[[2]][4,1],result.ismb[[3]][4,1],result.ismb[[4]][4,1],result.ismb[[5]][4,1])
width.cohort4 = 1.96*c(result.ismb[[1]][4,2],result.ismb[[2]][4,2],result.ismb[[3]][4,2],result.ismb[[4]][4,2],result.ismb[[5]][4,2])
plotCI(center.cohort4, y=NULL, uiw = width.cohort4, liw = width.cohort4, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort4",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort4 = center.cohort4 + width.cohort4
lower.cohort4 = center.cohort4 - width.cohort4
points(upper.cohort4, type="b", col="blue", lwd=2)
points(lower.cohort4, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 4 : cohort5
center.cohort5 = c(result.ismb[[1]][5,1],result.ismb[[2]][5,1],result.ismb[[3]][5,1],result.ismb[[4]][5,1],result.ismb[[5]][5,1])
width.cohort5 = 1.96*c(result.ismb[[1]][5,2],result.ismb[[2]][5,2],result.ismb[[3]][5,2],result.ismb[[4]][5,2],result.ismb[[5]][5,2])
plotCI(center.cohort5, y=NULL, uiw = width.cohort5, liw = width.cohort5, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "cohort5",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.cohort5 = center.cohort5 + width.cohort5
lower.cohort5 = center.cohort5 - width.cohort5
points(upper.cohort5, type="b", col="blue", lwd=2)
points(lower.cohort5, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 5 : predoc
center.predoc = c(result.ismb[[1]][6,1],result.ismb[[2]][6,1],result.ismb[[3]][6,1],result.ismb[[4]][6,1],result.ismb[[5]][6,1])
width.predoc = 1.96*c(result.ismb[[1]][6,2],result.ismb[[2]][6,2],result.ismb[[3]][6,2],result.ismb[[4]][6,2],result.ismb[[5]][6,2])
plotCI(center.predoc, y=NULL, uiw = width.predoc, liw = width.predoc, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "predoc",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.predoc = center.predoc + width.predoc
lower.predoc = center.predoc - width.predoc
points(upper.predoc, type="b", col="blue", lwd=2)
points(lower.predoc, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 6 : private
center.private = c(result.ismb[[1]][7,1],result.ismb[[2]][7,1],result.ismb[[3]][7,1],result.ismb[[4]][7,1],result.ismb[[5]][7,1])
width.private = 1.96*c(result.ismb[[1]][7,2],result.ismb[[2]][7,2],result.ismb[[3]][7,2],result.ismb[[4]][7,2],result.ismb[[5]][7,2])
plotCI(center.private, y=NULL, uiw = width.private, liw = width.private, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "private",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.private = center.private + width.private
lower.private = center.private - width.private
points(upper.private, type="b", col="blue", lwd=2)
points(lower.private, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 7 : XIX
center.XIX = c(result.ismb[[1]][8,1],result.ismb[[2]][8,1],result.ismb[[3]][8,1],result.ismb[[4]][8,1],result.ismb[[5]][8,1])
width.XIX = 1.96*c(result.ismb[[1]][8,2],result.ismb[[2]][8,2],result.ismb[[3]][8,2],result.ismb[[4]][8,2],result.ismb[[5]][8,2])
plotCI(center.XIX, y=NULL, uiw = width.XIX, liw = width.XIX, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "XIX",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.XIX = center.XIX + width.XIX
lower.XIX = center.XIX - width.XIX
points(upper.XIX, type="b", col="blue", lwd=2)
points(lower.XIX, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 8 : molar
center.molar = c(result.ismb[[1]][9,1],result.ismb[[2]][9,1],result.ismb[[3]][9,1],result.ismb[[4]][9,1],result.ismb[[5]][9,1])
width.molar = 1.96*c(result.ismb[[1]][9,2],result.ismb[[2]][9,2],result.ismb[[3]][9,2],result.ismb[[4]][9,2],result.ismb[[5]][9,2])
plotCI(center.molar, y=NULL, uiw = width.molar, liw = width.molar, axes=FALSE, xlab = "t", ylab = "beta", main = "molar",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.molar = center.molar + width.molar
lower.molar = center.molar - width.molar
points(upper.molar, type="b", col="blue", lwd=2)
points(lower.molar, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 9 : premolar
center.premolar = c(result.ismb[[1]][10,1],result.ismb[[2]][10,1],result.ismb[[3]][10,1],result.ismb[[4]][10,1],result.ismb[[5]][10,1])
width.premolar = 1.96*c(result.ismb[[1]][10,2],result.ismb[[2]][10,2],result.ismb[[3]][10,2],result.ismb[[4]][10,2],result.ismb[[5]][10,2])
plotCI(center.premolar, y=NULL, uiw = width.premolar, liw = width.premolar, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "premolar",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.premolar = center.premolar + width.premolar
lower.premolar = center.premolar - width.premolar
points(upper.premolar, type="b", col="blue", lwd=2)
points(lower.premolar, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 10 : amalgam
center.amalgam = c(result.ismb[[1]][11,1],result.ismb[[2]][11,1],result.ismb[[3]][11,1],result.ismb[[4]][11,1],result.ismb[[5]][11,1])
width.amalgam = 1.96*c(result.ismb[[1]][11,2],result.ismb[[2]][11,2],result.ismb[[3]][11,2],result.ismb[[4]][11,2],result.ismb[[5]][11,2])
plotCI(center.amalgam, y=NULL, uiw = width.amalgam, liw = width.amalgam, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "amalgam",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.amalgam = center.amalgam + width.amalgam
lower.amalgam = center.amalgam - width.amalgam
points(upper.amalgam, type="b", col="blue", lwd=2)
points(lower.amalgam, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 11 : composite
center.composite = c(result.ismb[[1]][12,1],result.ismb[[2]][12,1],result.ismb[[3]][12,1],result.ismb[[4]][12,1],result.ismb[[5]][12,1])
width.composite = 1.96*c(result.ismb[[1]][12,2],result.ismb[[2]][12,2],result.ismb[[3]][12,2],result.ismb[[4]][12,2],result.ismb[[5]][12,2])
plotCI(center.composite, y=NULL, uiw = width.composite, liw = width.composite, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "composite",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.composite = center.composite + width.composite
lower.composite = center.composite - width.composite
points(upper.composite, type="b", col="blue", lwd=2)
points(lower.composite, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 12 : GIC
center.GIC = c(result.ismb[[1]][13,1],result.ismb[[2]][13,1],result.ismb[[3]][13,1],result.ismb[[4]][13,1],result.ismb[[5]][13,1])
width.GIC = 1.96*c(result.ismb[[1]][13,2],result.ismb[[2]][13,2],result.ismb[[3]][13,2],result.ismb[[4]][13,2],result.ismb[[5]][13,2])
plotCI(center.GIC, y=NULL, uiw = width.GIC, liw = width.GIC, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "GIC",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
upper.GIC = center.GIC + width.GIC
lower.GIC = center.GIC - width.GIC
points(upper.GIC, type="b", col="blue", lwd=2)
points(lower.GIC, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# 
# # Covariate 2 : age^2
# center.age2 = c(result.ismb[[1]][4,1],result.ismb[[2]][4,1],result.ismb[[3]][4,1],result.ismb[[4]][4,1],result.ismb[[5]][4,1])
# width.age2 = 1.96*c(result.ismb[[1]][4,2],result.ismb[[2]][4,2],result.ismb[[3]][4,2],result.ismb[[4]][4,2],result.ismb[[5]][4,2])
# plotCI(center.age2, y=NULL, uiw = width.age2, liw = width.age2, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "age^2",type = "o")
# axis(side=2)
# axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
# upper.age2 = center.age2 + width.age2
# lower.age2 = center.age2 - width.age2
# points(upper.age2, type="b", col="blue", lwd=2)
# points(lower.age2, type="b", col="blue", lwd=2)
# abline(h=0, col="red",lty=2)

#### Draw CI with different t_0 ####
# Covariate 1 : gender
center.gender = c(tresult.ismb[[1]][2,1],tresult.ismb[[2]][2,1],tresult.ismb[[3]][2,1],tresult.ismb[[4]][2,1],tresult.ismb[[5]][2,1])
width.gender = 1.96*c(tresult.ismb[[1]][2,2],tresult.ismb[[2]][2,2],tresult.ismb[[3]][2,2],tresult.ismb[[4]][2,2],tresult.ismb[[5]][2,2])
plotCI(center.gender, y=NULL, uiw = width.gender, liw = width.gender, axes=FALSE, xlab = "t_0", ylab = "beta", main = "gender",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.gender = center.gender + width.gender
lower.gender = center.gender - width.gender
points(upper.gender, type="b", col="blue", lwd=2)
points(lower.gender, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 2 : age
center.age = c(tresult.ismb[[1]][3,1],tresult.ismb[[2]][3,1],tresult.ismb[[3]][3,1],tresult.ismb[[4]][3,1],tresult.ismb[[5]][3,1])
width.age = 1.96*c(tresult.ismb[[1]][3,2],tresult.ismb[[2]][3,2],tresult.ismb[[3]][3,2],tresult.ismb[[4]][3,2],tresult.ismb[[5]][3,2])
plotCI(center.age, y=NULL, uiw = width.age, liw = width.age, axes=FALSE, xlab = "t_0", ylab = "beta", main = "age",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.age = center.age + width.age
lower.age = center.age - width.age
points(upper.age, type="b", col="blue", lwd=2)
points(lower.age, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 3 : cohort 4
center.cohort4 = c(tresult.ismb[[1]][4,1],tresult.ismb[[2]][4,1],tresult.ismb[[3]][4,1],tresult.ismb[[4]][4,1],tresult.ismb[[5]][4,1])
width.cohort4 = 1.96*c(tresult.ismb[[1]][4,2],tresult.ismb[[2]][4,2],tresult.ismb[[3]][4,2],tresult.ismb[[4]][4,2],tresult.ismb[[5]][4,2])
plotCI(center.cohort4, y=NULL, uiw = width.cohort4, liw = width.cohort4, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort4",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.cohort4 = center.cohort4 + width.cohort4
lower.cohort4 = center.cohort4 - width.cohort4
points(upper.cohort4, type="b", col="blue", lwd=2)
points(lower.cohort4, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 4 : cohort5
center.cohort5 = c(tresult.ismb[[1]][5,1],tresult.ismb[[2]][5,1],tresult.ismb[[3]][5,1],tresult.ismb[[4]][5,1],tresult.ismb[[5]][5,1])
width.cohort5 = 1.96*c(tresult.ismb[[1]][5,2],tresult.ismb[[2]][5,2],tresult.ismb[[3]][5,2],tresult.ismb[[4]][5,2],tresult.ismb[[5]][5,2])
plotCI(center.cohort5, y=NULL, uiw = width.cohort5, liw = width.cohort5, axes=FALSE, xlab = "t_0", ylab = "beta", main = "cohort5",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.cohort5 = center.cohort5 + width.cohort5
lower.cohort5 = center.cohort5 - width.cohort5
points(upper.cohort5, type="b", col="blue", lwd=2)
points(lower.cohort5, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 5 : predoc
center.predoc = c(tresult.ismb[[1]][6,1],tresult.ismb[[2]][6,1],tresult.ismb[[3]][6,1],tresult.ismb[[4]][6,1],tresult.ismb[[5]][6,1])
width.predoc = 1.96*c(tresult.ismb[[1]][6,2],tresult.ismb[[2]][6,2],tresult.ismb[[3]][6,2],tresult.ismb[[4]][6,2],tresult.ismb[[5]][6,2])
plotCI(center.predoc, y=NULL, uiw = width.predoc, liw = width.predoc, axes=FALSE, xlab = "t_0", ylab = "beta", main = "predoc",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.predoc = center.predoc + width.predoc
lower.predoc = center.predoc - width.predoc
points(upper.predoc, type="b", col="blue", lwd=2)
points(lower.predoc, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 6 : private
center.private = c(tresult.ismb[[1]][7,1],tresult.ismb[[2]][7,1],tresult.ismb[[3]][7,1],tresult.ismb[[4]][7,1],tresult.ismb[[5]][7,1])
width.private = 1.96*c(tresult.ismb[[1]][7,2],tresult.ismb[[2]][7,2],tresult.ismb[[3]][7,2],tresult.ismb[[4]][7,2],tresult.ismb[[5]][7,2])
plotCI(center.private, y=NULL, uiw = width.private, liw = width.private, axes=FALSE, xlab = "t_0", ylab = "beta", main = "private",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.private = center.private + width.private
lower.private = center.private - width.private
points(upper.private, type="b", col="blue", lwd=2)
points(lower.private, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 7 : XIX
center.XIX = c(tresult.ismb[[1]][8,1],tresult.ismb[[2]][8,1],tresult.ismb[[3]][8,1],tresult.ismb[[4]][8,1],tresult.ismb[[5]][8,1])
width.XIX = 1.96*c(tresult.ismb[[1]][8,2],tresult.ismb[[2]][8,2],tresult.ismb[[3]][8,2],tresult.ismb[[4]][8,2],tresult.ismb[[5]][8,2])
plotCI(center.XIX, y=NULL, uiw = width.XIX, liw = width.XIX, axes=FALSE, xlab = "t_0", ylab = "beta", main = "XIX",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.XIX = center.XIX + width.XIX
lower.XIX = center.XIX - width.XIX
points(upper.XIX, type="b", col="blue", lwd=2)
points(lower.XIX, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 8 : molar
center.molar = c(tresult.ismb[[1]][9,1],tresult.ismb[[2]][9,1],tresult.ismb[[3]][9,1],tresult.ismb[[4]][9,1],tresult.ismb[[5]][9,1])
width.molar = 1.96*c(tresult.ismb[[1]][9,2],tresult.ismb[[2]][9,2],tresult.ismb[[3]][9,2],tresult.ismb[[4]][9,2],tresult.ismb[[5]][9,2])
plotCI(center.molar, y=NULL, uiw = width.molar, liw = width.molar, axes=FALSE, xlab = "t_0", ylab = "beta", main = "molar",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.molar = center.molar + width.molar
lower.molar = center.molar - width.molar
points(upper.molar, type="b", col="blue", lwd=2)
points(lower.molar, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 9 : premolar
center.premolar = c(tresult.ismb[[1]][10,1],tresult.ismb[[2]][10,1],tresult.ismb[[3]][10,1],tresult.ismb[[4]][10,1],tresult.ismb[[5]][10,1])
width.premolar = 1.96*c(tresult.ismb[[1]][10,2],tresult.ismb[[2]][10,2],tresult.ismb[[3]][10,2],tresult.ismb[[4]][10,2],tresult.ismb[[5]][10,2])
plotCI(center.premolar, y=NULL, uiw = width.premolar, liw = width.premolar, axes=FALSE, xlab = "t_0", ylab = "beta", main = "premolar",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.premolar = center.premolar + width.premolar
lower.premolar = center.premolar - width.premolar
points(upper.premolar, type="b", col="blue", lwd=2)
points(lower.premolar, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 10 : amalgam
center.amalgam = c(tresult.ismb[[1]][11,1],tresult.ismb[[2]][11,1],tresult.ismb[[3]][11,1],tresult.ismb[[4]][11,1],tresult.ismb[[5]][11,1])
width.amalgam = 1.96*c(tresult.ismb[[1]][11,2],tresult.ismb[[2]][11,2],tresult.ismb[[3]][11,2],tresult.ismb[[4]][11,2],tresult.ismb[[5]][11,2])
plotCI(center.amalgam, y=NULL, uiw = width.amalgam, liw = width.amalgam, axes=FALSE, xlab = "t_0", ylab = "beta", main = "amalgam",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.amalgam = center.amalgam + width.amalgam
lower.amalgam = center.amalgam - width.amalgam
points(upper.amalgam, type="b", col="blue", lwd=2)
points(lower.amalgam, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 11 : composite
center.composite = c(tresult.ismb[[1]][12,1],tresult.ismb[[2]][12,1],tresult.ismb[[3]][12,1],tresult.ismb[[4]][12,1],tresult.ismb[[5]][12,1])
width.composite = 1.96*c(tresult.ismb[[1]][12,2],tresult.ismb[[2]][12,2],tresult.ismb[[3]][12,2],tresult.ismb[[4]][12,2],tresult.ismb[[5]][12,2])
plotCI(center.composite, y=NULL, uiw = width.composite, liw = width.composite, axes=FALSE, xlab = "t_0", ylab = "beta", main = "composite",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.composite = center.composite + width.composite
lower.composite = center.composite - width.composite
points(upper.composite, type="b", col="blue", lwd=2)
points(lower.composite, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# Covariate 12 : GIC
center.GIC = c(tresult.ismb[[1]][13,1],tresult.ismb[[2]][13,1],tresult.ismb[[3]][13,1],tresult.ismb[[4]][13,1],tresult.ismb[[5]][13,1])
width.GIC = 1.96*c(tresult.ismb[[1]][13,2],tresult.ismb[[2]][13,2],tresult.ismb[[3]][13,2],tresult.ismb[[4]][13,2],tresult.ismb[[5]][13,2])
plotCI(center.GIC, y=NULL, uiw = width.GIC, liw = width.GIC, axes=FALSE, xlab = "t_0", ylab = "beta", main = "GIC",type = "o")
axis(side=2)
axis(side=1,at=1:5,label=c("0","0.5","1","1.5","2"))
upper.GIC = center.GIC + width.GIC
lower.GIC = center.GIC - width.GIC
points(upper.GIC, type="b", col="blue", lwd=2)
points(lower.GIC, type="b", col="blue", lwd=2)
abline(h=0, col="red",lty=2)

# 
# # Covariate 2 : age^2
# center.age2 = c(result.ismb[[1]][4,1],result.ismb[[2]][4,1],result.ismb[[3]][4,1],result.ismb[[4]][4,1],result.ismb[[5]][4,1])
# width.age2 = 1.96*c(result.ismb[[1]][4,2],result.ismb[[2]][4,2],result.ismb[[3]][4,2],result.ismb[[4]][4,2],result.ismb[[5]][4,2])
# plotCI(center.age2, y=NULL, uiw = width.age2, liw = width.age2, axes=FALSE, xlab = "Quantile", ylab = "beta", main = "age^2",type = "o")
# axis(side=2)
# axis(side=1,at=1:5,label=c("5%","10%","15%","20%","25%"))
# upper.age2 = center.age2 + width.age2
# lower.age2 = center.age2 - width.age2
# points(upper.age2, type="b", col="blue", lwd=2)
# points(lower.age2, type="b", col="blue", lwd=2)
# abline(h=0, col="red",lty=2)
# 
# cox.fit = coxph(Surv(survTime, event) ~ ., data=dummy)
# median.crq = crq.est(dummy$survTime, 13, covariate, dummy$event, 0, 0.5)
# median.ismb = ismb.est(dummy$survTime, 13, covariate, dummy$event, 0, 0.4 ,500)
