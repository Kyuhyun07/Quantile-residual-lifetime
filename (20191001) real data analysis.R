library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)
library(xtable)
library(dplyr)
#LG
load("C:/Users/com/Desktop/대학원/1. 박사 논문/2. 논문준비 (Quantile residual life)/4. Real data/spec.RData")
#Mac
#load("~/Desktop/Q/Paper/Real data/spec.RData")
# Select amalgam restoration
subDat = subset(spec, restorationType==1)
# Select molar
subDat = subset(subDat, toothType=="molar1"|toothType=="molar23")

#### Subset setting ####

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
data1002


#p=cbind(allevent$restorationSize,allevent$restorationType)
#Y=allevent$survTime
#nc=2
#covariate=p
#D=allevent$event
#t_0=0
#Q=0.1
#ne=100

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
# 1. covariate = restorationSize
for (t in 1:5){
  crq.est(data1002$survTime, 1, data1002$restorationSize, data1002$event, 0, t*0.05)
  ismb.est(data1002$survTime, 1, data1002$restorationSize, data1002$event, 0, t*0.05 ,500)
}

# 2. covariate = age
for (t in 1:5){
  crq.est(data1002$survTime, 1, data1002$age, data1002$event, 0, t*0.05)
  ismb.est(data1002$survTime, 1, data1002$age, data1002$event, 0, t*0.05 ,500)
}