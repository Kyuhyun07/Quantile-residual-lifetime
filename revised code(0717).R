library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)
library(xtable)

load("C:/Users/com/Desktop/real data/spec.RData")
## Paper data
subDat <- subset(spec, cohort %in% paste(4:6) & age >= 65)
subDat = subset(spec, first==1 & event==1)

#### Data generation function ####
data.gen<-function(censor, t_0){
  unif = runif(n=1000,min = 0,max = 1)
  sim=matrix(NA,1000,7)
  # Generate C_i
  sim[,2] = runif(1000,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,6] = rbinom(1000,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=25))
  for (q in 1:1000){
    sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
  }
  # Generate Y_i (min(T,C))
  sim[,3] = apply(sim[,1:2], 1, FUN=min)
  sim[,4] = sim[,3]-(t_0)
  sim[,5] = log(sim[,4])
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,7]=I(sim[,1]<sim[,2])
  # Delete data that smaller than t_0
  sim = na.omit(sim)
  # Ordering
  sim = sim[order(sim[,4]),]
  n = nrow(sim)
  sim = as.data.frame(sim)
}
#### Crq method ####
# Y = min(T, C), X = covariate, D = censoring indicator, t_0 = time base
# multi-covariate use cbind()
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
      } else
      {
        data[i,(nc+8)] = data[i-1,(nc+8)]}
    } else {
      data[i,(nc+8)] = data[i,(nc+3)]/data[i,(nc+7)]
    }
  }
  colnames(data)[1:2]=c("Y-t_0", "log(y-t_0)")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[3:(nc+2)]=covar
  #for (j in 3:(nc+2)){
  #  colnames(data)[j]="covariate"
  #}
  colnames(data)[(nc+3):(nc+8)] = c("delta","# at risk","# event","s/d","G_KM","Weight")
  
  cov = as.matrix(data[,3:(nc+2)])
  crq.fit = crq(Surv(data[,2],data[,(nc+3)]) ~ cov,method='Portnoy')
  a = summary(crq.fit,taus = c(0.2, Q))
  beta.sd = matrix(NA, nc+1, 2)
  for (j in 1:(nc+1)){
    beta.sd[j,] = a[2][[1]]$coefficient[j,c(1,4)]
  }
  print(beta.sd)
}

#### ISMB method ####
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
  colnames(data)[1:2]=c("Y-t_0", "log(Y-t_0)")
  for (j in 3:(nc+2)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+3):(nc+8)] = c("delta","# at risk","# event","s/d","G_KM","Weight")
  
  #### Define functions ####
  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(W*X) %*% (pnorm((logT-X%*%beta)/sqrt(diag(X %*% G %*% t(X)))) -1+Q)
  }
  
  #### revised object equation (with eta) ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*W*X) %*% (pnorm((logT-X%*%beta)/sqrt(diag(X %*% G %*% t(X)))) -1+Q)
  }
  
  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,m)),data[,3:(nc+2)]))
  W = data[,(nc+8)]
  logT = data[,2]
  G = diag(1/n, nc+1, nc+1)
  
  betastart = c(1,rep(1,nc))
  is.fit = nleqslv(betastart,objectF, control=list(ftol=1e-5, allowSingular=TRUE))
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(m,1)
      result = t(eta*W*X)%*%(pnorm((logT-X%*%solbeta)/sqrt(diag(X %*% G %*% t(X)))) -1+Q)
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% G %*% t(X)))))*X)%*%(X/sqrt(diag(X %*% G %*% t(X))))
    inva.beta = ginv(a.beta)
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
a=cbind(subDat$restorationSize,subDat$restorationType)
tic()
crq.est(subDat$survTime, 2, a, subDat$event, 0, 0.5)
toc()

tic()
ismb.est(subDat$survTime, 2, a, subDat$event, 0, 0.5,100)
toc()