#### ndata6(200518) : Condition ####
# data size = 200
# Only beta0 effective
# Quantile 75%
# simulation 2000
# eta = 100
# WKM Weight outside, G(Z_i)
# Crq vs rq vs ismb

library(quantreg)
library(survival)
library(nleqslv)
library(xtable)

#### True Beta ####
#beta_0    beta_1
#t_0=0 1.609438 0
#t_0=1 1.410748 0
#t_0=2 1.219403 0
#t_0=3 1.040613 0

#### Find censoring point when t_0=0####
d.exp.beta.initial.0=5
d.k=2
d.r.initial.0=(log(10/2.5)^(1/d.k))/d.exp.beta.initial.0
d.u=runif(n=100000,min = 0,max = 1)
d.x=rbinom(100000,size=1,p=0.5)
d.T<-c()
for (q in 1:100000){
    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.0
}
##Find c which is derermined to achieve 0,10,30,50,70% censoring rate in at t_0=0
i=9
censor.t0=c()
for (m in 1:4){
  while(1){
    d.C<-runif(1000000,0,i)
    dummy = matrix(NA, 1000000, 3)
    dummy[,1] = d.C
    dummy[,2] = d.T
    dummy[,3] = apply(dummy[,1:2], 1, FUN=min)
    newdummy = dummy[dummy[,3]>=0,]
    i=i+0.01
    if(sum(newdummy[,1]<newdummy[,2])/nrow(newdummy)<0.1*(9-2*m)) break
  }
  print(i)
  table(newdummy[,1]<newdummy[,2])[2]/nrow(newdummy)
  censor.t0[m]=i
}

##Find c which is derermined to achieve 0,10,30,50,70% censoring rate in at t_0=1
i=9
censor.t1=c()
for (m in 1:4){
  while(1){
    d.C<-runif(1000000,0,i)
    dummy = matrix(NA, 1000000, 3)
    dummy[,1] = d.C
    dummy[,2] = d.T
    dummy[,3] = apply(dummy[,1:2], 1, FUN=min)
    newdummy = dummy[dummy[,3]>=1,]
    i=i+0.01
    if(sum(newdummy[,1]<newdummy[,2])/nrow(newdummy)<0.1*(9-2*m)) break
  }
  print(i)
  table(newdummy[,1]<newdummy[,2])[2]/nrow(newdummy)
  censor.t1[m]=i
}

##Find c which is derermined to achieve 0,10,30,50,70% censoring rate in at t_0=2
i=9
censor.t2=c()
for (m in 1:4){
  while(1){
    d.C<-runif(1000000,0,i)
    dummy = matrix(NA, 1000000, 3)
    dummy[,1] = d.C
    dummy[,2] = d.T
    dummy[,3] = apply(dummy[,1:2], 1, FUN=min)
    newdummy = dummy[dummy[,3]>=2,]
    i=i+0.01
    if(sum(newdummy[,1]<newdummy[,2])/nrow(newdummy)<0.1*(9-2*m)) break
  }
  print(i)
  table(newdummy[,1]<newdummy[,2])[2]/nrow(newdummy)
  censor.t2[m]=i
}

##Find c which is derermined to achieve 0,10,30,50,70% censoring rate in at t_0=3
i=9
censor.t3=c()
for (m in 1:4){
  while(1){
    d.C<-runif(1000000,0,i)
    dummy = matrix(NA, 1000000, 3)
    dummy[,1] = d.C
    dummy[,2] = d.T
    dummy[,3] = apply(dummy[,1:2], 1, FUN=min)
    newdummy = dummy[dummy[,3]>=3,]
    i=i+0.01
    if(sum(newdummy[,1]<newdummy[,2])/nrow(newdummy)<0.1*(9-2*m)) break
  }
  print(i)
  table(newdummy[,1]<newdummy[,2])[2]/nrow(newdummy)
  censor.t3[m]=i
}

#### Data Generation function ####
data.gen<-function(samplesize, censor){
  sim=matrix(NA,samplesize,5)
  colnames(sim) = c("T","C","Z","X","censored")
  # Generate C_i
  sim[,2] = runif(samplesize,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,4] = rbinom(samplesize,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
  unif = runif(n=samplesize ,min = 0,max = 1)
  for (q in 1:samplesize){
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
  }
  # Generate Y_i (min(T,C))
  sim[,3] = apply(sim[,1:2], 1, FUN=min)
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,5]=I(sim[,1]<sim[,2])
  # Ordering
  sim = sim[order(sim[,3]),]
  n = nrow(sim)
  sim = as.data.frame(sim)
  return(sim)
}

#### Given information(My assumption) ####
exp.beta.initial.0=5
k=2
r.initial.0=(log(10/2.5))^(1/k)/exp.beta.initial.0

#### Example parameter ####
# c.0=5000000
# c.1=123.69
# c.3=41.27
# c.5=23.55
# c.7=14.2
# a<-data.gen(200,c.3)
# Z=a[,3]
# nc=1
# covariate=a[,4]
# D=a[,5]
# t_0=3
# Q=0.25
# ne=100

#### Estimation functions ####
# 1. Crq
crq.est = function(Z, nc, covariate, D, t_0, Q){
  n = length(Z)
  data = matrix(NA, n, nc+4)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data[is.nan(data[,2]),2]=NA
  data = as.data.frame(data)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4)] = c("delta")
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
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.nan(data[,2]),2]=NA
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data)

  # Weight
  # Kaplan-Meier estimator for censoring
  # survfit weight
  # fit = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # WKM surv weight
  # fit = WKM(data[,1],  1-data[,(nc+4)], zc = rep(1,n), w = rep(1,n))
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$times[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # weight using WKM jump
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  data[,(nc+6)] = fit$jump
  m = nrow(data)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")

  cov = as.matrix(data[,4:(nc+3)])
  rq.fit = rq(data[,2] ~ cov, tau=Q, weight = (data[,(nc+6)]*data[,3]))
  rst = summary.rq(rq.fit, se="boot")
  print(rst$coefficients[,c(1,2)])
}

# 3. ISMB & weight-out
ismbwo.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data)

  # Weight
  # Kaplan-Meier estimator for censoring
  # survfit weight
  # fit = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # WKM surv weight
  # fit = WKM(data[,1],  1-data[,(nc+4)], zc = rep(1,n), w = rep(1,n))
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$times[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # weight using WKM jump
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  data[,(nc+6)] = fit$jump
  m = nrow(data)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")

  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(X*I*W) %*% {(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
  }

  #### revised object equation ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*X*I*W) %*% {(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
  }

  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+6)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)

  # Change betastart when real data analysis c(1,rep(1,nc))
  betastart = crq.fit[,1]
  is.fit = nleqslv(betastart, objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I*W) %*% {(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(X*I*W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
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

# 4. ISMB & weight-in
ismbwi.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data)

  # Weight
  # Kaplan-Meier estimator for censoring
  # survfit weight
  # fit = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # WKM surv weight
  # fit = WKM(data[,1],  1-data[,(nc+4)], zc = rep(1,n), w = rep(1,n))
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$times[i],(nc+5)] = fit$surv[i]
  # }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  # }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # weight using WKM jump
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  data[,(nc+6)] = fit$jump*n
  m = nrow(data)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")

  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(X*I) %*% {W*(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
  }

  #### revised object equation ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*X*I) %*% {W*(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
  }

  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+6)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)

  # Change betastart when real data analysis c(1,rep(1,nc))
  betastart = crq.fit[,1]
  is.fit = nleqslv(betastart, objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I) %*% {W*(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(X*I*W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
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

#### Indicator function ####
ind=function(a,b,c){
  if (a>=b&a<=c) {
    result=1
  } else {
    result=0
  }
  print(result)
}

#### Make table for Beta estimation and variance estimation and Coverage ####
table0.crq<-matrix(NA,5,8)
rownames(table0.crq)<-c(0,10,30,50,70)
colnames(table0.crq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.crq<-matrix(NA,5,8)
rownames(table1.crq)<-c(0,10,30,50,70)
colnames(table1.crq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.crq<-matrix(NA,5,8)
rownames(table2.crq)<-c(0,10,30,50,70)
colnames(table2.crq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.crq<-matrix(NA,5,8)
rownames(table3.crq)<-c(0,10,30,50,70)
colnames(table3.crq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.rq<-matrix(NA,5,8)
rownames(table0.rq)<-c(0,10,30,50,70)
colnames(table0.rq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.rq<-matrix(NA,5,8)
rownames(table1.rq)<-c(0,10,30,50,70)
colnames(table1.rq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.rq<-matrix(NA,5,8)
rownames(table2.rq)<-c(0,10,30,50,70)
colnames(table2.rq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.rq<-matrix(NA,5,8)
rownames(table3.rq)<-c(0,10,30,50,70)
colnames(table3.rq)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.is<-matrix(NA,5,8)
rownames(table0.is)<-c(0,10,30,50,70)
colnames(table0.is)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.is<-matrix(NA,5,8)
rownames(table1.is)<-c(0,10,30,50,70)
colnames(table1.is)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.is<-matrix(NA,5,8)
rownames(table2.is)<-c(0,10,30,50,70)
colnames(table2.is)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.is<-matrix(NA,5,8)
rownames(table3.is)<-c(0,10,30,50,70)
colnames(table3.is)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

#### censoring point at t_0=0 ####
c.0=5000000
c.1=36.98
c.3=12.53
c.5=7.47
c.7=4.81
#### t_0=0 & c=0% ####
b0.crq.00 = c()
b0.crq.sd.00 = c()
b1.crq.00 = c()
b1.crq.sd.00 = c()
b0.rq.00 = c()
b0.rq.sd.00 = c()
b1.rq.00 = c()
b1.rq.sd.00 = c()
b0.is.00 = c()
b0.is.sd.00 = c()
b1.is.00 = c()
b1.is.sd.00 = c()
cover.crq.00=matrix(NA,2000,8)
colnames(cover.crq.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.00=matrix(NA,2000,8)
colnames(cover.rq.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.00=matrix(NA,2000,8)
colnames(cover.is.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.crq.00[i] = crq.fit[1,1]
    b0.crq.sd.00[i] = crq.fit[1,2]
    b1.crq.00[i] = crq.fit[2,1]
    b1.crq.sd.00[i] = crq.fit[2,2]
    cover.crq.00[i,1] = crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.00[i,2] = crq.fit[1,1]
    cover.crq.00[i,3] = crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.00[i,4] = ind(1.609438, cover.crq.00[i,1], cover.crq.00[i,3])
    cover.crq.00[i,5] = crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.00[i,6] = crq.fit[2,1]
    cover.crq.00[i,7] = crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.00[i,8] = ind(0, cover.crq.00[i,5], cover.crq.00[i,7])}
    ,error=function(e){
      b0.crq.00[i] = NA
      b0.crq.sd.00[i] = NA
      b1.crq.00[i] = NA
      b1.crq.sd.00[i] = NA
      cover.crq.00[i,1] = NA
      cover.crq.00[i,2] = NA
      cover.crq.00[i,3] = NA
      cover.crq.00[i,4] = NA
      cover.crq.00[i,5] = NA
      cover.crq.00[i,6] = NA
      cover.crq.00[i,7] = NA
      cover.crq.00[i,8] = NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.rq.00[i] = rq.fit[1,1]
    b0.rq.sd.00[i] = rq.fit[1,2]
    b1.rq.00[i] = rq.fit[2,1]
    b1.rq.sd.00[i] = rq.fit[2,2]
    cover.rq.00[i,1] = rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.00[i,2] = rq.fit[1,1]
    cover.rq.00[i,3] = rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.00[i,4] = ind(1.609438, cover.rq.00[i,1], cover.rq.00[i,3])
    cover.rq.00[i,5] = rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.00[i,6] = rq.fit[2,1]
    cover.rq.00[i,7] = rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.00[i,8] = ind(0, cover.rq.00[i,5], cover.rq.00[i,7])}
    ,error=function(e){
      b0.rq.00[i] = NA
      b0.rq.sd.00[i] = NA
      b1.rq.00[i] = NA
      b1.rq.sd.00[i] = NA
      cover.rq.00[i,1] = NA
      cover.rq.00[i,2] = NA
      cover.rq.00[i,3] = NA
      cover.rq.00[i,4] = NA
      cover.rq.00[i,5] = NA
      cover.rq.00[i,6] = NA
      cover.rq.00[i,7] = NA
      cover.rq.00[i,8] = NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 0, 0.75, 100)
    b0.is.00[i] = ismb.fit[1,1]
    b0.is.sd.00[i] = ismb.fit[1,2]
    b1.is.00[i] = ismb.fit[2,1]
    b1.is.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.is.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.00[i,2] = ismb.fit[1,1]
    cover.is.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.00[i,4] = ind(1.609438, cover.is.00[i,1], cover.is.00[i,3])
    cover.is.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.00[i,6] = ismb.fit[2,1]
    cover.is.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.00[i,8] = ind(0., cover.is.00[i,5], cover.is.00[i,7])}
    , error=function(e){
      b0.is.00[i] = NA
      b0.is.sd.00[i] = NA
      b1.is.00[i] = NA
      b1.is.sd.00[i] = NA
      # Coverage
      cover.is.00[i,1] = NA
      cover.is.00[i,2] = NA
      cover.is.00[i,3] = NA
      cover.is.00[i,4] = NA
      cover.is.00[i,5] = NA
      cover.is.00[i,6] = NA
      cover.is.00[i,7] = NA
      cover.is.00[i,8] = NA
    })
}
# Crq beta table
table0.crq[1,1]<-mean(b0.crq.00,na.rm=TRUE)
table0.crq[1,2]<-mean(b0.crq.sd.00,na.rm=TRUE)
table0.crq[1,3]<-sd(b0.crq.00,na.rm=TRUE)
table0.crq[1,4]<-mean(cover.crq.00[,4],na.rm=TRUE)
table0.crq[1,5]<-mean(b1.crq.00,na.rm=TRUE)
table0.crq[1,6]<-mean(b1.crq.sd.00,na.rm=TRUE)
table0.crq[1,7]<-sd(b1.crq.00,na.rm=TRUE)
table0.crq[1,8]<-mean(cover.crq.00[,8],na.rm=TRUE)

# Rq beta table
table0.rq[1,1]<-mean(b0.rq.00,na.rm=TRUE)
table0.rq[1,2]<-mean(b0.rq.sd.00,na.rm=TRUE)
table0.rq[1,3]<-sd(b0.rq.00,na.rm=TRUE)
table0.rq[1,4]<-mean(cover.rq.00[,4],na.rm=TRUE)
table0.rq[1,5]<-mean(b1.rq.00,na.rm=TRUE)
table0.rq[1,6]<-mean(b1.rq.sd.00,na.rm=TRUE)
table0.rq[1,7]<-sd(b1.rq.00,na.rm=TRUE)
table0.rq[1,8]<-mean(cover.rq.00[,8],na.rm=TRUE)

# IS beta table
table0.is[1,1]<-mean(b0.is.00,na.rm=TRUE)
table0.is[1,2]<-mean(b0.is.sd.00,na.rm=TRUE)
table0.is[1,3]<-sd(b0.is.00,na.rm=TRUE)
table0.is[1,4]<-mean(cover.is.00[,4],na.rm=TRUE)
table0.is[1,5]<-mean(b1.is.00,na.rm=TRUE)
table0.is[1,6]<-mean(b1.is.sd.00,na.rm=TRUE)
table0.is[1,7]<-sd(b1.is.00,na.rm=TRUE)
table0.is[1,8]<-mean(cover.is.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.crq.01<-c()
b0.crq.sd.01<-c()
b1.crq.01<-c()
b1.crq.sd.01<-c()
b0.rq.01<-c()
b0.rq.sd.01<-c()
b1.rq.01<-c()
b1.rq.sd.01<-c()
b0.is.01<-c()
b0.is.sd.01<-c()
b1.is.01<-c()
b1.is.sd.01<-c()
cover.crq.01=matrix(NA,2000,8)
colnames(cover.crq.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.01=matrix(NA,2000,8)
colnames(cover.rq.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.01=matrix(NA,2000,8)
colnames(cover.is.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.crq.01[i] = crq.fit[1,1]
    b0.crq.sd.01[i] = crq.fit[1,2]
    b1.crq.01[i] = crq.fit[2,1]
    b1.crq.sd.01[i] = crq.fit[2,2]
    cover.crq.01[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.01[i,2]=crq.fit[1,1]
    cover.crq.01[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.01[i,4]=ind(1.609438,cover.crq.01[i,1],cover.crq.01[i,3])
    cover.crq.01[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.01[i,6]=crq.fit[2,1]
    cover.crq.01[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.01[i,8]=ind(0,cover.crq.01[i,5],cover.crq.01[i,7])}
    ,error=function(e){
      b0.crq.01[i] = NA
      b0.crq.sd.01[i] = NA
      b1.crq.01[i] = NA
      b1.crq.sd.01[i] = NA
      cover.crq.01[i,1]=NA
      cover.crq.01[i,2]=NA
      cover.crq.01[i,3]=NA
      cover.crq.01[i,4]=NA
      cover.crq.01[i,5]=NA
      cover.crq.01[i,6]=NA
      cover.crq.01[i,7]=NA
      cover.crq.01[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.rq.01[i] = rq.fit[1,1]
    b0.rq.sd.01[i] = rq.fit[1,2]
    b1.rq.01[i] = rq.fit[2,1]
    b1.rq.sd.01[i] = rq.fit[2,2]
    cover.rq.01[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.01[i,2]=rq.fit[1,1]
    cover.rq.01[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.01[i,4]=ind(1.609438,cover.rq.01[i,1],cover.rq.01[i,3])
    cover.rq.01[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.01[i,6]=rq.fit[2,1]
    cover.rq.01[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.01[i,8]=ind(0,cover.rq.01[i,5],cover.rq.01[i,7])}
    ,error=function(e){
      b0.rq.01[i] = NA
      b0.rq.sd.01[i] = NA
      b1.rq.01[i] = NA
      b1.rq.sd.01[i] = NA
      cover.rq.01[i,1]=NA
      cover.rq.01[i,2]=NA
      cover.rq.01[i,3]=NA
      cover.rq.01[i,4]=NA
      cover.rq.01[i,5]=NA
      cover.rq.01[i,6]=NA
      cover.rq.01[i,7]=NA
      cover.rq.01[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 0, 0.75, 100)
    b0.is.01[i] = ismb.fit[1,1]
    b0.is.sd.01[i] = ismb.fit[1,2]
    b1.is.01[i] = ismb.fit[2,1]
    b1.is.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.is.01[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.01[i,2]=ismb.fit[1,1]
    cover.is.01[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.01[i,4]=ind(1.609438,cover.is.01[i,1],cover.is.01[i,3])
    cover.is.01[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.01[i,6]=ismb.fit[2,1]
    cover.is.01[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.01[i,8]=ind(0,cover.is.01[i,5],cover.is.01[i,7])}
    , error=function(e){
      b0.is.01[i] = NA
      b0.is.sd.01[i] = NA
      b1.is.01[i] = NA
      b1.is.sd.01[i] = NA
      # Coverage
      cover.is.01[i,1]=NA
      cover.is.01[i,2]=NA
      cover.is.01[i,3]=NA
      cover.is.01[i,4]=NA
      cover.is.01[i,5]=NA
      cover.is.01[i,6]=NA
      cover.is.01[i,7]=NA
      cover.is.01[i,8]=NA
    })
}

# Crq beta table
table0.crq[2,1]<-mean(b0.crq.01,na.rm=TRUE)
table0.crq[2,2]<-mean(b0.crq.sd.01,na.rm=TRUE)
table0.crq[2,3]<-sd(b0.crq.01,na.rm=TRUE)
table0.crq[2,4]<-mean(cover.crq.01[,4],na.rm=TRUE)
table0.crq[2,5]<-mean(b1.crq.01,na.rm=TRUE)
table0.crq[2,6]<-mean(b1.crq.sd.01,na.rm=TRUE)
table0.crq[2,7]<-sd(b1.crq.01,na.rm=TRUE)
table0.crq[2,8]<-mean(cover.crq.01[,8],na.rm=TRUE)

# rq beta table
table0.rq[2,1]<-mean(b0.rq.01,na.rm=TRUE)
table0.rq[2,2]<-mean(b0.rq.sd.01,na.rm=TRUE)
table0.rq[2,3]<-sd(b0.rq.01,na.rm=TRUE)
table0.rq[2,4]<-mean(cover.rq.01[,4],na.rm=TRUE)
table0.rq[2,5]<-mean(b1.rq.01,na.rm=TRUE)
table0.rq[2,6]<-mean(b1.rq.sd.01,na.rm=TRUE)
table0.rq[2,7]<-sd(b1.rq.01,na.rm=TRUE)
table0.rq[2,8]<-mean(cover.rq.01[,8],na.rm=TRUE)

# IS beta table
table0.is[2,1]<-mean(b0.is.01,na.rm=TRUE)
table0.is[2,2]<-mean(b0.is.sd.01,na.rm=TRUE)
table0.is[2,3]<-sd(b0.is.01,na.rm=TRUE)
table0.is[2,4]<-mean(cover.is.01[,4],na.rm=TRUE)
table0.is[2,5]<-mean(b1.is.01,na.rm=TRUE)
table0.is[2,6]<-mean(b1.is.sd.01,na.rm=TRUE)
table0.is[2,7]<-sd(b1.is.01,na.rm=TRUE)
table0.is[2,8]<-mean(cover.is.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.crq.03<-c()
b0.crq.sd.03<-c()
b1.crq.03<-c()
b1.crq.sd.03<-c()
b0.rq.03<-c()
b0.rq.sd.03<-c()
b1.rq.03<-c()
b1.rq.sd.03<-c()
b0.is.03<-c()
b0.is.sd.03<-c()
b1.is.03<-c()
b1.is.sd.03<-c()
cover.crq.03=matrix(NA,2000,8)
colnames(cover.crq.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.03=matrix(NA,2000,8)
colnames(cover.rq.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.03=matrix(NA,2000,8)
colnames(cover.is.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.crq.03[i] = crq.fit[1,1]
    b0.crq.sd.03[i] = crq.fit[1,2]
    b1.crq.03[i] = crq.fit[2,1]
    b1.crq.sd.03[i] = crq.fit[2,2]
    cover.crq.03[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.03[i,2]=crq.fit[1,1]
    cover.crq.03[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.03[i,4]=ind(1.609438,cover.crq.03[i,1],cover.crq.03[i,3])
    cover.crq.03[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.03[i,6]=crq.fit[2,1]
    cover.crq.03[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.03[i,8]=ind(0,cover.crq.03[i,5],cover.crq.03[i,7])}
    ,error=function(e){
      b0.crq.03[i] = NA
      b0.crq.sd.03[i] = NA
      b1.crq.03[i] = NA
      b1.crq.sd.03[i] = NA
      cover.crq.03[i,1]=NA
      cover.crq.03[i,2]=NA
      cover.crq.03[i,3]=NA
      cover.crq.03[i,4]=NA
      cover.crq.03[i,5]=NA
      cover.crq.03[i,6]=NA
      cover.crq.03[i,7]=NA
      cover.crq.03[i,8]=NA
    })
  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.rq.03[i] = rq.fit[1,1]
    b0.rq.sd.03[i] = rq.fit[1,2]
    b1.rq.03[i] = rq.fit[2,1]
    b1.rq.sd.03[i] = rq.fit[2,2]
    cover.rq.03[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.03[i,2]=rq.fit[1,1]
    cover.rq.03[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.03[i,4]=ind(1.609438,cover.rq.03[i,1],cover.rq.03[i,3])
    cover.rq.03[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.03[i,6]=rq.fit[2,1]
    cover.rq.03[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.03[i,8]=ind(0,cover.rq.03[i,5],cover.rq.03[i,7])}
    ,error=function(e){
      b0.rq.03[i] = NA
      b0.rq.sd.03[i] = NA
      b1.rq.03[i] = NA
      b1.rq.sd.03[i] = NA
      cover.rq.03[i,1]=NA
      cover.rq.03[i,2]=NA
      cover.rq.03[i,3]=NA
      cover.rq.03[i,4]=NA
      cover.rq.03[i,5]=NA
      cover.rq.03[i,6]=NA
      cover.rq.03[i,7]=NA
      cover.rq.03[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 0, 0.75, 100)
    b0.is.03[i] = ismb.fit[1,1]
    b0.is.sd.03[i] = ismb.fit[1,2]
    b1.is.03[i] = ismb.fit[2,1]
    b1.is.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.is.03[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.03[i,2]=ismb.fit[1,1]
    cover.is.03[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.03[i,4]=ind(1.609438,cover.is.03[i,1],cover.is.03[i,3])
    cover.is.03[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.03[i,6]=ismb.fit[2,1]
    cover.is.03[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.03[i,8]=ind(0,cover.is.03[i,5],cover.is.03[i,7])}
    , error=function(e){
      b0.is.03[i] = NA
      b0.is.sd.03[i] = NA
      b1.is.03[i] = NA
      b1.is.sd.03[i] = NA
      # Coverage
      cover.is.03[i,1]=NA
      cover.is.03[i,2]=NA
      cover.is.03[i,3]=NA
      cover.is.03[i,4]=NA
      cover.is.03[i,5]=NA
      cover.is.03[i,6]=NA
      cover.is.03[i,7]=NA
      cover.is.03[i,8]=NA
    })
}
# Crq beta table
table0.crq[3,1]<-mean(b0.crq.03,na.rm=TRUE)
table0.crq[3,2]<-mean(b0.crq.sd.03,na.rm=TRUE)
table0.crq[3,3]<-sd(b0.crq.03,na.rm=TRUE)
table0.crq[3,4]<-mean(cover.crq.03[,4],na.rm=TRUE)
table0.crq[3,5]<-mean(b1.crq.03,na.rm=TRUE)
table0.crq[3,6]<-mean(b1.crq.sd.03,na.rm=TRUE)
table0.crq[3,7]<-sd(b1.crq.03,na.rm=TRUE)
table0.crq[3,8]<-mean(cover.crq.03[,8],na.rm=TRUE)

# rq beta table
table0.rq[3,1]<-mean(b0.rq.03,na.rm=TRUE)
table0.rq[3,2]<-mean(b0.rq.sd.03,na.rm=TRUE)
table0.rq[3,3]<-sd(b0.rq.03,na.rm=TRUE)
table0.rq[3,4]<-mean(cover.rq.03[,4],na.rm=TRUE)
table0.rq[3,5]<-mean(b1.rq.03,na.rm=TRUE)
table0.rq[3,6]<-mean(b1.rq.sd.03,na.rm=TRUE)
table0.rq[3,7]<-sd(b1.rq.03,na.rm=TRUE)
table0.rq[3,8]<-mean(cover.rq.03[,8],na.rm=TRUE)

# IS beta table
table0.is[3,1]<-mean(b0.is.03,na.rm=TRUE)
table0.is[3,2]<-mean(b0.is.sd.03,na.rm=TRUE)
table0.is[3,3]<-sd(b0.is.03,na.rm=TRUE)
table0.is[3,4]<-mean(cover.is.03[,4],na.rm=TRUE)
table0.is[3,5]<-mean(b1.is.03,na.rm=TRUE)
table0.is[3,6]<-mean(b1.is.sd.03,na.rm=TRUE)
table0.is[3,7]<-sd(b1.is.03,na.rm=TRUE)
table0.is[3,8]<-mean(cover.is.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.crq.05<-c()
b0.crq.sd.05<-c()
b1.crq.05<-c()
b1.crq.sd.05<-c()
b0.rq.05<-c()
b0.rq.sd.05<-c()
b1.rq.05<-c()
b1.rq.sd.05<-c()
b0.is.05<-c()
b0.is.sd.05<-c()
b1.is.05<-c()
b1.is.sd.05<-c()
cover.crq.05=matrix(NA,2000,8)
colnames(cover.crq.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.05=matrix(NA,2000,8)
colnames(cover.rq.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.05=matrix(NA,2000,8)
colnames(cover.is.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.crq.05[i] = crq.fit[1,1]
    b0.crq.sd.05[i] = crq.fit[1,2]
    b1.crq.05[i] = crq.fit[2,1]
    b1.crq.sd.05[i] = crq.fit[2,2]
    cover.crq.05[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.05[i,2]=crq.fit[1,1]
    cover.crq.05[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.05[i,4]=ind(1.609438, cover.crq.05[i,1], cover.crq.05[i,3])
    cover.crq.05[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.05[i,6]=crq.fit[2,1]
    cover.crq.05[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.05[i,8]=ind(0, cover.crq.05[i,5], cover.crq.05[i,7])}
    ,error=function(e){
      b0.crq.05[i] = NA
      b0.crq.sd.05[i] = NA
      b1.crq.05[i] = NA
      b1.crq.sd.05[i] = NA
      cover.crq.05[i,1]=NA
      cover.crq.05[i,2]=NA
      cover.crq.05[i,3]=NA
      cover.crq.05[i,4]=NA
      cover.crq.05[i,5]=NA
      cover.crq.05[i,6]=NA
      cover.crq.05[i,7]=NA
      cover.crq.05[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.rq.05[i] = rq.fit[1,1]
    b0.rq.sd.05[i] = rq.fit[1,2]
    b1.rq.05[i] = rq.fit[2,1]
    b1.rq.sd.05[i] = rq.fit[2,2]
    cover.rq.05[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.05[i,2]=rq.fit[1,1]
    cover.rq.05[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.05[i,4]=ind(1.609438, cover.rq.05[i,1], cover.rq.05[i,3])
    cover.rq.05[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.05[i,6]=rq.fit[2,1]
    cover.rq.05[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.05[i,8]=ind(0, cover.rq.05[i,5], cover.rq.05[i,7])}
    ,error=function(e){
      b0.rq.05[i] = NA
      b0.rq.sd.05[i] = NA
      b1.rq.05[i] = NA
      b1.rq.sd.05[i] = NA
      cover.rq.05[i,1]=NA
      cover.rq.05[i,2]=NA
      cover.rq.05[i,3]=NA
      cover.rq.05[i,4]=NA
      cover.rq.05[i,5]=NA
      cover.rq.05[i,6]=NA
      cover.rq.05[i,7]=NA
      cover.rq.05[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 0, 0.75, 100)
    b0.is.05[i] = ismb.fit[1,1]
    b0.is.sd.05[i] = ismb.fit[1,2]
    b1.is.05[i] = ismb.fit[2,1]
    b1.is.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.is.05[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.05[i,2]=ismb.fit[1,1]
    cover.is.05[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.05[i,4]=ind(1.609438,cover.is.05[i,1],cover.is.05[i,3])
    cover.is.05[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.05[i,6]=ismb.fit[2,1]
    cover.is.05[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.05[i,8]=ind(0,cover.is.05[i,5],cover.is.05[i,7])}
    , error=function(e){
      b0.is.05[i] = NA
      b0.is.sd.05[i] = NA
      b1.is.05[i] = NA
      b1.is.sd.05[i] = NA
      # Coverage
      cover.is.05[i,1]=NA
      cover.is.05[i,2]=NA
      cover.is.05[i,3]=NA
      cover.is.05[i,4]=NA
      cover.is.05[i,5]=NA
      cover.is.05[i,6]=NA
      cover.is.05[i,7]=NA
      cover.is.05[i,8]=NA
    })
}

# Crq beta table
table0.crq[4,1]<-mean(b0.crq.05,na.rm=TRUE)
table0.crq[4,2]<-mean(b0.crq.sd.05,na.rm=TRUE)
table0.crq[4,3]<-sd(b0.crq.05,na.rm=TRUE)
table0.crq[4,4]<-mean(cover.crq.05[,4],na.rm=TRUE)
table0.crq[4,5]<-mean(b1.crq.05,na.rm=TRUE)
table0.crq[4,6]<-mean(b1.crq.sd.05,na.rm=TRUE)
table0.crq[4,7]<-sd(b1.crq.05,na.rm=TRUE)
table0.crq[4,8]<-mean(cover.crq.05[,8],na.rm=TRUE)

# rq beta table
table0.rq[4,1]<-mean(b0.rq.05,na.rm=TRUE)
table0.rq[4,2]<-mean(b0.rq.sd.05,na.rm=TRUE)
table0.rq[4,3]<-sd(b0.rq.05,na.rm=TRUE)
table0.rq[4,4]<-mean(cover.rq.05[,4],na.rm=TRUE)
table0.rq[4,5]<-mean(b1.rq.05,na.rm=TRUE)
table0.rq[4,6]<-mean(b1.rq.sd.05,na.rm=TRUE)
table0.rq[4,7]<-sd(b1.rq.05,na.rm=TRUE)
table0.rq[4,8]<-mean(cover.rq.05[,8],na.rm=TRUE)

# IS beta table
table0.is[4,1]<-mean(b0.is.05,na.rm=TRUE)
table0.is[4,2]<-mean(b0.is.sd.05,na.rm=TRUE)
table0.is[4,3]<-sd(b0.is.05,na.rm=TRUE)
table0.is[4,4]<-mean(cover.is.05[,4],na.rm=TRUE)
table0.is[4,5]<-mean(b1.is.05,na.rm=TRUE)
table0.is[4,6]<-mean(b1.is.sd.05,na.rm=TRUE)
table0.is[4,7]<-sd(b1.is.05,na.rm=TRUE)
table0.is[4,8]<-mean(cover.is.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.rq.07<-c()
b0.rq.sd.07<-c()
b1.rq.07<-c()
b1.rq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
cover.crq.07=matrix(NA,2000,8)
colnames(cover.crq.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.07=matrix(NA,2000,8)
colnames(cover.rq.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.07=matrix(NA,2000,8)
colnames(cover.is.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.crq.07[i] = crq.fit[1,1]
    b0.crq.sd.07[i] = crq.fit[1,2]
    b1.crq.07[i] = crq.fit[2,1]
    b1.crq.sd.07[i] = crq.fit[2,2]
    cover.crq.07[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.07[i,2]=crq.fit[1,1]
    cover.crq.07[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.07[i,4]=ind(1.609438,cover.crq.07[i,1],cover.crq.07[i,3])
    cover.crq.07[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.07[i,6]=crq.fit[2,1]
    cover.crq.07[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.07[i,8]=ind(0,cover.crq.07[i,5],cover.crq.07[i,7])}
    ,error=function(e){
      b0.crq.07[i] = NA
      b0.crq.sd.07[i] = NA
      b1.crq.07[i] = NA
      b1.crq.sd.07[i] = NA
      cover.crq.07[i,1]=NA
      cover.crq.07[i,2]=NA
      cover.crq.07[i,3]=NA
      cover.crq.07[i,4]=NA
      cover.crq.07[i,5]=NA
      cover.crq.07[i,6]=NA
      cover.crq.07[i,7]=NA
      cover.crq.07[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 0, 0.75)
    b0.rq.07[i] = rq.fit[1,1]
    b0.rq.sd.07[i] = rq.fit[1,2]
    b1.rq.07[i] = rq.fit[2,1]
    b1.rq.sd.07[i] = rq.fit[2,2]
    cover.rq.07[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.07[i,2]=rq.fit[1,1]
    cover.rq.07[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.07[i,4]=ind(1.609438,cover.rq.07[i,1],cover.rq.07[i,3])
    cover.rq.07[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.07[i,6]=rq.fit[2,1]
    cover.rq.07[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.07[i,8]=ind(0,cover.rq.07[i,5],cover.rq.07[i,7])}
    ,error=function(e){
      b0.rq.07[i] = NA
      b0.rq.sd.07[i] = NA
      b1.rq.07[i] = NA
      b1.rq.sd.07[i] = NA
      cover.rq.07[i,1]=NA
      cover.rq.07[i,2]=NA
      cover.rq.07[i,3]=NA
      cover.rq.07[i,4]=NA
      cover.rq.07[i,5]=NA
      cover.rq.07[i,6]=NA
      cover.rq.07[i,7]=NA
      cover.rq.07[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 0, 0.75, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.is.07[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.07[i,2]=ismb.fit[1,1]
    cover.is.07[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.07[i,4]=ind(1.609438,cover.is.07[i,1],cover.is.07[i,3])
    cover.is.07[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.07[i,6]=ismb.fit[2,1]
    cover.is.07[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.07[i,8]=ind(0,cover.is.07[i,5],cover.is.07[i,7])}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
      # Coverage
      cover.is.07[i,1]=NA
      cover.is.07[i,2]=NA
      cover.is.07[i,3]=NA
      cover.is.07[i,4]=NA
      cover.is.07[i,5]=NA
      cover.is.07[i,6]=NA
      cover.is.07[i,7]=NA
      cover.is.07[i,8]=NA
    })
}
# Crq beta table
table0.crq[5,1]<-mean(b0.crq.07,na.rm=TRUE)
table0.crq[5,2]<-mean(b0.crq.sd.07,na.rm=TRUE)
table0.crq[5,3]<-sd(b0.crq.07,na.rm=TRUE)
table0.crq[5,4]<-mean(cover.crq.07[,4],na.rm=TRUE)
table0.crq[5,5]<-mean(b1.crq.07,na.rm=TRUE)
table0.crq[5,6]<-mean(b1.crq.sd.07,na.rm=TRUE)
table0.crq[5,7]<-sd(b1.crq.07,na.rm=TRUE)
table0.crq[5,8]<-mean(cover.crq.07[,8],na.rm=TRUE)

# rq beta table
table0.rq[5,1]<-mean(b0.rq.07,na.rm=TRUE)
table0.rq[5,2]<-mean(b0.rq.sd.07,na.rm=TRUE)
table0.rq[5,3]<-sd(b0.rq.07,na.rm=TRUE)
table0.rq[5,4]<-mean(cover.rq.07[,4],na.rm=TRUE)
table0.rq[5,5]<-mean(b1.rq.07,na.rm=TRUE)
table0.rq[5,6]<-mean(b1.rq.sd.07,na.rm=TRUE)
table0.rq[5,7]<-sd(b1.rq.07,na.rm=TRUE)
table0.rq[5,8]<-mean(cover.rq.07[,8],na.rm=TRUE)

# IS beta table
table0.is[5,1]<-mean(b0.is.07,na.rm=TRUE)
table0.is[5,2]<-mean(b0.is.sd.07,na.rm=TRUE)
table0.is[5,3]<-sd(b0.is.07,na.rm=TRUE)
table0.is[5,4]<-mean(cover.is.07[,4],na.rm=TRUE)
table0.is[5,5]<-mean(b1.is.07,na.rm=TRUE)
table0.is[5,6]<-mean(b1.is.sd.07,na.rm=TRUE)
table0.is[5,7]<-sd(b1.is.07,na.rm=TRUE)
table0.is[5,8]<-mean(cover.is.07[,8],na.rm=TRUE)

#### censoring point at t_0=1 ####
c.0=5000000
c.1=30.1
c.3=10.82
c.5=6.68
c.7=4.39
#### t_0=1 & c=0% ####
b0.crq.10 = c()
b0.crq.sd.10 = c()
b1.crq.10 = c()
b1.crq.sd.10 = c()
b0.rq.10 = c()
b0.rq.sd.10 = c()
b1.rq.10 = c()
b1.rq.sd.10 = c()
b0.is.10 = c()
b0.is.sd.10 = c()
b1.is.10 = c()
b1.is.sd.10 = c()
delnumt_1=c()
cover.crq.10=matrix(NA,2000,8)
colnames(cover.crq.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.10=matrix(NA,2000,8)
colnames(cover.rq.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.10=matrix(NA,2000,8)
colnames(cover.is.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  delnumt_1 = c(delnumt_1, sum(a$Z<1))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.crq.10[i] = crq.fit[1,1]
    b0.crq.sd.10[i] = crq.fit[1,2]
    b1.crq.10[i] = crq.fit[2,1]
    b1.crq.sd.10[i] = crq.fit[2,2]
    cover.crq.10[i,1] = crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.10[i,2] = crq.fit[1,1]
    cover.crq.10[i,3] = crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.10[i,4] = ind(1.410748, cover.crq.10[i,1], cover.crq.10[i,3])
    cover.crq.10[i,5] = crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.10[i,6] = crq.fit[2,1]
    cover.crq.10[i,7] = crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.10[i,8] = ind(0, cover.crq.10[i,5], cover.crq.10[i,7])}
    ,error=function(e){
      b0.crq.10[i] = NA
      b0.crq.sd.10[i] = NA
      b1.crq.10[i] = NA
      b1.crq.sd.10[i] = NA
      cover.crq.10[i,1] = NA
      cover.crq.10[i,2] = NA
      cover.crq.10[i,3] = NA
      cover.crq.10[i,4] = NA
      cover.crq.10[i,5] = NA
      cover.crq.10[i,6] = NA
      cover.crq.10[i,7] = NA
      cover.crq.10[i,8] = NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.rq.10[i] = rq.fit[1,1]
    b0.rq.sd.10[i] = rq.fit[1,2]
    b1.rq.10[i] = rq.fit[2,1]
    b1.rq.sd.10[i] = rq.fit[2,2]
    cover.rq.10[i,1] = rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.10[i,2] = rq.fit[1,1]
    cover.rq.10[i,3] = rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.10[i,4] = ind(1.410748, cover.rq.10[i,1], cover.rq.10[i,3])
    cover.rq.10[i,5] = rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.10[i,6] = rq.fit[2,1]
    cover.rq.10[i,7] = rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.10[i,8] = ind(0, cover.rq.10[i,5], cover.rq.10[i,7])}
    ,error=function(e){
      b0.crq.10[i] = NA
      b0.crq.sd.10[i] = NA
      b1.crq.10[i] = NA
      b1.crq.sd.10[i] = NA
      cover.crq.10[i,1] = NA
      cover.crq.10[i,2] = NA
      cover.crq.10[i,3] = NA
      cover.crq.10[i,4] = NA
      cover.crq.10[i,5] = NA
      cover.crq.10[i,6] = NA
      cover.crq.10[i,7] = NA
      cover.crq.10[i,8] = NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 1, 0.75, 100)
    b0.is.10[i] = ismb.fit[1,1]
    b0.is.sd.10[i] = ismb.fit[1,2]
    b1.is.10[i] = ismb.fit[2,1]
    b1.is.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.is.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.10[i,2] = ismb.fit[1,1]
    cover.is.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.10[i,4] = ind(1.410748, cover.is.10[i,1], cover.is.10[i,3])
    cover.is.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.10[i,6] = ismb.fit[2,1]
    cover.is.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.10[i,8] = ind(0, cover.is.10[i,5], cover.is.10[i,7])}
    , error=function(e){
      b0.is.10[i] = NA
      b0.is.sd.10[i] = NA
      b1.is.10[i] = NA
      b1.is.sd.10[i] = NA
      # Coverage
      cover.is.10[i,1] = NA
      cover.is.10[i,2] = NA
      cover.is.10[i,3] = NA
      cover.is.10[i,4] = NA
      cover.is.10[i,5] = NA
      cover.is.10[i,6] = NA
      cover.is.10[i,7] = NA
      cover.is.10[i,8] = NA
    })
}
# Crq beta table
table1.crq[1,1]<-mean(b0.crq.10,na.rm=TRUE)
table1.crq[1,2]<-mean(b0.crq.sd.10,na.rm=TRUE)
table1.crq[1,3]<-sd(b0.crq.10,na.rm=TRUE)
table1.crq[1,4]<-mean(cover.crq.10[,4],na.rm=TRUE)
table1.crq[1,5]<-mean(b1.crq.10,na.rm=TRUE)
table1.crq[1,6]<-mean(b1.crq.sd.10,na.rm=TRUE)
table1.crq[1,7]<-sd(b1.crq.10,na.rm=TRUE)
table1.crq[1,8]<-mean(cover.crq.10[,8],na.rm=TRUE)

# Rq beta table
table1.rq[1,1]<-mean(b0.rq.10,na.rm=TRUE)
table1.rq[1,2]<-mean(b0.rq.sd.10,na.rm=TRUE)
table1.rq[1,3]<-sd(b0.rq.10,na.rm=TRUE)
table1.rq[1,4]<-mean(cover.rq.10[,4],na.rm=TRUE)
table1.rq[1,5]<-mean(b1.rq.10,na.rm=TRUE)
table1.rq[1,6]<-mean(b1.rq.sd.10,na.rm=TRUE)
table1.rq[1,7]<-sd(b1.rq.10,na.rm=TRUE)
table1.rq[1,8]<-mean(cover.rq.10[,8],na.rm=TRUE)

# IS beta table
table1.is[1,1]<-mean(b0.is.10,na.rm=TRUE)
table1.is[1,2]<-mean(b0.is.sd.10,na.rm=TRUE)
table1.is[1,3]<-sd(b0.is.10,na.rm=TRUE)
table1.is[1,4]<-mean(cover.is.10[,4],na.rm=TRUE)
table1.is[1,5]<-mean(b1.is.10,na.rm=TRUE)
table1.is[1,6]<-mean(b1.is.sd.10,na.rm=TRUE)
table1.is[1,7]<-sd(b1.is.10,na.rm=TRUE)
table1.is[1,8]<-mean(cover.is.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.crq.11<-c()
b0.crq.sd.11<-c()
b1.crq.11<-c()
b1.crq.sd.11<-c()
b0.rq.11<-c()
b0.rq.sd.11<-c()
b1.rq.11<-c()
b1.rq.sd.11<-c()
b0.is.11<-c()
b0.is.sd.11<-c()
b1.is.11<-c()
b1.is.sd.11<-c()
cover.crq.11=matrix(NA,2000,8)
colnames(cover.crq.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.11=matrix(NA,2000,8)
colnames(cover.rq.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.11=matrix(NA,2000,8)
colnames(cover.is.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  delnumt_1 = c(delnumt_1, sum(a$Z<1))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.crq.11[i] = crq.fit[1,1]
    b0.crq.sd.11[i] = crq.fit[1,2]
    b1.crq.11[i] = crq.fit[2,1]
    b1.crq.sd.11[i] = crq.fit[2,2]
    cover.crq.11[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.11[i,2]=crq.fit[1,1]
    cover.crq.11[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.11[i,4]=ind(1.410748,cover.crq.11[i,1],cover.crq.11[i,3])
    cover.crq.11[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.11[i,6]=crq.fit[2,1]
    cover.crq.11[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.11[i,8]=ind(0,cover.crq.11[i,5],cover.crq.11[i,7])}
    ,error=function(e){
      b0.crq.11[i] = NA
      b0.crq.sd.11[i] = NA
      b1.crq.11[i] = NA
      b1.crq.sd.11[i] = NA
      cover.crq.11[i,1]=NA
      cover.crq.11[i,2]=NA
      cover.crq.11[i,3]=NA
      cover.crq.11[i,4]=NA
      cover.crq.11[i,5]=NA
      cover.crq.11[i,6]=NA
      cover.crq.11[i,7]=NA
      cover.crq.11[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.rq.11[i] = rq.fit[1,1]
    b0.rq.sd.11[i] = rq.fit[1,2]
    b1.rq.11[i] = rq.fit[2,1]
    b1.rq.sd.11[i] = rq.fit[2,2]
    cover.rq.11[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.11[i,2]=rq.fit[1,1]
    cover.rq.11[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.11[i,4]=ind(1.410748,cover.rq.11[i,1],cover.rq.11[i,3])
    cover.rq.11[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.11[i,6]=rq.fit[2,1]
    cover.rq.11[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.11[i,8]=ind(0,cover.rq.11[i,5],cover.rq.11[i,7])}
    ,error=function(e){
      b0.rq.11[i] = NA
      b0.rq.sd.11[i] = NA
      b1.rq.11[i] = NA
      b1.rq.sd.11[i] = NA
      cover.rq.11[i,1]=NA
      cover.rq.11[i,2]=NA
      cover.rq.11[i,3]=NA
      cover.rq.11[i,4]=NA
      cover.rq.11[i,5]=NA
      cover.rq.11[i,6]=NA
      cover.rq.11[i,7]=NA
      cover.rq.11[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 1, 0.75, 100)
    b0.is.11[i] = ismb.fit[1,1]
    b0.is.sd.11[i] = ismb.fit[1,2]
    b1.is.11[i] = ismb.fit[2,1]
    b1.is.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.is.11[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.11[i,2]=ismb.fit[1,1]
    cover.is.11[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.11[i,4]=ind(1.410748,cover.is.11[i,1],cover.is.11[i,3])
    cover.is.11[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.11[i,6]=ismb.fit[2,1]
    cover.is.11[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.11[i,8]=ind(0,cover.is.11[i,5],cover.is.11[i,7])}
    , error=function(e){
      b0.is.11[i] = NA
      b0.is.sd.11[i] = NA
      b1.is.11[i] = NA
      b1.is.sd.11[i] = NA
      # Coverage
      cover.is.11[i,1]=NA
      cover.is.11[i,2]=NA
      cover.is.11[i,3]=NA
      cover.is.11[i,4]=NA
      cover.is.11[i,5]=NA
      cover.is.11[i,6]=NA
      cover.is.11[i,7]=NA
      cover.is.11[i,8]=NA
    })
}

# Crq beta table
table1.crq[2,1]<-mean(b0.crq.11,na.rm=TRUE)
table1.crq[2,2]<-mean(b0.crq.sd.11,na.rm=TRUE)
table1.crq[2,3]<-sd(b0.crq.11,na.rm=TRUE)
table1.crq[2,4]<-mean(cover.crq.11[,4],na.rm=TRUE)
table1.crq[2,5]<-mean(b1.crq.11,na.rm=TRUE)
table1.crq[2,6]<-mean(b1.crq.sd.11,na.rm=TRUE)
table1.crq[2,7]<-sd(b1.crq.11,na.rm=TRUE)
table1.crq[2,8]<-mean(cover.crq.11[,8],na.rm=TRUE)

# rq beta table
table1.rq[2,1]<-mean(b0.rq.11,na.rm=TRUE)
table1.rq[2,2]<-mean(b0.rq.sd.11,na.rm=TRUE)
table1.rq[2,3]<-sd(b0.rq.11,na.rm=TRUE)
table1.rq[2,4]<-mean(cover.rq.11[,4],na.rm=TRUE)
table1.rq[2,5]<-mean(b1.rq.11,na.rm=TRUE)
table1.rq[2,6]<-mean(b1.rq.sd.11,na.rm=TRUE)
table1.rq[2,7]<-sd(b1.rq.11,na.rm=TRUE)
table1.rq[2,8]<-mean(cover.rq.11[,8],na.rm=TRUE)

# IS beta table
table1.is[2,1]<-mean(b0.is.11,na.rm=TRUE)
table1.is[2,2]<-mean(b0.is.sd.11,na.rm=TRUE)
table1.is[2,3]<-sd(b0.is.11,na.rm=TRUE)
table1.is[2,4]<-mean(cover.is.11[,4],na.rm=TRUE)
table1.is[2,5]<-mean(b1.is.11,na.rm=TRUE)
table1.is[2,6]<-mean(b1.is.sd.11,na.rm=TRUE)
table1.is[2,7]<-sd(b1.is.11,na.rm=TRUE)
table1.is[2,8]<-mean(cover.is.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.crq.13<-c()
b0.crq.sd.13<-c()
b1.crq.13<-c()
b1.crq.sd.13<-c()
b0.rq.13<-c()
b0.rq.sd.13<-c()
b1.rq.13<-c()
b1.rq.sd.13<-c()
b0.is.13<-c()
b0.is.sd.13<-c()
b1.is.13<-c()
b1.is.sd.13<-c()
cover.crq.13=matrix(NA,2000,8)
colnames(cover.crq.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.13=matrix(NA,2000,8)
colnames(cover.rq.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.13=matrix(NA,2000,8)
colnames(cover.is.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  delnumt_1 = c(delnumt_1, sum(a$Z<1))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.crq.13[i] = crq.fit[1,1]
    b0.crq.sd.13[i] = crq.fit[1,2]
    b1.crq.13[i] = crq.fit[2,1]
    b1.crq.sd.13[i] = crq.fit[2,2]
    cover.crq.13[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.13[i,2]=crq.fit[1,1]
    cover.crq.13[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.13[i,4]=ind(1.410748,cover.crq.13[i,1],cover.crq.13[i,3])
    cover.crq.13[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.13[i,6]=crq.fit[2,1]
    cover.crq.13[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.13[i,8]=ind(0,cover.crq.13[i,5],cover.crq.13[i,7])}
    ,error=function(e){
      b0.crq.13[i] = NA
      b0.crq.sd.13[i] = NA
      b1.crq.13[i] = NA
      b1.crq.sd.13[i] = NA
      cover.crq.13[i,1]=NA
      cover.crq.13[i,2]=NA
      cover.crq.13[i,3]=NA
      cover.crq.13[i,4]=NA
      cover.crq.13[i,5]=NA
      cover.crq.13[i,6]=NA
      cover.crq.13[i,7]=NA
      cover.crq.13[i,8]=NA
    })
  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.rq.13[i] = rq.fit[1,1]
    b0.rq.sd.13[i] = rq.fit[1,2]
    b1.rq.13[i] = rq.fit[2,1]
    b1.rq.sd.13[i] = rq.fit[2,2]
    cover.rq.13[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.13[i,2]=rq.fit[1,1]
    cover.rq.13[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.13[i,4]=ind(1.410748,cover.rq.13[i,1],cover.rq.13[i,3])
    cover.rq.13[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.13[i,6]=rq.fit[2,1]
    cover.rq.13[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.13[i,8]=ind(0,cover.rq.13[i,5],cover.rq.13[i,7])}
    ,error=function(e){
      b0.rq.13[i] = NA
      b0.rq.sd.13[i] = NA
      b1.rq.13[i] = NA
      b1.rq.sd.13[i] = NA
      cover.rq.13[i,1]=NA
      cover.rq.13[i,2]=NA
      cover.rq.13[i,3]=NA
      cover.rq.13[i,4]=NA
      cover.rq.13[i,5]=NA
      cover.rq.13[i,6]=NA
      cover.rq.13[i,7]=NA
      cover.rq.13[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 1, 0.75, 100)
    b0.is.13[i] = ismb.fit[1,1]
    b0.is.sd.13[i] = ismb.fit[1,2]
    b1.is.13[i] = ismb.fit[2,1]
    b1.is.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.is.13[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.13[i,2]=ismb.fit[1,1]
    cover.is.13[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.13[i,4]=ind(1.410748,cover.is.13[i,1],cover.is.13[i,3])
    cover.is.13[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.13[i,6]=ismb.fit[2,1]
    cover.is.13[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.13[i,8]=ind(0,cover.is.13[i,5],cover.is.13[i,7])}
    , error=function(e){
      b0.is.13[i] = NA
      b0.is.sd.13[i] = NA
      b1.is.13[i] = NA
      b1.is.sd.13[i] = NA
      # Coverage
      cover.is.13[i,1]=NA
      cover.is.13[i,2]=NA
      cover.is.13[i,3]=NA
      cover.is.13[i,4]=NA
      cover.is.13[i,5]=NA
      cover.is.13[i,6]=NA
      cover.is.13[i,7]=NA
      cover.is.13[i,8]=NA
    })
}
# Crq beta table
table1.crq[3,1]<-mean(b0.crq.13,na.rm=TRUE)
table1.crq[3,2]<-mean(b0.crq.sd.13,na.rm=TRUE)
table1.crq[3,3]<-sd(b0.crq.13,na.rm=TRUE)
table1.crq[3,4]<-mean(cover.crq.13[,4],na.rm=TRUE)
table1.crq[3,5]<-mean(b1.crq.13,na.rm=TRUE)
table1.crq[3,6]<-mean(b1.crq.sd.13,na.rm=TRUE)
table1.crq[3,7]<-sd(b1.crq.13,na.rm=TRUE)
table1.crq[3,8]<-mean(cover.crq.13[,8],na.rm=TRUE)

# rq beta table
table1.rq[3,1]<-mean(b0.rq.13,na.rm=TRUE)
table1.rq[3,2]<-mean(b0.rq.sd.13,na.rm=TRUE)
table1.rq[3,3]<-sd(b0.rq.13,na.rm=TRUE)
table1.rq[3,4]<-mean(cover.rq.13[,4],na.rm=TRUE)
table1.rq[3,5]<-mean(b1.rq.13,na.rm=TRUE)
table1.rq[3,6]<-mean(b1.rq.sd.13,na.rm=TRUE)
table1.rq[3,7]<-sd(b1.rq.13,na.rm=TRUE)
table1.rq[3,8]<-mean(cover.rq.13[,8],na.rm=TRUE)

# IS beta table
table1.is[3,1]<-mean(b0.is.13,na.rm=TRUE)
table1.is[3,2]<-mean(b0.is.sd.13,na.rm=TRUE)
table1.is[3,3]<-sd(b0.is.13,na.rm=TRUE)
table1.is[3,4]<-mean(cover.is.13[,4],na.rm=TRUE)
table1.is[3,5]<-mean(b1.is.13,na.rm=TRUE)
table1.is[3,6]<-mean(b1.is.sd.13,na.rm=TRUE)
table1.is[3,7]<-sd(b1.is.13,na.rm=TRUE)
table1.is[3,8]<-mean(cover.is.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.crq.15<-c()
b0.crq.sd.15<-c()
b1.crq.15<-c()
b1.crq.sd.15<-c()
b0.rq.15<-c()
b0.rq.sd.15<-c()
b1.rq.15<-c()
b1.rq.sd.15<-c()
b0.is.15<-c()
b0.is.sd.15<-c()
b1.is.15<-c()
b1.is.sd.15<-c()
cover.crq.15=matrix(NA,2000,8)
colnames(cover.crq.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.15=matrix(NA,2000,8)
colnames(cover.rq.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.15=matrix(NA,2000,8)
colnames(cover.is.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  delnumt_1 = c(delnumt_1, sum(a$Z<1))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.crq.15[i] = crq.fit[1,1]
    b0.crq.sd.15[i] = crq.fit[1,2]
    b1.crq.15[i] = crq.fit[2,1]
    b1.crq.sd.15[i] = crq.fit[2,2]
    cover.crq.15[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.15[i,2]=crq.fit[1,1]
    cover.crq.15[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.15[i,4]=ind(1.410748, cover.crq.15[i,1], cover.crq.15[i,3])
    cover.crq.15[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.15[i,6]=crq.fit[2,1]
    cover.crq.15[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.15[i,8]=ind(0, cover.crq.15[i,5], cover.crq.15[i,7])}
    ,error=function(e){
      b0.crq.15[i] = NA
      b0.crq.sd.15[i] = NA
      b1.crq.15[i] = NA
      b1.crq.sd.15[i] = NA
      cover.crq.15[i,1]=NA
      cover.crq.15[i,2]=NA
      cover.crq.15[i,3]=NA
      cover.crq.15[i,4]=NA
      cover.crq.15[i,5]=NA
      cover.crq.15[i,6]=NA
      cover.crq.15[i,7]=NA
      cover.crq.15[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.rq.15[i] = rq.fit[1,1]
    b0.rq.sd.15[i] = rq.fit[1,2]
    b1.rq.15[i] = rq.fit[2,1]
    b1.rq.sd.15[i] = rq.fit[2,2]
    cover.rq.15[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.15[i,2]=rq.fit[1,1]
    cover.rq.15[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.15[i,4]=ind(1.410748, cover.rq.15[i,1], cover.rq.15[i,3])
    cover.rq.15[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.15[i,6]=rq.fit[2,1]
    cover.rq.15[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.15[i,8]=ind(0, cover.rq.15[i,5], cover.rq.15[i,7])}
    ,error=function(e){
      b0.rq.15[i] = NA
      b0.rq.sd.15[i] = NA
      b1.rq.15[i] = NA
      b1.rq.sd.15[i] = NA
      cover.rq.15[i,1]=NA
      cover.rq.15[i,2]=NA
      cover.rq.15[i,3]=NA
      cover.rq.15[i,4]=NA
      cover.rq.15[i,5]=NA
      cover.rq.15[i,6]=NA
      cover.rq.15[i,7]=NA
      cover.rq.15[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 1, 0.75, 100)
    b0.is.15[i] = ismb.fit[1,1]
    b0.is.sd.15[i] = ismb.fit[1,2]
    b1.is.15[i] = ismb.fit[2,1]
    b1.is.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.is.15[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.15[i,2]=ismb.fit[1,1]
    cover.is.15[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.15[i,4]=ind(1.410748,cover.is.15[i,1],cover.is.15[i,3])
    cover.is.15[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.15[i,6]=ismb.fit[2,1]
    cover.is.15[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.15[i,8]=ind(0,cover.is.15[i,5],cover.is.15[i,7])}
    , error=function(e){
      b0.is.15[i] = NA
      b0.is.sd.15[i] = NA
      b1.is.15[i] = NA
      b1.is.sd.15[i] = NA
      # Coverage
      cover.is.15[i,1]=NA
      cover.is.15[i,2]=NA
      cover.is.15[i,3]=NA
      cover.is.15[i,4]=NA
      cover.is.15[i,5]=NA
      cover.is.15[i,6]=NA
      cover.is.15[i,7]=NA
      cover.is.15[i,8]=NA
    })
}

# Crq beta table
table1.crq[4,1]<-mean(b0.crq.15,na.rm=TRUE)
table1.crq[4,2]<-mean(b0.crq.sd.15,na.rm=TRUE)
table1.crq[4,3]<-sd(b0.crq.15,na.rm=TRUE)
table1.crq[4,4]<-mean(cover.crq.15[,4],na.rm=TRUE)
table1.crq[4,5]<-mean(b1.crq.15,na.rm=TRUE)
table1.crq[4,6]<-mean(b1.crq.sd.15,na.rm=TRUE)
table1.crq[4,7]<-sd(b1.crq.15,na.rm=TRUE)
table1.crq[4,8]<-mean(cover.crq.15[,8],na.rm=TRUE)

# rq beta table
table1.rq[4,1]<-mean(b0.rq.15,na.rm=TRUE)
table1.rq[4,2]<-mean(b0.rq.sd.15,na.rm=TRUE)
table1.rq[4,3]<-sd(b0.rq.15,na.rm=TRUE)
table1.rq[4,4]<-mean(cover.rq.15[,4],na.rm=TRUE)
table1.rq[4,5]<-mean(b1.rq.15,na.rm=TRUE)
table1.rq[4,6]<-mean(b1.rq.sd.15,na.rm=TRUE)
table1.rq[4,7]<-sd(b1.rq.15,na.rm=TRUE)
table1.rq[4,8]<-mean(cover.rq.15[,8],na.rm=TRUE)

# IS beta table
table1.is[4,1]<-mean(b0.is.15,na.rm=TRUE)
table1.is[4,2]<-mean(b0.is.sd.15,na.rm=TRUE)
table1.is[4,3]<-sd(b0.is.15,na.rm=TRUE)
table1.is[4,4]<-mean(cover.is.15[,4],na.rm=TRUE)
table1.is[4,5]<-mean(b1.is.15,na.rm=TRUE)
table1.is[4,6]<-mean(b1.is.sd.15,na.rm=TRUE)
table1.is[4,7]<-sd(b1.is.15,na.rm=TRUE)
table1.is[4,8]<-mean(cover.is.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.rq.17<-c()
b0.rq.sd.17<-c()
b1.rq.17<-c()
b1.rq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
cover.crq.17=matrix(NA,2000,8)
colnames(cover.crq.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.17=matrix(NA,2000,8)
colnames(cover.rq.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.17=matrix(NA,2000,8)
colnames(cover.is.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  delnumt_1 = c(delnumt_1, sum(a$Z<1))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.crq.17[i] = crq.fit[1,1]
    b0.crq.sd.17[i] = crq.fit[1,2]
    b1.crq.17[i] = crq.fit[2,1]
    b1.crq.sd.17[i] = crq.fit[2,2]
    cover.crq.17[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.17[i,2]=crq.fit[1,1]
    cover.crq.17[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.17[i,4]=ind(1.410748,cover.crq.17[i,1],cover.crq.17[i,3])
    cover.crq.17[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.17[i,6]=crq.fit[2,1]
    cover.crq.17[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.17[i,8]=ind(0,cover.crq.17[i,5],cover.crq.17[i,7])}
    ,error=function(e){
      b0.crq.17[i] = NA
      b0.crq.sd.17[i] = NA
      b1.crq.17[i] = NA
      b1.crq.sd.17[i] = NA
      cover.crq.17[i,1]=NA
      cover.crq.17[i,2]=NA
      cover.crq.17[i,3]=NA
      cover.crq.17[i,4]=NA
      cover.crq.17[i,5]=NA
      cover.crq.17[i,6]=NA
      cover.crq.17[i,7]=NA
      cover.crq.17[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 1, 0.75)
    b0.rq.17[i] = rq.fit[1,1]
    b0.rq.sd.17[i] = rq.fit[1,2]
    b1.rq.17[i] = rq.fit[2,1]
    b1.rq.sd.17[i] = rq.fit[2,2]
    cover.rq.17[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.17[i,2]=rq.fit[1,1]
    cover.rq.17[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.17[i,4]=ind(1.410748,cover.rq.17[i,1],cover.rq.17[i,3])
    cover.rq.17[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.17[i,6]=rq.fit[2,1]
    cover.rq.17[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.17[i,8]=ind(0,cover.rq.17[i,5],cover.rq.17[i,7])}
    ,error=function(e){
      b0.rq.17[i] = NA
      b0.rq.sd.17[i] = NA
      b1.rq.17[i] = NA
      b1.rq.sd.17[i] = NA
      cover.rq.17[i,1]=NA
      cover.rq.17[i,2]=NA
      cover.rq.17[i,3]=NA
      cover.rq.17[i,4]=NA
      cover.rq.17[i,5]=NA
      cover.rq.17[i,6]=NA
      cover.rq.17[i,7]=NA
      cover.rq.17[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 1, 0.75, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.is.17[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.17[i,2]=ismb.fit[1,1]
    cover.is.17[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.17[i,4]=ind(1.410748,cover.is.17[i,1],cover.is.17[i,3])
    cover.is.17[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.17[i,6]=ismb.fit[2,1]
    cover.is.17[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.17[i,8]=ind(0,cover.is.17[i,5],cover.is.17[i,7])}
    , error=function(e){
      b0.is.17[i] = NA
      b0.is.sd.17[i] = NA
      b1.is.17[i] = NA
      b1.is.sd.17[i] = NA
      # Coverage
      cover.is.17[i,1]=NA
      cover.is.17[i,2]=NA
      cover.is.17[i,3]=NA
      cover.is.17[i,4]=NA
      cover.is.17[i,5]=NA
      cover.is.17[i,6]=NA
      cover.is.17[i,7]=NA
      cover.is.17[i,8]=NA
    })
}
# Crq beta table
table1.crq[5,1]<-mean(b0.crq.17,na.rm=TRUE)
table1.crq[5,2]<-mean(b0.crq.sd.17,na.rm=TRUE)
table1.crq[5,3]<-sd(b0.crq.17,na.rm=TRUE)
table1.crq[5,4]<-mean(cover.crq.17[,4],na.rm=TRUE)
table1.crq[5,5]<-mean(b1.crq.17,na.rm=TRUE)
table1.crq[5,6]<-mean(b1.crq.sd.17,na.rm=TRUE)
table1.crq[5,7]<-sd(b1.crq.17,na.rm=TRUE)
table1.crq[5,8]<-mean(cover.crq.17[,8],na.rm=TRUE)

# rq beta table
table1.rq[5,1]<-mean(b0.rq.17,na.rm=TRUE)
table1.rq[5,2]<-mean(b0.rq.sd.17,na.rm=TRUE)
table1.rq[5,3]<-sd(b0.rq.17,na.rm=TRUE)
table1.rq[5,4]<-mean(cover.rq.17[,4],na.rm=TRUE)
table1.rq[5,5]<-mean(b1.rq.17,na.rm=TRUE)
table1.rq[5,6]<-mean(b1.rq.sd.17,na.rm=TRUE)
table1.rq[5,7]<-sd(b1.rq.17,na.rm=TRUE)
table1.rq[5,8]<-mean(cover.rq.17[,8],na.rm=TRUE)

# IS beta table
table1.is[5,1]<-mean(b0.is.17,na.rm=TRUE)
table1.is[5,2]<-mean(b0.is.sd.17,na.rm=TRUE)
table1.is[5,3]<-sd(b0.is.17,na.rm=TRUE)
table1.is[5,4]<-mean(cover.is.17[,4],na.rm=TRUE)
table1.is[5,5]<-mean(b1.is.17,na.rm=TRUE)
table1.is[5,6]<-mean(b1.is.sd.17,na.rm=TRUE)
table1.is[5,7]<-sd(b1.is.17,na.rm=TRUE)
table1.is[5,8]<-mean(cover.is.17[,8],na.rm=TRUE)

#### censoring point at t_0=2 ####
c.0=5000000
c.1=25.65
c.3=9.89
c.5=6.44
c.7=4.5
#### t_0=2 & c=0% ####
b0.crq.20 = c()
b0.crq.sd.20 = c()
b1.crq.20 = c()
b1.crq.sd.20 = c()
b0.rq.20 = c()
b0.rq.sd.20 = c()
b1.rq.20 = c()
b1.rq.sd.20 = c()
b0.is.20 = c()
b0.is.sd.20 = c()
b1.is.20 = c()
b1.is.sd.20 = c()
delnumt_2 = c()
cover.crq.20=matrix(NA,2000,8)
colnames(cover.crq.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.20=matrix(NA,2000,8)
colnames(cover.rq.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.20=matrix(NA,2000,8)
colnames(cover.is.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  delnumt_2 = c(delnumt_2, sum(a$Z<2))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.crq.20[i] = crq.fit[1,1]
    b0.crq.sd.20[i] = crq.fit[1,2]
    b1.crq.20[i] = crq.fit[2,1]
    b1.crq.sd.20[i] = crq.fit[2,2]
    cover.crq.20[i,1] = crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.20[i,2] = crq.fit[1,1]
    cover.crq.20[i,3] = crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.20[i,4] = ind(1.219403, cover.crq.20[i,1], cover.crq.20[i,3])
    cover.crq.20[i,5] = crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.20[i,6] = crq.fit[2,1]
    cover.crq.20[i,7] = crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.20[i,8] = ind(0, cover.crq.20[i,5], cover.crq.20[i,7])}
    ,error=function(e){
      b0.crq.20[i] = NA
      b0.crq.sd.20[i] = NA
      b1.crq.20[i] = NA
      b1.crq.sd.20[i] = NA
      cover.crq.20[i,1] = NA
      cover.crq.20[i,2] = NA
      cover.crq.20[i,3] = NA
      cover.crq.20[i,4] = NA
      cover.crq.20[i,5] = NA
      cover.crq.20[i,6] = NA
      cover.crq.20[i,7] = NA
      cover.crq.20[i,8] = NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.rq.20[i] = rq.fit[1,1]
    b0.rq.sd.20[i] = rq.fit[1,2]
    b1.rq.20[i] = rq.fit[2,1]
    b1.rq.sd.20[i] = rq.fit[2,2]
    cover.rq.20[i,1] = rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.20[i,2] = rq.fit[1,1]
    cover.rq.20[i,3] = rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.20[i,4] = ind(1.219403, cover.rq.20[i,1], cover.rq.20[i,3])
    cover.rq.20[i,5] = rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.20[i,6] = rq.fit[2,1]
    cover.rq.20[i,7] = rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.20[i,8] = ind(0, cover.rq.20[i,5], cover.rq.20[i,7])}
    ,error=function(e){
      b0.crq.20[i] = NA
      b0.crq.sd.20[i] = NA
      b1.crq.20[i] = NA
      b1.crq.sd.20[i] = NA
      cover.crq.20[i,1] = NA
      cover.crq.20[i,2] = NA
      cover.crq.20[i,3] = NA
      cover.crq.20[i,4] = NA
      cover.crq.20[i,5] = NA
      cover.crq.20[i,6] = NA
      cover.crq.20[i,7] = NA
      cover.crq.20[i,8] = NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 2, 0.75, 100)
    b0.is.20[i] = ismb.fit[1,1]
    b0.is.sd.20[i] = ismb.fit[1,2]
    b1.is.20[i] = ismb.fit[2,1]
    b1.is.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.is.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.20[i,2] = ismb.fit[1,1]
    cover.is.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.20[i,4] = ind(1.219403, cover.is.20[i,1], cover.is.20[i,3])
    cover.is.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.20[i,6] = ismb.fit[2,1]
    cover.is.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.20[i,8] = ind(0, cover.is.20[i,5], cover.is.20[i,7])}
    , error=function(e){
      b0.is.20[i] = NA
      b0.is.sd.20[i] = NA
      b1.is.20[i] = NA
      b1.is.sd.20[i] = NA
      # Coverage
      cover.is.20[i,1] = NA
      cover.is.20[i,2] = NA
      cover.is.20[i,3] = NA
      cover.is.20[i,4] = NA
      cover.is.20[i,5] = NA
      cover.is.20[i,6] = NA
      cover.is.20[i,7] = NA
      cover.is.20[i,8] = NA
    })
}
# Crq beta table
table2.crq[1,1]<-mean(b0.crq.20,na.rm=TRUE)
table2.crq[1,2]<-mean(b0.crq.sd.20,na.rm=TRUE)
table2.crq[1,3]<-sd(b0.crq.20,na.rm=TRUE)
table2.crq[1,4]<-mean(cover.crq.20[,4],na.rm=TRUE)
table2.crq[1,5]<-mean(b1.crq.20,na.rm=TRUE)
table2.crq[1,6]<-mean(b1.crq.sd.20,na.rm=TRUE)
table2.crq[1,7]<-sd(b1.crq.20,na.rm=TRUE)
table2.crq[1,8]<-mean(cover.crq.20[,8],na.rm=TRUE)

# Rq beta table
table2.rq[1,1]<-mean(b0.rq.20,na.rm=TRUE)
table2.rq[1,2]<-mean(b0.rq.sd.20,na.rm=TRUE)
table2.rq[1,3]<-sd(b0.rq.20,na.rm=TRUE)
table2.rq[1,4]<-mean(cover.rq.20[,4],na.rm=TRUE)
table2.rq[1,5]<-mean(b1.rq.20,na.rm=TRUE)
table2.rq[1,6]<-mean(b1.rq.sd.20,na.rm=TRUE)
table2.rq[1,7]<-sd(b1.rq.20,na.rm=TRUE)
table2.rq[1,8]<-mean(cover.rq.20[,8],na.rm=TRUE)

# IS beta table
table2.is[1,1]<-mean(b0.is.20,na.rm=TRUE)
table2.is[1,2]<-mean(b0.is.sd.20,na.rm=TRUE)
table2.is[1,3]<-sd(b0.is.20,na.rm=TRUE)
table2.is[1,4]<-mean(cover.is.20[,4],na.rm=TRUE)
table2.is[1,5]<-mean(b1.is.20,na.rm=TRUE)
table2.is[1,6]<-mean(b1.is.sd.20,na.rm=TRUE)
table2.is[1,7]<-sd(b1.is.20,na.rm=TRUE)
table2.is[1,8]<-mean(cover.is.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.crq.21<-c()
b0.crq.sd.21<-c()
b1.crq.21<-c()
b1.crq.sd.21<-c()
b0.rq.21<-c()
b0.rq.sd.21<-c()
b1.rq.21<-c()
b1.rq.sd.21<-c()
b0.is.21<-c()
b0.is.sd.21<-c()
b1.is.21<-c()
b1.is.sd.21<-c()
cover.crq.21=matrix(NA,2000,8)
colnames(cover.crq.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.21=matrix(NA,2000,8)
colnames(cover.rq.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.21=matrix(NA,2000,8)
colnames(cover.is.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  delnumt_2 = c(delnumt_2, sum(a$Z<2))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.crq.21[i] = crq.fit[1,1]
    b0.crq.sd.21[i] = crq.fit[1,2]
    b1.crq.21[i] = crq.fit[2,1]
    b1.crq.sd.21[i] = crq.fit[2,2]
    cover.crq.21[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.21[i,2]=crq.fit[1,1]
    cover.crq.21[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.21[i,4]=ind(1.219403,cover.crq.21[i,1],cover.crq.21[i,3])
    cover.crq.21[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.21[i,6]=crq.fit[2,1]
    cover.crq.21[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.21[i,8]=ind(0,cover.crq.21[i,5],cover.crq.21[i,7])}
    ,error=function(e){
      b0.crq.21[i] = NA
      b0.crq.sd.21[i] = NA
      b1.crq.21[i] = NA
      b1.crq.sd.21[i] = NA
      cover.crq.21[i,1]=NA
      cover.crq.21[i,2]=NA
      cover.crq.21[i,3]=NA
      cover.crq.21[i,4]=NA
      cover.crq.21[i,5]=NA
      cover.crq.21[i,6]=NA
      cover.crq.21[i,7]=NA
      cover.crq.21[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.rq.21[i] = rq.fit[1,1]
    b0.rq.sd.21[i] = rq.fit[1,2]
    b1.rq.21[i] = rq.fit[2,1]
    b1.rq.sd.21[i] = rq.fit[2,2]
    cover.rq.21[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.21[i,2]=rq.fit[1,1]
    cover.rq.21[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.21[i,4]=ind(1.219403,cover.rq.21[i,1],cover.rq.21[i,3])
    cover.rq.21[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.21[i,6]=rq.fit[2,1]
    cover.rq.21[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.21[i,8]=ind(0,cover.rq.21[i,5],cover.rq.21[i,7])}
    ,error=function(e){
      b0.rq.21[i] = NA
      b0.rq.sd.21[i] = NA
      b1.rq.21[i] = NA
      b1.rq.sd.21[i] = NA
      cover.rq.21[i,1]=NA
      cover.rq.21[i,2]=NA
      cover.rq.21[i,3]=NA
      cover.rq.21[i,4]=NA
      cover.rq.21[i,5]=NA
      cover.rq.21[i,6]=NA
      cover.rq.21[i,7]=NA
      cover.rq.21[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 2, 0.75, 100)
    b0.is.21[i] = ismb.fit[1,1]
    b0.is.sd.21[i] = ismb.fit[1,2]
    b1.is.21[i] = ismb.fit[2,1]
    b1.is.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.is.21[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.21[i,2]=ismb.fit[1,1]
    cover.is.21[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.21[i,4]=ind(1.219403,cover.is.21[i,1],cover.is.21[i,3])
    cover.is.21[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.21[i,6]=ismb.fit[2,1]
    cover.is.21[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.21[i,8]=ind(0,cover.is.21[i,5],cover.is.21[i,7])}
    , error=function(e){
      b0.is.21[i] = NA
      b0.is.sd.21[i] = NA
      b1.is.21[i] = NA
      b1.is.sd.21[i] = NA
      # Coverage
      cover.is.21[i,1]=NA
      cover.is.21[i,2]=NA
      cover.is.21[i,3]=NA
      cover.is.21[i,4]=NA
      cover.is.21[i,5]=NA
      cover.is.21[i,6]=NA
      cover.is.21[i,7]=NA
      cover.is.21[i,8]=NA
    })
}

# Crq beta table
table2.crq[2,1]<-mean(b0.crq.21,na.rm=TRUE)
table2.crq[2,2]<-mean(b0.crq.sd.21,na.rm=TRUE)
table2.crq[2,3]<-sd(b0.crq.21,na.rm=TRUE)
table2.crq[2,4]<-mean(cover.crq.21[,4],na.rm=TRUE)
table2.crq[2,5]<-mean(b1.crq.21,na.rm=TRUE)
table2.crq[2,6]<-mean(b1.crq.sd.21,na.rm=TRUE)
table2.crq[2,7]<-sd(b1.crq.21,na.rm=TRUE)
table2.crq[2,8]<-mean(cover.crq.21[,8],na.rm=TRUE)

# rq beta table
table2.rq[2,1]<-mean(b0.rq.21,na.rm=TRUE)
table2.rq[2,2]<-mean(b0.rq.sd.21,na.rm=TRUE)
table2.rq[2,3]<-sd(b0.rq.21,na.rm=TRUE)
table2.rq[2,4]<-mean(cover.rq.21[,4],na.rm=TRUE)
table2.rq[2,5]<-mean(b1.rq.21,na.rm=TRUE)
table2.rq[2,6]<-mean(b1.rq.sd.21,na.rm=TRUE)
table2.rq[2,7]<-sd(b1.rq.21,na.rm=TRUE)
table2.rq[2,8]<-mean(cover.rq.21[,8],na.rm=TRUE)

# IS beta table
table2.is[2,1]<-mean(b0.is.21,na.rm=TRUE)
table2.is[2,2]<-mean(b0.is.sd.21,na.rm=TRUE)
table2.is[2,3]<-sd(b0.is.21,na.rm=TRUE)
table2.is[2,4]<-mean(cover.is.21[,4],na.rm=TRUE)
table2.is[2,5]<-mean(b1.is.21,na.rm=TRUE)
table2.is[2,6]<-mean(b1.is.sd.21,na.rm=TRUE)
table2.is[2,7]<-sd(b1.is.21,na.rm=TRUE)
table2.is[2,8]<-mean(cover.is.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.crq.23<-c()
b0.crq.sd.23<-c()
b1.crq.23<-c()
b1.crq.sd.23<-c()
b0.rq.23<-c()
b0.rq.sd.23<-c()
b1.rq.23<-c()
b1.rq.sd.23<-c()
b0.is.23<-c()
b0.is.sd.23<-c()
b1.is.23<-c()
b1.is.sd.23<-c()
cover.crq.23=matrix(NA,2000,8)
colnames(cover.crq.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.23=matrix(NA,2000,8)
colnames(cover.rq.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.23=matrix(NA,2000,8)
colnames(cover.is.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  delnumt_2 = c(delnumt_2, sum(a$Z<2))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.crq.23[i] = crq.fit[1,1]
    b0.crq.sd.23[i] = crq.fit[1,2]
    b1.crq.23[i] = crq.fit[2,1]
    b1.crq.sd.23[i] = crq.fit[2,2]
    cover.crq.23[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.23[i,2]=crq.fit[1,1]
    cover.crq.23[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.23[i,4]=ind(1.219403,cover.crq.23[i,1],cover.crq.23[i,3])
    cover.crq.23[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.23[i,6]=crq.fit[2,1]
    cover.crq.23[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.23[i,8]=ind(0,cover.crq.23[i,5],cover.crq.23[i,7])}
    ,error=function(e){
      b0.crq.23[i] = NA
      b0.crq.sd.23[i] = NA
      b1.crq.23[i] = NA
      b1.crq.sd.23[i] = NA
      cover.crq.23[i,1]=NA
      cover.crq.23[i,2]=NA
      cover.crq.23[i,3]=NA
      cover.crq.23[i,4]=NA
      cover.crq.23[i,5]=NA
      cover.crq.23[i,6]=NA
      cover.crq.23[i,7]=NA
      cover.crq.23[i,8]=NA
    })
  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.rq.23[i] = rq.fit[1,1]
    b0.rq.sd.23[i] = rq.fit[1,2]
    b1.rq.23[i] = rq.fit[2,1]
    b1.rq.sd.23[i] = rq.fit[2,2]
    cover.rq.23[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.23[i,2]=rq.fit[1,1]
    cover.rq.23[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.23[i,4]=ind(1.219403,cover.rq.23[i,1],cover.rq.23[i,3])
    cover.rq.23[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.23[i,6]=rq.fit[2,1]
    cover.rq.23[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.23[i,8]=ind(0,cover.rq.23[i,5],cover.rq.23[i,7])}
    ,error=function(e){
      b0.rq.23[i] = NA
      b0.rq.sd.23[i] = NA
      b1.rq.23[i] = NA
      b1.rq.sd.23[i] = NA
      cover.rq.23[i,1]=NA
      cover.rq.23[i,2]=NA
      cover.rq.23[i,3]=NA
      cover.rq.23[i,4]=NA
      cover.rq.23[i,5]=NA
      cover.rq.23[i,6]=NA
      cover.rq.23[i,7]=NA
      cover.rq.23[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 2, 0.75, 100)
    b0.is.23[i] = ismb.fit[1,1]
    b0.is.sd.23[i] = ismb.fit[1,2]
    b1.is.23[i] = ismb.fit[2,1]
    b1.is.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.is.23[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.23[i,2]=ismb.fit[1,1]
    cover.is.23[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.23[i,4]=ind(1.219403,cover.is.23[i,1],cover.is.23[i,3])
    cover.is.23[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.23[i,6]=ismb.fit[2,1]
    cover.is.23[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.23[i,8]=ind(0,cover.is.23[i,5],cover.is.23[i,7])}
    , error=function(e){
      b0.is.23[i] = NA
      b0.is.sd.23[i] = NA
      b1.is.23[i] = NA
      b1.is.sd.23[i] = NA
      # Coverage
      cover.is.23[i,1]=NA
      cover.is.23[i,2]=NA
      cover.is.23[i,3]=NA
      cover.is.23[i,4]=NA
      cover.is.23[i,5]=NA
      cover.is.23[i,6]=NA
      cover.is.23[i,7]=NA
      cover.is.23[i,8]=NA
    })
}
# Crq beta table
table2.crq[3,1]<-mean(b0.crq.23,na.rm=TRUE)
table2.crq[3,2]<-mean(b0.crq.sd.23,na.rm=TRUE)
table2.crq[3,3]<-sd(b0.crq.23,na.rm=TRUE)
table2.crq[3,4]<-mean(cover.crq.23[,4],na.rm=TRUE)
table2.crq[3,5]<-mean(b1.crq.23,na.rm=TRUE)
table2.crq[3,6]<-mean(b1.crq.sd.23,na.rm=TRUE)
table2.crq[3,7]<-sd(b1.crq.23,na.rm=TRUE)
table2.crq[3,8]<-mean(cover.crq.23[,8],na.rm=TRUE)

# rq beta table
table2.rq[3,1]<-mean(b0.rq.23,na.rm=TRUE)
table2.rq[3,2]<-mean(b0.rq.sd.23,na.rm=TRUE)
table2.rq[3,3]<-sd(b0.rq.23,na.rm=TRUE)
table2.rq[3,4]<-mean(cover.rq.23[,4],na.rm=TRUE)
table2.rq[3,5]<-mean(b1.rq.23,na.rm=TRUE)
table2.rq[3,6]<-mean(b1.rq.sd.23,na.rm=TRUE)
table2.rq[3,7]<-sd(b1.rq.23,na.rm=TRUE)
table2.rq[3,8]<-mean(cover.rq.23[,8],na.rm=TRUE)

# IS beta table
table2.is[3,1]<-mean(b0.is.23,na.rm=TRUE)
table2.is[3,2]<-mean(b0.is.sd.23,na.rm=TRUE)
table2.is[3,3]<-sd(b0.is.23,na.rm=TRUE)
table2.is[3,4]<-mean(cover.is.23[,4],na.rm=TRUE)
table2.is[3,5]<-mean(b1.is.23,na.rm=TRUE)
table2.is[3,6]<-mean(b1.is.sd.23,na.rm=TRUE)
table2.is[3,7]<-sd(b1.is.23,na.rm=TRUE)
table2.is[3,8]<-mean(cover.is.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.crq.25<-c()
b0.crq.sd.25<-c()
b1.crq.25<-c()
b1.crq.sd.25<-c()
b0.rq.25<-c()
b0.rq.sd.25<-c()
b1.rq.25<-c()
b1.rq.sd.25<-c()
b0.is.25<-c()
b0.is.sd.25<-c()
b1.is.25<-c()
b1.is.sd.25<-c()
cover.crq.25=matrix(NA,2000,8)
colnames(cover.crq.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.25=matrix(NA,2000,8)
colnames(cover.rq.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.25=matrix(NA,2000,8)
colnames(cover.is.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  delnumt_2 = c(delnumt_2, sum(a$Z<2))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.crq.25[i] = crq.fit[1,1]
    b0.crq.sd.25[i] = crq.fit[1,2]
    b1.crq.25[i] = crq.fit[2,1]
    b1.crq.sd.25[i] = crq.fit[2,2]
    cover.crq.25[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.25[i,2]=crq.fit[1,1]
    cover.crq.25[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.25[i,4]=ind(1.219403, cover.crq.25[i,1], cover.crq.25[i,3])
    cover.crq.25[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.25[i,6]=crq.fit[2,1]
    cover.crq.25[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.25[i,8]=ind(0, cover.crq.25[i,5], cover.crq.25[i,7])}
    ,error=function(e){
      b0.crq.25[i] = NA
      b0.crq.sd.25[i] = NA
      b1.crq.25[i] = NA
      b1.crq.sd.25[i] = NA
      cover.crq.25[i,1]=NA
      cover.crq.25[i,2]=NA
      cover.crq.25[i,3]=NA
      cover.crq.25[i,4]=NA
      cover.crq.25[i,5]=NA
      cover.crq.25[i,6]=NA
      cover.crq.25[i,7]=NA
      cover.crq.25[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.rq.25[i] = rq.fit[1,1]
    b0.rq.sd.25[i] = rq.fit[1,2]
    b1.rq.25[i] = rq.fit[2,1]
    b1.rq.sd.25[i] = rq.fit[2,2]
    cover.rq.25[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.25[i,2]=rq.fit[1,1]
    cover.rq.25[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.25[i,4]=ind(1.219403, cover.rq.25[i,1], cover.rq.25[i,3])
    cover.rq.25[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.25[i,6]=rq.fit[2,1]
    cover.rq.25[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.25[i,8]=ind(0, cover.rq.25[i,5], cover.rq.25[i,7])}
    ,error=function(e){
      b0.rq.25[i] = NA
      b0.rq.sd.25[i] = NA
      b1.rq.25[i] = NA
      b1.rq.sd.25[i] = NA
      cover.rq.25[i,1]=NA
      cover.rq.25[i,2]=NA
      cover.rq.25[i,3]=NA
      cover.rq.25[i,4]=NA
      cover.rq.25[i,5]=NA
      cover.rq.25[i,6]=NA
      cover.rq.25[i,7]=NA
      cover.rq.25[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 2, 0.75, 100)
    b0.is.25[i] = ismb.fit[1,1]
    b0.is.sd.25[i] = ismb.fit[1,2]
    b1.is.25[i] = ismb.fit[2,1]
    b1.is.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.is.25[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.25[i,2]=ismb.fit[1,1]
    cover.is.25[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.25[i,4]=ind(1.219403,cover.is.25[i,1],cover.is.25[i,3])
    cover.is.25[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.25[i,6]=ismb.fit[2,1]
    cover.is.25[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.25[i,8]=ind(0,cover.is.25[i,5],cover.is.25[i,7])}
    , error=function(e){
      b0.is.25[i] = NA
      b0.is.sd.25[i] = NA
      b1.is.25[i] = NA
      b1.is.sd.25[i] = NA
      # Coverage
      cover.is.25[i,1]=NA
      cover.is.25[i,2]=NA
      cover.is.25[i,3]=NA
      cover.is.25[i,4]=NA
      cover.is.25[i,5]=NA
      cover.is.25[i,6]=NA
      cover.is.25[i,7]=NA
      cover.is.25[i,8]=NA
    })
}

# Crq beta table
table2.crq[4,1]<-mean(b0.crq.25,na.rm=TRUE)
table2.crq[4,2]<-mean(b0.crq.sd.25,na.rm=TRUE)
table2.crq[4,3]<-sd(b0.crq.25,na.rm=TRUE)
table2.crq[4,4]<-mean(cover.crq.25[,4],na.rm=TRUE)
table2.crq[4,5]<-mean(b1.crq.25,na.rm=TRUE)
table2.crq[4,6]<-mean(b1.crq.sd.25,na.rm=TRUE)
table2.crq[4,7]<-sd(b1.crq.25,na.rm=TRUE)
table2.crq[4,8]<-mean(cover.crq.25[,8],na.rm=TRUE)

# rq beta table
table2.rq[4,1]<-mean(b0.rq.25,na.rm=TRUE)
table2.rq[4,2]<-mean(b0.rq.sd.25,na.rm=TRUE)
table2.rq[4,3]<-sd(b0.rq.25,na.rm=TRUE)
table2.rq[4,4]<-mean(cover.rq.25[,4],na.rm=TRUE)
table2.rq[4,5]<-mean(b1.rq.25,na.rm=TRUE)
table2.rq[4,6]<-mean(b1.rq.sd.25,na.rm=TRUE)
table2.rq[4,7]<-sd(b1.rq.25,na.rm=TRUE)
table2.rq[4,8]<-mean(cover.rq.25[,8],na.rm=TRUE)

# IS beta table
table2.is[4,1]<-mean(b0.is.25,na.rm=TRUE)
table2.is[4,2]<-mean(b0.is.sd.25,na.rm=TRUE)
table2.is[4,3]<-sd(b0.is.25,na.rm=TRUE)
table2.is[4,4]<-mean(cover.is.25[,4],na.rm=TRUE)
table2.is[4,5]<-mean(b1.is.25,na.rm=TRUE)
table2.is[4,6]<-mean(b1.is.sd.25,na.rm=TRUE)
table2.is[4,7]<-sd(b1.is.25,na.rm=TRUE)
table2.is[4,8]<-mean(cover.is.25[,8],na.rm=TRUE)

#### t_0=2 & c=70% ####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.rq.27<-c()
b0.rq.sd.27<-c()
b1.rq.27<-c()
b1.rq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
cover.crq.27=matrix(NA,2000,8)
colnames(cover.crq.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.27=matrix(NA,2000,8)
colnames(cover.rq.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.27=matrix(NA,2000,8)
colnames(cover.is.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  delnumt_2 = c(delnumt_2, sum(a$Z<2))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.crq.27[i] = crq.fit[1,1]
    b0.crq.sd.27[i] = crq.fit[1,2]
    b1.crq.27[i] = crq.fit[2,1]
    b1.crq.sd.27[i] = crq.fit[2,2]
    cover.crq.27[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.27[i,2]=crq.fit[1,1]
    cover.crq.27[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.27[i,4]=ind(1.219403,cover.crq.27[i,1],cover.crq.27[i,3])
    cover.crq.27[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.27[i,6]=crq.fit[2,1]
    cover.crq.27[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.27[i,8]=ind(0,cover.crq.27[i,5],cover.crq.27[i,7])}
    ,error=function(e){
      b0.crq.27[i] = NA
      b0.crq.sd.27[i] = NA
      b1.crq.27[i] = NA
      b1.crq.sd.27[i] = NA
      cover.crq.27[i,1]=NA
      cover.crq.27[i,2]=NA
      cover.crq.27[i,3]=NA
      cover.crq.27[i,4]=NA
      cover.crq.27[i,5]=NA
      cover.crq.27[i,6]=NA
      cover.crq.27[i,7]=NA
      cover.crq.27[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 2, 0.75)
    b0.rq.27[i] = rq.fit[1,1]
    b0.rq.sd.27[i] = rq.fit[1,2]
    b1.rq.27[i] = rq.fit[2,1]
    b1.rq.sd.27[i] = rq.fit[2,2]
    cover.rq.27[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.27[i,2]=rq.fit[1,1]
    cover.rq.27[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.27[i,4]=ind(1.219403,cover.rq.27[i,1],cover.rq.27[i,3])
    cover.rq.27[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.27[i,6]=rq.fit[2,1]
    cover.rq.27[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.27[i,8]=ind(0,cover.rq.27[i,5],cover.rq.27[i,7])}
    ,error=function(e){
      b0.rq.27[i] = NA
      b0.rq.sd.27[i] = NA
      b1.rq.27[i] = NA
      b1.rq.sd.27[i] = NA
      cover.rq.27[i,1]=NA
      cover.rq.27[i,2]=NA
      cover.rq.27[i,3]=NA
      cover.rq.27[i,4]=NA
      cover.rq.27[i,5]=NA
      cover.rq.27[i,6]=NA
      cover.rq.27[i,7]=NA
      cover.rq.27[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 2, 0.75, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.is.27[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.27[i,2]=ismb.fit[1,1]
    cover.is.27[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.27[i,4]=ind(1.219403,cover.is.27[i,1],cover.is.27[i,3])
    cover.is.27[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.27[i,6]=ismb.fit[2,1]
    cover.is.27[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.27[i,8]=ind(0,cover.is.27[i,5],cover.is.27[i,7])}
    , error=function(e){
      b0.is.27[i] = NA
      b0.is.sd.27[i] = NA
      b1.is.27[i] = NA
      b1.is.sd.27[i] = NA
      # Coverage
      cover.is.27[i,1]=NA
      cover.is.27[i,2]=NA
      cover.is.27[i,3]=NA
      cover.is.27[i,4]=NA
      cover.is.27[i,5]=NA
      cover.is.27[i,6]=NA
      cover.is.27[i,7]=NA
      cover.is.27[i,8]=NA
    })
}
# Crq beta table
table2.crq[5,1]<-mean(b0.crq.27,na.rm=TRUE)
table2.crq[5,2]<-mean(b0.crq.sd.27,na.rm=TRUE)
table2.crq[5,3]<-sd(b0.crq.27,na.rm=TRUE)
table2.crq[5,4]<-mean(cover.crq.27[,4],na.rm=TRUE)
table2.crq[5,5]<-mean(b1.crq.27,na.rm=TRUE)
table2.crq[5,6]<-mean(b1.crq.sd.27,na.rm=TRUE)
table2.crq[5,7]<-sd(b1.crq.27,na.rm=TRUE)
table2.crq[5,8]<-mean(cover.crq.27[,8],na.rm=TRUE)

# rq beta table
table2.rq[5,1]<-mean(b0.rq.27,na.rm=TRUE)
table2.rq[5,2]<-mean(b0.rq.sd.27,na.rm=TRUE)
table2.rq[5,3]<-sd(b0.rq.27,na.rm=TRUE)
table2.rq[5,4]<-mean(cover.rq.27[,4],na.rm=TRUE)
table2.rq[5,5]<-mean(b1.rq.27,na.rm=TRUE)
table2.rq[5,6]<-mean(b1.rq.sd.27,na.rm=TRUE)
table2.rq[5,7]<-sd(b1.rq.27,na.rm=TRUE)
table2.rq[5,8]<-mean(cover.rq.27[,8],na.rm=TRUE)

# IS beta table
table2.is[5,1]<-mean(b0.is.27,na.rm=TRUE)
table2.is[5,2]<-mean(b0.is.sd.27,na.rm=TRUE)
table2.is[5,3]<-sd(b0.is.27,na.rm=TRUE)
table2.is[5,4]<-mean(cover.is.27[,4],na.rm=TRUE)
table2.is[5,5]<-mean(b1.is.27,na.rm=TRUE)
table2.is[5,6]<-mean(b1.is.sd.27,na.rm=TRUE)
table2.is[5,7]<-sd(b1.is.27,na.rm=TRUE)
table2.is[5,8]<-mean(cover.is.27[,8],na.rm=TRUE)

#### censoring point at t_0=3 ####
c.0=5000000
c.1=21.92
c.3=9.46
c.5=6.59
c.7=4.91
#### t_0=3 & c=0% ####
b0.crq.30 = c()
b0.crq.sd.30 = c()
b1.crq.30 = c()
b1.crq.sd.30 = c()
b0.rq.30 = c()
b0.rq.sd.30 = c()
b1.rq.30 = c()
b1.rq.sd.30 = c()
b0.is.30 = c()
b0.is.sd.30 = c()
b1.is.30 = c()
b1.is.sd.30 = c()
delnumt_3 = c()
cover.crq.30=matrix(NA,2000,8)
colnames(cover.crq.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.30=matrix(NA,2000,8)
colnames(cover.rq.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.30=matrix(NA,2000,8)
colnames(cover.is.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  delnumt_3 = c(delnumt_3, sum(a$Z<3))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.crq.30[i] = crq.fit[1,1]
    b0.crq.sd.30[i] = crq.fit[1,2]
    b1.crq.30[i] = crq.fit[2,1]
    b1.crq.sd.30[i] = crq.fit[2,2]
    cover.crq.30[i,1] = crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.30[i,2] = crq.fit[1,1]
    cover.crq.30[i,3] = crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.30[i,4] = ind(1.040613, cover.crq.30[i,1], cover.crq.30[i,3])
    cover.crq.30[i,5] = crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.30[i,6] = crq.fit[2,1]
    cover.crq.30[i,7] = crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.30[i,8] = ind(0, cover.crq.30[i,5], cover.crq.30[i,7])}
    ,error=function(e){
      b0.crq.30[i] = NA
      b0.crq.sd.30[i] = NA
      b1.crq.30[i] = NA
      b1.crq.sd.30[i] = NA
      cover.crq.30[i,1] = NA
      cover.crq.30[i,2] = NA
      cover.crq.30[i,3] = NA
      cover.crq.30[i,4] = NA
      cover.crq.30[i,5] = NA
      cover.crq.30[i,6] = NA
      cover.crq.30[i,7] = NA
      cover.crq.30[i,8] = NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.rq.30[i] = rq.fit[1,1]
    b0.rq.sd.30[i] = rq.fit[1,2]
    b1.rq.30[i] = rq.fit[2,1]
    b1.rq.sd.30[i] = rq.fit[2,2]
    cover.rq.30[i,1] = rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.30[i,2] = rq.fit[1,1]
    cover.rq.30[i,3] = rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.30[i,4] = ind(1.040613, cover.rq.30[i,1], cover.rq.30[i,3])
    cover.rq.30[i,5] = rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.30[i,6] = rq.fit[2,1]
    cover.rq.30[i,7] = rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.30[i,8] = ind(0, cover.rq.30[i,5], cover.rq.30[i,7])}
    ,error=function(e){
      b0.crq.30[i] = NA
      b0.crq.sd.30[i] = NA
      b1.crq.30[i] = NA
      b1.crq.sd.30[i] = NA
      cover.crq.30[i,1] = NA
      cover.crq.30[i,2] = NA
      cover.crq.30[i,3] = NA
      cover.crq.30[i,4] = NA
      cover.crq.30[i,5] = NA
      cover.crq.30[i,6] = NA
      cover.crq.30[i,7] = NA
      cover.crq.30[i,8] = NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 3, 0.75, 100)
    b0.is.30[i] = ismb.fit[1,1]
    b0.is.sd.30[i] = ismb.fit[1,2]
    b1.is.30[i] = ismb.fit[2,1]
    b1.is.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.is.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.30[i,2] = ismb.fit[1,1]
    cover.is.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.30[i,4] = ind(1.040613, cover.is.30[i,1], cover.is.30[i,3])
    cover.is.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.30[i,6] = ismb.fit[2,1]
    cover.is.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.30[i,8] = ind(0, cover.is.30[i,5], cover.is.30[i,7])}
    , error=function(e){
      b0.is.30[i] = NA
      b0.is.sd.30[i] = NA
      b1.is.30[i] = NA
      b1.is.sd.30[i] = NA
      # Coverage
      cover.is.30[i,1] = NA
      cover.is.30[i,2] = NA
      cover.is.30[i,3] = NA
      cover.is.30[i,4] = NA
      cover.is.30[i,5] = NA
      cover.is.30[i,6] = NA
      cover.is.30[i,7] = NA
      cover.is.30[i,8] = NA
    })
}
# Crq beta table
table3.crq[1,1]<-mean(b0.crq.30,na.rm=TRUE)
table3.crq[1,2]<-mean(b0.crq.sd.30,na.rm=TRUE)
table3.crq[1,3]<-sd(b0.crq.30,na.rm=TRUE)
table3.crq[1,4]<-mean(cover.crq.30[,4],na.rm=TRUE)
table3.crq[1,5]<-mean(b1.crq.30,na.rm=TRUE)
table3.crq[1,6]<-mean(b1.crq.sd.30,na.rm=TRUE)
table3.crq[1,7]<-sd(b1.crq.30,na.rm=TRUE)
table3.crq[1,8]<-mean(cover.crq.30[,8],na.rm=TRUE)

# Rq beta table
table3.rq[1,1]<-mean(b0.rq.30,na.rm=TRUE)
table3.rq[1,2]<-mean(b0.rq.sd.30,na.rm=TRUE)
table3.rq[1,3]<-sd(b0.rq.30,na.rm=TRUE)
table3.rq[1,4]<-mean(cover.rq.30[,4],na.rm=TRUE)
table3.rq[1,5]<-mean(b1.rq.30,na.rm=TRUE)
table3.rq[1,6]<-mean(b1.rq.sd.30,na.rm=TRUE)
table3.rq[1,7]<-sd(b1.rq.30,na.rm=TRUE)
table3.rq[1,8]<-mean(cover.rq.30[,8],na.rm=TRUE)

# IS beta table
table3.is[1,1]<-mean(b0.is.30,na.rm=TRUE)
table3.is[1,2]<-mean(b0.is.sd.30,na.rm=TRUE)
table3.is[1,3]<-sd(b0.is.30,na.rm=TRUE)
table3.is[1,4]<-mean(cover.is.30[,4],na.rm=TRUE)
table3.is[1,5]<-mean(b1.is.30,na.rm=TRUE)
table3.is[1,6]<-mean(b1.is.sd.30,na.rm=TRUE)
table3.is[1,7]<-sd(b1.is.30,na.rm=TRUE)
table3.is[1,8]<-mean(cover.is.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.crq.31<-c()
b0.crq.sd.31<-c()
b1.crq.31<-c()
b1.crq.sd.31<-c()
b0.rq.31<-c()
b0.rq.sd.31<-c()
b1.rq.31<-c()
b1.rq.sd.31<-c()
b0.is.31<-c()
b0.is.sd.31<-c()
b1.is.31<-c()
b1.is.sd.31<-c()
cover.crq.31=matrix(NA,2000,8)
colnames(cover.crq.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.31=matrix(NA,2000,8)
colnames(cover.rq.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.31=matrix(NA,2000,8)
colnames(cover.is.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  delnumt_3 = c(delnumt_3, sum(a$Z<3))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.crq.31[i] = crq.fit[1,1]
    b0.crq.sd.31[i] = crq.fit[1,2]
    b1.crq.31[i] = crq.fit[2,1]
    b1.crq.sd.31[i] = crq.fit[2,2]
    cover.crq.31[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.31[i,2]=crq.fit[1,1]
    cover.crq.31[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.31[i,4]=ind(1.040613,cover.crq.31[i,1],cover.crq.31[i,3])
    cover.crq.31[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.31[i,6]=crq.fit[2,1]
    cover.crq.31[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.31[i,8]=ind(0,cover.crq.31[i,5],cover.crq.31[i,7])}
    ,error=function(e){
      b0.crq.31[i] = NA
      b0.crq.sd.31[i] = NA
      b1.crq.31[i] = NA
      b1.crq.sd.31[i] = NA
      cover.crq.31[i,1]=NA
      cover.crq.31[i,2]=NA
      cover.crq.31[i,3]=NA
      cover.crq.31[i,4]=NA
      cover.crq.31[i,5]=NA
      cover.crq.31[i,6]=NA
      cover.crq.31[i,7]=NA
      cover.crq.31[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.rq.31[i] = rq.fit[1,1]
    b0.rq.sd.31[i] = rq.fit[1,2]
    b1.rq.31[i] = rq.fit[2,1]
    b1.rq.sd.31[i] = rq.fit[2,2]
    cover.rq.31[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.31[i,2]=rq.fit[1,1]
    cover.rq.31[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.31[i,4]=ind(1.040613,cover.rq.31[i,1],cover.rq.31[i,3])
    cover.rq.31[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.31[i,6]=rq.fit[2,1]
    cover.rq.31[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.31[i,8]=ind(0,cover.rq.31[i,5],cover.rq.31[i,7])}
    ,error=function(e){
      b0.rq.31[i] = NA
      b0.rq.sd.31[i] = NA
      b1.rq.31[i] = NA
      b1.rq.sd.31[i] = NA
      cover.rq.31[i,1]=NA
      cover.rq.31[i,2]=NA
      cover.rq.31[i,3]=NA
      cover.rq.31[i,4]=NA
      cover.rq.31[i,5]=NA
      cover.rq.31[i,6]=NA
      cover.rq.31[i,7]=NA
      cover.rq.31[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 3, 0.75, 100)
    b0.is.31[i] = ismb.fit[1,1]
    b0.is.sd.31[i] = ismb.fit[1,2]
    b1.is.31[i] = ismb.fit[2,1]
    b1.is.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.is.31[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.31[i,2]=ismb.fit[1,1]
    cover.is.31[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.31[i,4]=ind(1.040613,cover.is.31[i,1],cover.is.31[i,3])
    cover.is.31[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.31[i,6]=ismb.fit[2,1]
    cover.is.31[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.31[i,8]=ind(0,cover.is.31[i,5],cover.is.31[i,7])}
    , error=function(e){
      b0.is.31[i] = NA
      b0.is.sd.31[i] = NA
      b1.is.31[i] = NA
      b1.is.sd.31[i] = NA
      # Coverage
      cover.is.31[i,1]=NA
      cover.is.31[i,2]=NA
      cover.is.31[i,3]=NA
      cover.is.31[i,4]=NA
      cover.is.31[i,5]=NA
      cover.is.31[i,6]=NA
      cover.is.31[i,7]=NA
      cover.is.31[i,8]=NA
    })
}

# Crq beta table
table3.crq[2,1]<-mean(b0.crq.31,na.rm=TRUE)
table3.crq[2,2]<-mean(b0.crq.sd.31,na.rm=TRUE)
table3.crq[2,3]<-sd(b0.crq.31,na.rm=TRUE)
table3.crq[2,4]<-mean(cover.crq.31[,4],na.rm=TRUE)
table3.crq[2,5]<-mean(b1.crq.31,na.rm=TRUE)
table3.crq[2,6]<-mean(b1.crq.sd.31,na.rm=TRUE)
table3.crq[2,7]<-sd(b1.crq.31,na.rm=TRUE)
table3.crq[2,8]<-mean(cover.crq.31[,8],na.rm=TRUE)

# rq beta table
table3.rq[2,1]<-mean(b0.rq.31,na.rm=TRUE)
table3.rq[2,2]<-mean(b0.rq.sd.31,na.rm=TRUE)
table3.rq[2,3]<-sd(b0.rq.31,na.rm=TRUE)
table3.rq[2,4]<-mean(cover.rq.31[,4],na.rm=TRUE)
table3.rq[2,5]<-mean(b1.rq.31,na.rm=TRUE)
table3.rq[2,6]<-mean(b1.rq.sd.31,na.rm=TRUE)
table3.rq[2,7]<-sd(b1.rq.31,na.rm=TRUE)
table3.rq[2,8]<-mean(cover.rq.31[,8],na.rm=TRUE)

# IS beta table
table3.is[2,1]<-mean(b0.is.31,na.rm=TRUE)
table3.is[2,2]<-mean(b0.is.sd.31,na.rm=TRUE)
table3.is[2,3]<-sd(b0.is.31,na.rm=TRUE)
table3.is[2,4]<-mean(cover.is.31[,4],na.rm=TRUE)
table3.is[2,5]<-mean(b1.is.31,na.rm=TRUE)
table3.is[2,6]<-mean(b1.is.sd.31,na.rm=TRUE)
table3.is[2,7]<-sd(b1.is.31,na.rm=TRUE)
table3.is[2,8]<-mean(cover.is.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
b0.crq.33<-c()
b0.crq.sd.33<-c()
b1.crq.33<-c()
b1.crq.sd.33<-c()
b0.rq.33<-c()
b0.rq.sd.33<-c()
b1.rq.33<-c()
b1.rq.sd.33<-c()
b0.is.33<-c()
b0.is.sd.33<-c()
b1.is.33<-c()
b1.is.sd.33<-c()
cover.crq.33=matrix(NA,2000,8)
colnames(cover.crq.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.33=matrix(NA,2000,8)
colnames(cover.rq.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.33=matrix(NA,2000,8)
colnames(cover.is.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  delnumt_3 = c(delnumt_3, sum(a$Z<3))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.crq.33[i] = crq.fit[1,1]
    b0.crq.sd.33[i] = crq.fit[1,2]
    b1.crq.33[i] = crq.fit[2,1]
    b1.crq.sd.33[i] = crq.fit[2,2]
    cover.crq.33[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.33[i,2]=crq.fit[1,1]
    cover.crq.33[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.33[i,4]=ind(1.040613,cover.crq.33[i,1],cover.crq.33[i,3])
    cover.crq.33[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.33[i,6]=crq.fit[2,1]
    cover.crq.33[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.33[i,8]=ind(0,cover.crq.33[i,5],cover.crq.33[i,7])}
    ,error=function(e){
      b0.crq.33[i] = NA
      b0.crq.sd.33[i] = NA
      b1.crq.33[i] = NA
      b1.crq.sd.33[i] = NA
      cover.crq.33[i,1]=NA
      cover.crq.33[i,2]=NA
      cover.crq.33[i,3]=NA
      cover.crq.33[i,4]=NA
      cover.crq.33[i,5]=NA
      cover.crq.33[i,6]=NA
      cover.crq.33[i,7]=NA
      cover.crq.33[i,8]=NA
    })
  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.rq.33[i] = rq.fit[1,1]
    b0.rq.sd.33[i] = rq.fit[1,2]
    b1.rq.33[i] = rq.fit[2,1]
    b1.rq.sd.33[i] = rq.fit[2,2]
    cover.rq.33[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.33[i,2]=rq.fit[1,1]
    cover.rq.33[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.33[i,4]=ind(1.040613,cover.rq.33[i,1],cover.rq.33[i,3])
    cover.rq.33[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.33[i,6]=rq.fit[2,1]
    cover.rq.33[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.33[i,8]=ind(0,cover.rq.33[i,5],cover.rq.33[i,7])}
    ,error=function(e){
      b0.rq.33[i] = NA
      b0.rq.sd.33[i] = NA
      b1.rq.33[i] = NA
      b1.rq.sd.33[i] = NA
      cover.rq.33[i,1]=NA
      cover.rq.33[i,2]=NA
      cover.rq.33[i,3]=NA
      cover.rq.33[i,4]=NA
      cover.rq.33[i,5]=NA
      cover.rq.33[i,6]=NA
      cover.rq.33[i,7]=NA
      cover.rq.33[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 3, 0.75, 100)
    b0.is.33[i] = ismb.fit[1,1]
    b0.is.sd.33[i] = ismb.fit[1,2]
    b1.is.33[i] = ismb.fit[2,1]
    b1.is.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.is.33[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.33[i,2]=ismb.fit[1,1]
    cover.is.33[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.33[i,4]=ind(1.040613,cover.is.33[i,1],cover.is.33[i,3])
    cover.is.33[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.33[i,6]=ismb.fit[2,1]
    cover.is.33[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.33[i,8]=ind(0,cover.is.33[i,5],cover.is.33[i,7])}
    , error=function(e){
      b0.is.33[i] = NA
      b0.is.sd.33[i] = NA
      b1.is.33[i] = NA
      b1.is.sd.33[i] = NA
      # Coverage
      cover.is.33[i,1]=NA
      cover.is.33[i,2]=NA
      cover.is.33[i,3]=NA
      cover.is.33[i,4]=NA
      cover.is.33[i,5]=NA
      cover.is.33[i,6]=NA
      cover.is.33[i,7]=NA
      cover.is.33[i,8]=NA
    })
}
# Crq beta table
table3.crq[3,1]<-mean(b0.crq.33,na.rm=TRUE)
table3.crq[3,2]<-mean(b0.crq.sd.33,na.rm=TRUE)
table3.crq[3,3]<-sd(b0.crq.33,na.rm=TRUE)
table3.crq[3,4]<-mean(cover.crq.33[,4],na.rm=TRUE)
table3.crq[3,5]<-mean(b1.crq.33,na.rm=TRUE)
table3.crq[3,6]<-mean(b1.crq.sd.33,na.rm=TRUE)
table3.crq[3,7]<-sd(b1.crq.33,na.rm=TRUE)
table3.crq[3,8]<-mean(cover.crq.33[,8],na.rm=TRUE)

# rq beta table
table3.rq[3,1]<-mean(b0.rq.33,na.rm=TRUE)
table3.rq[3,2]<-mean(b0.rq.sd.33,na.rm=TRUE)
table3.rq[3,3]<-sd(b0.rq.33,na.rm=TRUE)
table3.rq[3,4]<-mean(cover.rq.33[,4],na.rm=TRUE)
table3.rq[3,5]<-mean(b1.rq.33,na.rm=TRUE)
table3.rq[3,6]<-mean(b1.rq.sd.33,na.rm=TRUE)
table3.rq[3,7]<-sd(b1.rq.33,na.rm=TRUE)
table3.rq[3,8]<-mean(cover.rq.33[,8],na.rm=TRUE)

# IS beta table
table3.is[3,1]<-mean(b0.is.33,na.rm=TRUE)
table3.is[3,2]<-mean(b0.is.sd.33,na.rm=TRUE)
table3.is[3,3]<-sd(b0.is.33,na.rm=TRUE)
table3.is[3,4]<-mean(cover.is.33[,4],na.rm=TRUE)
table3.is[3,5]<-mean(b1.is.33,na.rm=TRUE)
table3.is[3,6]<-mean(b1.is.sd.33,na.rm=TRUE)
table3.is[3,7]<-sd(b1.is.33,na.rm=TRUE)
table3.is[3,8]<-mean(cover.is.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.crq.35<-c()
b0.crq.sd.35<-c()
b1.crq.35<-c()
b1.crq.sd.35<-c()
b0.rq.35<-c()
b0.rq.sd.35<-c()
b1.rq.35<-c()
b1.rq.sd.35<-c()
b0.is.35<-c()
b0.is.sd.35<-c()
b1.is.35<-c()
b1.is.sd.35<-c()
cover.crq.35=matrix(NA,2000,8)
colnames(cover.crq.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.35=matrix(NA,2000,8)
colnames(cover.rq.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.35=matrix(NA,2000,8)
colnames(cover.is.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  delnumt_3 = c(delnumt_3, sum(a$Z<3))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.crq.35[i] = crq.fit[1,1]
    b0.crq.sd.35[i] = crq.fit[1,2]
    b1.crq.35[i] = crq.fit[2,1]
    b1.crq.sd.35[i] = crq.fit[2,2]
    cover.crq.35[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.35[i,2]=crq.fit[1,1]
    cover.crq.35[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.35[i,4]=ind(1.040613, cover.crq.35[i,1], cover.crq.35[i,3])
    cover.crq.35[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.35[i,6]=crq.fit[2,1]
    cover.crq.35[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.35[i,8]=ind(0, cover.crq.35[i,5], cover.crq.35[i,7])}
    ,error=function(e){
      b0.crq.35[i] = NA
      b0.crq.sd.35[i] = NA
      b1.crq.35[i] = NA
      b1.crq.sd.35[i] = NA
      cover.crq.35[i,1]=NA
      cover.crq.35[i,2]=NA
      cover.crq.35[i,3]=NA
      cover.crq.35[i,4]=NA
      cover.crq.35[i,5]=NA
      cover.crq.35[i,6]=NA
      cover.crq.35[i,7]=NA
      cover.crq.35[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.rq.35[i] = rq.fit[1,1]
    b0.rq.sd.35[i] = rq.fit[1,2]
    b1.rq.35[i] = rq.fit[2,1]
    b1.rq.sd.35[i] = rq.fit[2,2]
    cover.rq.35[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.35[i,2]=rq.fit[1,1]
    cover.rq.35[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.35[i,4]=ind(1.040613, cover.rq.35[i,1], cover.rq.35[i,3])
    cover.rq.35[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.35[i,6]=rq.fit[2,1]
    cover.rq.35[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.35[i,8]=ind(0, cover.rq.35[i,5], cover.rq.35[i,7])}
    ,error=function(e){
      b0.rq.35[i] = NA
      b0.rq.sd.35[i] = NA
      b1.rq.35[i] = NA
      b1.rq.sd.35[i] = NA
      cover.rq.35[i,1]=NA
      cover.rq.35[i,2]=NA
      cover.rq.35[i,3]=NA
      cover.rq.35[i,4]=NA
      cover.rq.35[i,5]=NA
      cover.rq.35[i,6]=NA
      cover.rq.35[i,7]=NA
      cover.rq.35[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 3, 0.75, 100)
    b0.is.35[i] = ismb.fit[1,1]
    b0.is.sd.35[i] = ismb.fit[1,2]
    b1.is.35[i] = ismb.fit[2,1]
    b1.is.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.is.35[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.35[i,2]=ismb.fit[1,1]
    cover.is.35[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.35[i,4]=ind(1.040613,cover.is.35[i,1],cover.is.35[i,3])
    cover.is.35[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.35[i,6]=ismb.fit[2,1]
    cover.is.35[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.35[i,8]=ind(0,cover.is.35[i,5],cover.is.35[i,7])}
    , error=function(e){
      b0.is.35[i] = NA
      b0.is.sd.35[i] = NA
      b1.is.35[i] = NA
      b1.is.sd.35[i] = NA
      # Coverage
      cover.is.35[i,1]=NA
      cover.is.35[i,2]=NA
      cover.is.35[i,3]=NA
      cover.is.35[i,4]=NA
      cover.is.35[i,5]=NA
      cover.is.35[i,6]=NA
      cover.is.35[i,7]=NA
      cover.is.35[i,8]=NA
    })
}

# Crq beta table
table3.crq[4,1]<-mean(b0.crq.35,na.rm=TRUE)
table3.crq[4,2]<-mean(b0.crq.sd.35,na.rm=TRUE)
table3.crq[4,3]<-sd(b0.crq.35,na.rm=TRUE)
table3.crq[4,4]<-mean(cover.crq.35[,4],na.rm=TRUE)
table3.crq[4,5]<-mean(b1.crq.35,na.rm=TRUE)
table3.crq[4,6]<-mean(b1.crq.sd.35,na.rm=TRUE)
table3.crq[4,7]<-sd(b1.crq.35,na.rm=TRUE)
table3.crq[4,8]<-mean(cover.crq.35[,8],na.rm=TRUE)

# rq beta table
table3.rq[4,1]<-mean(b0.rq.35,na.rm=TRUE)
table3.rq[4,2]<-mean(b0.rq.sd.35,na.rm=TRUE)
table3.rq[4,3]<-sd(b0.rq.35,na.rm=TRUE)
table3.rq[4,4]<-mean(cover.rq.35[,4],na.rm=TRUE)
table3.rq[4,5]<-mean(b1.rq.35,na.rm=TRUE)
table3.rq[4,6]<-mean(b1.rq.sd.35,na.rm=TRUE)
table3.rq[4,7]<-sd(b1.rq.35,na.rm=TRUE)
table3.rq[4,8]<-mean(cover.rq.35[,8],na.rm=TRUE)

# IS beta table
table3.is[4,1]<-mean(b0.is.35,na.rm=TRUE)
table3.is[4,2]<-mean(b0.is.sd.35,na.rm=TRUE)
table3.is[4,3]<-sd(b0.is.35,na.rm=TRUE)
table3.is[4,4]<-mean(cover.is.35[,4],na.rm=TRUE)
table3.is[4,5]<-mean(b1.is.35,na.rm=TRUE)
table3.is[4,6]<-mean(b1.is.sd.35,na.rm=TRUE)
table3.is[4,7]<-sd(b1.is.35,na.rm=TRUE)
table3.is[4,8]<-mean(cover.is.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.rq.37<-c()
b0.rq.sd.37<-c()
b1.rq.37<-c()
b1.rq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
cover.crq.37=matrix(NA,2000,8)
colnames(cover.crq.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.rq.37=matrix(NA,2000,8)
colnames(cover.rq.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.37=matrix(NA,2000,8)
colnames(cover.is.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  delnumt_3 = c(delnumt_3, sum(a$Z<3))
  # Beta estimation method 1 : Crq package
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.crq.37[i] = crq.fit[1,1]
    b0.crq.sd.37[i] = crq.fit[1,2]
    b1.crq.37[i] = crq.fit[2,1]
    b1.crq.sd.37[i] = crq.fit[2,2]
    cover.crq.37[i,1]=crq.fit[1,1]-1.96*crq.fit[1,2]
    cover.crq.37[i,2]=crq.fit[1,1]
    cover.crq.37[i,3]=crq.fit[1,1]+1.96*crq.fit[1,2]
    cover.crq.37[i,4]=ind(1.040613,cover.crq.37[i,1],cover.crq.37[i,3])
    cover.crq.37[i,5]=crq.fit[2,1]-1.96*crq.fit[2,2]
    cover.crq.37[i,6]=crq.fit[2,1]
    cover.crq.37[i,7]=crq.fit[2,1]+1.96*crq.fit[2,2]
    cover.crq.37[i,8]=ind(0,cover.crq.37[i,5],cover.crq.37[i,7])}
    ,error=function(e){
      b0.crq.37[i] = NA
      b0.crq.sd.37[i] = NA
      b1.crq.37[i] = NA
      b1.crq.sd.37[i] = NA
      cover.crq.37[i,1]=NA
      cover.crq.37[i,2]=NA
      cover.crq.37[i,3]=NA
      cover.crq.37[i,4]=NA
      cover.crq.37[i,5]=NA
      cover.crq.37[i,6]=NA
      cover.crq.37[i,7]=NA
      cover.crq.37[i,8]=NA
    })

  # Beta estimation method 2 : rq with weight
  tryCatch({
    rq.fit = rq.est(a[,3], 1, a[,4], a[,5], 3, 0.75)
    b0.rq.37[i] = rq.fit[1,1]
    b0.rq.sd.37[i] = rq.fit[1,2]
    b1.rq.37[i] = rq.fit[2,1]
    b1.rq.sd.37[i] = rq.fit[2,2]
    cover.rq.37[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.37[i,2]=rq.fit[1,1]
    cover.rq.37[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.37[i,4]=ind(1.040613,cover.rq.37[i,1],cover.rq.37[i,3])
    cover.rq.37[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.37[i,6]=rq.fit[2,1]
    cover.rq.37[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.37[i,8]=ind(0,cover.rq.37[i,5],cover.rq.37[i,7])}
    ,error=function(e){
      b0.rq.37[i] = NA
      b0.rq.sd.37[i] = NA
      b1.rq.37[i] = NA
      b1.rq.sd.37[i] = NA
      cover.rq.37[i,1]=NA
      cover.rq.37[i,2]=NA
      cover.rq.37[i,3]=NA
      cover.rq.37[i,4]=NA
      cover.rq.37[i,5]=NA
      cover.rq.37[i,6]=NA
      cover.rq.37[i,7]=NA
      cover.rq.37[i,8]=NA
    })

  # Beta estimation method 3: ISMB
  tryCatch({
    ismb.fit = ismbwo.est(a[,3], 1, a[,4], a[,5], 3, 0.75, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.is.37[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.37[i,2]=ismb.fit[1,1]
    cover.is.37[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.37[i,4]=ind(1.040613,cover.is.37[i,1],cover.is.37[i,3])
    cover.is.37[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.37[i,6]=ismb.fit[2,1]
    cover.is.37[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.37[i,8]=ind(0,cover.is.37[i,5],cover.is.37[i,7])}
    , error=function(e){
      b0.is.37[i] = NA
      b0.is.sd.37[i] = NA
      b1.is.37[i] = NA
      b1.is.sd.37[i] = NA
      # Coverage
      cover.is.37[i,1]=NA
      cover.is.37[i,2]=NA
      cover.is.37[i,3]=NA
      cover.is.37[i,4]=NA
      cover.is.37[i,5]=NA
      cover.is.37[i,6]=NA
      cover.is.37[i,7]=NA
      cover.is.37[i,8]=NA
    })
}
# Crq beta table
table3.crq[5,1]<-mean(b0.crq.37,na.rm=TRUE)
table3.crq[5,2]<-mean(b0.crq.sd.37,na.rm=TRUE)
table3.crq[5,3]<-sd(b0.crq.37,na.rm=TRUE)
table3.crq[5,4]<-mean(cover.crq.37[,4],na.rm=TRUE)
table3.crq[5,5]<-mean(b1.crq.37,na.rm=TRUE)
table3.crq[5,6]<-mean(b1.crq.sd.37,na.rm=TRUE)
table3.crq[5,7]<-sd(b1.crq.37,na.rm=TRUE)
table3.crq[5,8]<-mean(cover.crq.37[,8],na.rm=TRUE)

# rq beta table
table3.rq[5,1]<-mean(b0.rq.37,na.rm=TRUE)
table3.rq[5,2]<-mean(b0.rq.sd.37,na.rm=TRUE)
table3.rq[5,3]<-sd(b0.rq.37,na.rm=TRUE)
table3.rq[5,4]<-mean(cover.rq.37[,4],na.rm=TRUE)
table3.rq[5,5]<-mean(b1.rq.37,na.rm=TRUE)
table3.rq[5,6]<-mean(b1.rq.sd.37,na.rm=TRUE)
table3.rq[5,7]<-sd(b1.rq.37,na.rm=TRUE)
table3.rq[5,8]<-mean(cover.rq.37[,8],na.rm=TRUE)

# IS beta table
table3.is[5,1]<-mean(b0.is.37,na.rm=TRUE)
table3.is[5,2]<-mean(b0.is.sd.37,na.rm=TRUE)
table3.is[5,3]<-sd(b0.is.37,na.rm=TRUE)
table3.is[5,4]<-mean(cover.is.37[,4],na.rm=TRUE)
table3.is[5,5]<-mean(b1.is.37,na.rm=TRUE)
table3.is[5,6]<-mean(b1.is.sd.37,na.rm=TRUE)
table3.is[5,7]<-sd(b1.is.37,na.rm=TRUE)
table3.is[5,8]<-mean(cover.is.37[,8],na.rm=TRUE)

sum(is.na(b0.is.05))
sum(is.na(b0.is.07))
sum(is.na(b0.is.15))
sum(is.na(b0.is.17))
sum(is.na(b0.is.25))
sum(is.na(b0.is.27))
sum(is.na(b0.is.35))
sum(is.na(b0.is.37))
mean(delnumt_1)
mean(delnumt_2)
mean(delnumt_3)
