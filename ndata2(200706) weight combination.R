#### ndata2(200706) : Condition ####
# data size = 200
# beta0, beta1 effective
# Quantile 50%
# simulation 2000
# eta = 100
# Use WKM$surv as weight
# Weight combination 1 : weight-in & G(t)
# Weight combination 2 : weight-in & G(t|t>t_0)
# Weight combination 3 : weight-out & G(t)
# Weight combination 4 : weight-out & G(t|t>t_0)

library(survival)
library(nleqslv)
library(xtable)
library(emplik)

#### True Beta ####
#beta_0    beta_1
#t_0=0 1.609438 0.6931472
#t_0=1 1.410748 0.7974189
#t_0=2 1.219403 0.9070615
#t_0=3 1.040613 1.0174711


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
    if (sim[q,4]==0){
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
    } else {
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.1
    }
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
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/5))^(1/k)/exp.beta.initial.1

#### Example parameter ####
# c.0=5000000
# c.1=61.17
# c.3=22.53
# c.5=13.43
# c.7=8.5
# set.seed(1)
# a<-data.gen(200,c.3)
# Z=a[,3]
# nc=1
# covariate=a[,4]
# D=a[,5]
# t_0=3
# Q=0.5
# ne=100

#### Estimation functions ####
# 1. ISMB with Weight combination 1 : weight-in & G(t)
ismbw1.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data) #head(data,20)

  # G(t) weight using WKM jump
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  data[,(nc+6)] = fit$jump*n
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
  betastart = c(1,1)
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

# 2. ISMB with Weight combination 2 : weight-in & G(T|T>t_0)
ismbw2.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data) #head(data,20)
  
  # using valid observations only (Z > t0)
  data2 <- data[data[,2]!=-10,] #head(data2,20)
  n2 <- nrow(data2)
  
  # weight using WKM jump
  fit = WKM(data2[,1], (data2[,(nc+4)]), zc = rep(1,n2), w = rep(1,n2))
  data2[,(nc+6)] = fit$jump*n2
  m = nrow(data2)
  colnames(data2)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data2)[j]="covariate"
  }
  colnames(data2)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  
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
  X = as.matrix(cbind(c(rep(1,n2)),data2[,4:(nc+3)]))
  W = data2[,(nc+6)]
  logT = data2[,2]
  I = data2[,3]
  H = diag(1/n, nc+1, nc+1)
  
  # Change betastart when real data analysis c(1,rep(1,nc))
  betastart = c(1,1)
  is.fit = nleqslv(betastart, objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n2,1)
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

# 3. ISMB with Weight combination 3 : weight-out & G(t)
ismbw3.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[200,(nc+4)] = 1
  data = as.data.frame(data)

  # weight using WKM jump
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  data[,(nc+6)] = fit$jump
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
  betastart = c(1,1)
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

# 4. ISMB with Weight combination 4 : weight-out & G(t|t>t_0)
ismbw4.est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data) #head(data,20)
  
  # using valid observations only (Z > t0)
  data2 <- data[data[,2]!=-10,] #head(data2,20)
  n2 <- nrow(data2)
  
  # weight using WKM jump
  fit = WKM(data2[,1], (data2[,(nc+4)]), zc = rep(1,n2), w = rep(1,n2))
  data2[,(nc+6)] = fit$jump*n2
  m = nrow(data2)
  colnames(data2)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data2)[j]="covariate"
  }
  colnames(data2)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  
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
  X = as.matrix(cbind(c(rep(1,n2)),data2[,4:(nc+3)]))
  W = data2[,(nc+6)]
  logT = data2[,2]
  I = data2[,3]
  H = diag(1/n, nc+1, nc+1)
  
  # Change betastart when real data analysis c(1,rep(1,nc))
  betastart = c(1,1)
  is.fit = nleqslv(betastart, objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n2,1)
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
table0.isw1<-matrix(NA,5,8)
rownames(table0.isw1)<-c(0,10,30,50,70)
colnames(table0.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw1<-matrix(NA,5,8)
rownames(table1.isw1)<-c(0,10,30,50,70)
colnames(table1.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw1<-matrix(NA,5,8)
rownames(table2.isw1)<-c(0,10,30,50,70)
colnames(table2.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw1<-matrix(NA,5,8)
rownames(table3.isw1)<-c(0,10,30,50,70)
colnames(table3.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.isw2<-matrix(NA,5,8)
rownames(table0.isw2)<-c(0,10,30,50,70)
colnames(table0.isw2)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw2<-matrix(NA,5,8)
rownames(table1.isw2)<-c(0,10,30,50,70)
colnames(table1.isw2)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw2<-matrix(NA,5,8)
rownames(table2.isw2)<-c(0,10,30,50,70)
colnames(table2.isw2)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw2<-matrix(NA,5,8)
rownames(table3.isw2)<-c(0,10,30,50,70)
colnames(table3.isw2)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.isw3<-matrix(NA,5,8)
rownames(table0.isw3)<-c(0,10,30,50,70)
colnames(table0.isw3)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw3<-matrix(NA,5,8)
rownames(table1.isw3)<-c(0,10,30,50,70)
colnames(table1.isw3)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw3<-matrix(NA,5,8)
rownames(table2.isw3)<-c(0,10,30,50,70)
colnames(table2.isw3)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw3<-matrix(NA,5,8)
rownames(table3.isw3)<-c(0,10,30,50,70)
colnames(table3.isw3)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.isw4<-matrix(NA,5,8)
rownames(table0.isw4)<-c(0,10,30,50,70)
colnames(table0.isw4)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw4<-matrix(NA,5,8)
rownames(table1.isw4)<-c(0,10,30,50,70)
colnames(table1.isw4)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw4<-matrix(NA,5,8)
rownames(table2.isw4)<-c(0,10,30,50,70)
colnames(table2.isw4)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw4<-matrix(NA,5,8)
rownames(table3.isw4)<-c(0,10,30,50,70)
colnames(table3.isw4)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

#### censoring point at t_0=0 ####
c.0=5000000
c.1=78.11
c.3=26.36
c.5=15.08
c.7=9.09
#### t_0=0 & c=0% ####
b0.isw1.00 = c()
b0.isw1.sd.00 = c()
b1.isw1.00 = c()
b1.isw1.sd.00 = c()
cover.isw1.00=matrix(NA,2000,8)
colnames(cover.isw1.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.00 = c()
b0.isw2.sd.00 = c()
b1.isw2.00 = c()
b1.isw2.sd.00 = c()
cover.isw2.00=matrix(NA,2000,8)
colnames(cover.isw2.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.00 = c()
b0.isw3.sd.00 = c()
b1.isw3.00 = c()
b1.isw3.sd.00 = c()
cover.isw3.00=matrix(NA,2000,8)
colnames(cover.isw3.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.00 = c()
b0.isw4.sd.00 = c()
b1.isw4.00 = c()
b1.isw4.sd.00 = c()
cover.isw4.00=matrix(NA,2000,8)
colnames(cover.isw4.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.00[i] = ismb.fit[1,1]
    b0.isw1.sd.00[i] = ismb.fit[1,2]
    b1.isw1.00[i] = ismb.fit[2,1]
    b1.isw1.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.00[i,2] = ismb.fit[1,1]
    cover.isw1.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.00[i,4] = ind(1.609438, cover.isw1.00[i,1], cover.isw1.00[i,3])
    cover.isw1.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.00[i,6] = ismb.fit[2,1]
    cover.isw1.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.00[i,8] = ind(0.6931472, cover.isw1.00[i,5], cover.isw1.00[i,7])}
    , error=function(e){
      b0.isw1.00[i] = NA
      b0.isw1.sd.00[i] = NA
      b1.isw1.00[i] = NA
      b1.isw1.sd.00[i] = NA
      # Coverage
      cover.isw1.00[i,1] = NA
      cover.isw1.00[i,2] = NA
      cover.isw1.00[i,3] = NA
      cover.isw1.00[i,4] = NA
      cover.isw1.00[i,5] = NA
      cover.isw1.00[i,6] = NA
      cover.isw1.00[i,7] = NA
      cover.isw1.00[i,8] = NA
    })

  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw2.00[i] = ismb.fit[1,1]
    b0.isw2.sd.00[i] = ismb.fit[1,2]
    b1.isw2.00[i] = ismb.fit[2,1]
    b1.isw2.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.00[i,2] = ismb.fit[1,1]
    cover.isw2.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.00[i,4] = ind(1.609438, cover.isw2.00[i,1], cover.isw2.00[i,3])
    cover.isw2.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.00[i,6] = ismb.fit[2,1]
    cover.isw2.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.00[i,8] = ind(0.6931472, cover.isw2.00[i,5], cover.isw2.00[i,7])}
    , error=function(e){
      b0.isw2.00[i] = NA
      b0.isw2.sd.00[i] = NA
      b1.isw2.00[i] = NA
      b1.isw2.sd.00[i] = NA
      # Coverage
      cover.isw2.00[i,1] = NA
      cover.isw2.00[i,2] = NA
      cover.isw2.00[i,3] = NA
      cover.isw2.00[i,4] = NA
      cover.isw2.00[i,5] = NA
      cover.isw2.00[i,6] = NA
      cover.isw2.00[i,7] = NA
      cover.isw2.00[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw3.00[i] = ismb.fit[1,1]
    b0.isw3.sd.00[i] = ismb.fit[1,2]
    b1.isw3.00[i] = ismb.fit[2,1]
    b1.isw3.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.00[i,2] = ismb.fit[1,1]
    cover.isw3.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.00[i,4] = ind(1.609438, cover.isw3.00[i,1], cover.isw3.00[i,3])
    cover.isw3.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.00[i,6] = ismb.fit[2,1]
    cover.isw3.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.00[i,8] = ind(0.6931472, cover.isw3.00[i,5], cover.isw3.00[i,7])}
    , error=function(e){
      b0.isw3.00[i] = NA
      b0.isw3.sd.00[i] = NA
      b1.isw3.00[i] = NA
      b1.isw3.sd.00[i] = NA
      # Coverage
      cover.isw3.00[i,1] = NA
      cover.isw3.00[i,2] = NA
      cover.isw3.00[i,3] = NA
      cover.isw3.00[i,4] = NA
      cover.isw3.00[i,5] = NA
      cover.isw3.00[i,6] = NA
      cover.isw3.00[i,7] = NA
      cover.isw3.00[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw4.00[i] = ismb.fit[1,1]
    b0.isw4.sd.00[i] = ismb.fit[1,2]
    b1.isw4.00[i] = ismb.fit[2,1]
    b1.isw4.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.00[i,2] = ismb.fit[1,1]
    cover.isw4.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.00[i,4] = ind(1.609438, cover.isw4.00[i,1], cover.isw4.00[i,3])
    cover.isw4.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.00[i,6] = ismb.fit[2,1]
    cover.isw4.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.00[i,8] = ind(0.6931472, cover.isw4.00[i,5], cover.isw4.00[i,7])}
    , error=function(e){
      b0.isw4.00[i] = NA
      b0.isw4.sd.00[i] = NA
      b1.isw4.00[i] = NA
      b1.isw4.sd.00[i] = NA
      # Coverage
      cover.isw4.00[i,1] = NA
      cover.isw4.00[i,2] = NA
      cover.isw4.00[i,3] = NA
      cover.isw4.00[i,4] = NA
      cover.isw4.00[i,5] = NA
      cover.isw4.00[i,6] = NA
      cover.isw4.00[i,7] = NA
      cover.isw4.00[i,8] = NA
    })
}

# IS beta table
table0.isw1[1,1]<-mean(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,2]<-mean(b0.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,3]<-sd(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,4]<-mean(cover.isw1.00[,4],na.rm=TRUE)
table0.isw1[1,5]<-mean(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,6]<-mean(b1.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,7]<-sd(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,8]<-mean(cover.isw1.00[,8],na.rm=TRUE)

table0.isw2[1,1]<-mean(b0.isw2.00,na.rm=TRUE)
table0.isw2[1,2]<-mean(b0.isw2.sd.00,na.rm=TRUE)
table0.isw2[1,3]<-sd(b0.isw2.00,na.rm=TRUE)
table0.isw2[1,4]<-mean(cover.isw2.00[,4],na.rm=TRUE)
table0.isw2[1,5]<-mean(b1.isw2.00,na.rm=TRUE)
table0.isw2[1,6]<-mean(b1.isw2.sd.00,na.rm=TRUE)
table0.isw2[1,7]<-sd(b1.isw2.00,na.rm=TRUE)
table0.isw2[1,8]<-mean(cover.isw2.00[,8],na.rm=TRUE)

table0.isw3[1,1]<-mean(b0.isw3.00,na.rm=TRUE)
table0.isw3[1,2]<-mean(b0.isw3.sd.00,na.rm=TRUE)
table0.isw3[1,3]<-sd(b0.isw3.00,na.rm=TRUE)
table0.isw3[1,4]<-mean(cover.isw3.00[,4],na.rm=TRUE)
table0.isw3[1,5]<-mean(b1.isw3.00,na.rm=TRUE)
table0.isw3[1,6]<-mean(b1.isw3.sd.00,na.rm=TRUE)
table0.isw3[1,7]<-sd(b1.isw3.00,na.rm=TRUE)
table0.isw3[1,8]<-mean(cover.isw3.00[,8],na.rm=TRUE)

table0.isw4[1,1]<-mean(b0.isw4.00,na.rm=TRUE)
table0.isw4[1,2]<-mean(b0.isw4.sd.00,na.rm=TRUE)
table0.isw4[1,3]<-sd(b0.isw4.00,na.rm=TRUE)
table0.isw4[1,4]<-mean(cover.isw4.00[,4],na.rm=TRUE)
table0.isw4[1,5]<-mean(b1.isw4.00,na.rm=TRUE)
table0.isw4[1,6]<-mean(b1.isw4.sd.00,na.rm=TRUE)
table0.isw4[1,7]<-sd(b1.isw4.00,na.rm=TRUE)
table0.isw4[1,8]<-mean(cover.isw4.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.isw1.01 = c()
b0.isw1.sd.01 = c()
b1.isw1.01 = c()
b1.isw1.sd.01 = c()
cover.isw1.01=matrix(NA,2000,8)
colnames(cover.isw1.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.01 = c()
b0.isw2.sd.01 = c()
b1.isw2.01 = c()
b1.isw2.sd.01 = c()
cover.isw2.01=matrix(NA,2000,8)
colnames(cover.isw2.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.01 = c()
b0.isw3.sd.01 = c()
b1.isw3.01 = c()
b1.isw3.sd.01 = c()
cover.isw3.01=matrix(NA,2000,8)
colnames(cover.isw3.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.01 = c()
b0.isw4.sd.01 = c()
b1.isw4.01 = c()
b1.isw4.sd.01 = c()
cover.isw4.01=matrix(NA,2000,8)
colnames(cover.isw4.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.01[i] = ismb.fit[1,1]
    b0.isw1.sd.01[i] = ismb.fit[1,2]
    b1.isw1.01[i] = ismb.fit[2,1]
    b1.isw1.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.01[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.01[i,2] = ismb.fit[1,1]
    cover.isw1.01[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.01[i,4] = ind(1.609438, cover.isw1.01[i,1], cover.isw1.01[i,3])
    cover.isw1.01[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.01[i,6] = ismb.fit[2,1]
    cover.isw1.01[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.01[i,8] = ind(0.6931472, cover.isw1.01[i,5], cover.isw1.01[i,7])}
    , error=function(e){
      b0.isw1.01[i] = NA
      b0.isw1.sd.01[i] = NA
      b1.isw1.01[i] = NA
      b1.isw1.sd.01[i] = NA
      # Coverage
      cover.isw1.01[i,1] = NA
      cover.isw1.01[i,2] = NA
      cover.isw1.01[i,3] = NA
      cover.isw1.01[i,4] = NA
      cover.isw1.01[i,5] = NA
      cover.isw1.01[i,6] = NA
      cover.isw1.01[i,7] = NA
      cover.isw1.01[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw2.01[i] = ismb.fit[1,1]
    b0.isw2.sd.01[i] = ismb.fit[1,2]
    b1.isw2.01[i] = ismb.fit[2,1]
    b1.isw2.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.01[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.01[i,2] = ismb.fit[1,1]
    cover.isw2.01[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.01[i,4] = ind(1.609438, cover.isw2.01[i,1], cover.isw2.01[i,3])
    cover.isw2.01[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.01[i,6] = ismb.fit[2,1]
    cover.isw2.01[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.01[i,8] = ind(0.6931472, cover.isw2.01[i,5], cover.isw2.01[i,7])}
    , error=function(e){
      b0.isw2.01[i] = NA
      b0.isw2.sd.01[i] = NA
      b1.isw2.01[i] = NA
      b1.isw2.sd.01[i] = NA
      # Coverage
      cover.isw2.01[i,1] = NA
      cover.isw2.01[i,2] = NA
      cover.isw2.01[i,3] = NA
      cover.isw2.01[i,4] = NA
      cover.isw2.01[i,5] = NA
      cover.isw2.01[i,6] = NA
      cover.isw2.01[i,7] = NA
      cover.isw2.01[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw3.01[i] = ismb.fit[1,1]
    b0.isw3.sd.01[i] = ismb.fit[1,2]
    b1.isw3.01[i] = ismb.fit[2,1]
    b1.isw3.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.01[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.01[i,2] = ismb.fit[1,1]
    cover.isw3.01[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.01[i,4] = ind(1.609438, cover.isw3.01[i,1], cover.isw3.01[i,3])
    cover.isw3.01[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.01[i,6] = ismb.fit[2,1]
    cover.isw3.01[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.01[i,8] = ind(0.6931472, cover.isw3.01[i,5], cover.isw3.01[i,7])}
    , error=function(e){
      b0.isw3.01[i] = NA
      b0.isw3.sd.01[i] = NA
      b1.isw3.01[i] = NA
      b1.isw3.sd.01[i] = NA
      # Coverage
      cover.isw3.01[i,1] = NA
      cover.isw3.01[i,2] = NA
      cover.isw3.01[i,3] = NA
      cover.isw3.01[i,4] = NA
      cover.isw3.01[i,5] = NA
      cover.isw3.01[i,6] = NA
      cover.isw3.01[i,7] = NA
      cover.isw3.01[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw4.01[i] = ismb.fit[1,1]
    b0.isw4.sd.01[i] = ismb.fit[1,2]
    b1.isw4.01[i] = ismb.fit[2,1]
    b1.isw4.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.01[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.01[i,2] = ismb.fit[1,1]
    cover.isw4.01[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.01[i,4] = ind(1.609438, cover.isw4.01[i,1], cover.isw4.01[i,3])
    cover.isw4.01[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.01[i,6] = ismb.fit[2,1]
    cover.isw4.01[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.01[i,8] = ind(0.6931472, cover.isw4.01[i,5], cover.isw4.01[i,7])}
    , error=function(e){
      b0.isw4.01[i] = NA
      b0.isw4.sd.01[i] = NA
      b1.isw4.01[i] = NA
      b1.isw4.sd.01[i] = NA
      # Coverage
      cover.isw4.01[i,1] = NA
      cover.isw4.01[i,2] = NA
      cover.isw4.01[i,3] = NA
      cover.isw4.01[i,4] = NA
      cover.isw4.01[i,5] = NA
      cover.isw4.01[i,6] = NA
      cover.isw4.01[i,7] = NA
      cover.isw4.01[i,8] = NA
    })
}

# IS beta table
table0.isw1[2,1]<-mean(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,2]<-mean(b0.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,3]<-sd(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,4]<-mean(cover.isw1.01[,4],na.rm=TRUE)
table0.isw1[2,5]<-mean(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,6]<-mean(b1.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,7]<-sd(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,8]<-mean(cover.isw1.01[,8],na.rm=TRUE)

table0.isw2[2,1]<-mean(b0.isw2.01,na.rm=TRUE)
table0.isw2[2,2]<-mean(b0.isw2.sd.01,na.rm=TRUE)
table0.isw2[2,3]<-sd(b0.isw2.01,na.rm=TRUE)
table0.isw2[2,4]<-mean(cover.isw2.01[,4],na.rm=TRUE)
table0.isw2[2,5]<-mean(b1.isw2.01,na.rm=TRUE)
table0.isw2[2,6]<-mean(b1.isw2.sd.01,na.rm=TRUE)
table0.isw2[2,7]<-sd(b1.isw2.01,na.rm=TRUE)
table0.isw2[2,8]<-mean(cover.isw2.01[,8],na.rm=TRUE)

table0.isw3[2,1]<-mean(b0.isw3.01,na.rm=TRUE)
table0.isw3[2,2]<-mean(b0.isw3.sd.01,na.rm=TRUE)
table0.isw3[2,3]<-sd(b0.isw3.01,na.rm=TRUE)
table0.isw3[2,4]<-mean(cover.isw3.01[,4],na.rm=TRUE)
table0.isw3[2,5]<-mean(b1.isw3.01,na.rm=TRUE)
table0.isw3[2,6]<-mean(b1.isw3.sd.01,na.rm=TRUE)
table0.isw3[2,7]<-sd(b1.isw3.01,na.rm=TRUE)
table0.isw3[2,8]<-mean(cover.isw3.01[,8],na.rm=TRUE)

table0.isw4[2,1]<-mean(b0.isw4.01,na.rm=TRUE)
table0.isw4[2,2]<-mean(b0.isw4.sd.01,na.rm=TRUE)
table0.isw4[2,3]<-sd(b0.isw4.01,na.rm=TRUE)
table0.isw4[2,4]<-mean(cover.isw4.01[,4],na.rm=TRUE)
table0.isw4[2,5]<-mean(b1.isw4.01,na.rm=TRUE)
table0.isw4[2,6]<-mean(b1.isw4.sd.01,na.rm=TRUE)
table0.isw4[2,7]<-sd(b1.isw4.01,na.rm=TRUE)
table0.isw4[2,8]<-mean(cover.isw4.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.isw1.03 = c()
b0.isw1.sd.03 = c()
b1.isw1.03 = c()
b1.isw1.sd.03 = c()
cover.isw1.03=matrix(NA,2000,8)
colnames(cover.isw1.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.03 = c()
b0.isw2.sd.03 = c()
b1.isw2.03 = c()
b1.isw2.sd.03 = c()
cover.isw2.03=matrix(NA,2000,8)
colnames(cover.isw2.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.03 = c()
b0.isw3.sd.03 = c()
b1.isw3.03 = c()
b1.isw3.sd.03 = c()
cover.isw3.03=matrix(NA,2000,8)
colnames(cover.isw3.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.03 = c()
b0.isw4.sd.03 = c()
b1.isw4.03 = c()
b1.isw4.sd.03 = c()
cover.isw4.03=matrix(NA,2000,8)
colnames(cover.isw4.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.3)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.03[i] = ismb.fit[1,1]
    b0.isw1.sd.03[i] = ismb.fit[1,2]
    b1.isw1.03[i] = ismb.fit[2,1]
    b1.isw1.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.03[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.03[i,2] = ismb.fit[1,1]
    cover.isw1.03[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.03[i,4] = ind(1.609438, cover.isw1.03[i,1], cover.isw1.03[i,3])
    cover.isw1.03[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.03[i,6] = ismb.fit[2,1]
    cover.isw1.03[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.03[i,8] = ind(0.6931472, cover.isw1.03[i,5], cover.isw1.03[i,7])}
    , error=function(e){
      b0.isw1.03[i] = NA
      b0.isw1.sd.03[i] = NA
      b1.isw1.03[i] = NA
      b1.isw1.sd.03[i] = NA
      # Coverage
      cover.isw1.03[i,1] = NA
      cover.isw1.03[i,2] = NA
      cover.isw1.03[i,3] = NA
      cover.isw1.03[i,4] = NA
      cover.isw1.03[i,5] = NA
      cover.isw1.03[i,6] = NA
      cover.isw1.03[i,7] = NA
      cover.isw1.03[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw2.03[i] = ismb.fit[1,1]
    b0.isw2.sd.03[i] = ismb.fit[1,2]
    b1.isw2.03[i] = ismb.fit[2,1]
    b1.isw2.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.03[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.03[i,2] = ismb.fit[1,1]
    cover.isw2.03[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.03[i,4] = ind(1.609438, cover.isw2.03[i,1], cover.isw2.03[i,3])
    cover.isw2.03[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.03[i,6] = ismb.fit[2,1]
    cover.isw2.03[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.03[i,8] = ind(0.6931472, cover.isw2.03[i,5], cover.isw2.03[i,7])}
    , error=function(e){
      b0.isw2.03[i] = NA
      b0.isw2.sd.03[i] = NA
      b1.isw2.03[i] = NA
      b1.isw2.sd.03[i] = NA
      # Coverage
      cover.isw2.03[i,1] = NA
      cover.isw2.03[i,2] = NA
      cover.isw2.03[i,3] = NA
      cover.isw2.03[i,4] = NA
      cover.isw2.03[i,5] = NA
      cover.isw2.03[i,6] = NA
      cover.isw2.03[i,7] = NA
      cover.isw2.03[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw3.03[i] = ismb.fit[1,1]
    b0.isw3.sd.03[i] = ismb.fit[1,2]
    b1.isw3.03[i] = ismb.fit[2,1]
    b1.isw3.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.03[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.03[i,2] = ismb.fit[1,1]
    cover.isw3.03[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.03[i,4] = ind(1.609438, cover.isw3.03[i,1], cover.isw3.03[i,3])
    cover.isw3.03[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.03[i,6] = ismb.fit[2,1]
    cover.isw3.03[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.03[i,8] = ind(0.6931472, cover.isw3.03[i,5], cover.isw3.03[i,7])}
    , error=function(e){
      b0.isw3.03[i] = NA
      b0.isw3.sd.03[i] = NA
      b1.isw3.03[i] = NA
      b1.isw3.sd.03[i] = NA
      # Coverage
      cover.isw3.03[i,1] = NA
      cover.isw3.03[i,2] = NA
      cover.isw3.03[i,3] = NA
      cover.isw3.03[i,4] = NA
      cover.isw3.03[i,5] = NA
      cover.isw3.03[i,6] = NA
      cover.isw3.03[i,7] = NA
      cover.isw3.03[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw4.03[i] = ismb.fit[1,1]
    b0.isw4.sd.03[i] = ismb.fit[1,2]
    b1.isw4.03[i] = ismb.fit[2,1]
    b1.isw4.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.03[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.03[i,2] = ismb.fit[1,1]
    cover.isw4.03[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.03[i,4] = ind(1.609438, cover.isw4.03[i,1], cover.isw4.03[i,3])
    cover.isw4.03[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.03[i,6] = ismb.fit[2,1]
    cover.isw4.03[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.03[i,8] = ind(0.6931472, cover.isw4.03[i,5], cover.isw4.03[i,7])}
    , error=function(e){
      b0.isw4.03[i] = NA
      b0.isw4.sd.03[i] = NA
      b1.isw4.03[i] = NA
      b1.isw4.sd.03[i] = NA
      # Coverage
      cover.isw4.03[i,1] = NA
      cover.isw4.03[i,2] = NA
      cover.isw4.03[i,3] = NA
      cover.isw4.03[i,4] = NA
      cover.isw4.03[i,5] = NA
      cover.isw4.03[i,6] = NA
      cover.isw4.03[i,7] = NA
      cover.isw4.03[i,8] = NA
    })
}

# IS beta table
table0.isw1[3,1]<-mean(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,2]<-mean(b0.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,3]<-sd(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,4]<-mean(cover.isw1.03[,4],na.rm=TRUE)
table0.isw1[3,5]<-mean(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,6]<-mean(b1.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,7]<-sd(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,8]<-mean(cover.isw1.03[,8],na.rm=TRUE)

table0.isw2[3,1]<-mean(b0.isw2.03,na.rm=TRUE)
table0.isw2[3,2]<-mean(b0.isw2.sd.03,na.rm=TRUE)
table0.isw2[3,3]<-sd(b0.isw2.03,na.rm=TRUE)
table0.isw2[3,4]<-mean(cover.isw2.03[,4],na.rm=TRUE)
table0.isw2[3,5]<-mean(b1.isw2.03,na.rm=TRUE)
table0.isw2[3,6]<-mean(b1.isw2.sd.03,na.rm=TRUE)
table0.isw2[3,7]<-sd(b1.isw2.03,na.rm=TRUE)
table0.isw2[3,8]<-mean(cover.isw2.03[,8],na.rm=TRUE)

table0.isw3[3,1]<-mean(b0.isw3.03,na.rm=TRUE)
table0.isw3[3,2]<-mean(b0.isw3.sd.03,na.rm=TRUE)
table0.isw3[3,3]<-sd(b0.isw3.03,na.rm=TRUE)
table0.isw3[3,4]<-mean(cover.isw3.03[,4],na.rm=TRUE)
table0.isw3[3,5]<-mean(b1.isw3.03,na.rm=TRUE)
table0.isw3[3,6]<-mean(b1.isw3.sd.03,na.rm=TRUE)
table0.isw3[3,7]<-sd(b1.isw3.03,na.rm=TRUE)
table0.isw3[3,8]<-mean(cover.isw3.03[,8],na.rm=TRUE)

table0.isw4[3,1]<-mean(b0.isw4.03,na.rm=TRUE)
table0.isw4[3,2]<-mean(b0.isw4.sd.03,na.rm=TRUE)
table0.isw4[3,3]<-sd(b0.isw4.03,na.rm=TRUE)
table0.isw4[3,4]<-mean(cover.isw4.03[,4],na.rm=TRUE)
table0.isw4[3,5]<-mean(b1.isw4.03,na.rm=TRUE)
table0.isw4[3,6]<-mean(b1.isw4.sd.03,na.rm=TRUE)
table0.isw4[3,7]<-sd(b1.isw4.03,na.rm=TRUE)
table0.isw4[3,8]<-mean(cover.isw4.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.isw1.05 = c()
b0.isw1.sd.05 = c()
b1.isw1.05 = c()
b1.isw1.sd.05 = c()
cover.isw1.05=matrix(NA,2000,8)
colnames(cover.isw1.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.05 = c()
b0.isw2.sd.05 = c()
b1.isw2.05 = c()
b1.isw2.sd.05 = c()
cover.isw2.05=matrix(NA,2000,8)
colnames(cover.isw2.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.05 = c()
b0.isw3.sd.05 = c()
b1.isw3.05 = c()
b1.isw3.sd.05 = c()
cover.isw3.05=matrix(NA,2000,8)
colnames(cover.isw3.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.05 = c()
b0.isw4.sd.05 = c()
b1.isw4.05 = c()
b1.isw4.sd.05 = c()
cover.isw4.05=matrix(NA,2000,8)
colnames(cover.isw4.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.5)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.05[i] = ismb.fit[1,1]
    b0.isw1.sd.05[i] = ismb.fit[1,2]
    b1.isw1.05[i] = ismb.fit[2,1]
    b1.isw1.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.05[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.05[i,2] = ismb.fit[1,1]
    cover.isw1.05[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.05[i,4] = ind(1.609438, cover.isw1.05[i,1], cover.isw1.05[i,3])
    cover.isw1.05[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.05[i,6] = ismb.fit[2,1]
    cover.isw1.05[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.05[i,8] = ind(0.6931472, cover.isw1.05[i,5], cover.isw1.05[i,7])}
    , error=function(e){
      b0.isw1.05[i] = NA
      b0.isw1.sd.05[i] = NA
      b1.isw1.05[i] = NA
      b1.isw1.sd.05[i] = NA
      # Coverage
      cover.isw1.05[i,1] = NA
      cover.isw1.05[i,2] = NA
      cover.isw1.05[i,3] = NA
      cover.isw1.05[i,4] = NA
      cover.isw1.05[i,5] = NA
      cover.isw1.05[i,6] = NA
      cover.isw1.05[i,7] = NA
      cover.isw1.05[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw2.05[i] = ismb.fit[1,1]
    b0.isw2.sd.05[i] = ismb.fit[1,2]
    b1.isw2.05[i] = ismb.fit[2,1]
    b1.isw2.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.05[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.05[i,2] = ismb.fit[1,1]
    cover.isw2.05[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.05[i,4] = ind(1.609438, cover.isw2.05[i,1], cover.isw2.05[i,3])
    cover.isw2.05[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.05[i,6] = ismb.fit[2,1]
    cover.isw2.05[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.05[i,8] = ind(0.6931472, cover.isw2.05[i,5], cover.isw2.05[i,7])}
    , error=function(e){
      b0.isw2.05[i] = NA
      b0.isw2.sd.05[i] = NA
      b1.isw2.05[i] = NA
      b1.isw2.sd.05[i] = NA
      # Coverage
      cover.isw2.05[i,1] = NA
      cover.isw2.05[i,2] = NA
      cover.isw2.05[i,3] = NA
      cover.isw2.05[i,4] = NA
      cover.isw2.05[i,5] = NA
      cover.isw2.05[i,6] = NA
      cover.isw2.05[i,7] = NA
      cover.isw2.05[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw3.05[i] = ismb.fit[1,1]
    b0.isw3.sd.05[i] = ismb.fit[1,2]
    b1.isw3.05[i] = ismb.fit[2,1]
    b1.isw3.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.05[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.05[i,2] = ismb.fit[1,1]
    cover.isw3.05[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.05[i,4] = ind(1.609438, cover.isw3.05[i,1], cover.isw3.05[i,3])
    cover.isw3.05[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.05[i,6] = ismb.fit[2,1]
    cover.isw3.05[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.05[i,8] = ind(0.6931472, cover.isw3.05[i,5], cover.isw3.05[i,7])}
    , error=function(e){
      b0.isw3.05[i] = NA
      b0.isw3.sd.05[i] = NA
      b1.isw3.05[i] = NA
      b1.isw3.sd.05[i] = NA
      # Coverage
      cover.isw3.05[i,1] = NA
      cover.isw3.05[i,2] = NA
      cover.isw3.05[i,3] = NA
      cover.isw3.05[i,4] = NA
      cover.isw3.05[i,5] = NA
      cover.isw3.05[i,6] = NA
      cover.isw3.05[i,7] = NA
      cover.isw3.05[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw4.05[i] = ismb.fit[1,1]
    b0.isw4.sd.05[i] = ismb.fit[1,2]
    b1.isw4.05[i] = ismb.fit[2,1]
    b1.isw4.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.05[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.05[i,2] = ismb.fit[1,1]
    cover.isw4.05[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.05[i,4] = ind(1.609438, cover.isw4.05[i,1], cover.isw4.05[i,3])
    cover.isw4.05[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.05[i,6] = ismb.fit[2,1]
    cover.isw4.05[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.05[i,8] = ind(0.6931472, cover.isw4.05[i,5], cover.isw4.05[i,7])}
    , error=function(e){
      b0.isw4.05[i] = NA
      b0.isw4.sd.05[i] = NA
      b1.isw4.05[i] = NA
      b1.isw4.sd.05[i] = NA
      # Coverage
      cover.isw4.05[i,1] = NA
      cover.isw4.05[i,2] = NA
      cover.isw4.05[i,3] = NA
      cover.isw4.05[i,4] = NA
      cover.isw4.05[i,5] = NA
      cover.isw4.05[i,6] = NA
      cover.isw4.05[i,7] = NA
      cover.isw4.05[i,8] = NA
    })
}

# IS beta table
table0.isw1[4,1]<-mean(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,2]<-mean(b0.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,3]<-sd(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,4]<-mean(cover.isw1.05[,4],na.rm=TRUE)
table0.isw1[4,5]<-mean(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,6]<-mean(b1.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,7]<-sd(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,8]<-mean(cover.isw1.05[,8],na.rm=TRUE)

table0.isw2[4,1]<-mean(b0.isw2.05,na.rm=TRUE)
table0.isw2[4,2]<-mean(b0.isw2.sd.05,na.rm=TRUE)
table0.isw2[4,3]<-sd(b0.isw2.05,na.rm=TRUE)
table0.isw2[4,4]<-mean(cover.isw2.05[,4],na.rm=TRUE)
table0.isw2[4,5]<-mean(b1.isw2.05,na.rm=TRUE)
table0.isw2[4,6]<-mean(b1.isw2.sd.05,na.rm=TRUE)
table0.isw2[4,7]<-sd(b1.isw2.05,na.rm=TRUE)
table0.isw2[4,8]<-mean(cover.isw2.05[,8],na.rm=TRUE)

table0.isw3[4,1]<-mean(b0.isw3.05,na.rm=TRUE)
table0.isw3[4,2]<-mean(b0.isw3.sd.05,na.rm=TRUE)
table0.isw3[4,3]<-sd(b0.isw3.05,na.rm=TRUE)
table0.isw3[4,4]<-mean(cover.isw3.05[,4],na.rm=TRUE)
table0.isw3[4,5]<-mean(b1.isw3.05,na.rm=TRUE)
table0.isw3[4,6]<-mean(b1.isw3.sd.05,na.rm=TRUE)
table0.isw3[4,7]<-sd(b1.isw3.05,na.rm=TRUE)
table0.isw3[4,8]<-mean(cover.isw3.05[,8],na.rm=TRUE)

table0.isw4[4,1]<-mean(b0.isw4.05,na.rm=TRUE)
table0.isw4[4,2]<-mean(b0.isw4.sd.05,na.rm=TRUE)
table0.isw4[4,3]<-sd(b0.isw4.05,na.rm=TRUE)
table0.isw4[4,4]<-mean(cover.isw4.05[,4],na.rm=TRUE)
table0.isw4[4,5]<-mean(b1.isw4.05,na.rm=TRUE)
table0.isw4[4,6]<-mean(b1.isw4.sd.05,na.rm=TRUE)
table0.isw4[4,7]<-sd(b1.isw4.05,na.rm=TRUE)
table0.isw4[4,8]<-mean(cover.isw4.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.isw1.07 = c()
b0.isw1.sd.07 = c()
b1.isw1.07 = c()
b1.isw1.sd.07 = c()
cover.isw1.07=matrix(NA,2000,8)
colnames(cover.isw1.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.07 = c()
b0.isw2.sd.07 = c()
b1.isw2.07 = c()
b1.isw2.sd.07 = c()
cover.isw2.07=matrix(NA,2000,8)
colnames(cover.isw2.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.07 = c()
b0.isw3.sd.07 = c()
b1.isw3.07 = c()
b1.isw3.sd.07 = c()
cover.isw3.07=matrix(NA,2000,8)
colnames(cover.isw3.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.07 = c()
b0.isw4.sd.07 = c()
b1.isw4.07 = c()
b1.isw4.sd.07 = c()
cover.isw4.07=matrix(NA,2000,8)
colnames(cover.isw4.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.7)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.07[i] = ismb.fit[1,1]
    b0.isw1.sd.07[i] = ismb.fit[1,2]
    b1.isw1.07[i] = ismb.fit[2,1]
    b1.isw1.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.07[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.07[i,2] = ismb.fit[1,1]
    cover.isw1.07[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.07[i,4] = ind(1.609438, cover.isw1.07[i,1], cover.isw1.07[i,3])
    cover.isw1.07[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.07[i,6] = ismb.fit[2,1]
    cover.isw1.07[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.07[i,8] = ind(0.6931472, cover.isw1.07[i,5], cover.isw1.07[i,7])}
    , error=function(e){
      b0.isw1.07[i] = NA
      b0.isw1.sd.07[i] = NA
      b1.isw1.07[i] = NA
      b1.isw1.sd.07[i] = NA
      # Coverage
      cover.isw1.07[i,1] = NA
      cover.isw1.07[i,2] = NA
      cover.isw1.07[i,3] = NA
      cover.isw1.07[i,4] = NA
      cover.isw1.07[i,5] = NA
      cover.isw1.07[i,6] = NA
      cover.isw1.07[i,7] = NA
      cover.isw1.07[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw2.07[i] = ismb.fit[1,1]
    b0.isw2.sd.07[i] = ismb.fit[1,2]
    b1.isw2.07[i] = ismb.fit[2,1]
    b1.isw2.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.07[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.07[i,2] = ismb.fit[1,1]
    cover.isw2.07[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.07[i,4] = ind(1.609438, cover.isw2.07[i,1], cover.isw2.07[i,3])
    cover.isw2.07[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.07[i,6] = ismb.fit[2,1]
    cover.isw2.07[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.07[i,8] = ind(0.6931472, cover.isw2.07[i,5], cover.isw2.07[i,7])}
    , error=function(e){
      b0.isw2.07[i] = NA
      b0.isw2.sd.07[i] = NA
      b1.isw2.07[i] = NA
      b1.isw2.sd.07[i] = NA
      # Coverage
      cover.isw2.07[i,1] = NA
      cover.isw2.07[i,2] = NA
      cover.isw2.07[i,3] = NA
      cover.isw2.07[i,4] = NA
      cover.isw2.07[i,5] = NA
      cover.isw2.07[i,6] = NA
      cover.isw2.07[i,7] = NA
      cover.isw2.07[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw3.07[i] = ismb.fit[1,1]
    b0.isw3.sd.07[i] = ismb.fit[1,2]
    b1.isw3.07[i] = ismb.fit[2,1]
    b1.isw3.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.07[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.07[i,2] = ismb.fit[1,1]
    cover.isw3.07[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.07[i,4] = ind(1.609438, cover.isw3.07[i,1], cover.isw3.07[i,3])
    cover.isw3.07[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.07[i,6] = ismb.fit[2,1]
    cover.isw3.07[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.07[i,8] = ind(0.6931472, cover.isw3.07[i,5], cover.isw3.07[i,7])}
    , error=function(e){
      b0.isw3.07[i] = NA
      b0.isw3.sd.07[i] = NA
      b1.isw3.07[i] = NA
      b1.isw3.sd.07[i] = NA
      # Coverage
      cover.isw3.07[i,1] = NA
      cover.isw3.07[i,2] = NA
      cover.isw3.07[i,3] = NA
      cover.isw3.07[i,4] = NA
      cover.isw3.07[i,5] = NA
      cover.isw3.07[i,6] = NA
      cover.isw3.07[i,7] = NA
      cover.isw3.07[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw4.07[i] = ismb.fit[1,1]
    b0.isw4.sd.07[i] = ismb.fit[1,2]
    b1.isw4.07[i] = ismb.fit[2,1]
    b1.isw4.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.07[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.07[i,2] = ismb.fit[1,1]
    cover.isw4.07[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.07[i,4] = ind(1.609438, cover.isw4.07[i,1], cover.isw4.07[i,3])
    cover.isw4.07[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.07[i,6] = ismb.fit[2,1]
    cover.isw4.07[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.07[i,8] = ind(0.6931472, cover.isw4.07[i,5], cover.isw4.07[i,7])}
    , error=function(e){
      b0.isw4.07[i] = NA
      b0.isw4.sd.07[i] = NA
      b1.isw4.07[i] = NA
      b1.isw4.sd.07[i] = NA
      # Coverage
      cover.isw4.07[i,1] = NA
      cover.isw4.07[i,2] = NA
      cover.isw4.07[i,3] = NA
      cover.isw4.07[i,4] = NA
      cover.isw4.07[i,5] = NA
      cover.isw4.07[i,6] = NA
      cover.isw4.07[i,7] = NA
      cover.isw4.07[i,8] = NA
    })
}

# IS beta table
table0.isw1[5,1]<-mean(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,2]<-mean(b0.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,3]<-sd(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,4]<-mean(cover.isw1.07[,4],na.rm=TRUE)
table0.isw1[5,5]<-mean(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,6]<-mean(b1.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,7]<-sd(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,8]<-mean(cover.isw1.07[,8],na.rm=TRUE)

table0.isw2[5,1]<-mean(b0.isw2.07,na.rm=TRUE)
table0.isw2[5,2]<-mean(b0.isw2.sd.07,na.rm=TRUE)
table0.isw2[5,3]<-sd(b0.isw2.07,na.rm=TRUE)
table0.isw2[5,4]<-mean(cover.isw2.07[,4],na.rm=TRUE)
table0.isw2[5,5]<-mean(b1.isw2.07,na.rm=TRUE)
table0.isw2[5,6]<-mean(b1.isw2.sd.07,na.rm=TRUE)
table0.isw2[5,7]<-sd(b1.isw2.07,na.rm=TRUE)
table0.isw2[5,8]<-mean(cover.isw2.07[,8],na.rm=TRUE)

table0.isw3[5,1]<-mean(b0.isw3.07,na.rm=TRUE)
table0.isw3[5,2]<-mean(b0.isw3.sd.07,na.rm=TRUE)
table0.isw3[5,3]<-sd(b0.isw3.07,na.rm=TRUE)
table0.isw3[5,4]<-mean(cover.isw3.07[,4],na.rm=TRUE)
table0.isw3[5,5]<-mean(b1.isw3.07,na.rm=TRUE)
table0.isw3[5,6]<-mean(b1.isw3.sd.07,na.rm=TRUE)
table0.isw3[5,7]<-sd(b1.isw3.07,na.rm=TRUE)
table0.isw3[5,8]<-mean(cover.isw3.07[,8],na.rm=TRUE)

table0.isw4[5,1]<-mean(b0.isw4.07,na.rm=TRUE)
table0.isw4[5,2]<-mean(b0.isw4.sd.07,na.rm=TRUE)
table0.isw4[5,3]<-sd(b0.isw4.07,na.rm=TRUE)
table0.isw4[5,4]<-mean(cover.isw4.07[,4],na.rm=TRUE)
table0.isw4[5,5]<-mean(b1.isw4.07,na.rm=TRUE)
table0.isw4[5,6]<-mean(b1.isw4.sd.07,na.rm=TRUE)
table0.isw4[5,7]<-sd(b1.isw4.07,na.rm=TRUE)
table0.isw4[5,8]<-mean(cover.isw4.07[,8],na.rm=TRUE)

#### censoring point at t_0=1 ####
c.0=5000000
c.1=70.39
c.3=24.35
c.5=14.07
c.7=8.49
#### t_0=1 & c=0% ####
b0.isw1.10 = c()
b0.isw1.sd.10 = c()
b1.isw1.10 = c()
b1.isw1.sd.10 = c()
cover.isw1.10=matrix(NA,2000,8)
colnames(cover.isw1.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.10 = c()
b0.isw2.sd.10 = c()
b1.isw2.10 = c()
b1.isw2.sd.10 = c()
cover.isw2.10=matrix(NA,2000,8)
colnames(cover.isw2.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.10 = c()
b0.isw3.sd.10 = c()
b1.isw3.10 = c()
b1.isw3.sd.10 = c()
cover.isw3.10=matrix(NA,2000,8)
colnames(cover.isw3.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.10 = c()
b0.isw4.sd.10 = c()
b1.isw4.10 = c()
b1.isw4.sd.10 = c()
cover.isw4.10=matrix(NA,2000,8)
colnames(cover.isw4.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.10[i] = ismb.fit[1,1]
    b0.isw1.sd.10[i] = ismb.fit[1,2]
    b1.isw1.10[i] = ismb.fit[2,1]
    b1.isw1.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.10[i,2] = ismb.fit[1,1]
    cover.isw1.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.10[i,4] = ind(1.410748, cover.isw1.10[i,1], cover.isw1.10[i,3])
    cover.isw1.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.10[i,6] = ismb.fit[2,1]
    cover.isw1.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.10[i,8] = ind(0.7974189, cover.isw1.10[i,5], cover.isw1.10[i,7])}
    , error=function(e){
      b0.isw1.10[i] = NA
      b0.isw1.sd.10[i] = NA
      b1.isw1.10[i] = NA
      b1.isw1.sd.10[i] = NA
      # Coverage
      cover.isw1.10[i,1] = NA
      cover.isw1.10[i,2] = NA
      cover.isw1.10[i,3] = NA
      cover.isw1.10[i,4] = NA
      cover.isw1.10[i,5] = NA
      cover.isw1.10[i,6] = NA
      cover.isw1.10[i,7] = NA
      cover.isw1.10[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw2.10[i] = ismb.fit[1,1]
    b0.isw2.sd.10[i] = ismb.fit[1,2]
    b1.isw2.10[i] = ismb.fit[2,1]
    b1.isw2.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.10[i,2] = ismb.fit[1,1]
    cover.isw2.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.10[i,4] = ind(1.410748, cover.isw2.10[i,1], cover.isw2.10[i,3])
    cover.isw2.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.10[i,6] = ismb.fit[2,1]
    cover.isw2.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.10[i,8] = ind(0.7974189, cover.isw2.10[i,5], cover.isw2.10[i,7])}
    , error=function(e){
      b0.isw2.10[i] = NA
      b0.isw2.sd.10[i] = NA
      b1.isw2.10[i] = NA
      b1.isw2.sd.10[i] = NA
      # Coverage
      cover.isw2.10[i,1] = NA
      cover.isw2.10[i,2] = NA
      cover.isw2.10[i,3] = NA
      cover.isw2.10[i,4] = NA
      cover.isw2.10[i,5] = NA
      cover.isw2.10[i,6] = NA
      cover.isw2.10[i,7] = NA
      cover.isw2.10[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw3.10[i] = ismb.fit[1,1]
    b0.isw3.sd.10[i] = ismb.fit[1,2]
    b1.isw3.10[i] = ismb.fit[2,1]
    b1.isw3.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.10[i,2] = ismb.fit[1,1]
    cover.isw3.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.10[i,4] = ind(1.410748, cover.isw3.10[i,1], cover.isw3.10[i,3])
    cover.isw3.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.10[i,6] = ismb.fit[2,1]
    cover.isw3.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.10[i,8] = ind(0.7974189, cover.isw3.10[i,5], cover.isw3.10[i,7])}
    , error=function(e){
      b0.isw3.10[i] = NA
      b0.isw3.sd.10[i] = NA
      b1.isw3.10[i] = NA
      b1.isw3.sd.10[i] = NA
      # Coverage
      cover.isw3.10[i,1] = NA
      cover.isw3.10[i,2] = NA
      cover.isw3.10[i,3] = NA
      cover.isw3.10[i,4] = NA
      cover.isw3.10[i,5] = NA
      cover.isw3.10[i,6] = NA
      cover.isw3.10[i,7] = NA
      cover.isw3.10[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw4.10[i] = ismb.fit[1,1]
    b0.isw4.sd.10[i] = ismb.fit[1,2]
    b1.isw4.10[i] = ismb.fit[2,1]
    b1.isw4.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.10[i,2] = ismb.fit[1,1]
    cover.isw4.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.10[i,4] = ind(1.410748, cover.isw4.10[i,1], cover.isw4.10[i,3])
    cover.isw4.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.10[i,6] = ismb.fit[2,1]
    cover.isw4.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.10[i,8] = ind(0.7974189, cover.isw4.10[i,5], cover.isw4.10[i,7])}
    , error=function(e){
      b0.isw4.10[i] = NA
      b0.isw4.sd.10[i] = NA
      b1.isw4.10[i] = NA
      b1.isw4.sd.10[i] = NA
      # Coverage
      cover.isw4.10[i,1] = NA
      cover.isw4.10[i,2] = NA
      cover.isw4.10[i,3] = NA
      cover.isw4.10[i,4] = NA
      cover.isw4.10[i,5] = NA
      cover.isw4.10[i,6] = NA
      cover.isw4.10[i,7] = NA
      cover.isw4.10[i,8] = NA
    })
}

# IS beta table
table1.isw1[1,1]<-mean(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,2]<-mean(b0.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,3]<-sd(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,4]<-mean(cover.isw1.10[,4],na.rm=TRUE)
table1.isw1[1,5]<-mean(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,6]<-mean(b1.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,7]<-sd(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,8]<-mean(cover.isw1.10[,8],na.rm=TRUE)

table1.isw2[1,1]<-mean(b0.isw2.10,na.rm=TRUE)
table1.isw2[1,2]<-mean(b0.isw2.sd.10,na.rm=TRUE)
table1.isw2[1,3]<-sd(b0.isw2.10,na.rm=TRUE)
table1.isw2[1,4]<-mean(cover.isw2.10[,4],na.rm=TRUE)
table1.isw2[1,5]<-mean(b1.isw2.10,na.rm=TRUE)
table1.isw2[1,6]<-mean(b1.isw2.sd.10,na.rm=TRUE)
table1.isw2[1,7]<-sd(b1.isw2.10,na.rm=TRUE)
table1.isw2[1,8]<-mean(cover.isw2.10[,8],na.rm=TRUE)

table1.isw3[1,1]<-mean(b0.isw3.10,na.rm=TRUE)
table1.isw3[1,2]<-mean(b0.isw3.sd.10,na.rm=TRUE)
table1.isw3[1,3]<-sd(b0.isw3.10,na.rm=TRUE)
table1.isw3[1,4]<-mean(cover.isw3.10[,4],na.rm=TRUE)
table1.isw3[1,5]<-mean(b1.isw3.10,na.rm=TRUE)
table1.isw3[1,6]<-mean(b1.isw3.sd.10,na.rm=TRUE)
table1.isw3[1,7]<-sd(b1.isw3.10,na.rm=TRUE)
table1.isw3[1,8]<-mean(cover.isw3.10[,8],na.rm=TRUE)

table1.isw4[1,1]<-mean(b0.isw4.10,na.rm=TRUE)
table1.isw4[1,2]<-mean(b0.isw4.sd.10,na.rm=TRUE)
table1.isw4[1,3]<-sd(b0.isw4.10,na.rm=TRUE)
table1.isw4[1,4]<-mean(cover.isw4.10[,4],na.rm=TRUE)
table1.isw4[1,5]<-mean(b1.isw4.10,na.rm=TRUE)
table1.isw4[1,6]<-mean(b1.isw4.sd.10,na.rm=TRUE)
table1.isw4[1,7]<-sd(b1.isw4.10,na.rm=TRUE)
table1.isw4[1,8]<-mean(cover.isw4.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.isw1.11 = c()
b0.isw1.sd.11 = c()
b1.isw1.11 = c()
b1.isw1.sd.11 = c()
cover.isw1.11=matrix(NA,2000,8)
colnames(cover.isw1.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.11 = c()
b0.isw2.sd.11 = c()
b1.isw2.11 = c()
b1.isw2.sd.11 = c()
cover.isw2.11=matrix(NA,2000,8)
colnames(cover.isw2.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.11 = c()
b0.isw3.sd.11 = c()
b1.isw3.11 = c()
b1.isw3.sd.11 = c()
cover.isw3.11=matrix(NA,2000,8)
colnames(cover.isw3.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.11 = c()
b0.isw4.sd.11 = c()
b1.isw4.11 = c()
b1.isw4.sd.11 = c()
cover.isw4.11=matrix(NA,2000,8)
colnames(cover.isw4.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.1)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.11[i] = ismb.fit[1,1]
    b0.isw1.sd.11[i] = ismb.fit[1,2]
    b1.isw1.11[i] = ismb.fit[2,1]
    b1.isw1.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.11[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.11[i,2] = ismb.fit[1,1]
    cover.isw1.11[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.11[i,4] = ind(1.410748, cover.isw1.11[i,1], cover.isw1.11[i,3])
    cover.isw1.11[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.11[i,6] = ismb.fit[2,1]
    cover.isw1.11[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.11[i,8] = ind(0.7974189, cover.isw1.11[i,5], cover.isw1.11[i,7])}
    , error=function(e){
      b0.isw1.11[i] = NA
      b0.isw1.sd.11[i] = NA
      b1.isw1.11[i] = NA
      b1.isw1.sd.11[i] = NA
      # Coverage
      cover.isw1.11[i,1] = NA
      cover.isw1.11[i,2] = NA
      cover.isw1.11[i,3] = NA
      cover.isw1.11[i,4] = NA
      cover.isw1.11[i,5] = NA
      cover.isw1.11[i,6] = NA
      cover.isw1.11[i,7] = NA
      cover.isw1.11[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw2.11[i] = ismb.fit[1,1]
    b0.isw2.sd.11[i] = ismb.fit[1,2]
    b1.isw2.11[i] = ismb.fit[2,1]
    b1.isw2.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.11[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.11[i,2] = ismb.fit[1,1]
    cover.isw2.11[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.11[i,4] = ind(1.410748, cover.isw2.11[i,1], cover.isw2.11[i,3])
    cover.isw2.11[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.11[i,6] = ismb.fit[2,1]
    cover.isw2.11[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.11[i,8] = ind(0.7974189, cover.isw2.11[i,5], cover.isw2.11[i,7])}
    , error=function(e){
      b0.isw2.11[i] = NA
      b0.isw2.sd.11[i] = NA
      b1.isw2.11[i] = NA
      b1.isw2.sd.11[i] = NA
      # Coverage
      cover.isw2.11[i,1] = NA
      cover.isw2.11[i,2] = NA
      cover.isw2.11[i,3] = NA
      cover.isw2.11[i,4] = NA
      cover.isw2.11[i,5] = NA
      cover.isw2.11[i,6] = NA
      cover.isw2.11[i,7] = NA
      cover.isw2.11[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw3.11[i] = ismb.fit[1,1]
    b0.isw3.sd.11[i] = ismb.fit[1,2]
    b1.isw3.11[i] = ismb.fit[2,1]
    b1.isw3.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.11[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.11[i,2] = ismb.fit[1,1]
    cover.isw3.11[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.11[i,4] = ind(1.410748, cover.isw3.11[i,1], cover.isw3.11[i,3])
    cover.isw3.11[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.11[i,6] = ismb.fit[2,1]
    cover.isw3.11[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.11[i,8] = ind(0.7974189, cover.isw3.11[i,5], cover.isw3.11[i,7])}
    , error=function(e){
      b0.isw3.11[i] = NA
      b0.isw3.sd.11[i] = NA
      b1.isw3.11[i] = NA
      b1.isw3.sd.11[i] = NA
      # Coverage
      cover.isw3.11[i,1] = NA
      cover.isw3.11[i,2] = NA
      cover.isw3.11[i,3] = NA
      cover.isw3.11[i,4] = NA
      cover.isw3.11[i,5] = NA
      cover.isw3.11[i,6] = NA
      cover.isw3.11[i,7] = NA
      cover.isw3.11[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw4.11[i] = ismb.fit[1,1]
    b0.isw4.sd.11[i] = ismb.fit[1,2]
    b1.isw4.11[i] = ismb.fit[2,1]
    b1.isw4.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.11[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.11[i,2] = ismb.fit[1,1]
    cover.isw4.11[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.11[i,4] = ind(1.410748, cover.isw4.11[i,1], cover.isw4.11[i,3])
    cover.isw4.11[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.11[i,6] = ismb.fit[2,1]
    cover.isw4.11[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.11[i,8] = ind(0.7974189, cover.isw4.11[i,5], cover.isw4.11[i,7])}
    , error=function(e){
      b0.isw4.11[i] = NA
      b0.isw4.sd.11[i] = NA
      b1.isw4.11[i] = NA
      b1.isw4.sd.11[i] = NA
      # Coverage
      cover.isw4.11[i,1] = NA
      cover.isw4.11[i,2] = NA
      cover.isw4.11[i,3] = NA
      cover.isw4.11[i,4] = NA
      cover.isw4.11[i,5] = NA
      cover.isw4.11[i,6] = NA
      cover.isw4.11[i,7] = NA
      cover.isw4.11[i,8] = NA
    })
}

# IS beta table
table1.isw1[2,1]<-mean(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,2]<-mean(b0.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,3]<-sd(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,4]<-mean(cover.isw1.11[,4],na.rm=TRUE)
table1.isw1[2,5]<-mean(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,6]<-mean(b1.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,7]<-sd(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,8]<-mean(cover.isw1.11[,8],na.rm=TRUE)

table1.isw2[2,1]<-mean(b0.isw2.11,na.rm=TRUE)
table1.isw2[2,2]<-mean(b0.isw2.sd.11,na.rm=TRUE)
table1.isw2[2,3]<-sd(b0.isw2.11,na.rm=TRUE)
table1.isw2[2,4]<-mean(cover.isw2.11[,4],na.rm=TRUE)
table1.isw2[2,5]<-mean(b1.isw2.11,na.rm=TRUE)
table1.isw2[2,6]<-mean(b1.isw2.sd.11,na.rm=TRUE)
table1.isw2[2,7]<-sd(b1.isw2.11,na.rm=TRUE)
table1.isw2[2,8]<-mean(cover.isw2.11[,8],na.rm=TRUE)

table1.isw3[2,1]<-mean(b0.isw3.11,na.rm=TRUE)
table1.isw3[2,2]<-mean(b0.isw3.sd.11,na.rm=TRUE)
table1.isw3[2,3]<-sd(b0.isw3.11,na.rm=TRUE)
table1.isw3[2,4]<-mean(cover.isw3.11[,4],na.rm=TRUE)
table1.isw3[2,5]<-mean(b1.isw3.11,na.rm=TRUE)
table1.isw3[2,6]<-mean(b1.isw3.sd.11,na.rm=TRUE)
table1.isw3[2,7]<-sd(b1.isw3.11,na.rm=TRUE)
table1.isw3[2,8]<-mean(cover.isw3.11[,8],na.rm=TRUE)

table1.isw4[2,1]<-mean(b0.isw4.11,na.rm=TRUE)
table1.isw4[2,2]<-mean(b0.isw4.sd.11,na.rm=TRUE)
table1.isw4[2,3]<-sd(b0.isw4.11,na.rm=TRUE)
table1.isw4[2,4]<-mean(cover.isw4.11[,4],na.rm=TRUE)
table1.isw4[2,5]<-mean(b1.isw4.11,na.rm=TRUE)
table1.isw4[2,6]<-mean(b1.isw4.sd.11,na.rm=TRUE)
table1.isw4[2,7]<-sd(b1.isw4.11,na.rm=TRUE)
table1.isw4[2,8]<-mean(cover.isw4.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.isw1.13 = c()
b0.isw1.sd.13 = c()
b1.isw1.13 = c()
b1.isw1.sd.13 = c()
cover.isw1.13=matrix(NA,2000,8)
colnames(cover.isw1.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.13 = c()
b0.isw2.sd.13 = c()
b1.isw2.13 = c()
b1.isw2.sd.13 = c()
cover.isw2.13=matrix(NA,2000,8)
colnames(cover.isw2.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.13 = c()
b0.isw3.sd.13 = c()
b1.isw3.13 = c()
b1.isw3.sd.13 = c()
cover.isw3.13=matrix(NA,2000,8)
colnames(cover.isw3.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.13 = c()
b0.isw4.sd.13 = c()
b1.isw4.13 = c()
b1.isw4.sd.13 = c()
cover.isw4.13=matrix(NA,2000,8)
colnames(cover.isw4.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.3)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.13[i] = ismb.fit[1,1]
    b0.isw1.sd.13[i] = ismb.fit[1,2]
    b1.isw1.13[i] = ismb.fit[2,1]
    b1.isw1.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.13[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.13[i,2] = ismb.fit[1,1]
    cover.isw1.13[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.13[i,4] = ind(1.410748, cover.isw1.13[i,1], cover.isw1.13[i,3])
    cover.isw1.13[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.13[i,6] = ismb.fit[2,1]
    cover.isw1.13[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.13[i,8] = ind(0.7974189, cover.isw1.13[i,5], cover.isw1.13[i,7])}
    , error=function(e){
      b0.isw1.13[i] = NA
      b0.isw1.sd.13[i] = NA
      b1.isw1.13[i] = NA
      b1.isw1.sd.13[i] = NA
      # Coverage
      cover.isw1.13[i,1] = NA
      cover.isw1.13[i,2] = NA
      cover.isw1.13[i,3] = NA
      cover.isw1.13[i,4] = NA
      cover.isw1.13[i,5] = NA
      cover.isw1.13[i,6] = NA
      cover.isw1.13[i,7] = NA
      cover.isw1.13[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw2.13[i] = ismb.fit[1,1]
    b0.isw2.sd.13[i] = ismb.fit[1,2]
    b1.isw2.13[i] = ismb.fit[2,1]
    b1.isw2.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.13[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.13[i,2] = ismb.fit[1,1]
    cover.isw2.13[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.13[i,4] = ind(1.410748, cover.isw2.13[i,1], cover.isw2.13[i,3])
    cover.isw2.13[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.13[i,6] = ismb.fit[2,1]
    cover.isw2.13[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.13[i,8] = ind(0.7974189, cover.isw2.13[i,5], cover.isw2.13[i,7])}
    , error=function(e){
      b0.isw2.13[i] = NA
      b0.isw2.sd.13[i] = NA
      b1.isw2.13[i] = NA
      b1.isw2.sd.13[i] = NA
      # Coverage
      cover.isw2.13[i,1] = NA
      cover.isw2.13[i,2] = NA
      cover.isw2.13[i,3] = NA
      cover.isw2.13[i,4] = NA
      cover.isw2.13[i,5] = NA
      cover.isw2.13[i,6] = NA
      cover.isw2.13[i,7] = NA
      cover.isw2.13[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw3.13[i] = ismb.fit[1,1]
    b0.isw3.sd.13[i] = ismb.fit[1,2]
    b1.isw3.13[i] = ismb.fit[2,1]
    b1.isw3.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.13[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.13[i,2] = ismb.fit[1,1]
    cover.isw3.13[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.13[i,4] = ind(1.410748, cover.isw3.13[i,1], cover.isw3.13[i,3])
    cover.isw3.13[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.13[i,6] = ismb.fit[2,1]
    cover.isw3.13[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.13[i,8] = ind(0.7974189, cover.isw3.13[i,5], cover.isw3.13[i,7])}
    , error=function(e){
      b0.isw3.13[i] = NA
      b0.isw3.sd.13[i] = NA
      b1.isw3.13[i] = NA
      b1.isw3.sd.13[i] = NA
      # Coverage
      cover.isw3.13[i,1] = NA
      cover.isw3.13[i,2] = NA
      cover.isw3.13[i,3] = NA
      cover.isw3.13[i,4] = NA
      cover.isw3.13[i,5] = NA
      cover.isw3.13[i,6] = NA
      cover.isw3.13[i,7] = NA
      cover.isw3.13[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw4.13[i] = ismb.fit[1,1]
    b0.isw4.sd.13[i] = ismb.fit[1,2]
    b1.isw4.13[i] = ismb.fit[2,1]
    b1.isw4.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.13[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.13[i,2] = ismb.fit[1,1]
    cover.isw4.13[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.13[i,4] = ind(1.410748, cover.isw4.13[i,1], cover.isw4.13[i,3])
    cover.isw4.13[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.13[i,6] = ismb.fit[2,1]
    cover.isw4.13[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.13[i,8] = ind(0.7974189, cover.isw4.13[i,5], cover.isw4.13[i,7])}
    , error=function(e){
      b0.isw4.13[i] = NA
      b0.isw4.sd.13[i] = NA
      b1.isw4.13[i] = NA
      b1.isw4.sd.13[i] = NA
      # Coverage
      cover.isw4.13[i,1] = NA
      cover.isw4.13[i,2] = NA
      cover.isw4.13[i,3] = NA
      cover.isw4.13[i,4] = NA
      cover.isw4.13[i,5] = NA
      cover.isw4.13[i,6] = NA
      cover.isw4.13[i,7] = NA
      cover.isw4.13[i,8] = NA
    })
}

# IS beta table
table1.isw1[3,1]<-mean(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,2]<-mean(b0.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,3]<-sd(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,4]<-mean(cover.isw1.13[,4],na.rm=TRUE)
table1.isw1[3,5]<-mean(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,6]<-mean(b1.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,7]<-sd(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,8]<-mean(cover.isw1.13[,8],na.rm=TRUE)

table1.isw2[3,1]<-mean(b0.isw2.13,na.rm=TRUE)
table1.isw2[3,2]<-mean(b0.isw2.sd.13,na.rm=TRUE)
table1.isw2[3,3]<-sd(b0.isw2.13,na.rm=TRUE)
table1.isw2[3,4]<-mean(cover.isw2.13[,4],na.rm=TRUE)
table1.isw2[3,5]<-mean(b1.isw2.13,na.rm=TRUE)
table1.isw2[3,6]<-mean(b1.isw2.sd.13,na.rm=TRUE)
table1.isw2[3,7]<-sd(b1.isw2.13,na.rm=TRUE)
table1.isw2[3,8]<-mean(cover.isw2.13[,8],na.rm=TRUE)

table1.isw3[3,1]<-mean(b0.isw3.13,na.rm=TRUE)
table1.isw3[3,2]<-mean(b0.isw3.sd.13,na.rm=TRUE)
table1.isw3[3,3]<-sd(b0.isw3.13,na.rm=TRUE)
table1.isw3[3,4]<-mean(cover.isw3.13[,4],na.rm=TRUE)
table1.isw3[3,5]<-mean(b1.isw3.13,na.rm=TRUE)
table1.isw3[3,6]<-mean(b1.isw3.sd.13,na.rm=TRUE)
table1.isw3[3,7]<-sd(b1.isw3.13,na.rm=TRUE)
table1.isw3[3,8]<-mean(cover.isw3.13[,8],na.rm=TRUE)

table1.isw4[3,1]<-mean(b0.isw4.13,na.rm=TRUE)
table1.isw4[3,2]<-mean(b0.isw4.sd.13,na.rm=TRUE)
table1.isw4[3,3]<-sd(b0.isw4.13,na.rm=TRUE)
table1.isw4[3,4]<-mean(cover.isw4.13[,4],na.rm=TRUE)
table1.isw4[3,5]<-mean(b1.isw4.13,na.rm=TRUE)
table1.isw4[3,6]<-mean(b1.isw4.sd.13,na.rm=TRUE)
table1.isw4[3,7]<-sd(b1.isw4.13,na.rm=TRUE)
table1.isw4[3,8]<-mean(cover.isw4.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.isw1.15 = c()
b0.isw1.sd.15 = c()
b1.isw1.15 = c()
b1.isw1.sd.15 = c()
cover.isw1.15=matrix(NA,2000,8)
colnames(cover.isw1.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.15 = c()
b0.isw2.sd.15 = c()
b1.isw2.15 = c()
b1.isw2.sd.15 = c()
cover.isw2.15=matrix(NA,2000,8)
colnames(cover.isw2.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.15 = c()
b0.isw3.sd.15 = c()
b1.isw3.15 = c()
b1.isw3.sd.15 = c()
cover.isw3.15=matrix(NA,2000,8)
colnames(cover.isw3.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.15 = c()
b0.isw4.sd.15 = c()
b1.isw4.15 = c()
b1.isw4.sd.15 = c()
cover.isw4.15=matrix(NA,2000,8)
colnames(cover.isw4.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.5)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.15[i] = ismb.fit[1,1]
    b0.isw1.sd.15[i] = ismb.fit[1,2]
    b1.isw1.15[i] = ismb.fit[2,1]
    b1.isw1.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.15[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.15[i,2] = ismb.fit[1,1]
    cover.isw1.15[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.15[i,4] = ind(1.410748, cover.isw1.15[i,1], cover.isw1.15[i,3])
    cover.isw1.15[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.15[i,6] = ismb.fit[2,1]
    cover.isw1.15[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.15[i,8] = ind(0.7974189, cover.isw1.15[i,5], cover.isw1.15[i,7])}
    , error=function(e){
      b0.isw1.15[i] = NA
      b0.isw1.sd.15[i] = NA
      b1.isw1.15[i] = NA
      b1.isw1.sd.15[i] = NA
      # Coverage
      cover.isw1.15[i,1] = NA
      cover.isw1.15[i,2] = NA
      cover.isw1.15[i,3] = NA
      cover.isw1.15[i,4] = NA
      cover.isw1.15[i,5] = NA
      cover.isw1.15[i,6] = NA
      cover.isw1.15[i,7] = NA
      cover.isw1.15[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw2.15[i] = ismb.fit[1,1]
    b0.isw2.sd.15[i] = ismb.fit[1,2]
    b1.isw2.15[i] = ismb.fit[2,1]
    b1.isw2.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.15[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.15[i,2] = ismb.fit[1,1]
    cover.isw2.15[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.15[i,4] = ind(1.410748, cover.isw2.15[i,1], cover.isw2.15[i,3])
    cover.isw2.15[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.15[i,6] = ismb.fit[2,1]
    cover.isw2.15[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.15[i,8] = ind(0.7974189, cover.isw2.15[i,5], cover.isw2.15[i,7])}
    , error=function(e){
      b0.isw2.15[i] = NA
      b0.isw2.sd.15[i] = NA
      b1.isw2.15[i] = NA
      b1.isw2.sd.15[i] = NA
      # Coverage
      cover.isw2.15[i,1] = NA
      cover.isw2.15[i,2] = NA
      cover.isw2.15[i,3] = NA
      cover.isw2.15[i,4] = NA
      cover.isw2.15[i,5] = NA
      cover.isw2.15[i,6] = NA
      cover.isw2.15[i,7] = NA
      cover.isw2.15[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw3.15[i] = ismb.fit[1,1]
    b0.isw3.sd.15[i] = ismb.fit[1,2]
    b1.isw3.15[i] = ismb.fit[2,1]
    b1.isw3.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.15[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.15[i,2] = ismb.fit[1,1]
    cover.isw3.15[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.15[i,4] = ind(1.410748, cover.isw3.15[i,1], cover.isw3.15[i,3])
    cover.isw3.15[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.15[i,6] = ismb.fit[2,1]
    cover.isw3.15[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.15[i,8] = ind(0.7974189, cover.isw3.15[i,5], cover.isw3.15[i,7])}
    , error=function(e){
      b0.isw3.15[i] = NA
      b0.isw3.sd.15[i] = NA
      b1.isw3.15[i] = NA
      b1.isw3.sd.15[i] = NA
      # Coverage
      cover.isw3.15[i,1] = NA
      cover.isw3.15[i,2] = NA
      cover.isw3.15[i,3] = NA
      cover.isw3.15[i,4] = NA
      cover.isw3.15[i,5] = NA
      cover.isw3.15[i,6] = NA
      cover.isw3.15[i,7] = NA
      cover.isw3.15[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw4.15[i] = ismb.fit[1,1]
    b0.isw4.sd.15[i] = ismb.fit[1,2]
    b1.isw4.15[i] = ismb.fit[2,1]
    b1.isw4.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.15[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.15[i,2] = ismb.fit[1,1]
    cover.isw4.15[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.15[i,4] = ind(1.410748, cover.isw4.15[i,1], cover.isw4.15[i,3])
    cover.isw4.15[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.15[i,6] = ismb.fit[2,1]
    cover.isw4.15[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.15[i,8] = ind(0.7974189, cover.isw4.15[i,5], cover.isw4.15[i,7])}
    , error=function(e){
      b0.isw4.15[i] = NA
      b0.isw4.sd.15[i] = NA
      b1.isw4.15[i] = NA
      b1.isw4.sd.15[i] = NA
      # Coverage
      cover.isw4.15[i,1] = NA
      cover.isw4.15[i,2] = NA
      cover.isw4.15[i,3] = NA
      cover.isw4.15[i,4] = NA
      cover.isw4.15[i,5] = NA
      cover.isw4.15[i,6] = NA
      cover.isw4.15[i,7] = NA
      cover.isw4.15[i,8] = NA
    })
}

# IS beta table
table1.isw1[4,1]<-mean(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,2]<-mean(b0.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,3]<-sd(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,4]<-mean(cover.isw1.15[,4],na.rm=TRUE)
table1.isw1[4,5]<-mean(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,6]<-mean(b1.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,7]<-sd(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,8]<-mean(cover.isw1.15[,8],na.rm=TRUE)

table1.isw2[4,1]<-mean(b0.isw2.15,na.rm=TRUE)
table1.isw2[4,2]<-mean(b0.isw2.sd.15,na.rm=TRUE)
table1.isw2[4,3]<-sd(b0.isw2.15,na.rm=TRUE)
table1.isw2[4,4]<-mean(cover.isw2.15[,4],na.rm=TRUE)
table1.isw2[4,5]<-mean(b1.isw2.15,na.rm=TRUE)
table1.isw2[4,6]<-mean(b1.isw2.sd.15,na.rm=TRUE)
table1.isw2[4,7]<-sd(b1.isw2.15,na.rm=TRUE)
table1.isw2[4,8]<-mean(cover.isw2.15[,8],na.rm=TRUE)

table1.isw3[4,1]<-mean(b0.isw3.15,na.rm=TRUE)
table1.isw3[4,2]<-mean(b0.isw3.sd.15,na.rm=TRUE)
table1.isw3[4,3]<-sd(b0.isw3.15,na.rm=TRUE)
table1.isw3[4,4]<-mean(cover.isw3.15[,4],na.rm=TRUE)
table1.isw3[4,5]<-mean(b1.isw3.15,na.rm=TRUE)
table1.isw3[4,6]<-mean(b1.isw3.sd.15,na.rm=TRUE)
table1.isw3[4,7]<-sd(b1.isw3.15,na.rm=TRUE)
table1.isw3[4,8]<-mean(cover.isw3.15[,8],na.rm=TRUE)

table1.isw4[4,1]<-mean(b0.isw4.15,na.rm=TRUE)
table1.isw4[4,2]<-mean(b0.isw4.sd.15,na.rm=TRUE)
table1.isw4[4,3]<-sd(b0.isw4.15,na.rm=TRUE)
table1.isw4[4,4]<-mean(cover.isw4.15[,4],na.rm=TRUE)
table1.isw4[4,5]<-mean(b1.isw4.15,na.rm=TRUE)
table1.isw4[4,6]<-mean(b1.isw4.sd.15,na.rm=TRUE)
table1.isw4[4,7]<-sd(b1.isw4.15,na.rm=TRUE)
table1.isw4[4,8]<-mean(cover.isw4.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.isw1.17 = c()
b0.isw1.sd.17 = c()
b1.isw1.17 = c()
b1.isw1.sd.17 = c()
cover.isw1.17=matrix(NA,2000,8)
colnames(cover.isw1.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.17 = c()
b0.isw2.sd.17 = c()
b1.isw2.17 = c()
b1.isw2.sd.17 = c()
cover.isw2.17=matrix(NA,2000,8)
colnames(cover.isw2.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.17 = c()
b0.isw3.sd.17 = c()
b1.isw3.17 = c()
b1.isw3.sd.17 = c()
cover.isw3.17=matrix(NA,2000,8)
colnames(cover.isw3.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.17 = c()
b0.isw4.sd.17 = c()
b1.isw4.17 = c()
b1.isw4.sd.17 = c()
cover.isw4.17=matrix(NA,2000,8)
colnames(cover.isw4.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.7)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.17[i] = ismb.fit[1,1]
    b0.isw1.sd.17[i] = ismb.fit[1,2]
    b1.isw1.17[i] = ismb.fit[2,1]
    b1.isw1.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.17[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.17[i,2] = ismb.fit[1,1]
    cover.isw1.17[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.17[i,4] = ind(1.410748, cover.isw1.17[i,1], cover.isw1.17[i,3])
    cover.isw1.17[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.17[i,6] = ismb.fit[2,1]
    cover.isw1.17[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.17[i,8] = ind(0.7974189, cover.isw1.17[i,5], cover.isw1.17[i,7])}
    , error=function(e){
      b0.isw1.17[i] = NA
      b0.isw1.sd.17[i] = NA
      b1.isw1.17[i] = NA
      b1.isw1.sd.17[i] = NA
      # Coverage
      cover.isw1.17[i,1] = NA
      cover.isw1.17[i,2] = NA
      cover.isw1.17[i,3] = NA
      cover.isw1.17[i,4] = NA
      cover.isw1.17[i,5] = NA
      cover.isw1.17[i,6] = NA
      cover.isw1.17[i,7] = NA
      cover.isw1.17[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw2.17[i] = ismb.fit[1,1]
    b0.isw2.sd.17[i] = ismb.fit[1,2]
    b1.isw2.17[i] = ismb.fit[2,1]
    b1.isw2.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.17[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.17[i,2] = ismb.fit[1,1]
    cover.isw2.17[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.17[i,4] = ind(1.410748, cover.isw2.17[i,1], cover.isw2.17[i,3])
    cover.isw2.17[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.17[i,6] = ismb.fit[2,1]
    cover.isw2.17[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.17[i,8] = ind(0.7974189, cover.isw2.17[i,5], cover.isw2.17[i,7])}
    , error=function(e){
      b0.isw2.17[i] = NA
      b0.isw2.sd.17[i] = NA
      b1.isw2.17[i] = NA
      b1.isw2.sd.17[i] = NA
      # Coverage
      cover.isw2.17[i,1] = NA
      cover.isw2.17[i,2] = NA
      cover.isw2.17[i,3] = NA
      cover.isw2.17[i,4] = NA
      cover.isw2.17[i,5] = NA
      cover.isw2.17[i,6] = NA
      cover.isw2.17[i,7] = NA
      cover.isw2.17[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw3.17[i] = ismb.fit[1,1]
    b0.isw3.sd.17[i] = ismb.fit[1,2]
    b1.isw3.17[i] = ismb.fit[2,1]
    b1.isw3.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.17[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.17[i,2] = ismb.fit[1,1]
    cover.isw3.17[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.17[i,4] = ind(1.410748, cover.isw3.17[i,1], cover.isw3.17[i,3])
    cover.isw3.17[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.17[i,6] = ismb.fit[2,1]
    cover.isw3.17[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.17[i,8] = ind(0.7974189, cover.isw3.17[i,5], cover.isw3.17[i,7])}
    , error=function(e){
      b0.isw3.17[i] = NA
      b0.isw3.sd.17[i] = NA
      b1.isw3.17[i] = NA
      b1.isw3.sd.17[i] = NA
      # Coverage
      cover.isw3.17[i,1] = NA
      cover.isw3.17[i,2] = NA
      cover.isw3.17[i,3] = NA
      cover.isw3.17[i,4] = NA
      cover.isw3.17[i,5] = NA
      cover.isw3.17[i,6] = NA
      cover.isw3.17[i,7] = NA
      cover.isw3.17[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw4.17[i] = ismb.fit[1,1]
    b0.isw4.sd.17[i] = ismb.fit[1,2]
    b1.isw4.17[i] = ismb.fit[2,1]
    b1.isw4.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.17[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.17[i,2] = ismb.fit[1,1]
    cover.isw4.17[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.17[i,4] = ind(1.410748, cover.isw4.17[i,1], cover.isw4.17[i,3])
    cover.isw4.17[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.17[i,6] = ismb.fit[2,1]
    cover.isw4.17[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.17[i,8] = ind(0.7974189, cover.isw4.17[i,5], cover.isw4.17[i,7])}
    , error=function(e){
      b0.isw4.17[i] = NA
      b0.isw4.sd.17[i] = NA
      b1.isw4.17[i] = NA
      b1.isw4.sd.17[i] = NA
      # Coverage
      cover.isw4.17[i,1] = NA
      cover.isw4.17[i,2] = NA
      cover.isw4.17[i,3] = NA
      cover.isw4.17[i,4] = NA
      cover.isw4.17[i,5] = NA
      cover.isw4.17[i,6] = NA
      cover.isw4.17[i,7] = NA
      cover.isw4.17[i,8] = NA
    })
}

# IS beta table
table1.isw1[5,1]<-mean(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,2]<-mean(b0.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,3]<-sd(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,4]<-mean(cover.isw1.17[,4],na.rm=TRUE)
table1.isw1[5,5]<-mean(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,6]<-mean(b1.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,7]<-sd(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,8]<-mean(cover.isw1.17[,8],na.rm=TRUE)

table1.isw2[5,1]<-mean(b0.isw2.17,na.rm=TRUE)
table1.isw2[5,2]<-mean(b0.isw2.sd.17,na.rm=TRUE)
table1.isw2[5,3]<-sd(b0.isw2.17,na.rm=TRUE)
table1.isw2[5,4]<-mean(cover.isw2.17[,4],na.rm=TRUE)
table1.isw2[5,5]<-mean(b1.isw2.17,na.rm=TRUE)
table1.isw2[5,6]<-mean(b1.isw2.sd.17,na.rm=TRUE)
table1.isw2[5,7]<-sd(b1.isw2.17,na.rm=TRUE)
table1.isw2[5,8]<-mean(cover.isw2.17[,8],na.rm=TRUE)

table1.isw3[5,1]<-mean(b0.isw3.17,na.rm=TRUE)
table1.isw3[5,2]<-mean(b0.isw3.sd.17,na.rm=TRUE)
table1.isw3[5,3]<-sd(b0.isw3.17,na.rm=TRUE)
table1.isw3[5,4]<-mean(cover.isw3.17[,4],na.rm=TRUE)
table1.isw3[5,5]<-mean(b1.isw3.17,na.rm=TRUE)
table1.isw3[5,6]<-mean(b1.isw3.sd.17,na.rm=TRUE)
table1.isw3[5,7]<-sd(b1.isw3.17,na.rm=TRUE)
table1.isw3[5,8]<-mean(cover.isw3.17[,8],na.rm=TRUE)

table1.isw4[5,1]<-mean(b0.isw4.17,na.rm=TRUE)
table1.isw4[5,2]<-mean(b0.isw4.sd.17,na.rm=TRUE)
table1.isw4[5,3]<-sd(b0.isw4.17,na.rm=TRUE)
table1.isw4[5,4]<-mean(cover.isw4.17[,4],na.rm=TRUE)
table1.isw4[5,5]<-mean(b1.isw4.17,na.rm=TRUE)
table1.isw4[5,6]<-mean(b1.isw4.sd.17,na.rm=TRUE)
table1.isw4[5,7]<-sd(b1.isw4.17,na.rm=TRUE)
table1.isw4[5,8]<-mean(cover.isw4.17[,8],na.rm=TRUE)

#### censoring point at t_0=2 ####
c.0=5000000
c.1=64.86
c.3=23.34
c.5=13.62
c.7=8.36
#### t_0=2 & c=0% ####
b0.isw1.20 = c()
b0.isw1.sd.20 = c()
b1.isw1.20 = c()
b1.isw1.sd.20 = c()
cover.isw1.20=matrix(NA,2000,8)
colnames(cover.isw1.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.20 = c()
b0.isw2.sd.20 = c()
b1.isw2.20 = c()
b1.isw2.sd.20 = c()
cover.isw2.20=matrix(NA,2000,8)
colnames(cover.isw2.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.20 = c()
b0.isw3.sd.20 = c()
b1.isw3.20 = c()
b1.isw3.sd.20 = c()
cover.isw3.20=matrix(NA,2000,8)
colnames(cover.isw3.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.20 = c()
b0.isw4.sd.20 = c()
b1.isw4.20 = c()
b1.isw4.sd.20 = c()
cover.isw4.20=matrix(NA,2000,8)
colnames(cover.isw4.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.20[i] = ismb.fit[1,1]
    b0.isw1.sd.20[i] = ismb.fit[1,2]
    b1.isw1.20[i] = ismb.fit[2,1]
    b1.isw1.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.20[i,2] = ismb.fit[1,1]
    cover.isw1.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.20[i,4] = ind(1.219403, cover.isw1.20[i,1], cover.isw1.20[i,3])
    cover.isw1.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.20[i,6] = ismb.fit[2,1]
    cover.isw1.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.20[i,8] = ind(0.9070615, cover.isw1.20[i,5], cover.isw1.20[i,7])}
    , error=function(e){
      b0.isw1.20[i] = NA
      b0.isw1.sd.20[i] = NA
      b1.isw1.20[i] = NA
      b1.isw1.sd.20[i] = NA
      # Coverage
      cover.isw1.20[i,1] = NA
      cover.isw1.20[i,2] = NA
      cover.isw1.20[i,3] = NA
      cover.isw1.20[i,4] = NA
      cover.isw1.20[i,5] = NA
      cover.isw1.20[i,6] = NA
      cover.isw1.20[i,7] = NA
      cover.isw1.20[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw2.20[i] = ismb.fit[1,1]
    b0.isw2.sd.20[i] = ismb.fit[1,2]
    b1.isw2.20[i] = ismb.fit[2,1]
    b1.isw2.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.20[i,2] = ismb.fit[1,1]
    cover.isw2.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.20[i,4] = ind(1.219403, cover.isw2.20[i,1], cover.isw2.20[i,3])
    cover.isw2.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.20[i,6] = ismb.fit[2,1]
    cover.isw2.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.20[i,8] = ind(0.9070615, cover.isw2.20[i,5], cover.isw2.20[i,7])}
    , error=function(e){
      b0.isw2.20[i] = NA
      b0.isw2.sd.20[i] = NA
      b1.isw2.20[i] = NA
      b1.isw2.sd.20[i] = NA
      # Coverage
      cover.isw2.20[i,1] = NA
      cover.isw2.20[i,2] = NA
      cover.isw2.20[i,3] = NA
      cover.isw2.20[i,4] = NA
      cover.isw2.20[i,5] = NA
      cover.isw2.20[i,6] = NA
      cover.isw2.20[i,7] = NA
      cover.isw2.20[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw3.20[i] = ismb.fit[1,1]
    b0.isw3.sd.20[i] = ismb.fit[1,2]
    b1.isw3.20[i] = ismb.fit[2,1]
    b1.isw3.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.20[i,2] = ismb.fit[1,1]
    cover.isw3.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.20[i,4] = ind(1.219403, cover.isw3.20[i,1], cover.isw3.20[i,3])
    cover.isw3.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.20[i,6] = ismb.fit[2,1]
    cover.isw3.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.20[i,8] = ind(0.9070615, cover.isw3.20[i,5], cover.isw3.20[i,7])}
    , error=function(e){
      b0.isw3.20[i] = NA
      b0.isw3.sd.20[i] = NA
      b1.isw3.20[i] = NA
      b1.isw3.sd.20[i] = NA
      # Coverage
      cover.isw3.20[i,1] = NA
      cover.isw3.20[i,2] = NA
      cover.isw3.20[i,3] = NA
      cover.isw3.20[i,4] = NA
      cover.isw3.20[i,5] = NA
      cover.isw3.20[i,6] = NA
      cover.isw3.20[i,7] = NA
      cover.isw3.20[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw4.20[i] = ismb.fit[1,1]
    b0.isw4.sd.20[i] = ismb.fit[1,2]
    b1.isw4.20[i] = ismb.fit[2,1]
    b1.isw4.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.20[i,2] = ismb.fit[1,1]
    cover.isw4.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.20[i,4] = ind(1.219403, cover.isw4.20[i,1], cover.isw4.20[i,3])
    cover.isw4.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.20[i,6] = ismb.fit[2,1]
    cover.isw4.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.20[i,8] = ind(0.9070615, cover.isw4.20[i,5], cover.isw4.20[i,7])}
    , error=function(e){
      b0.isw4.20[i] = NA
      b0.isw4.sd.20[i] = NA
      b1.isw4.20[i] = NA
      b1.isw4.sd.20[i] = NA
      # Coverage
      cover.isw4.20[i,1] = NA
      cover.isw4.20[i,2] = NA
      cover.isw4.20[i,3] = NA
      cover.isw4.20[i,4] = NA
      cover.isw4.20[i,5] = NA
      cover.isw4.20[i,6] = NA
      cover.isw4.20[i,7] = NA
      cover.isw4.20[i,8] = NA
    })
}

# IS beta table
table2.isw1[1,1]<-mean(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,2]<-mean(b0.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,3]<-sd(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,4]<-mean(cover.isw1.20[,4],na.rm=TRUE)
table2.isw1[1,5]<-mean(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,6]<-mean(b1.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,7]<-sd(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,8]<-mean(cover.isw1.20[,8],na.rm=TRUE)

table2.isw2[1,1]<-mean(b0.isw2.20,na.rm=TRUE)
table2.isw2[1,2]<-mean(b0.isw2.sd.20,na.rm=TRUE)
table2.isw2[1,3]<-sd(b0.isw2.20,na.rm=TRUE)
table2.isw2[1,4]<-mean(cover.isw2.20[,4],na.rm=TRUE)
table2.isw2[1,5]<-mean(b1.isw2.20,na.rm=TRUE)
table2.isw2[1,6]<-mean(b1.isw2.sd.20,na.rm=TRUE)
table2.isw2[1,7]<-sd(b1.isw2.20,na.rm=TRUE)
table2.isw2[1,8]<-mean(cover.isw2.20[,8],na.rm=TRUE)

table2.isw3[1,1]<-mean(b0.isw3.20,na.rm=TRUE)
table2.isw3[1,2]<-mean(b0.isw3.sd.20,na.rm=TRUE)
table2.isw3[1,3]<-sd(b0.isw3.20,na.rm=TRUE)
table2.isw3[1,4]<-mean(cover.isw3.20[,4],na.rm=TRUE)
table2.isw3[1,5]<-mean(b1.isw3.20,na.rm=TRUE)
table2.isw3[1,6]<-mean(b1.isw3.sd.20,na.rm=TRUE)
table2.isw3[1,7]<-sd(b1.isw3.20,na.rm=TRUE)
table2.isw3[1,8]<-mean(cover.isw3.20[,8],na.rm=TRUE)

table2.isw4[1,1]<-mean(b0.isw4.20,na.rm=TRUE)
table2.isw4[1,2]<-mean(b0.isw4.sd.20,na.rm=TRUE)
table2.isw4[1,3]<-sd(b0.isw4.20,na.rm=TRUE)
table2.isw4[1,4]<-mean(cover.isw4.20[,4],na.rm=TRUE)
table2.isw4[1,5]<-mean(b1.isw4.20,na.rm=TRUE)
table2.isw4[1,6]<-mean(b1.isw4.sd.20,na.rm=TRUE)
table2.isw4[1,7]<-sd(b1.isw4.20,na.rm=TRUE)
table2.isw4[1,8]<-mean(cover.isw4.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.isw1.21 = c()
b0.isw1.sd.21 = c()
b1.isw1.21 = c()
b1.isw1.sd.21 = c()
cover.isw1.21=matrix(NA,2000,8)
colnames(cover.isw1.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.21 = c()
b0.isw2.sd.21 = c()
b1.isw2.21 = c()
b1.isw2.sd.21 = c()
cover.isw2.21=matrix(NA,2000,8)
colnames(cover.isw2.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.21 = c()
b0.isw3.sd.21 = c()
b1.isw3.21 = c()
b1.isw3.sd.21 = c()
cover.isw3.21=matrix(NA,2000,8)
colnames(cover.isw3.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.21 = c()
b0.isw4.sd.21 = c()
b1.isw4.21 = c()
b1.isw4.sd.21 = c()
cover.isw4.21=matrix(NA,2000,8)
colnames(cover.isw4.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.1)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.21[i] = ismb.fit[1,1]
    b0.isw1.sd.21[i] = ismb.fit[1,2]
    b1.isw1.21[i] = ismb.fit[2,1]
    b1.isw1.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.21[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.21[i,2] = ismb.fit[1,1]
    cover.isw1.21[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.21[i,4] = ind(1.219403, cover.isw1.21[i,1], cover.isw1.21[i,3])
    cover.isw1.21[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.21[i,6] = ismb.fit[2,1]
    cover.isw1.21[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.21[i,8] = ind(0.9070615, cover.isw1.21[i,5], cover.isw1.21[i,7])}
    , error=function(e){
      b0.isw1.21[i] = NA
      b0.isw1.sd.21[i] = NA
      b1.isw1.21[i] = NA
      b1.isw1.sd.21[i] = NA
      # Coverage
      cover.isw1.21[i,1] = NA
      cover.isw1.21[i,2] = NA
      cover.isw1.21[i,3] = NA
      cover.isw1.21[i,4] = NA
      cover.isw1.21[i,5] = NA
      cover.isw1.21[i,6] = NA
      cover.isw1.21[i,7] = NA
      cover.isw1.21[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw2.21[i] = ismb.fit[1,1]
    b0.isw2.sd.21[i] = ismb.fit[1,2]
    b1.isw2.21[i] = ismb.fit[2,1]
    b1.isw2.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.21[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.21[i,2] = ismb.fit[1,1]
    cover.isw2.21[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.21[i,4] = ind(1.219403, cover.isw2.21[i,1], cover.isw2.21[i,3])
    cover.isw2.21[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.21[i,6] = ismb.fit[2,1]
    cover.isw2.21[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.21[i,8] = ind(0.9070615, cover.isw2.21[i,5], cover.isw2.21[i,7])}
    , error=function(e){
      b0.isw2.21[i] = NA
      b0.isw2.sd.21[i] = NA
      b1.isw2.21[i] = NA
      b1.isw2.sd.21[i] = NA
      # Coverage
      cover.isw2.21[i,1] = NA
      cover.isw2.21[i,2] = NA
      cover.isw2.21[i,3] = NA
      cover.isw2.21[i,4] = NA
      cover.isw2.21[i,5] = NA
      cover.isw2.21[i,6] = NA
      cover.isw2.21[i,7] = NA
      cover.isw2.21[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw3.21[i] = ismb.fit[1,1]
    b0.isw3.sd.21[i] = ismb.fit[1,2]
    b1.isw3.21[i] = ismb.fit[2,1]
    b1.isw3.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.21[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.21[i,2] = ismb.fit[1,1]
    cover.isw3.21[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.21[i,4] = ind(1.219403, cover.isw3.21[i,1], cover.isw3.21[i,3])
    cover.isw3.21[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.21[i,6] = ismb.fit[2,1]
    cover.isw3.21[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.21[i,8] = ind(0.9070615, cover.isw3.21[i,5], cover.isw3.21[i,7])}
    , error=function(e){
      b0.isw3.21[i] = NA
      b0.isw3.sd.21[i] = NA
      b1.isw3.21[i] = NA
      b1.isw3.sd.21[i] = NA
      # Coverage
      cover.isw3.21[i,1] = NA
      cover.isw3.21[i,2] = NA
      cover.isw3.21[i,3] = NA
      cover.isw3.21[i,4] = NA
      cover.isw3.21[i,5] = NA
      cover.isw3.21[i,6] = NA
      cover.isw3.21[i,7] = NA
      cover.isw3.21[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw4.21[i] = ismb.fit[1,1]
    b0.isw4.sd.21[i] = ismb.fit[1,2]
    b1.isw4.21[i] = ismb.fit[2,1]
    b1.isw4.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.21[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.21[i,2] = ismb.fit[1,1]
    cover.isw4.21[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.21[i,4] = ind(1.219403, cover.isw4.21[i,1], cover.isw4.21[i,3])
    cover.isw4.21[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.21[i,6] = ismb.fit[2,1]
    cover.isw4.21[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.21[i,8] = ind(0.9070615, cover.isw4.21[i,5], cover.isw4.21[i,7])}
    , error=function(e){
      b0.isw4.21[i] = NA
      b0.isw4.sd.21[i] = NA
      b1.isw4.21[i] = NA
      b1.isw4.sd.21[i] = NA
      # Coverage
      cover.isw4.21[i,1] = NA
      cover.isw4.21[i,2] = NA
      cover.isw4.21[i,3] = NA
      cover.isw4.21[i,4] = NA
      cover.isw4.21[i,5] = NA
      cover.isw4.21[i,6] = NA
      cover.isw4.21[i,7] = NA
      cover.isw4.21[i,8] = NA
    })
}

# IS beta table
table2.isw1[2,1]<-mean(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,2]<-mean(b0.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,3]<-sd(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,4]<-mean(cover.isw1.21[,4],na.rm=TRUE)
table2.isw1[2,5]<-mean(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,6]<-mean(b1.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,7]<-sd(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,8]<-mean(cover.isw1.21[,8],na.rm=TRUE)

table2.isw2[2,1]<-mean(b0.isw2.21,na.rm=TRUE)
table2.isw2[2,2]<-mean(b0.isw2.sd.21,na.rm=TRUE)
table2.isw2[2,3]<-sd(b0.isw2.21,na.rm=TRUE)
table2.isw2[2,4]<-mean(cover.isw2.21[,4],na.rm=TRUE)
table2.isw2[2,5]<-mean(b1.isw2.21,na.rm=TRUE)
table2.isw2[2,6]<-mean(b1.isw2.sd.21,na.rm=TRUE)
table2.isw2[2,7]<-sd(b1.isw2.21,na.rm=TRUE)
table2.isw2[2,8]<-mean(cover.isw2.21[,8],na.rm=TRUE)

table2.isw3[2,1]<-mean(b0.isw3.21,na.rm=TRUE)
table2.isw3[2,2]<-mean(b0.isw3.sd.21,na.rm=TRUE)
table2.isw3[2,3]<-sd(b0.isw3.21,na.rm=TRUE)
table2.isw3[2,4]<-mean(cover.isw3.21[,4],na.rm=TRUE)
table2.isw3[2,5]<-mean(b1.isw3.21,na.rm=TRUE)
table2.isw3[2,6]<-mean(b1.isw3.sd.21,na.rm=TRUE)
table2.isw3[2,7]<-sd(b1.isw3.21,na.rm=TRUE)
table2.isw3[2,8]<-mean(cover.isw3.21[,8],na.rm=TRUE)

table2.isw4[2,1]<-mean(b0.isw4.21,na.rm=TRUE)
table2.isw4[2,2]<-mean(b0.isw4.sd.21,na.rm=TRUE)
table2.isw4[2,3]<-sd(b0.isw4.21,na.rm=TRUE)
table2.isw4[2,4]<-mean(cover.isw4.21[,4],na.rm=TRUE)
table2.isw4[2,5]<-mean(b1.isw4.21,na.rm=TRUE)
table2.isw4[2,6]<-mean(b1.isw4.sd.21,na.rm=TRUE)
table2.isw4[2,7]<-sd(b1.isw4.21,na.rm=TRUE)
table2.isw4[2,8]<-mean(cover.isw4.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.isw1.23 = c()
b0.isw1.sd.23 = c()
b1.isw1.23 = c()
b1.isw1.sd.23 = c()
cover.isw1.23=matrix(NA,2000,8)
colnames(cover.isw1.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.23 = c()
b0.isw2.sd.23 = c()
b1.isw2.23 = c()
b1.isw2.sd.23 = c()
cover.isw2.23=matrix(NA,2000,8)
colnames(cover.isw2.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.23 = c()
b0.isw3.sd.23 = c()
b1.isw3.23 = c()
b1.isw3.sd.23 = c()
cover.isw3.23=matrix(NA,2000,8)
colnames(cover.isw3.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.23 = c()
b0.isw4.sd.23 = c()
b1.isw4.23 = c()
b1.isw4.sd.23 = c()
cover.isw4.23=matrix(NA,2000,8)
colnames(cover.isw4.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.3)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.23[i] = ismb.fit[1,1]
    b0.isw1.sd.23[i] = ismb.fit[1,2]
    b1.isw1.23[i] = ismb.fit[2,1]
    b1.isw1.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.23[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.23[i,2] = ismb.fit[1,1]
    cover.isw1.23[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.23[i,4] = ind(1.219403, cover.isw1.23[i,1], cover.isw1.23[i,3])
    cover.isw1.23[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.23[i,6] = ismb.fit[2,1]
    cover.isw1.23[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.23[i,8] = ind(0.9070615, cover.isw1.23[i,5], cover.isw1.23[i,7])}
    , error=function(e){
      b0.isw1.23[i] = NA
      b0.isw1.sd.23[i] = NA
      b1.isw1.23[i] = NA
      b1.isw1.sd.23[i] = NA
      # Coverage
      cover.isw1.23[i,1] = NA
      cover.isw1.23[i,2] = NA
      cover.isw1.23[i,3] = NA
      cover.isw1.23[i,4] = NA
      cover.isw1.23[i,5] = NA
      cover.isw1.23[i,6] = NA
      cover.isw1.23[i,7] = NA
      cover.isw1.23[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw2.23[i] = ismb.fit[1,1]
    b0.isw2.sd.23[i] = ismb.fit[1,2]
    b1.isw2.23[i] = ismb.fit[2,1]
    b1.isw2.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.23[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.23[i,2] = ismb.fit[1,1]
    cover.isw2.23[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.23[i,4] = ind(1.219403, cover.isw2.23[i,1], cover.isw2.23[i,3])
    cover.isw2.23[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.23[i,6] = ismb.fit[2,1]
    cover.isw2.23[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.23[i,8] = ind(0.9070615, cover.isw2.23[i,5], cover.isw2.23[i,7])}
    , error=function(e){
      b0.isw2.23[i] = NA
      b0.isw2.sd.23[i] = NA
      b1.isw2.23[i] = NA
      b1.isw2.sd.23[i] = NA
      # Coverage
      cover.isw2.23[i,1] = NA
      cover.isw2.23[i,2] = NA
      cover.isw2.23[i,3] = NA
      cover.isw2.23[i,4] = NA
      cover.isw2.23[i,5] = NA
      cover.isw2.23[i,6] = NA
      cover.isw2.23[i,7] = NA
      cover.isw2.23[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw3.23[i] = ismb.fit[1,1]
    b0.isw3.sd.23[i] = ismb.fit[1,2]
    b1.isw3.23[i] = ismb.fit[2,1]
    b1.isw3.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.23[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.23[i,2] = ismb.fit[1,1]
    cover.isw3.23[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.23[i,4] = ind(1.219403, cover.isw3.23[i,1], cover.isw3.23[i,3])
    cover.isw3.23[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.23[i,6] = ismb.fit[2,1]
    cover.isw3.23[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.23[i,8] = ind(0.9070615, cover.isw3.23[i,5], cover.isw3.23[i,7])}
    , error=function(e){
      b0.isw3.23[i] = NA
      b0.isw3.sd.23[i] = NA
      b1.isw3.23[i] = NA
      b1.isw3.sd.23[i] = NA
      # Coverage
      cover.isw3.23[i,1] = NA
      cover.isw3.23[i,2] = NA
      cover.isw3.23[i,3] = NA
      cover.isw3.23[i,4] = NA
      cover.isw3.23[i,5] = NA
      cover.isw3.23[i,6] = NA
      cover.isw3.23[i,7] = NA
      cover.isw3.23[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw4.23[i] = ismb.fit[1,1]
    b0.isw4.sd.23[i] = ismb.fit[1,2]
    b1.isw4.23[i] = ismb.fit[2,1]
    b1.isw4.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.23[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.23[i,2] = ismb.fit[1,1]
    cover.isw4.23[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.23[i,4] = ind(1.219403, cover.isw4.23[i,1], cover.isw4.23[i,3])
    cover.isw4.23[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.23[i,6] = ismb.fit[2,1]
    cover.isw4.23[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.23[i,8] = ind(0.9070615, cover.isw4.23[i,5], cover.isw4.23[i,7])}
    , error=function(e){
      b0.isw4.23[i] = NA
      b0.isw4.sd.23[i] = NA
      b1.isw4.23[i] = NA
      b1.isw4.sd.23[i] = NA
      # Coverage
      cover.isw4.23[i,1] = NA
      cover.isw4.23[i,2] = NA
      cover.isw4.23[i,3] = NA
      cover.isw4.23[i,4] = NA
      cover.isw4.23[i,5] = NA
      cover.isw4.23[i,6] = NA
      cover.isw4.23[i,7] = NA
      cover.isw4.23[i,8] = NA
    })
}

# IS beta table
table2.isw1[3,1]<-mean(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,2]<-mean(b0.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,3]<-sd(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,4]<-mean(cover.isw1.23[,4],na.rm=TRUE)
table2.isw1[3,5]<-mean(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,6]<-mean(b1.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,7]<-sd(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,8]<-mean(cover.isw1.23[,8],na.rm=TRUE)

table2.isw2[3,1]<-mean(b0.isw2.23,na.rm=TRUE)
table2.isw2[3,2]<-mean(b0.isw2.sd.23,na.rm=TRUE)
table2.isw2[3,3]<-sd(b0.isw2.23,na.rm=TRUE)
table2.isw2[3,4]<-mean(cover.isw2.23[,4],na.rm=TRUE)
table2.isw2[3,5]<-mean(b1.isw2.23,na.rm=TRUE)
table2.isw2[3,6]<-mean(b1.isw2.sd.23,na.rm=TRUE)
table2.isw2[3,7]<-sd(b1.isw2.23,na.rm=TRUE)
table2.isw2[3,8]<-mean(cover.isw2.23[,8],na.rm=TRUE)

table2.isw3[3,1]<-mean(b0.isw3.23,na.rm=TRUE)
table2.isw3[3,2]<-mean(b0.isw3.sd.23,na.rm=TRUE)
table2.isw3[3,3]<-sd(b0.isw3.23,na.rm=TRUE)
table2.isw3[3,4]<-mean(cover.isw3.23[,4],na.rm=TRUE)
table2.isw3[3,5]<-mean(b1.isw3.23,na.rm=TRUE)
table2.isw3[3,6]<-mean(b1.isw3.sd.23,na.rm=TRUE)
table2.isw3[3,7]<-sd(b1.isw3.23,na.rm=TRUE)
table2.isw3[3,8]<-mean(cover.isw3.23[,8],na.rm=TRUE)

table2.isw4[3,1]<-mean(b0.isw4.23,na.rm=TRUE)
table2.isw4[3,2]<-mean(b0.isw4.sd.23,na.rm=TRUE)
table2.isw4[3,3]<-sd(b0.isw4.23,na.rm=TRUE)
table2.isw4[3,4]<-mean(cover.isw4.23[,4],na.rm=TRUE)
table2.isw4[3,5]<-mean(b1.isw4.23,na.rm=TRUE)
table2.isw4[3,6]<-mean(b1.isw4.sd.23,na.rm=TRUE)
table2.isw4[3,7]<-sd(b1.isw4.23,na.rm=TRUE)
table2.isw4[3,8]<-mean(cover.isw4.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.isw1.25 = c()
b0.isw1.sd.25 = c()
b1.isw1.25 = c()
b1.isw1.sd.25 = c()
cover.isw1.25=matrix(NA,2000,8)
colnames(cover.isw1.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.25 = c()
b0.isw2.sd.25 = c()
b1.isw2.25 = c()
b1.isw2.sd.25 = c()
cover.isw2.25=matrix(NA,2000,8)
colnames(cover.isw2.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.25 = c()
b0.isw3.sd.25 = c()
b1.isw3.25 = c()
b1.isw3.sd.25 = c()
cover.isw3.25=matrix(NA,2000,8)
colnames(cover.isw3.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.25 = c()
b0.isw4.sd.25 = c()
b1.isw4.25 = c()
b1.isw4.sd.25 = c()
cover.isw4.25=matrix(NA,2000,8)
colnames(cover.isw4.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.5)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.25[i] = ismb.fit[1,1]
    b0.isw1.sd.25[i] = ismb.fit[1,2]
    b1.isw1.25[i] = ismb.fit[2,1]
    b1.isw1.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.25[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.25[i,2] = ismb.fit[1,1]
    cover.isw1.25[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.25[i,4] = ind(1.219403, cover.isw1.25[i,1], cover.isw1.25[i,3])
    cover.isw1.25[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.25[i,6] = ismb.fit[2,1]
    cover.isw1.25[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.25[i,8] = ind(0.9070615, cover.isw1.25[i,5], cover.isw1.25[i,7])}
    , error=function(e){
      b0.isw1.25[i] = NA
      b0.isw1.sd.25[i] = NA
      b1.isw1.25[i] = NA
      b1.isw1.sd.25[i] = NA
      # Coverage
      cover.isw1.25[i,1] = NA
      cover.isw1.25[i,2] = NA
      cover.isw1.25[i,3] = NA
      cover.isw1.25[i,4] = NA
      cover.isw1.25[i,5] = NA
      cover.isw1.25[i,6] = NA
      cover.isw1.25[i,7] = NA
      cover.isw1.25[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw2.25[i] = ismb.fit[1,1]
    b0.isw2.sd.25[i] = ismb.fit[1,2]
    b1.isw2.25[i] = ismb.fit[2,1]
    b1.isw2.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.25[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.25[i,2] = ismb.fit[1,1]
    cover.isw2.25[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.25[i,4] = ind(1.219403, cover.isw2.25[i,1], cover.isw2.25[i,3])
    cover.isw2.25[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.25[i,6] = ismb.fit[2,1]
    cover.isw2.25[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.25[i,8] = ind(0.9070615, cover.isw2.25[i,5], cover.isw2.25[i,7])}
    , error=function(e){
      b0.isw2.25[i] = NA
      b0.isw2.sd.25[i] = NA
      b1.isw2.25[i] = NA
      b1.isw2.sd.25[i] = NA
      # Coverage
      cover.isw2.25[i,1] = NA
      cover.isw2.25[i,2] = NA
      cover.isw2.25[i,3] = NA
      cover.isw2.25[i,4] = NA
      cover.isw2.25[i,5] = NA
      cover.isw2.25[i,6] = NA
      cover.isw2.25[i,7] = NA
      cover.isw2.25[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw3.25[i] = ismb.fit[1,1]
    b0.isw3.sd.25[i] = ismb.fit[1,2]
    b1.isw3.25[i] = ismb.fit[2,1]
    b1.isw3.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.25[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.25[i,2] = ismb.fit[1,1]
    cover.isw3.25[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.25[i,4] = ind(1.219403, cover.isw3.25[i,1], cover.isw3.25[i,3])
    cover.isw3.25[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.25[i,6] = ismb.fit[2,1]
    cover.isw3.25[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.25[i,8] = ind(0.9070615, cover.isw3.25[i,5], cover.isw3.25[i,7])}
    , error=function(e){
      b0.isw3.25[i] = NA
      b0.isw3.sd.25[i] = NA
      b1.isw3.25[i] = NA
      b1.isw3.sd.25[i] = NA
      # Coverage
      cover.isw3.25[i,1] = NA
      cover.isw3.25[i,2] = NA
      cover.isw3.25[i,3] = NA
      cover.isw3.25[i,4] = NA
      cover.isw3.25[i,5] = NA
      cover.isw3.25[i,6] = NA
      cover.isw3.25[i,7] = NA
      cover.isw3.25[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw4.25[i] = ismb.fit[1,1]
    b0.isw4.sd.25[i] = ismb.fit[1,2]
    b1.isw4.25[i] = ismb.fit[2,1]
    b1.isw4.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.25[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.25[i,2] = ismb.fit[1,1]
    cover.isw4.25[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.25[i,4] = ind(1.219403, cover.isw4.25[i,1], cover.isw4.25[i,3])
    cover.isw4.25[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.25[i,6] = ismb.fit[2,1]
    cover.isw4.25[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.25[i,8] = ind(0.9070615, cover.isw4.25[i,5], cover.isw4.25[i,7])}
    , error=function(e){
      b0.isw4.25[i] = NA
      b0.isw4.sd.25[i] = NA
      b1.isw4.25[i] = NA
      b1.isw4.sd.25[i] = NA
      # Coverage
      cover.isw4.25[i,1] = NA
      cover.isw4.25[i,2] = NA
      cover.isw4.25[i,3] = NA
      cover.isw4.25[i,4] = NA
      cover.isw4.25[i,5] = NA
      cover.isw4.25[i,6] = NA
      cover.isw4.25[i,7] = NA
      cover.isw4.25[i,8] = NA
    })
}

# IS beta table
table2.isw1[4,1]<-mean(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,2]<-mean(b0.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,3]<-sd(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,4]<-mean(cover.isw1.25[,4],na.rm=TRUE)
table2.isw1[4,5]<-mean(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,6]<-mean(b1.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,7]<-sd(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,8]<-mean(cover.isw1.25[,8],na.rm=TRUE)

table2.isw2[4,1]<-mean(b0.isw2.25,na.rm=TRUE)
table2.isw2[4,2]<-mean(b0.isw2.sd.25,na.rm=TRUE)
table2.isw2[4,3]<-sd(b0.isw2.25,na.rm=TRUE)
table2.isw2[4,4]<-mean(cover.isw2.25[,4],na.rm=TRUE)
table2.isw2[4,5]<-mean(b1.isw2.25,na.rm=TRUE)
table2.isw2[4,6]<-mean(b1.isw2.sd.25,na.rm=TRUE)
table2.isw2[4,7]<-sd(b1.isw2.25,na.rm=TRUE)
table2.isw2[4,8]<-mean(cover.isw2.25[,8],na.rm=TRUE)

table2.isw3[4,1]<-mean(b0.isw3.25,na.rm=TRUE)
table2.isw3[4,2]<-mean(b0.isw3.sd.25,na.rm=TRUE)
table2.isw3[4,3]<-sd(b0.isw3.25,na.rm=TRUE)
table2.isw3[4,4]<-mean(cover.isw3.25[,4],na.rm=TRUE)
table2.isw3[4,5]<-mean(b1.isw3.25,na.rm=TRUE)
table2.isw3[4,6]<-mean(b1.isw3.sd.25,na.rm=TRUE)
table2.isw3[4,7]<-sd(b1.isw3.25,na.rm=TRUE)
table2.isw3[4,8]<-mean(cover.isw3.25[,8],na.rm=TRUE)

table2.isw4[4,1]<-mean(b0.isw4.25,na.rm=TRUE)
table2.isw4[4,2]<-mean(b0.isw4.sd.25,na.rm=TRUE)
table2.isw4[4,3]<-sd(b0.isw4.25,na.rm=TRUE)
table2.isw4[4,4]<-mean(cover.isw4.25[,4],na.rm=TRUE)
table2.isw4[4,5]<-mean(b1.isw4.25,na.rm=TRUE)
table2.isw4[4,6]<-mean(b1.isw4.sd.25,na.rm=TRUE)
table2.isw4[4,7]<-sd(b1.isw4.25,na.rm=TRUE)
table2.isw4[4,8]<-mean(cover.isw4.25[,8],na.rm=TRUE)

#### t_0=2 & c=70% ####
b0.isw1.27 = c()
b0.isw1.sd.27 = c()
b1.isw1.27 = c()
b1.isw1.sd.27 = c()
cover.isw1.27=matrix(NA,2000,8)
colnames(cover.isw1.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.27 = c()
b0.isw2.sd.27 = c()
b1.isw2.27 = c()
b1.isw2.sd.27 = c()
cover.isw2.27=matrix(NA,2000,8)
colnames(cover.isw2.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.27 = c()
b0.isw3.sd.27 = c()
b1.isw3.27 = c()
b1.isw3.sd.27 = c()
cover.isw3.27=matrix(NA,2000,8)
colnames(cover.isw3.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.27 = c()
b0.isw4.sd.27 = c()
b1.isw4.27 = c()
b1.isw4.sd.27 = c()
cover.isw4.27=matrix(NA,2000,8)
colnames(cover.isw4.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.7)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.27[i] = ismb.fit[1,1]
    b0.isw1.sd.27[i] = ismb.fit[1,2]
    b1.isw1.27[i] = ismb.fit[2,1]
    b1.isw1.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.27[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.27[i,2] = ismb.fit[1,1]
    cover.isw1.27[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.27[i,4] = ind(1.219403, cover.isw1.27[i,1], cover.isw1.27[i,3])
    cover.isw1.27[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.27[i,6] = ismb.fit[2,1]
    cover.isw1.27[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.27[i,8] = ind(0.9070615, cover.isw1.27[i,5], cover.isw1.27[i,7])}
    , error=function(e){
      b0.isw1.27[i] = NA
      b0.isw1.sd.27[i] = NA
      b1.isw1.27[i] = NA
      b1.isw1.sd.27[i] = NA
      # Coverage
      cover.isw1.27[i,1] = NA
      cover.isw1.27[i,2] = NA
      cover.isw1.27[i,3] = NA
      cover.isw1.27[i,4] = NA
      cover.isw1.27[i,5] = NA
      cover.isw1.27[i,6] = NA
      cover.isw1.27[i,7] = NA
      cover.isw1.27[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw2.27[i] = ismb.fit[1,1]
    b0.isw2.sd.27[i] = ismb.fit[1,2]
    b1.isw2.27[i] = ismb.fit[2,1]
    b1.isw2.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.27[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.27[i,2] = ismb.fit[1,1]
    cover.isw2.27[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.27[i,4] = ind(1.219403, cover.isw2.27[i,1], cover.isw2.27[i,3])
    cover.isw2.27[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.27[i,6] = ismb.fit[2,1]
    cover.isw2.27[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.27[i,8] = ind(0.9070615, cover.isw2.27[i,5], cover.isw2.27[i,7])}
    , error=function(e){
      b0.isw2.27[i] = NA
      b0.isw2.sd.27[i] = NA
      b1.isw2.27[i] = NA
      b1.isw2.sd.27[i] = NA
      # Coverage
      cover.isw2.27[i,1] = NA
      cover.isw2.27[i,2] = NA
      cover.isw2.27[i,3] = NA
      cover.isw2.27[i,4] = NA
      cover.isw2.27[i,5] = NA
      cover.isw2.27[i,6] = NA
      cover.isw2.27[i,7] = NA
      cover.isw2.27[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw3.27[i] = ismb.fit[1,1]
    b0.isw3.sd.27[i] = ismb.fit[1,2]
    b1.isw3.27[i] = ismb.fit[2,1]
    b1.isw3.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.27[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.27[i,2] = ismb.fit[1,1]
    cover.isw3.27[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.27[i,4] = ind(1.219403, cover.isw3.27[i,1], cover.isw3.27[i,3])
    cover.isw3.27[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.27[i,6] = ismb.fit[2,1]
    cover.isw3.27[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.27[i,8] = ind(0.9070615, cover.isw3.27[i,5], cover.isw3.27[i,7])}
    , error=function(e){
      b0.isw3.27[i] = NA
      b0.isw3.sd.27[i] = NA
      b1.isw3.27[i] = NA
      b1.isw3.sd.27[i] = NA
      # Coverage
      cover.isw3.27[i,1] = NA
      cover.isw3.27[i,2] = NA
      cover.isw3.27[i,3] = NA
      cover.isw3.27[i,4] = NA
      cover.isw3.27[i,5] = NA
      cover.isw3.27[i,6] = NA
      cover.isw3.27[i,7] = NA
      cover.isw3.27[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw4.27[i] = ismb.fit[1,1]
    b0.isw4.sd.27[i] = ismb.fit[1,2]
    b1.isw4.27[i] = ismb.fit[2,1]
    b1.isw4.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.27[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.27[i,2] = ismb.fit[1,1]
    cover.isw4.27[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.27[i,4] = ind(1.219403, cover.isw4.27[i,1], cover.isw4.27[i,3])
    cover.isw4.27[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.27[i,6] = ismb.fit[2,1]
    cover.isw4.27[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.27[i,8] = ind(0.9070615, cover.isw4.27[i,5], cover.isw4.27[i,7])}
    , error=function(e){
      b0.isw4.27[i] = NA
      b0.isw4.sd.27[i] = NA
      b1.isw4.27[i] = NA
      b1.isw4.sd.27[i] = NA
      # Coverage
      cover.isw4.27[i,1] = NA
      cover.isw4.27[i,2] = NA
      cover.isw4.27[i,3] = NA
      cover.isw4.27[i,4] = NA
      cover.isw4.27[i,5] = NA
      cover.isw4.27[i,6] = NA
      cover.isw4.27[i,7] = NA
      cover.isw4.27[i,8] = NA
    })
}

# IS beta table
table2.isw1[5,1]<-mean(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,2]<-mean(b0.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,3]<-sd(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,4]<-mean(cover.isw1.27[,4],na.rm=TRUE)
table2.isw1[5,5]<-mean(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,6]<-mean(b1.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,7]<-sd(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,8]<-mean(cover.isw1.27[,8],na.rm=TRUE)

table2.isw2[5,1]<-mean(b0.isw2.27,na.rm=TRUE)
table2.isw2[5,2]<-mean(b0.isw2.sd.27,na.rm=TRUE)
table2.isw2[5,3]<-sd(b0.isw2.27,na.rm=TRUE)
table2.isw2[5,4]<-mean(cover.isw2.27[,4],na.rm=TRUE)
table2.isw2[5,5]<-mean(b1.isw2.27,na.rm=TRUE)
table2.isw2[5,6]<-mean(b1.isw2.sd.27,na.rm=TRUE)
table2.isw2[5,7]<-sd(b1.isw2.27,na.rm=TRUE)
table2.isw2[5,8]<-mean(cover.isw2.27[,8],na.rm=TRUE)

table2.isw3[5,1]<-mean(b0.isw3.27,na.rm=TRUE)
table2.isw3[5,2]<-mean(b0.isw3.sd.27,na.rm=TRUE)
table2.isw3[5,3]<-sd(b0.isw3.27,na.rm=TRUE)
table2.isw3[5,4]<-mean(cover.isw3.27[,4],na.rm=TRUE)
table2.isw3[5,5]<-mean(b1.isw3.27,na.rm=TRUE)
table2.isw3[5,6]<-mean(b1.isw3.sd.27,na.rm=TRUE)
table2.isw3[5,7]<-sd(b1.isw3.27,na.rm=TRUE)
table2.isw3[5,8]<-mean(cover.isw3.27[,8],na.rm=TRUE)

table2.isw4[5,1]<-mean(b0.isw4.27,na.rm=TRUE)
table2.isw4[5,2]<-mean(b0.isw4.sd.27,na.rm=TRUE)
table2.isw4[5,3]<-sd(b0.isw4.27,na.rm=TRUE)
table2.isw4[5,4]<-mean(cover.isw4.27[,4],na.rm=TRUE)
table2.isw4[5,5]<-mean(b1.isw4.27,na.rm=TRUE)
table2.isw4[5,6]<-mean(b1.isw4.sd.27,na.rm=TRUE)
table2.isw4[5,7]<-sd(b1.isw4.27,na.rm=TRUE)
table2.isw4[5,8]<-mean(cover.isw4.27[,8],na.rm=TRUE)

#### censoring point at t_0=3 ####
c.0=5000000
c.1=61.17
c.3=22.53
c.5=13.43
c.7=8.5
#### t_0=3 & c=0% ####
b0.isw1.30 = c()
b0.isw1.sd.30 = c()
b1.isw1.30 = c()
b1.isw1.sd.30 = c()
cover.isw1.30=matrix(NA,2000,8)
colnames(cover.isw1.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.30 = c()
b0.isw2.sd.30 = c()
b1.isw2.30 = c()
b1.isw2.sd.30 = c()
cover.isw2.30=matrix(NA,2000,8)
colnames(cover.isw2.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.30 = c()
b0.isw3.sd.30 = c()
b1.isw3.30 = c()
b1.isw3.sd.30 = c()
cover.isw3.30=matrix(NA,2000,8)
colnames(cover.isw3.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.30 = c()
b0.isw4.sd.30 = c()
b1.isw4.30 = c()
b1.isw4.sd.30 = c()
cover.isw4.30=matrix(NA,2000,8)
colnames(cover.isw4.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.30[i] = ismb.fit[1,1]
    b0.isw1.sd.30[i] = ismb.fit[1,2]
    b1.isw1.30[i] = ismb.fit[2,1]
    b1.isw1.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.30[i,2] = ismb.fit[1,1]
    cover.isw1.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.30[i,4] = ind(1.040613, cover.isw1.30[i,1], cover.isw1.30[i,3])
    cover.isw1.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.30[i,6] = ismb.fit[2,1]
    cover.isw1.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.30[i,8] = ind(1.0174711, cover.isw1.30[i,5], cover.isw1.30[i,7])}
    , error=function(e){
      b0.isw1.30[i] = NA
      b0.isw1.sd.30[i] = NA
      b1.isw1.30[i] = NA
      b1.isw1.sd.30[i] = NA
      # Coverage
      cover.isw1.30[i,1] = NA
      cover.isw1.30[i,2] = NA
      cover.isw1.30[i,3] = NA
      cover.isw1.30[i,4] = NA
      cover.isw1.30[i,5] = NA
      cover.isw1.30[i,6] = NA
      cover.isw1.30[i,7] = NA
      cover.isw1.30[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw2.30[i] = ismb.fit[1,1]
    b0.isw2.sd.30[i] = ismb.fit[1,2]
    b1.isw2.30[i] = ismb.fit[2,1]
    b1.isw2.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.30[i,2] = ismb.fit[1,1]
    cover.isw2.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.30[i,4] = ind(1.040613, cover.isw2.30[i,1], cover.isw2.30[i,3])
    cover.isw2.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.30[i,6] = ismb.fit[2,1]
    cover.isw2.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.30[i,8] = ind(1.0174711, cover.isw2.30[i,5], cover.isw2.30[i,7])}
    , error=function(e){
      b0.isw2.30[i] = NA
      b0.isw2.sd.30[i] = NA
      b1.isw2.30[i] = NA
      b1.isw2.sd.30[i] = NA
      # Coverage
      cover.isw2.30[i,1] = NA
      cover.isw2.30[i,2] = NA
      cover.isw2.30[i,3] = NA
      cover.isw2.30[i,4] = NA
      cover.isw2.30[i,5] = NA
      cover.isw2.30[i,6] = NA
      cover.isw2.30[i,7] = NA
      cover.isw2.30[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw3.30[i] = ismb.fit[1,1]
    b0.isw3.sd.30[i] = ismb.fit[1,2]
    b1.isw3.30[i] = ismb.fit[2,1]
    b1.isw3.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.30[i,2] = ismb.fit[1,1]
    cover.isw3.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.30[i,4] = ind(1.040613, cover.isw3.30[i,1], cover.isw3.30[i,3])
    cover.isw3.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.30[i,6] = ismb.fit[2,1]
    cover.isw3.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.30[i,8] = ind(1.0174711, cover.isw3.30[i,5], cover.isw3.30[i,7])}
    , error=function(e){
      b0.isw3.30[i] = NA
      b0.isw3.sd.30[i] = NA
      b1.isw3.30[i] = NA
      b1.isw3.sd.30[i] = NA
      # Coverage
      cover.isw3.30[i,1] = NA
      cover.isw3.30[i,2] = NA
      cover.isw3.30[i,3] = NA
      cover.isw3.30[i,4] = NA
      cover.isw3.30[i,5] = NA
      cover.isw3.30[i,6] = NA
      cover.isw3.30[i,7] = NA
      cover.isw3.30[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw4.30[i] = ismb.fit[1,1]
    b0.isw4.sd.30[i] = ismb.fit[1,2]
    b1.isw4.30[i] = ismb.fit[2,1]
    b1.isw4.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.30[i,2] = ismb.fit[1,1]
    cover.isw4.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.30[i,4] = ind(1.040613, cover.isw4.30[i,1], cover.isw4.30[i,3])
    cover.isw4.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.30[i,6] = ismb.fit[2,1]
    cover.isw4.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.30[i,8] = ind(1.0174711, cover.isw4.30[i,5], cover.isw4.30[i,7])}
    , error=function(e){
      b0.isw4.30[i] = NA
      b0.isw4.sd.30[i] = NA
      b1.isw4.30[i] = NA
      b1.isw4.sd.30[i] = NA
      # Coverage
      cover.isw4.30[i,1] = NA
      cover.isw4.30[i,2] = NA
      cover.isw4.30[i,3] = NA
      cover.isw4.30[i,4] = NA
      cover.isw4.30[i,5] = NA
      cover.isw4.30[i,6] = NA
      cover.isw4.30[i,7] = NA
      cover.isw4.30[i,8] = NA
    })
}

# IS beta table
table3.isw1[1,1]<-mean(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,2]<-mean(b0.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,3]<-sd(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,4]<-mean(cover.isw1.30[,4],na.rm=TRUE)
table3.isw1[1,5]<-mean(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,6]<-mean(b1.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,7]<-sd(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,8]<-mean(cover.isw1.30[,8],na.rm=TRUE)

table3.isw2[1,1]<-mean(b0.isw2.30,na.rm=TRUE)
table3.isw2[1,2]<-mean(b0.isw2.sd.30,na.rm=TRUE)
table3.isw2[1,3]<-sd(b0.isw2.30,na.rm=TRUE)
table3.isw2[1,4]<-mean(cover.isw2.30[,4],na.rm=TRUE)
table3.isw2[1,5]<-mean(b1.isw2.30,na.rm=TRUE)
table3.isw2[1,6]<-mean(b1.isw2.sd.30,na.rm=TRUE)
table3.isw2[1,7]<-sd(b1.isw2.30,na.rm=TRUE)
table3.isw2[1,8]<-mean(cover.isw2.30[,8],na.rm=TRUE)

table3.isw3[1,1]<-mean(b0.isw3.30,na.rm=TRUE)
table3.isw3[1,2]<-mean(b0.isw3.sd.30,na.rm=TRUE)
table3.isw3[1,3]<-sd(b0.isw3.30,na.rm=TRUE)
table3.isw3[1,4]<-mean(cover.isw3.30[,4],na.rm=TRUE)
table3.isw3[1,5]<-mean(b1.isw3.30,na.rm=TRUE)
table3.isw3[1,6]<-mean(b1.isw3.sd.30,na.rm=TRUE)
table3.isw3[1,7]<-sd(b1.isw3.30,na.rm=TRUE)
table3.isw3[1,8]<-mean(cover.isw3.30[,8],na.rm=TRUE)

table3.isw4[1,1]<-mean(b0.isw4.30,na.rm=TRUE)
table3.isw4[1,2]<-mean(b0.isw4.sd.30,na.rm=TRUE)
table3.isw4[1,3]<-sd(b0.isw4.30,na.rm=TRUE)
table3.isw4[1,4]<-mean(cover.isw4.30[,4],na.rm=TRUE)
table3.isw4[1,5]<-mean(b1.isw4.30,na.rm=TRUE)
table3.isw4[1,6]<-mean(b1.isw4.sd.30,na.rm=TRUE)
table3.isw4[1,7]<-sd(b1.isw4.30,na.rm=TRUE)
table3.isw4[1,8]<-mean(cover.isw4.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.isw1.31 = c()
b0.isw1.sd.31 = c()
b1.isw1.31 = c()
b1.isw1.sd.31 = c()
cover.isw1.31=matrix(NA,2000,8)
colnames(cover.isw1.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.31 = c()
b0.isw2.sd.31 = c()
b1.isw2.31 = c()
b1.isw2.sd.31 = c()
cover.isw2.31=matrix(NA,2000,8)
colnames(cover.isw2.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.31 = c()
b0.isw3.sd.31 = c()
b1.isw3.31 = c()
b1.isw3.sd.31 = c()
cover.isw3.31=matrix(NA,2000,8)
colnames(cover.isw3.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.31 = c()
b0.isw4.sd.31 = c()
b1.isw4.31 = c()
b1.isw4.sd.31 = c()
cover.isw4.31=matrix(NA,2000,8)
colnames(cover.isw4.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.1)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.31[i] = ismb.fit[1,1]
    b0.isw1.sd.31[i] = ismb.fit[1,2]
    b1.isw1.31[i] = ismb.fit[2,1]
    b1.isw1.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.31[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.31[i,2] = ismb.fit[1,1]
    cover.isw1.31[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.31[i,4] = ind(1.040613, cover.isw1.31[i,1], cover.isw1.31[i,3])
    cover.isw1.31[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.31[i,6] = ismb.fit[2,1]
    cover.isw1.31[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.31[i,8] = ind(1.0174711, cover.isw1.31[i,5], cover.isw1.31[i,7])}
    , error=function(e){
      b0.isw1.31[i] = NA
      b0.isw1.sd.31[i] = NA
      b1.isw1.31[i] = NA
      b1.isw1.sd.31[i] = NA
      # Coverage
      cover.isw1.31[i,1] = NA
      cover.isw1.31[i,2] = NA
      cover.isw1.31[i,3] = NA
      cover.isw1.31[i,4] = NA
      cover.isw1.31[i,5] = NA
      cover.isw1.31[i,6] = NA
      cover.isw1.31[i,7] = NA
      cover.isw1.31[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw2.31[i] = ismb.fit[1,1]
    b0.isw2.sd.31[i] = ismb.fit[1,2]
    b1.isw2.31[i] = ismb.fit[2,1]
    b1.isw2.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.31[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.31[i,2] = ismb.fit[1,1]
    cover.isw2.31[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.31[i,4] = ind(1.040613, cover.isw2.31[i,1], cover.isw2.31[i,3])
    cover.isw2.31[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.31[i,6] = ismb.fit[2,1]
    cover.isw2.31[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.31[i,8] = ind(1.0174711, cover.isw2.31[i,5], cover.isw2.31[i,7])}
    , error=function(e){
      b0.isw2.31[i] = NA
      b0.isw2.sd.31[i] = NA
      b1.isw2.31[i] = NA
      b1.isw2.sd.31[i] = NA
      # Coverage
      cover.isw2.31[i,1] = NA
      cover.isw2.31[i,2] = NA
      cover.isw2.31[i,3] = NA
      cover.isw2.31[i,4] = NA
      cover.isw2.31[i,5] = NA
      cover.isw2.31[i,6] = NA
      cover.isw2.31[i,7] = NA
      cover.isw2.31[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw3.31[i] = ismb.fit[1,1]
    b0.isw3.sd.31[i] = ismb.fit[1,2]
    b1.isw3.31[i] = ismb.fit[2,1]
    b1.isw3.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.31[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.31[i,2] = ismb.fit[1,1]
    cover.isw3.31[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.31[i,4] = ind(1.040613, cover.isw3.31[i,1], cover.isw3.31[i,3])
    cover.isw3.31[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.31[i,6] = ismb.fit[2,1]
    cover.isw3.31[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.31[i,8] = ind(1.0174711, cover.isw3.31[i,5], cover.isw3.31[i,7])}
    , error=function(e){
      b0.isw3.31[i] = NA
      b0.isw3.sd.31[i] = NA
      b1.isw3.31[i] = NA
      b1.isw3.sd.31[i] = NA
      # Coverage
      cover.isw3.31[i,1] = NA
      cover.isw3.31[i,2] = NA
      cover.isw3.31[i,3] = NA
      cover.isw3.31[i,4] = NA
      cover.isw3.31[i,5] = NA
      cover.isw3.31[i,6] = NA
      cover.isw3.31[i,7] = NA
      cover.isw3.31[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw4.31[i] = ismb.fit[1,1]
    b0.isw4.sd.31[i] = ismb.fit[1,2]
    b1.isw4.31[i] = ismb.fit[2,1]
    b1.isw4.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.31[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.31[i,2] = ismb.fit[1,1]
    cover.isw4.31[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.31[i,4] = ind(1.040613, cover.isw4.31[i,1], cover.isw4.31[i,3])
    cover.isw4.31[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.31[i,6] = ismb.fit[2,1]
    cover.isw4.31[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.31[i,8] = ind(1.0174711, cover.isw4.31[i,5], cover.isw4.31[i,7])}
    , error=function(e){
      b0.isw4.31[i] = NA
      b0.isw4.sd.31[i] = NA
      b1.isw4.31[i] = NA
      b1.isw4.sd.31[i] = NA
      # Coverage
      cover.isw4.31[i,1] = NA
      cover.isw4.31[i,2] = NA
      cover.isw4.31[i,3] = NA
      cover.isw4.31[i,4] = NA
      cover.isw4.31[i,5] = NA
      cover.isw4.31[i,6] = NA
      cover.isw4.31[i,7] = NA
      cover.isw4.31[i,8] = NA
    })
}

# IS beta table
table3.isw1[2,1]<-mean(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,2]<-mean(b0.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,3]<-sd(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,4]<-mean(cover.isw1.31[,4],na.rm=TRUE)
table3.isw1[2,5]<-mean(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,6]<-mean(b1.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,7]<-sd(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,8]<-mean(cover.isw1.31[,8],na.rm=TRUE)

table3.isw2[2,1]<-mean(b0.isw2.31,na.rm=TRUE)
table3.isw2[2,2]<-mean(b0.isw2.sd.31,na.rm=TRUE)
table3.isw2[2,3]<-sd(b0.isw2.31,na.rm=TRUE)
table3.isw2[2,4]<-mean(cover.isw2.31[,4],na.rm=TRUE)
table3.isw2[2,5]<-mean(b1.isw2.31,na.rm=TRUE)
table3.isw2[2,6]<-mean(b1.isw2.sd.31,na.rm=TRUE)
table3.isw2[2,7]<-sd(b1.isw2.31,na.rm=TRUE)
table3.isw2[2,8]<-mean(cover.isw2.31[,8],na.rm=TRUE)

table3.isw3[2,1]<-mean(b0.isw3.31,na.rm=TRUE)
table3.isw3[2,2]<-mean(b0.isw3.sd.31,na.rm=TRUE)
table3.isw3[2,3]<-sd(b0.isw3.31,na.rm=TRUE)
table3.isw3[2,4]<-mean(cover.isw3.31[,4],na.rm=TRUE)
table3.isw3[2,5]<-mean(b1.isw3.31,na.rm=TRUE)
table3.isw3[2,6]<-mean(b1.isw3.sd.31,na.rm=TRUE)
table3.isw3[2,7]<-sd(b1.isw3.31,na.rm=TRUE)
table3.isw3[2,8]<-mean(cover.isw3.31[,8],na.rm=TRUE)

table3.isw4[2,1]<-mean(b0.isw4.31,na.rm=TRUE)
table3.isw4[2,2]<-mean(b0.isw4.sd.31,na.rm=TRUE)
table3.isw4[2,3]<-sd(b0.isw4.31,na.rm=TRUE)
table3.isw4[2,4]<-mean(cover.isw4.31[,4],na.rm=TRUE)
table3.isw4[2,5]<-mean(b1.isw4.31,na.rm=TRUE)
table3.isw4[2,6]<-mean(b1.isw4.sd.31,na.rm=TRUE)
table3.isw4[2,7]<-sd(b1.isw4.31,na.rm=TRUE)
table3.isw4[2,8]<-mean(cover.isw4.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
b0.isw1.33 = c()
b0.isw1.sd.33 = c()
b1.isw1.33 = c()
b1.isw1.sd.33 = c()
cover.isw1.33=matrix(NA,2000,8)
colnames(cover.isw1.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.33 = c()
b0.isw2.sd.33 = c()
b1.isw2.33 = c()
b1.isw2.sd.33 = c()
cover.isw2.33=matrix(NA,2000,8)
colnames(cover.isw2.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.33 = c()
b0.isw3.sd.33 = c()
b1.isw3.33 = c()
b1.isw3.sd.33 = c()
cover.isw3.33=matrix(NA,2000,8)
colnames(cover.isw3.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.33 = c()
b0.isw4.sd.33 = c()
b1.isw4.33 = c()
b1.isw4.sd.33 = c()
cover.isw4.33=matrix(NA,2000,8)
colnames(cover.isw4.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.3)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.33[i] = ismb.fit[1,1]
    b0.isw1.sd.33[i] = ismb.fit[1,2]
    b1.isw1.33[i] = ismb.fit[2,1]
    b1.isw1.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.33[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.33[i,2] = ismb.fit[1,1]
    cover.isw1.33[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.33[i,4] = ind(1.040613, cover.isw1.33[i,1], cover.isw1.33[i,3])
    cover.isw1.33[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.33[i,6] = ismb.fit[2,1]
    cover.isw1.33[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.33[i,8] = ind(1.0174711, cover.isw1.33[i,5], cover.isw1.33[i,7])}
    , error=function(e){
      b0.isw1.33[i] = NA
      b0.isw1.sd.33[i] = NA
      b1.isw1.33[i] = NA
      b1.isw1.sd.33[i] = NA
      # Coverage
      cover.isw1.33[i,1] = NA
      cover.isw1.33[i,2] = NA
      cover.isw1.33[i,3] = NA
      cover.isw1.33[i,4] = NA
      cover.isw1.33[i,5] = NA
      cover.isw1.33[i,6] = NA
      cover.isw1.33[i,7] = NA
      cover.isw1.33[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw2.33[i] = ismb.fit[1,1]
    b0.isw2.sd.33[i] = ismb.fit[1,2]
    b1.isw2.33[i] = ismb.fit[2,1]
    b1.isw2.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.33[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.33[i,2] = ismb.fit[1,1]
    cover.isw2.33[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.33[i,4] = ind(1.040613, cover.isw2.33[i,1], cover.isw2.33[i,3])
    cover.isw2.33[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.33[i,6] = ismb.fit[2,1]
    cover.isw2.33[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.33[i,8] = ind(1.0174711, cover.isw2.33[i,5], cover.isw2.33[i,7])}
    , error=function(e){
      b0.isw2.33[i] = NA
      b0.isw2.sd.33[i] = NA
      b1.isw2.33[i] = NA
      b1.isw2.sd.33[i] = NA
      # Coverage
      cover.isw2.33[i,1] = NA
      cover.isw2.33[i,2] = NA
      cover.isw2.33[i,3] = NA
      cover.isw2.33[i,4] = NA
      cover.isw2.33[i,5] = NA
      cover.isw2.33[i,6] = NA
      cover.isw2.33[i,7] = NA
      cover.isw2.33[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw3.33[i] = ismb.fit[1,1]
    b0.isw3.sd.33[i] = ismb.fit[1,2]
    b1.isw3.33[i] = ismb.fit[2,1]
    b1.isw3.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.33[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.33[i,2] = ismb.fit[1,1]
    cover.isw3.33[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.33[i,4] = ind(1.040613, cover.isw3.33[i,1], cover.isw3.33[i,3])
    cover.isw3.33[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.33[i,6] = ismb.fit[2,1]
    cover.isw3.33[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.33[i,8] = ind(1.0174711, cover.isw3.33[i,5], cover.isw3.33[i,7])}
    , error=function(e){
      b0.isw3.33[i] = NA
      b0.isw3.sd.33[i] = NA
      b1.isw3.33[i] = NA
      b1.isw3.sd.33[i] = NA
      # Coverage
      cover.isw3.33[i,1] = NA
      cover.isw3.33[i,2] = NA
      cover.isw3.33[i,3] = NA
      cover.isw3.33[i,4] = NA
      cover.isw3.33[i,5] = NA
      cover.isw3.33[i,6] = NA
      cover.isw3.33[i,7] = NA
      cover.isw3.33[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw4.33[i] = ismb.fit[1,1]
    b0.isw4.sd.33[i] = ismb.fit[1,2]
    b1.isw4.33[i] = ismb.fit[2,1]
    b1.isw4.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.33[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.33[i,2] = ismb.fit[1,1]
    cover.isw4.33[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.33[i,4] = ind(1.040613, cover.isw4.33[i,1], cover.isw4.33[i,3])
    cover.isw4.33[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.33[i,6] = ismb.fit[2,1]
    cover.isw4.33[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.33[i,8] = ind(1.0174711, cover.isw4.33[i,5], cover.isw4.33[i,7])}
    , error=function(e){
      b0.isw4.33[i] = NA
      b0.isw4.sd.33[i] = NA
      b1.isw4.33[i] = NA
      b1.isw4.sd.33[i] = NA
      # Coverage
      cover.isw4.33[i,1] = NA
      cover.isw4.33[i,2] = NA
      cover.isw4.33[i,3] = NA
      cover.isw4.33[i,4] = NA
      cover.isw4.33[i,5] = NA
      cover.isw4.33[i,6] = NA
      cover.isw4.33[i,7] = NA
      cover.isw4.33[i,8] = NA
    })
}

# IS beta table
table3.isw1[3,1]<-mean(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,2]<-mean(b0.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,3]<-sd(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,4]<-mean(cover.isw1.33[,4],na.rm=TRUE)
table3.isw1[3,5]<-mean(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,6]<-mean(b1.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,7]<-sd(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,8]<-mean(cover.isw1.33[,8],na.rm=TRUE)

table3.isw2[3,1]<-mean(b0.isw2.33,na.rm=TRUE)
table3.isw2[3,2]<-mean(b0.isw2.sd.33,na.rm=TRUE)
table3.isw2[3,3]<-sd(b0.isw2.33,na.rm=TRUE)
table3.isw2[3,4]<-mean(cover.isw2.33[,4],na.rm=TRUE)
table3.isw2[3,5]<-mean(b1.isw2.33,na.rm=TRUE)
table3.isw2[3,6]<-mean(b1.isw2.sd.33,na.rm=TRUE)
table3.isw2[3,7]<-sd(b1.isw2.33,na.rm=TRUE)
table3.isw2[3,8]<-mean(cover.isw2.33[,8],na.rm=TRUE)

table3.isw3[3,1]<-mean(b0.isw3.33,na.rm=TRUE)
table3.isw3[3,2]<-mean(b0.isw3.sd.33,na.rm=TRUE)
table3.isw3[3,3]<-sd(b0.isw3.33,na.rm=TRUE)
table3.isw3[3,4]<-mean(cover.isw3.33[,4],na.rm=TRUE)
table3.isw3[3,5]<-mean(b1.isw3.33,na.rm=TRUE)
table3.isw3[3,6]<-mean(b1.isw3.sd.33,na.rm=TRUE)
table3.isw3[3,7]<-sd(b1.isw3.33,na.rm=TRUE)
table3.isw3[3,8]<-mean(cover.isw3.33[,8],na.rm=TRUE)

table3.isw4[3,1]<-mean(b0.isw4.33,na.rm=TRUE)
table3.isw4[3,2]<-mean(b0.isw4.sd.33,na.rm=TRUE)
table3.isw4[3,3]<-sd(b0.isw4.33,na.rm=TRUE)
table3.isw4[3,4]<-mean(cover.isw4.33[,4],na.rm=TRUE)
table3.isw4[3,5]<-mean(b1.isw4.33,na.rm=TRUE)
table3.isw4[3,6]<-mean(b1.isw4.sd.33,na.rm=TRUE)
table3.isw4[3,7]<-sd(b1.isw4.33,na.rm=TRUE)
table3.isw4[3,8]<-mean(cover.isw4.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.isw1.35 = c()
b0.isw1.sd.35 = c()
b1.isw1.35 = c()
b1.isw1.sd.35 = c()
cover.isw1.35=matrix(NA,2000,8)
colnames(cover.isw1.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.35 = c()
b0.isw2.sd.35 = c()
b1.isw2.35 = c()
b1.isw2.sd.35 = c()
cover.isw2.35=matrix(NA,2000,8)
colnames(cover.isw2.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.35 = c()
b0.isw3.sd.35 = c()
b1.isw3.35 = c()
b1.isw3.sd.35 = c()
cover.isw3.35=matrix(NA,2000,8)
colnames(cover.isw3.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.35 = c()
b0.isw4.sd.35 = c()
b1.isw4.35 = c()
b1.isw4.sd.35 = c()
cover.isw4.35=matrix(NA,2000,8)
colnames(cover.isw4.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.5)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.35[i] = ismb.fit[1,1]
    b0.isw1.sd.35[i] = ismb.fit[1,2]
    b1.isw1.35[i] = ismb.fit[2,1]
    b1.isw1.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.35[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.35[i,2] = ismb.fit[1,1]
    cover.isw1.35[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.35[i,4] = ind(1.040613, cover.isw1.35[i,1], cover.isw1.35[i,3])
    cover.isw1.35[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.35[i,6] = ismb.fit[2,1]
    cover.isw1.35[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.35[i,8] = ind(1.0174711, cover.isw1.35[i,5], cover.isw1.35[i,7])}
    , error=function(e){
      b0.isw1.35[i] = NA
      b0.isw1.sd.35[i] = NA
      b1.isw1.35[i] = NA
      b1.isw1.sd.35[i] = NA
      # Coverage
      cover.isw1.35[i,1] = NA
      cover.isw1.35[i,2] = NA
      cover.isw1.35[i,3] = NA
      cover.isw1.35[i,4] = NA
      cover.isw1.35[i,5] = NA
      cover.isw1.35[i,6] = NA
      cover.isw1.35[i,7] = NA
      cover.isw1.35[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw2.35[i] = ismb.fit[1,1]
    b0.isw2.sd.35[i] = ismb.fit[1,2]
    b1.isw2.35[i] = ismb.fit[2,1]
    b1.isw2.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.35[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.35[i,2] = ismb.fit[1,1]
    cover.isw2.35[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.35[i,4] = ind(1.040613, cover.isw2.35[i,1], cover.isw2.35[i,3])
    cover.isw2.35[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.35[i,6] = ismb.fit[2,1]
    cover.isw2.35[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.35[i,8] = ind(1.0174711, cover.isw2.35[i,5], cover.isw2.35[i,7])}
    , error=function(e){
      b0.isw2.35[i] = NA
      b0.isw2.sd.35[i] = NA
      b1.isw2.35[i] = NA
      b1.isw2.sd.35[i] = NA
      # Coverage
      cover.isw2.35[i,1] = NA
      cover.isw2.35[i,2] = NA
      cover.isw2.35[i,3] = NA
      cover.isw2.35[i,4] = NA
      cover.isw2.35[i,5] = NA
      cover.isw2.35[i,6] = NA
      cover.isw2.35[i,7] = NA
      cover.isw2.35[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw3.35[i] = ismb.fit[1,1]
    b0.isw3.sd.35[i] = ismb.fit[1,2]
    b1.isw3.35[i] = ismb.fit[2,1]
    b1.isw3.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.35[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.35[i,2] = ismb.fit[1,1]
    cover.isw3.35[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.35[i,4] = ind(1.040613, cover.isw3.35[i,1], cover.isw3.35[i,3])
    cover.isw3.35[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.35[i,6] = ismb.fit[2,1]
    cover.isw3.35[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.35[i,8] = ind(1.0174711, cover.isw3.35[i,5], cover.isw3.35[i,7])}
    , error=function(e){
      b0.isw3.35[i] = NA
      b0.isw3.sd.35[i] = NA
      b1.isw3.35[i] = NA
      b1.isw3.sd.35[i] = NA
      # Coverage
      cover.isw3.35[i,1] = NA
      cover.isw3.35[i,2] = NA
      cover.isw3.35[i,3] = NA
      cover.isw3.35[i,4] = NA
      cover.isw3.35[i,5] = NA
      cover.isw3.35[i,6] = NA
      cover.isw3.35[i,7] = NA
      cover.isw3.35[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw4.35[i] = ismb.fit[1,1]
    b0.isw4.sd.35[i] = ismb.fit[1,2]
    b1.isw4.35[i] = ismb.fit[2,1]
    b1.isw4.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.35[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.35[i,2] = ismb.fit[1,1]
    cover.isw4.35[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.35[i,4] = ind(1.040613, cover.isw4.35[i,1], cover.isw4.35[i,3])
    cover.isw4.35[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.35[i,6] = ismb.fit[2,1]
    cover.isw4.35[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.35[i,8] = ind(1.0174711, cover.isw4.35[i,5], cover.isw4.35[i,7])}
    , error=function(e){
      b0.isw4.35[i] = NA
      b0.isw4.sd.35[i] = NA
      b1.isw4.35[i] = NA
      b1.isw4.sd.35[i] = NA
      # Coverage
      cover.isw4.35[i,1] = NA
      cover.isw4.35[i,2] = NA
      cover.isw4.35[i,3] = NA
      cover.isw4.35[i,4] = NA
      cover.isw4.35[i,5] = NA
      cover.isw4.35[i,6] = NA
      cover.isw4.35[i,7] = NA
      cover.isw4.35[i,8] = NA
    })
}

# IS beta table
table3.isw1[4,1]<-mean(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,2]<-mean(b0.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,3]<-sd(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,4]<-mean(cover.isw1.35[,4],na.rm=TRUE)
table3.isw1[4,5]<-mean(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,6]<-mean(b1.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,7]<-sd(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,8]<-mean(cover.isw1.35[,8],na.rm=TRUE)

table3.isw2[4,1]<-mean(b0.isw2.35,na.rm=TRUE)
table3.isw2[4,2]<-mean(b0.isw2.sd.35,na.rm=TRUE)
table3.isw2[4,3]<-sd(b0.isw2.35,na.rm=TRUE)
table3.isw2[4,4]<-mean(cover.isw2.35[,4],na.rm=TRUE)
table3.isw2[4,5]<-mean(b1.isw2.35,na.rm=TRUE)
table3.isw2[4,6]<-mean(b1.isw2.sd.35,na.rm=TRUE)
table3.isw2[4,7]<-sd(b1.isw2.35,na.rm=TRUE)
table3.isw2[4,8]<-mean(cover.isw2.35[,8],na.rm=TRUE)

table3.isw3[4,1]<-mean(b0.isw3.35,na.rm=TRUE)
table3.isw3[4,2]<-mean(b0.isw3.sd.35,na.rm=TRUE)
table3.isw3[4,3]<-sd(b0.isw3.35,na.rm=TRUE)
table3.isw3[4,4]<-mean(cover.isw3.35[,4],na.rm=TRUE)
table3.isw3[4,5]<-mean(b1.isw3.35,na.rm=TRUE)
table3.isw3[4,6]<-mean(b1.isw3.sd.35,na.rm=TRUE)
table3.isw3[4,7]<-sd(b1.isw3.35,na.rm=TRUE)
table3.isw3[4,8]<-mean(cover.isw3.35[,8],na.rm=TRUE)

table3.isw4[4,1]<-mean(b0.isw4.35,na.rm=TRUE)
table3.isw4[4,2]<-mean(b0.isw4.sd.35,na.rm=TRUE)
table3.isw4[4,3]<-sd(b0.isw4.35,na.rm=TRUE)
table3.isw4[4,4]<-mean(cover.isw4.35[,4],na.rm=TRUE)
table3.isw4[4,5]<-mean(b1.isw4.35,na.rm=TRUE)
table3.isw4[4,6]<-mean(b1.isw4.sd.35,na.rm=TRUE)
table3.isw4[4,7]<-sd(b1.isw4.35,na.rm=TRUE)
table3.isw4[4,8]<-mean(cover.isw4.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.isw1.37 = c()
b0.isw1.sd.37 = c()
b1.isw1.37 = c()
b1.isw1.sd.37 = c()
cover.isw1.37=matrix(NA,2000,8)
colnames(cover.isw1.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw2.37 = c()
b0.isw2.sd.37 = c()
b1.isw2.37 = c()
b1.isw2.sd.37 = c()
cover.isw2.37=matrix(NA,2000,8)
colnames(cover.isw2.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw3.37 = c()
b0.isw3.sd.37 = c()
b1.isw3.37 = c()
b1.isw3.sd.37 = c()
cover.isw3.37=matrix(NA,2000,8)
colnames(cover.isw3.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.isw4.37 = c()
b0.isw4.sd.37 = c()
b1.isw4.37 = c()
b1.isw4.sd.37 = c()
cover.isw4.37=matrix(NA,2000,8)
colnames(cover.isw4.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.7)
  # Beta estimation : ismb with weight combination1
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.37[i] = ismb.fit[1,1]
    b0.isw1.sd.37[i] = ismb.fit[1,2]
    b1.isw1.37[i] = ismb.fit[2,1]
    b1.isw1.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.37[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.37[i,2] = ismb.fit[1,1]
    cover.isw1.37[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.37[i,4] = ind(1.040613, cover.isw1.37[i,1], cover.isw1.37[i,3])
    cover.isw1.37[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.37[i,6] = ismb.fit[2,1]
    cover.isw1.37[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.37[i,8] = ind(1.0174711, cover.isw1.37[i,5], cover.isw1.37[i,7])}
    , error=function(e){
      b0.isw1.37[i] = NA
      b0.isw1.sd.37[i] = NA
      b1.isw1.37[i] = NA
      b1.isw1.sd.37[i] = NA
      # Coverage
      cover.isw1.37[i,1] = NA
      cover.isw1.37[i,2] = NA
      cover.isw1.37[i,3] = NA
      cover.isw1.37[i,4] = NA
      cover.isw1.37[i,5] = NA
      cover.isw1.37[i,6] = NA
      cover.isw1.37[i,7] = NA
      cover.isw1.37[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination2
  tryCatch({
    ismb.fit = ismbw2.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw2.37[i] = ismb.fit[1,1]
    b0.isw2.sd.37[i] = ismb.fit[1,2]
    b1.isw2.37[i] = ismb.fit[2,1]
    b1.isw2.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.isw2.37[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw2.37[i,2] = ismb.fit[1,1]
    cover.isw2.37[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw2.37[i,4] = ind(1.040613, cover.isw2.37[i,1], cover.isw2.37[i,3])
    cover.isw2.37[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw2.37[i,6] = ismb.fit[2,1]
    cover.isw2.37[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw2.37[i,8] = ind(1.0174711, cover.isw2.37[i,5], cover.isw2.37[i,7])}
    , error=function(e){
      b0.isw2.37[i] = NA
      b0.isw2.sd.37[i] = NA
      b1.isw2.37[i] = NA
      b1.isw2.sd.37[i] = NA
      # Coverage
      cover.isw2.37[i,1] = NA
      cover.isw2.37[i,2] = NA
      cover.isw2.37[i,3] = NA
      cover.isw2.37[i,4] = NA
      cover.isw2.37[i,5] = NA
      cover.isw2.37[i,6] = NA
      cover.isw2.37[i,7] = NA
      cover.isw2.37[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination3
  tryCatch({
    ismb.fit = ismbw3.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw3.37[i] = ismb.fit[1,1]
    b0.isw3.sd.37[i] = ismb.fit[1,2]
    b1.isw3.37[i] = ismb.fit[2,1]
    b1.isw3.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.isw3.37[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw3.37[i,2] = ismb.fit[1,1]
    cover.isw3.37[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw3.37[i,4] = ind(1.040613, cover.isw3.37[i,1], cover.isw3.37[i,3])
    cover.isw3.37[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw3.37[i,6] = ismb.fit[2,1]
    cover.isw3.37[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw3.37[i,8] = ind(1.0174711, cover.isw3.37[i,5], cover.isw3.37[i,7])}
    , error=function(e){
      b0.isw3.37[i] = NA
      b0.isw3.sd.37[i] = NA
      b1.isw3.37[i] = NA
      b1.isw3.sd.37[i] = NA
      # Coverage
      cover.isw3.37[i,1] = NA
      cover.isw3.37[i,2] = NA
      cover.isw3.37[i,3] = NA
      cover.isw3.37[i,4] = NA
      cover.isw3.37[i,5] = NA
      cover.isw3.37[i,6] = NA
      cover.isw3.37[i,7] = NA
      cover.isw3.37[i,8] = NA
    })
  
  #  Beta estimation : ismb with weight combination4
  tryCatch({
    ismb.fit = ismbw4.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw4.37[i] = ismb.fit[1,1]
    b0.isw4.sd.37[i] = ismb.fit[1,2]
    b1.isw4.37[i] = ismb.fit[2,1]
    b1.isw4.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.isw4.37[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw4.37[i,2] = ismb.fit[1,1]
    cover.isw4.37[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw4.37[i,4] = ind(1.040613, cover.isw4.37[i,1], cover.isw4.37[i,3])
    cover.isw4.37[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw4.37[i,6] = ismb.fit[2,1]
    cover.isw4.37[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw4.37[i,8] = ind(1.0174711, cover.isw4.37[i,5], cover.isw4.37[i,7])}
    , error=function(e){
      b0.isw4.37[i] = NA
      b0.isw4.sd.37[i] = NA
      b1.isw4.37[i] = NA
      b1.isw4.sd.37[i] = NA
      # Coverage
      cover.isw4.37[i,1] = NA
      cover.isw4.37[i,2] = NA
      cover.isw4.37[i,3] = NA
      cover.isw4.37[i,4] = NA
      cover.isw4.37[i,5] = NA
      cover.isw4.37[i,6] = NA
      cover.isw4.37[i,7] = NA
      cover.isw4.37[i,8] = NA
    })
}

# IS beta table
table3.isw1[5,1]<-mean(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,2]<-mean(b0.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,3]<-sd(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,4]<-mean(cover.isw1.37[,4],na.rm=TRUE)
table3.isw1[5,5]<-mean(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,6]<-mean(b1.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,7]<-sd(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,8]<-mean(cover.isw1.37[,8],na.rm=TRUE)

table3.isw2[5,1]<-mean(b0.isw2.37,na.rm=TRUE)
table3.isw2[5,2]<-mean(b0.isw2.sd.37,na.rm=TRUE)
table3.isw2[5,3]<-sd(b0.isw2.37,na.rm=TRUE)
table3.isw2[5,4]<-mean(cover.isw2.37[,4],na.rm=TRUE)
table3.isw2[5,5]<-mean(b1.isw2.37,na.rm=TRUE)
table3.isw2[5,6]<-mean(b1.isw2.sd.37,na.rm=TRUE)
table3.isw2[5,7]<-sd(b1.isw2.37,na.rm=TRUE)
table3.isw2[5,8]<-mean(cover.isw2.37[,8],na.rm=TRUE)

table3.isw3[5,1]<-mean(b0.isw3.37,na.rm=TRUE)
table3.isw3[5,2]<-mean(b0.isw3.sd.37,na.rm=TRUE)
table3.isw3[5,3]<-sd(b0.isw3.37,na.rm=TRUE)
table3.isw3[5,4]<-mean(cover.isw3.37[,4],na.rm=TRUE)
table3.isw3[5,5]<-mean(b1.isw3.37,na.rm=TRUE)
table3.isw3[5,6]<-mean(b1.isw3.sd.37,na.rm=TRUE)
table3.isw3[5,7]<-sd(b1.isw3.37,na.rm=TRUE)
table3.isw3[5,8]<-mean(cover.isw3.37[,8],na.rm=TRUE)

table3.isw4[5,1]<-mean(b0.isw4.37,na.rm=TRUE)
table3.isw4[5,2]<-mean(b0.isw4.sd.37,na.rm=TRUE)
table3.isw4[5,3]<-sd(b0.isw4.37,na.rm=TRUE)
table3.isw4[5,4]<-mean(cover.isw4.37[,4],na.rm=TRUE)
table3.isw4[5,5]<-mean(b1.isw4.37,na.rm=TRUE)
table3.isw4[5,6]<-mean(b1.isw4.sd.37,na.rm=TRUE)
table3.isw4[5,7]<-sd(b1.isw4.37,na.rm=TRUE)
table3.isw4[5,8]<-mean(cover.isw4.37[,8],na.rm=TRUE)