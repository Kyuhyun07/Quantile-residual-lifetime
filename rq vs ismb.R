library(quantreg)
library(survival)
library(nleqslv)
library(xtable)
library(tictoc)

#### Data Generation function ####
data.gen<-function(samplesize, censor, t_0){
  unif = runif(n=samplesize ,min = 0,max = 1)
  sim=matrix(NA,samplesize,7)
  # Generate C_i
  sim[,2] = runif(samplesize,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,6] = rbinom(samplesize,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
  for (q in 1:samplesize){
    if (sim[q,6]==0){
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
    } else {
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.1
    }
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
  # Column names
  colnames(sim) = c("T","C","Z","Z.diff","log(Z.diff)","X","censored")
  return(sim)
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
#### Given information (t_0,no censoring) ####
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/5))^(1/k)/exp.beta.initial.1
c.0=5000000
c.1=78.11
c.3=26.36
c.5=15.08
c.7=9.09

a<-data.gen(200,c.0,0)

n=200
Y=a[,3]
nc=1
covariate=a[,6]
D=a[,7]
t_0=0
Q=0.5
ne=100

#### Weight calculation ####
weight.cal = function(Y, nc, covariate, D, t_0){
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
  print(data)
}


#### Sample data after weight calculation ####
a.weight = weight.cal(a[,3], 1, a[,6], a[,7], 0)
#### Define covariate ####
cov = as.matrix(a.weight[,3:(nc+2)])

#### Rq method ####
# method = pfn, se=boot, 0.013 sec 
tic()
rq.fit = rq(a.weight[,2] ~ cov, tau=0.5, weight = a.weight[,9], method="pfn")
summary.rq(rq.fit, se="boot")
toc()

#### Crq method ####
# preprocessing data
a.weight[,1]=pmin(a.weight[,1],a.weight[,2])
crq.fit = crq(Curv(a[,1],a[,2],ctype="right") ~ cov, tau=0.5, method="Powell")
result = summary(crq.fit,taus = c(0.1,0.5))

#### ISMB method ####
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

tic()
# Covariate setting (1 covariate)
m = nrow(a.weight)
X = as.matrix(cbind(c(rep(1,m)),a.weight[,3:(nc+2)]))
W = a.weight[,(nc+8)]
logT = a.weight[,2]
G = diag(1/n, nc+1, nc+1)
  
betastart = c(1,rep(1,nc))

is.fit = nleqslv(betastart,objectF, control=list(ftol=1e-5))
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
toc()


