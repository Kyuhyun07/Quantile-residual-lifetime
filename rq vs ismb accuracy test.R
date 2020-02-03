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
  return(data)
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
#### Estimating function (Rq with weight and ISMB)#### 
# 1. Rq with weight
rq.est = function(logT, nc, cov, Q, W){
  rq.fit = rq(logT ~ cov, tau=Q, weight = W, method="pfn")
  a=summary.rq(rq.fit, se="boot")
  beta.sd = matrix(NA, nc+1, 2)
  for (j in 1:(nc+1)){
    beta.sd[j,]=a[3][[1]][j,c(1,2)]
  }
  print(beta.sd)
}

# 2. ISMB
ismb.est = function(n, logT, nc ,cov, Q, W, ne){
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
  m = length(logT)
  X = as.matrix(cbind(c(rep(1,m)),cov))
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

#### result table ####
result_table<-matrix(NA,2,8)
rownames(result_table)<-c("rq","ISMB")
colnames(result_table)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

#### t_0=2 & c=30% ####
b0.rq.00<-c()
b0.rq.sd.00<-c()
b1.rq.00<-c()
b1.rq.sd.00<-c()
b0.is.00<-c()
b0.is.sd.00<-c()
b1.is.00<-c()
b1.is.sd.00<-c()
cover.rq.00=matrix(NA,2000,8)
colnames(cover.rq.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
cover.is.00=matrix(NA,2000,8)
colnames(cover.is.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Rq package
  a<-data.gen(200,c.3,2)
  nc = 1
  a.weight=weight.cal(a[,3],nc,a[,6],a[,7],2)
  cov = as.matrix(a.weight[,3:(nc+2)])
  tryCatch({
    rq.fit = rq.est(a.weight[,2], 1, cov, 0.5, a.weight[,9])
    b0.rq.00[i] = rq.fit[1,1]
    b0.rq.sd.00[i] = rq.fit[1,2]
    b1.rq.00[i] = rq.fit[2,1]
    b1.rq.sd.00[i] = rq.fit[2,2]
    cover.rq.00[i,1]=rq.fit[1,1]-1.96*rq.fit[1,2]
    cover.rq.00[i,2]=rq.fit[1,1]
    cover.rq.00[i,3]=rq.fit[1,1]+1.96*rq.fit[1,2]
    cover.rq.00[i,4]=ind(1.609438,rq.fit[1,1]-1.96*rq.fit[1,2],rq.fit[1,1]+1.96*rq.fit[1,2])
    cover.rq.00[i,5]=rq.fit[2,1]-1.96*rq.fit[2,2]
    cover.rq.00[i,6]=rq.fit[2,1]
    cover.rq.00[i,7]=rq.fit[2,1]+1.96*rq.fit[2,2]
    cover.rq.00[i,8]=ind(0.6931472,rq.fit[2,1]-1.96*rq.fit[2,2],rq.fit[2,1]+1.96*rq.fit[2,2])}
    ,error=function(e){
      b0.rq.00[i] = NA
      b0.rq.sd.00[i] = NA
      b1.rq.00[i] = NA
      b1.rq.sd.00[i] = NA
      cover.rq.00[i,1]=NA
      cover.rq.00[i,2]=NA
      cover.rq.00[i,3]=NA
      cover.rq.00[i,4]=NA
      cover.rq.00[i,5]=NA
      cover.rq.00[i,6]=NA
      cover.rq.00[i,7]=NA
      cover.rq.00[i,8]=NA
    })
  # Beta estimation method 2: ISMB
  tryCatch({
    ismb.fit = ismb.est(200, a.weight[,2], 1, cov, 0.5, a.weight[,9], 100)
    b0.is.00[i] = ismb.fit[1,1]
    b0.is.sd.00[i] = ismb.fit[1,2]
    b1.is.00[i] = ismb.fit[2,1]
    b1.is.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.is.00[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.is.00[i,2]=ismb.fit[1,1]
    cover.is.00[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.is.00[i,4]=ind(1.609438,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.is.00[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.is.00[i,6]=ismb.fit[2,1]
    cover.is.00[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.is.00[i,8]=ind(0.6931472,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.00[i] = NA
      b0.is.sd.00[i] = NA
      b1.is.00[i] = NA
      b1.is.sd.00[i] = NA
      # Coverage
      cover.is.00[i,1]=NA
      cover.is.00[i,2]=NA
      cover.is.00[i,3]=NA
      cover.is.00[i,4]=NA
      cover.is.00[i,5]=NA
      cover.is.00[i,6]=NA
      cover.is.00[i,7]=NA
      cover.is.00[i,8]=NA
    })
}
# Crq beta table
result_table[1,1]<-mean(b0.rq.00,na.rm=TRUE)
result_table[1,2]<-mean(b0.rq.sd.00,na.rm=TRUE)
result_table[1,3]<-sd(b0.rq.00,na.rm=TRUE)
result_table[1,4]<-mean(cover.rq.00[,4],na.rm=TRUE)
result_table[1,5]<-mean(b1.rq.00,na.rm=TRUE)
result_table[1,6]<-mean(b1.rq.sd.00,na.rm=TRUE)
result_table[1,7]<-sd(b1.rq.00,na.rm=TRUE)
result_table[1,8]<-mean(cover.rq.00[,8],na.rm=TRUE)

# IS beta table
result_table[2,1]<-mean(b0.is.00,na.rm=TRUE)
result_table[2,2]<-mean(b0.is.sd.00,na.rm=TRUE)
result_table[2,3]<-sd(b0.is.00,na.rm=TRUE)
result_table[2,4]<-mean(cover.is.00[,4],na.rm=TRUE)
result_table[2,5]<-mean(b1.is.00,na.rm=TRUE)
result_table[2,6]<-mean(b1.is.sd.00,na.rm=TRUE)
result_table[2,7]<-sd(b1.is.00,na.rm=TRUE)
result_table[2,8]<-mean(cover.is.00[,8],na.rm=TRUE)