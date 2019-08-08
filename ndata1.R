library(quantreg)
library(survival)
library(nleqslv)
library(xtable)

#### Condition ####
# data size = 200
# beta0, beta1 effective
# Quantile 25%
# simulation 2000
# eta = 100
#### True Beta ####
# When T_0=0, beta_0=1.61, beta_1=1.61
# When T_0=1, beta_0=1.41, beta_1=1.77
# When T_0=2, beta_0=1.22, beta_1=1.92
# When T_0=3, beta_0=1.04, beta_1=2.06

#### Find censoring point ####
#d.exp.beta.initial.0=5
#d.exp.beta.initial.1=25
#d.k=2
#d.r.initial.0=(log(4/3)^(1/d.k))/d.exp.beta.initial.0
#d.r.initial.1=(log(4/3)^(1/d.k))/d.exp.beta.initial.1
#d.u=runif(n=1000000,min = 0,max = 1)
#d.x=rbinom(1000000,size=1,p=0.5)
#d.T<-c()
#for (q in 1:1000000){
#  if (d.x[q]==0){
#    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.0
#  } else {
#    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.1
#  }
#}
##Find c which is derermined to achieve 0,10,20,30% censoring rate in case 1
#i=123.66
#while(1){
#  d.C<-runif(1000000,0,i)
#  i=i+0.01
#  if(sum(d.C<d.T)<=100000) break # 0% = 0, 10% = 100000, 20% = 200000, 30% = 300000)
#}
#print(i)
#table(d.C<d.T)

#### Data Generation function ####
data.gen<-function(censor, t_0){
  unif = runif(n=200,min = 0,max = 1)
  sim=matrix(NA,200,7)
  # Generate C_i
  sim[,2] = runif(200,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,6] = rbinom(200,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=25))
  for (q in 1:200){
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
  ## Calcuation W_i
  for (w in 1:n){
    sim[w,8]<-sum(sim[,3]>=sim[w,3])
    sim[w,9]<-sum(sim[,3]==sim[w,3])
    if (sim[w,7]==0){
      sim[w,10]=(sim[w,8]-sim[w,9])/sim[w,8]
    } else {
      sim[w,10]=1
    }
    if (w==1){
      sim[w,11]<-sim[w,10]
    } else {
      sim[w,11]<-sim[w-1,11]*sim[w,10]
    }
    if (sim[w,11]==0){
      if (sim[w,7]==0){
        sim[w,12]=0
      } else {
        sim[w,12]<-sim[w-1,12]}
    }
    else {
      sim[w,12]<-sim[w,7]/(sim[w,11])}
  }
  # Column names
  colnames(sim) = c("T","C","Z","Z.diff","log(Z.diff)","X","censored","# at risk","# event","s/d","G_KM","Weight")
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
  G = diag(1/m, nc+1, nc+1)
  
  betastart = c(1,rep(1,nc))
  is.fit = nleqslv(betastart,objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
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


#### Given Information ####
c.0=5000000
c.1=246.31
c.2=123.66
c.3=81.63
c.5=27.37
c.7=11.8

## My assumption
exp.beta.initial.0=5
exp.beta.initial.1=25
k=2
r.initial.0=(log(4/3))^(1/k)/exp.beta.initial.0
r.initial.1=(log(4/3))^(1/k)/exp.beta.initial.1

#### Make table for Beta estimation and variance estimation and Coverage ####
# t_0=0
table0.crq<-matrix(NA,5,4)
rownames(table0.crq)<-c(0,10,30,50,70)
colnames(table0.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1.crq<-matrix(NA,5,4)
rownames(table1.crq)<-c(0,10,30,50,70)
colnames(table1.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table2.crq<-matrix(NA,5,4)
rownames(table2.crq)<-c(0,10,30,50,70)
colnames(table2.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3.crq<-matrix(NA,5,4)
rownames(table3.crq)<-c(0,10,30,50,70)
colnames(table3.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

table0.is<-matrix(NA,5,4)
rownames(table0.is)<-c(0,10,30,50,70)
colnames(table0.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1.is<-matrix(NA,5,4)
rownames(table1.is)<-c(0,10,30,50,70)
colnames(table1.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table2.is<-matrix(NA,5,4)
rownames(table2.is)<-c(0,10,30,50,70)
colnames(table2.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3.is<-matrix(NA,5,4)
rownames(table3.is)<-c(0,10,30,50,70)
colnames(table3.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

b0var.table0=matrix(NA,5,3)
colnames(b0var.table0)=c("true","Crq","ISMB")
rownames(b0var.table0)=c(0,10,30,50,70)
b0var.table1=matrix(NA,5,3)
colnames(b0var.table1)=c("true","Crq","ISMB")
rownames(b0var.table1)=c(0,10,30,50,70)
b0var.table2=matrix(NA,5,3)
colnames(b0var.table2)=c("true","Crq","ISMB")
rownames(b0var.table2)=c(0,10,30,50,70)
b0var.table3=matrix(NA,5,3)
colnames(b0var.table3)=c("true","Crq","ISMB")
rownames(b0var.table3)=c(0,10,30,50,70)

b1var.table0=matrix(NA,5,3)
colnames(b1var.table0)=c("true","Crq","ISMB")
rownames(b1var.table0)=c(0,10,30,50,70)
b1var.table1=matrix(NA,5,3)
colnames(b1var.table1)=c("true","Crq","ISMB")
rownames(b1var.table1)=c(0,10,30,50,70)
b1var.table2=matrix(NA,5,3)
colnames(b1var.table2)=c("true","Crq","ISMB")
rownames(b1var.table2)=c(0,10,30,50,70)
b1var.table3=matrix(NA,5,3)
colnames(b1var.table3)=c("true","Crq","ISMB")
rownames(b1var.table3)=c(0,10,30,50,70)

b0.coverage=matrix(NA,4,5)
colnames(b0.coverage)=c("0%","10%","30%","50%","70%")
rownames(b0.coverage)=c("t0=0","t0=1","t0=2","t0=3")

b1.coverage=matrix(NA,4,5)
colnames(b1.coverage)=c("0%","10%","30%","50%","70%")
rownames(b1.coverage)=c("t0=0","t0=1","t0=2","t0=3")

#### t_0=0 & c=0% ####
b0.crq.00<-c()
b0.crq.sd.00<-c()
b1.crq.00<-c()
b1.crq.sd.00<-c()
b0.is.00<-c()
b0.is.sd.00<-c()
b1.is.00<-c()
b1.is.sd.00<-c()
cover.00=matrix(NA,2000,8)
colnames(cover.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,0)
  a[,2]<-100
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.25)
    b0.crq.00[i] = crq.fit[1,1]
    b0.crq.sd.00[i] = crq.fit[1,2]
    b1.crq.00[i] = crq.fit[2,1]
    b1.crq.sd.00[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.00[i] = NA
      b0.crq.sd.00[i] = NA
      b1.crq.00[i] = NA
      b1.crq.sd.00[i] = NA
    })
}
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 2: ISMB
  a<-data.gen(c.0,0)
  a[,2]<-100
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.25, 100)
    b0.is.00[i] = ismb.fit[1,1]
    b0.is.sd.00[i] = ismb.fit[1,2]
    b1.is.00[i] = ismb.fit[2,1]
    b1.is.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.00[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.00[i,2]=ismb.fit[1,1]
    cover.00[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.00[i,4]=ind(1.61,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.00[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.00[i,6]=ismb.fit[2,1]
    cover.00[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.00[i,8]=ind(1.61,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.00[i] = NA
      b0.is.sd.00[i] = NA
      b1.is.00[i] = NA
      b1.is.sd.00[i] = NA
      # Coverage
      cover.00[i,1]=NA
      cover.00[i,2]=NA
      cover.00[i,3]=NA
      cover.00[i,4]=NA
      cover.00[i,5]=NA
      cover.00[i,6]=NA
      cover.00[i,7]=NA
      cover.00[i,8]=NA
    })
}
# Crq beta table
table0.crq[1,1]<-mean(b0.crq.00,na.rm=TRUE)
table0.crq[1,2]<-sd(b0.crq.00,na.rm=TRUE)
table0.crq[1,3]<-mean(b1.crq.00,na.rm=TRUE)
table0.crq[1,4]<-sd(b1.crq.00,na.rm=TRUE)

# IS beta table
table0.is[1,1]<-mean(b0.is.00,na.rm=TRUE)
table0.is[1,2]<-sd(b0.is.00,na.rm=TRUE)
table0.is[1,3]<-mean(b1.is.00,na.rm=TRUE)
table0.is[1,4]<-sd(b1.is.00,na.rm=TRUE)

# Variance table (Beta0)
b0var.table0[1,1]=sd(b0.is.00,na.rm=TRUE)
b0var.table0[1,2]=mean(b0.crq.sd.00,na.rm=TRUE)
b0var.table0[1,3]=mean(b0.is.sd.00,na.rm=TRUE)

# Variance table (Beta1)
b1var.table0[1,1]=sd(b1.is.00,na.rm=TRUE)
b1var.table0[1,2]=mean(b1.crq.sd.00,na.rm=TRUE)
b1var.table0[1,3]=mean(b1.is.sd.00,na.rm=TRUE)

# Coverage table
b0.coverage[1,1]=mean(cover.00[,4],na.rm=TRUE)
b1.coverage[1,1]=mean(cover.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.crq.01<-c()
b0.crq.sd.01<-c()
b1.crq.01<-c()
b1.crq.sd.01<-c()
b0.is.01<-c()
b0.is.sd.01<-c()
b1.is.01<-c()
b1.is.sd.01<-c()
cover.01=matrix(NA,2000,8)
colnames(cover.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,0)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.25)
    b0.crq.01[i] = crq.fit[1,1]
    b0.crq.sd.01[i] = crq.fit[1,2]
    b1.crq.01[i] = crq.fit[2,1]
    b1.crq.sd.01[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.01[i] = NA
      b0.crq.sd.01[i] = NA
      b1.crq.01[i] = NA
      b1.crq.sd.01[i] = NA
    })
}
set.seed(1)
for (i in 1:2000){
  a<-data.gen(c.1,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.25, 100)
    b0.is.01[i] = ismb.fit[1,1]
    b0.is.sd.01[i] = ismb.fit[1,2]
    b1.is.01[i] = ismb.fit[2,1]
    b1.is.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.01[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.01[i,2]=ismb.fit[1,1]
    cover.01[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.01[i,4]=ind(1.61,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.01[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.01[i,6]=ismb.fit[2,1]
    cover.01[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.01[i,8]=ind(1.61,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.01[i] = NA
      b0.is.sd.01[i] = NA
      b1.is.01[i] = NA
      b1.is.sd.01[i] = NA
      # Coverage
      cover.01[i,1]=NA
      cover.01[i,2]=NA
      cover.01[i,3]=NA
      cover.01[i,4]=NA
      cover.01[i,5]=NA
      cover.01[i,6]=NA
      cover.01[i,7]=NA
      cover.01[i,8]=NA
    })
}

# Crq beta table
table0.crq[2,1]<-mean(b0.crq.01,na.rm=TRUE)
table0.crq[2,2]<-sd(b0.crq.01,na.rm=TRUE)
table0.crq[2,3]<-mean(b1.crq.01,na.rm=TRUE)
table0.crq[2,4]<-sd(b1.crq.01,na.rm=TRUE)

# IS beta table
table0.is[2,1]<-mean(b0.is.01,na.rm=TRUE)
table0.is[2,2]<-sd(b0.is.01,na.rm=TRUE)
table0.is[2,3]<-mean(b1.is.01,na.rm=TRUE)
table0.is[2,4]<-sd(b1.is.01,na.rm=TRUE)

# Variance table (Beta0)
b0var.table0[2,1]=sd(b0.is.01,na.rm=TRUE)
b0var.table0[2,2]=mean(b0.crq.sd.01,na.rm=TRUE)
b0var.table0[2,3]=mean(b0.is.sd.01,na.rm=TRUE)

# Variance table (Beta1)
b1var.table0[2,1]=sd(b1.is.01,na.rm=TRUE)
b1var.table0[2,2]=mean(b1.crq.sd.01,na.rm=TRUE)
b1var.table0[2,3]=mean(b1.is.sd.01,na.rm=TRUE)

# Coverage table
b0.coverage[1,2]=mean(cover.01[,4],na.rm=TRUE)
b1.coverage[1,2]=mean(cover.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.crq.03<-c()
b0.crq.sd.03<-c()
b1.crq.03<-c()
b1.crq.sd.03<-c()
b0.is.03<-c()
b0.is.sd.03<-c()
b1.is.03<-c()
b1.is.sd.03<-c()
cover.03=matrix(NA,2000,8)
colnames(cover.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,0)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.25)
    b0.crq.03[i] = crq.fit[1,1]
    b0.crq.sd.03[i] = crq.fit[1,2]
    b1.crq.03[i] = crq.fit[2,1]
    b1.crq.sd.03[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.03[i] = NA
      b0.crq.sd.03[i] = NA
      b1.crq.03[i] = NA
      b1.crq.sd.03[i] = NA
    })
}
set.seed(1)
for (i in 1:2000){
  a<-data.gen(c.3,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.25, 100)
    b0.is.03[i] = ismb.fit[1,1]
    b0.is.sd.03[i] = ismb.fit[1,2]
    b1.is.03[i] = ismb.fit[2,1]
    b1.is.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.03[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.03[i,2]=ismb.fit[1,1]
    cover.03[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.03[i,4]=ind(1.61,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.03[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.03[i,6]=ismb.fit[2,1]
    cover.03[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.03[i,8]=ind(1.61,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.03[i] = NA
      b0.is.sd.03[i] = NA
      b1.is.03[i] = NA
      b1.is.sd.03[i] = NA
      # Coverage
      cover.03[i,1]=NA
      cover.03[i,2]=NA
      cover.03[i,3]=NA
      cover.03[i,4]=NA
      cover.03[i,5]=NA
      cover.03[i,6]=NA
      cover.03[i,7]=NA
      cover.03[i,8]=NA
    })
}
# Crq beta table
table0.crq[3,1]<-mean(b0.crq.03,na.rm=TRUE)
table0.crq[3,2]<-sd(b0.crq.03,na.rm=TRUE)
table0.crq[3,3]<-mean(b1.crq.03,na.rm=TRUE)
table0.crq[3,4]<-sd(b1.crq.03,na.rm=TRUE)

# IS beta table
table0.is[3,1]<-mean(b0.is.03,na.rm=TRUE)
table0.is[3,2]<-sd(b0.is.03,na.rm=TRUE)
table0.is[3,3]<-mean(b1.is.03,na.rm=TRUE)
table0.is[3,4]<-sd(b1.is.03,na.rm=TRUE)

# Variance table (Beta0)
b0var.table0[3,1]=sd(b0.is.03,na.rm=TRUE)
b0var.table0[3,2]=mean(b0.crq.sd.03,na.rm=TRUE)
b0var.table0[3,3]=mean(b0.is.sd.03,na.rm=TRUE)

# Variance table (Beta1)
b1var.table0[3,1]=sd(b1.is.03,na.rm=TRUE)
b1var.table0[3,2]=mean(b1.crq.sd.03,na.rm=TRUE)
b1var.table0[3,3]=mean(b1.is.sd.03,na.rm=TRUE)

# Coverage table
b0.coverage[1,3]=mean(cover.03[,4],na.rm=TRUE)
b1.coverage[1,3]=mean(cover.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.crq.05<-c()
b0.crq.sd.05<-c()
b1.crq.05<-c()
b1.crq.sd.05<-c()
b0.is.05<-c()
b0.is.sd.05<-c()
b1.is.05<-c()
b1.is.sd.05<-c()
cover.05=matrix(NA,2000,8)
colnames(cover.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.5,0)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.25)
    b0.crq.05[i] = crq.fit[1,1]
    b0.crq.sd.05[i] = crq.fit[1,2]
    b1.crq.05[i] = crq.fit[2,1]
    b1.crq.sd.05[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.05[i] = NA
      b0.crq.sd.05[i] = NA
      b1.crq.05[i] = NA
      b1.crq.sd.05[i] = NA
    })
}
set.seed(1)
for (i in 1:2000){
  a<-data.gen(c.5,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.25, 100)
    b0.is.05[i] = ismb.fit[1,1]
    b0.is.sd.05[i] = ismb.fit[1,2]
    b1.is.05[i] = ismb.fit[2,1]
    b1.is.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.05[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.05[i,2]=ismb.fit[1,1]
    cover.05[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.05[i,4]=ind(1.61,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.05[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.05[i,6]=ismb.fit[2,1]
    cover.05[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.05[i,8]=ind(1.61,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.05[i] = NA
      b0.is.sd.05[i] = NA
      b1.is.05[i] = NA
      b1.is.sd.05[i] = NA
      # Coverage
      cover.05[i,1]=NA
      cover.05[i,2]=NA
      cover.05[i,3]=NA
      cover.05[i,4]=NA
      cover.05[i,5]=NA
      cover.05[i,6]=NA
      cover.05[i,7]=NA
      cover.05[i,8]=NA
    })
}

# Crq beta table
table0.crq[4,1]<-mean(b0.crq.05,na.rm=TRUE)
table0.crq[4,2]<-sd(b0.crq.05,na.rm=TRUE)
table0.crq[4,3]<-mean(b1.crq.05,na.rm=TRUE)
table0.crq[4,4]<-sd(b1.crq.05,na.rm=TRUE)

# IS beta table
table0.is[4,1]<-mean(b0.is.05,na.rm=TRUE)
table0.is[4,2]<-sd(b0.is.05,na.rm=TRUE)
table0.is[4,3]<-mean(b1.is.05,na.rm=TRUE)
table0.is[4,4]<-sd(b1.is.05,na.rm=TRUE)

# Variance table (Beta0)
b0var.table0[4,1]=sd(b0.is.05,na.rm=TRUE)
b0var.table0[4,2]=mean(b0.crq.sd.05,na.rm=TRUE)
b0var.table0[4,3]=mean(b0.is.sd.05,na.rm=TRUE)

# Variance table (Beta1)
b1var.table0[4,1]=sd(b1.is.05,na.rm=TRUE)
b1var.table0[4,2]=mean(b1.crq.sd.05,na.rm=TRUE)
b1var.table0[4,3]=mean(b1.is.sd.05,na.rm=TRUE)

# Coverage table
b0.coverage[1,4]=mean(cover.05[,4],na.rm=TRUE)
b1.coverage[1,4]=mean(cover.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
cover.07=matrix(NA,2000,8)
colnames(cover.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.7,0)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.25)
    b0.crq.07[i] = crq.fit[1,1]
    b0.crq.sd.07[i] = crq.fit[1,2]
    b1.crq.07[i] = crq.fit[2,1]
    b1.crq.sd.07[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.07[i] = NA
      b0.crq.sd.07[i] = NA
      b1.crq.07[i] = NA
      b1.crq.sd.07[i] = NA
    })
}
set.seed(1)
for (i in 1:2000){
  a<-data.gen(c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.25, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.07[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.07[i,2]=ismb.fit[1,1]
    cover.07[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.07[i,4]=ind(1.61,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.07[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.07[i,6]=ismb.fit[2,1]
    cover.07[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.07[i,8]=ind(1.61,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
      # Coverage
      cover.07[i,1]=NA
      cover.07[i,2]=NA
      cover.07[i,3]=NA
      cover.07[i,4]=NA
      cover.07[i,5]=NA
      cover.07[i,6]=NA
      cover.07[i,7]=NA
      cover.07[i,8]=NA
    })
}
# Crq beta table
table0.crq[5,1]<-mean(b0.crq.07,na.rm=TRUE)
table0.crq[5,2]<-sd(b0.crq.07,na.rm=TRUE)
table0.crq[5,3]<-mean(b1.crq.07,na.rm=TRUE)
table0.crq[5,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table0.is[5,1]<-mean(b0.is.07,na.rm=TRUE)
table0.is[5,2]<-sd(b0.is.07,na.rm=TRUE)
table0.is[5,3]<-mean(b1.is.07,na.rm=TRUE)
table0.is[5,4]<-sd(b1.is.07,na.rm=TRUE)

# Variance table (Beta0)
b0var.table0[5,1]=sd(b0.is.07,na.rm=TRUE)
b0var.table0[5,2]=mean(b0.crq.sd.07,na.rm=TRUE)
b0var.table0[5,3]=mean(b0.is.sd.07,na.rm=TRUE)

# Variance table (Beta1)
b1var.table0[5,1]=sd(b1.is.07,na.rm=TRUE)
b1var.table0[5,2]=mean(b1.crq.sd.07,na.rm=TRUE)
b1var.table0[5,3]=mean(b1.is.sd.07,na.rm=TRUE)

# Coverage table
b0.coverage[1,5]=mean(cover.07[,4],na.rm=TRUE)
b1.coverage[1,5]=mean(cover.07[,8],na.rm=TRUE)


#### t_0=1 & c=0% ####
b0.crq.10<-c()
b0.crq.sd.10<-c()
b1.crq.10<-c()
b1.crq.sd.10<-c()
b0.is.10<-c()
b0.is.sd.10<-c()
b1.is.10<-c()
b1.is.sd.10<-c()
cover.10=matrix(NA,2000,8)
colnames(cover.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(11)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.25)
    b0.crq.10[i] = crq.fit[1,1]
    b0.crq.sd.10[i] = crq.fit[1,2]
    b1.crq.10[i] = crq.fit[2,1]
    b1.crq.sd.10[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.10[i] = NA
      b0.crq.sd.10[i] = NA
      b1.crq.10[i] = NA
      b1.crq.sd.10[i] = NA
    })
}
set.seed(11)
for (i in 1:2000){
  a<-data.gen(c.0,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.25, 100)
    b0.is.10[i] = ismb.fit[1,1]
    b0.is.sd.10[i] = ismb.fit[1,2]
    b1.is.10[i] = ismb.fit[2,1]
    b1.is.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.10[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.10[i,2]=ismb.fit[1,1]
    cover.10[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.10[i,4]=ind(1.41,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.10[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.10[i,6]=ismb.fit[2,1]
    cover.10[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.10[i,8]=ind(1.77,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.10[i] = NA
      b0.is.sd.10[i] = NA
      b1.is.10[i] = NA
      b1.is.sd.10[i] = NA
      # Coverage
      cover.10[i,1]=NA
      cover.10[i,2]=NA
      cover.10[i,3]=NA
      cover.10[i,4]=NA
      cover.10[i,5]=NA
      cover.10[i,6]=NA
      cover.10[i,7]=NA
      cover.10[i,8]=NA
    })
}
# Crq beta table
table1.crq[1,1]<-mean(b0.crq.10,na.rm=TRUE)
table1.crq[1,2]<-sd(b0.crq.10,na.rm=TRUE)
table1.crq[1,3]<-mean(b1.crq.10,na.rm=TRUE)
table1.crq[1,4]<-sd(b1.crq.10,na.rm=TRUE)

# IS beta table
table1.is[1,1]<-mean(b0.is.10,na.rm=TRUE)
table1.is[1,2]<-sd(b0.is.10,na.rm=TRUE)
table1.is[1,3]<-mean(b1.is.10,na.rm=TRUE)
table1.is[1,4]<-sd(b1.is.10,na.rm=TRUE)

# Variance table (Beta0)
b0var.table1[1,1]=sd(b0.is.10,na.rm=TRUE)
b0var.table1[1,2]=mean(b0.crq.sd.10,na.rm=TRUE)
b0var.table1[1,3]=mean(b0.is.sd.10,na.rm=TRUE)

# Variance table (Beta1)
b1var.table1[1,1]=sd(b1.is.10,na.rm=TRUE)
b1var.table1[1,2]=mean(b1.crq.sd.10,na.rm=TRUE)
b1var.table1[1,3]=mean(b1.is.sd.10,na.rm=TRUE)

# Coverage table
b0.coverage[2,1]=mean(cover.10[,4],na.rm=TRUE)
b1.coverage[2,1]=mean(cover.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.crq.11<-c()
b0.crq.sd.11<-c()
b1.crq.11<-c()
b1.crq.sd.11<-c()
b0.is.11<-c()
b0.is.sd.11<-c()
b1.is.11<-c()
b1.is.sd.11<-c()
cover.11=matrix(NA,2000,8)
colnames(cover.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(11)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.25)
    b0.crq.11[i] = crq.fit[1,1]
    b0.crq.sd.11[i] = crq.fit[1,2]
    b1.crq.11[i] = crq.fit[2,1]
    b1.crq.sd.11[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.11[i] = NA
      b0.crq.sd.11[i] = NA
      b1.crq.11[i] = NA
      b1.crq.sd.11[i] = NA
    })
}
set.seed(11)
for (i in 1:2000){
  a<-data.gen(c.1,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.25, 100)
    b0.is.11[i] = ismb.fit[1,1]
    b0.is.sd.11[i] = ismb.fit[1,2]
    b1.is.11[i] = ismb.fit[2,1]
    b1.is.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.11[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.11[i,2]=ismb.fit[1,1]
    cover.11[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.11[i,4]=ind(1.41,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.11[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.11[i,6]=ismb.fit[2,1]
    cover.11[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.11[i,8]=ind(1.77,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.11[i] = NA
      b0.is.sd.11[i] = NA
      b1.is.11[i] = NA
      b1.is.sd.11[i] = NA
      # Coverage
      cover.11[i,1]=NA
      cover.11[i,2]=NA
      cover.11[i,3]=NA
      cover.11[i,4]=NA
      cover.11[i,5]=NA
      cover.11[i,6]=NA
      cover.11[i,7]=NA
      cover.11[i,8]=NA
    })
}
# Crq beta table
table1.crq[2,1]<-mean(b0.crq.11,na.rm=TRUE)
table1.crq[2,2]<-sd(b0.crq.11,na.rm=TRUE)
table1.crq[2,3]<-mean(b1.crq.11,na.rm=TRUE)
table1.crq[2,4]<-sd(b1.crq.11,na.rm=TRUE)

# IS beta table
table1.is[2,1]<-mean(b0.is.11,na.rm=TRUE)
table1.is[2,2]<-sd(b0.is.11,na.rm=TRUE)
table1.is[2,3]<-mean(b1.is.11,na.rm=TRUE)
table1.is[2,4]<-sd(b1.is.11,na.rm=TRUE)

# Variance table (Beta0)
b0var.table1[2,1]=sd(b0.is.11,na.rm=TRUE)
b0var.table1[2,2]=mean(b0.crq.sd.11,na.rm=TRUE)
b0var.table1[2,3]=mean(b0.is.sd.11,na.rm=TRUE)

# Variance table (Beta1)
b1var.table1[2,1]=sd(b1.is.11,na.rm=TRUE)
b1var.table1[2,2]=mean(b1.crq.sd.11,na.rm=TRUE)
b1var.table1[2,3]=mean(b1.is.sd.11,na.rm=TRUE)

# Coverage table
b0.coverage[2,2]=mean(cover.11[,4],na.rm=TRUE)
b1.coverage[2,2]=mean(cover.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.crq.13<-c()
b0.crq.sd.13<-c()
b1.crq.13<-c()
b1.crq.sd.13<-c()
b0.is.13<-c()
b0.is.sd.13<-c()
b1.is.13<-c()
b1.is.sd.13<-c()
cover.13=matrix(NA,2000,8)
colnames(cover.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(11)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.25)
    b0.crq.13[i] = crq.fit[1,1]
    b0.crq.sd.13[i] = crq.fit[1,2]
    b1.crq.13[i] = crq.fit[2,1]
    b1.crq.sd.13[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.13[i] = NA
      b0.crq.sd.13[i] = NA
      b1.crq.13[i] = NA
      b1.crq.sd.13[i] = NA
    })
}
set.seed(11)
for (i in 1:2000){
  a<-data.gen(c.3,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.25, 100)
    b0.is.13[i] = ismb.fit[1,1]
    b0.is.sd.13[i] = ismb.fit[1,2]
    b1.is.13[i] = ismb.fit[2,1]
    b1.is.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.13[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.13[i,2]=ismb.fit[1,1]
    cover.13[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.13[i,4]=ind(1.41,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.13[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.13[i,6]=ismb.fit[2,1]
    cover.13[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.13[i,8]=ind(1.77,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.13[i] = NA
      b0.is.sd.13[i] = NA
      b1.is.13[i] = NA
      b1.is.sd.13[i] = NA
      # Coverage
      cover.13[i,1]=NA
      cover.13[i,2]=NA
      cover.13[i,3]=NA
      cover.13[i,4]=NA
      cover.13[i,5]=NA
      cover.13[i,6]=NA
      cover.13[i,7]=NA
      cover.13[i,8]=NA
    })
}

# Crq beta table
table1.crq[3,1]<-mean(b0.crq.13,na.rm=TRUE)
table1.crq[3,2]<-sd(b0.crq.13,na.rm=TRUE)
table1.crq[3,3]<-mean(b1.crq.13,na.rm=TRUE)
table1.crq[3,4]<-sd(b1.crq.13,na.rm=TRUE)

# IS beta table
table1.is[3,1]<-mean(b0.is.13,na.rm=TRUE)
table1.is[3,2]<-sd(b0.is.13,na.rm=TRUE)
table1.is[3,3]<-mean(b1.is.13,na.rm=TRUE)
table1.is[3,4]<-sd(b1.is.13,na.rm=TRUE)

# Variance table (Beta0)
b0var.table1[3,1]=sd(b0.is.13,na.rm=TRUE)
b0var.table1[3,2]=mean(b0.crq.sd.13,na.rm=TRUE)
b0var.table1[3,3]=mean(b0.is.sd.13,na.rm=TRUE)

# Variance table (Beta1)
b1var.table1[3,1]=sd(b1.is.13,na.rm=TRUE)
b1var.table1[3,2]=mean(b1.crq.sd.13,na.rm=TRUE)
b1var.table1[3,3]=mean(b1.is.sd.13,na.rm=TRUE)

# Coverage table
b0.coverage[2,3]=mean(cover.13[,4],na.rm=TRUE)
b1.coverage[2,3]=mean(cover.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.crq.15<-c()
b0.crq.sd.15<-c()
b1.crq.15<-c()
b1.crq.sd.15<-c()
b0.is.15<-c()
b0.is.sd.15<-c()
b1.is.15<-c()
b1.is.sd.15<-c()
cover.15=matrix(NA,2000,8)
colnames(cover.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(11)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.5,1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.25)
    b0.crq.15[i] = crq.fit[1,1]
    b0.crq.sd.15[i] = crq.fit[1,2]
    b1.crq.15[i] = crq.fit[2,1]
    b1.crq.sd.15[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.15[i] = NA
      b0.crq.sd.15[i] = NA
      b1.crq.15[i] = NA
      b1.crq.sd.15[i] = NA
    })
}
set.seed(11)
for (i in 1:2000){
  a<-data.gen(c.5,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.25, 100)
    b0.is.15[i] = ismb.fit[1,1]
    b0.is.sd.15[i] = ismb.fit[1,2]
    b1.is.15[i] = ismb.fit[2,1]
    b1.is.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.15[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.15[i,2]=ismb.fit[1,1]
    cover.15[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.15[i,4]=ind(1.41,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.15[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.15[i,6]=ismb.fit[2,1]
    cover.15[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.15[i,8]=ind(1.77,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.15[i] = NA
      b0.is.sd.15[i] = NA
      b1.is.15[i] = NA
      b1.is.sd.15[i] = NA
      # Coverage
      cover.15[i,1]=NA
      cover.15[i,2]=NA
      cover.15[i,3]=NA
      cover.15[i,4]=NA
      cover.15[i,5]=NA
      cover.15[i,6]=NA
      cover.15[i,7]=NA
      cover.15[i,8]=NA
    })
}

# Crq beta table
table1.crq[4,1]<-mean(b0.crq.15,na.rm=TRUE)
table1.crq[4,2]<-sd(b0.crq.15,na.rm=TRUE)
table1.crq[4,3]<-mean(b1.crq.15,na.rm=TRUE)
table1.crq[4,4]<-sd(b1.crq.15,na.rm=TRUE)

# IS beta table
table1.is[4,1]<-mean(b0.is.15,na.rm=TRUE)
table1.is[4,2]<-sd(b0.is.15,na.rm=TRUE)
table1.is[4,3]<-mean(b1.is.15,na.rm=TRUE)
table1.is[4,4]<-sd(b1.is.15,na.rm=TRUE)

# Variance table (Beta0)
b0var.table1[4,1]=sd(b0.is.15,na.rm=TRUE)
b0var.table1[4,2]=mean(b0.crq.sd.15,na.rm=TRUE)
b0var.table1[4,3]=mean(b0.is.sd.15,na.rm=TRUE)

# Variance table (Beta1)
b1var.table1[4,1]=sd(b1.is.15,na.rm=TRUE)
b1var.table1[4,2]=mean(b1.crq.sd.15,na.rm=TRUE)
b1var.table1[4,3]=mean(b1.is.sd.15,na.rm=TRUE)

# Coverage table
b0.coverage[2,4]=mean(cover.15[,4],na.rm=TRUE)
b1.coverage[2,4]=mean(cover.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
cover.17=matrix(NA,2000,8)
colnames(cover.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(11)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.7,1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.25)
    b0.crq.17[i] = crq.fit[1,1]
    b0.crq.sd.17[i] = crq.fit[1,2]
    b1.crq.17[i] = crq.fit[2,1]
    b1.crq.sd.17[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.17[i] = NA
      b0.crq.sd.17[i] = NA
      b1.crq.17[i] = NA
      b1.crq.sd.17[i] = NA
    })
}
set.seed(11)
for (i in 1:2000){
  a<-data.gen(c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.25, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.17[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.17[i,2]=ismb.fit[1,1]
    cover.17[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.17[i,4]=ind(1.41,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.17[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.17[i,6]=ismb.fit[2,1]
    cover.17[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.17[i,8]=ind(1.77,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.17[i] = NA
      b0.is.sd.17[i] = NA
      b1.is.17[i] = NA
      b1.is.sd.17[i] = NA
      # Coverage
      cover.17[i,1]=NA
      cover.17[i,2]=NA
      cover.17[i,3]=NA
      cover.17[i,4]=NA
      cover.17[i,5]=NA
      cover.17[i,6]=NA
      cover.17[i,7]=NA
      cover.17[i,8]=NA
    })
}

# Crq beta table
table1.crq[5,1]<-mean(b0.crq.17,na.rm=TRUE)
table1.crq[5,2]<-sd(b0.crq.17,na.rm=TRUE)
table1.crq[5,3]<-mean(b1.crq.17,na.rm=TRUE)
table1.crq[5,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table1.is[5,1]<-mean(b0.is.17,na.rm=TRUE)
table1.is[5,2]<-sd(b0.is.17,na.rm=TRUE)
table1.is[5,3]<-mean(b1.is.17,na.rm=TRUE)
table1.is[5,4]<-sd(b1.is.17,na.rm=TRUE)

# Variance table (Beta0)
b0var.table1[5,1]=sd(b0.is.17,na.rm=TRUE)
b0var.table1[5,2]=mean(b0.crq.sd.17,na.rm=TRUE)
b0var.table1[5,3]=mean(b0.is.sd.17,na.rm=TRUE)

# Variance table (Beta1)
b1var.table1[5,1]=sd(b1.is.17,na.rm=TRUE)
b1var.table1[5,2]=mean(b1.crq.sd.17,na.rm=TRUE)
b1var.table1[5,3]=mean(b1.is.sd.17,na.rm=TRUE)

# Coverage table
b0.coverage[2,5]=mean(cover.17[,4],na.rm=TRUE)
b1.coverage[2,5]=mean(cover.17[,8],na.rm=TRUE)


#### t_0=2 & c=0% ####
b0.crq.20<-c()
b0.crq.sd.20<-c()
b1.crq.20<-c()
b1.crq.sd.20<-c()
b0.is.20<-c()
b0.is.sd.20<-c()
b1.is.20<-c()
b1.is.sd.20<-c()
cover.20=matrix(NA,2000,8)
colnames(cover.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(21)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,2)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.25)
    b0.crq.20[i] = crq.fit[1,1]
    b0.crq.sd.20[i] = crq.fit[1,2]
    b1.crq.20[i] = crq.fit[2,1]
    b1.crq.sd.20[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.20[i] = NA
      b0.crq.sd.20[i] = NA
      b1.crq.20[i] = NA
      b1.crq.sd.20[i] = NA
    })
}
set.seed(21)
for (i in 1:2000){
  a<-data.gen(c.0,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.25, 100)
    b0.is.20[i] = ismb.fit[1,1]
    b0.is.sd.20[i] = ismb.fit[1,2]
    b1.is.20[i] = ismb.fit[2,1]
    b1.is.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.20[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.20[i,2]=ismb.fit[1,1]
    cover.20[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.20[i,4]=ind(1.22,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.20[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.20[i,6]=ismb.fit[2,1]
    cover.20[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.20[i,8]=ind(1.92,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.20[i] = NA
      b0.is.sd.20[i] = NA
      b1.is.20[i] = NA
      b1.is.sd.20[i] = NA
      # Coverage
      cover.20[i,1]=NA
      cover.20[i,2]=NA
      cover.20[i,3]=NA
      cover.20[i,4]=NA
      cover.20[i,5]=NA
      cover.20[i,6]=NA
      cover.20[i,7]=NA
      cover.20[i,8]=NA
    })
}
# Crq beta table
table2.crq[1,1]<-mean(b0.crq.20,na.rm=TRUE)
table2.crq[1,2]<-sd(b0.crq.20,na.rm=TRUE)
table2.crq[1,3]<-mean(b1.crq.20,na.rm=TRUE)
table2.crq[1,4]<-sd(b1.crq.20,na.rm=TRUE)

# IS beta table
table2.is[1,1]<-mean(b0.is.20,na.rm=TRUE)
table2.is[1,2]<-sd(b0.is.20,na.rm=TRUE)
table2.is[1,3]<-mean(b1.is.20,na.rm=TRUE)
table2.is[1,4]<-sd(b1.is.20,na.rm=TRUE)

# Variance table (Beta0)
b0var.table2[1,1]=sd(b0.is.20,na.rm=TRUE)
b0var.table2[1,2]=mean(b0.crq.sd.20,na.rm=TRUE)
b0var.table2[1,3]=mean(b0.is.sd.20,na.rm=TRUE)

# Variance table (Beta1)
b1var.table2[1,1]=sd(b1.is.20,na.rm=TRUE)
b1var.table2[1,2]=mean(b1.crq.sd.20,na.rm=TRUE)
b1var.table2[1,3]=mean(b1.is.sd.20,na.rm=TRUE)

# Coverage table
b0.coverage[3,1]=mean(cover.20[,4],na.rm=TRUE)
b1.coverage[3,1]=mean(cover.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.crq.21<-c()
b0.crq.sd.21<-c()
b1.crq.21<-c()
b1.crq.sd.21<-c()
b0.is.21<-c()
b0.is.sd.21<-c()
b1.is.21<-c()
b1.is.sd.21<-c()
cover.21=matrix(NA,2000,8)
colnames(cover.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(21)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,2)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.25)
    b0.crq.21[i] = crq.fit[1,1]
    b0.crq.sd.21[i] = crq.fit[1,2]
    b1.crq.21[i] = crq.fit[2,1]
    b1.crq.sd.21[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.21[i] = NA
      b0.crq.sd.21[i] = NA
      b1.crq.21[i] = NA
      b1.crq.sd.21[i] = NA
    })
}
set.seed(21)
for (i in 1:2000){
  a<-data.gen(c.1,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.25, 100)
    b0.is.21[i] = ismb.fit[1,1]
    b0.is.sd.21[i] = ismb.fit[1,2]
    b1.is.21[i] = ismb.fit[2,1]
    b1.is.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.21[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.21[i,2]=ismb.fit[1,1]
    cover.21[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.21[i,4]=ind(1.22,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.21[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.21[i,6]=ismb.fit[2,1]
    cover.21[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.21[i,8]=ind(1.92,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.21[i] = NA
      b0.is.sd.21[i] = NA
      b1.is.21[i] = NA
      b1.is.sd.21[i] = NA
      # Coverage
      cover.21[i,1]=NA
      cover.21[i,2]=NA
      cover.21[i,3]=NA
      cover.21[i,4]=NA
      cover.21[i,5]=NA
      cover.21[i,6]=NA
      cover.21[i,7]=NA
      cover.21[i,8]=NA
    })
}

# Crq beta table
table2.crq[2,1]<-mean(b0.crq.21,na.rm=TRUE)
table2.crq[2,2]<-sd(b0.crq.21,na.rm=TRUE)
table2.crq[2,3]<-mean(b1.crq.21,na.rm=TRUE)
table2.crq[2,4]<-sd(b1.crq.21,na.rm=TRUE)

# IS beta table
table2.is[2,1]<-mean(b0.is.21,na.rm=TRUE)
table2.is[2,2]<-sd(b0.is.21,na.rm=TRUE)
table2.is[2,3]<-mean(b1.is.21,na.rm=TRUE)
table2.is[2,4]<-sd(b1.is.21,na.rm=TRUE)

# Variance table (Beta0)
b0var.table2[2,1]=sd(b0.is.21,na.rm=TRUE)
b0var.table2[2,2]=mean(b0.crq.sd.21,na.rm=TRUE)
b0var.table2[2,3]=mean(b0.is.sd.21,na.rm=TRUE)

# Variance table (Beta1)
b1var.table2[2,1]=sd(b1.is.21,na.rm=TRUE)
b1var.table2[2,2]=mean(b1.crq.sd.21,na.rm=TRUE)
b1var.table2[2,3]=mean(b1.is.sd.21,na.rm=TRUE)

# Coverage table
b0.coverage[3,2]=mean(cover.21[,4],na.rm=TRUE)
b1.coverage[3,2]=mean(cover.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.crq.23<-c()
b0.crq.sd.23<-c()
b1.crq.23<-c()
b1.crq.sd.23<-c()
b0.is.23<-c()
b0.is.sd.23<-c()
b1.is.23<-c()
b1.is.sd.23<-c()
cover.23=matrix(NA,2000,8)
colnames(cover.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(21)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,2)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.25)
    b0.crq.23[i] = crq.fit[1,1]
    b0.crq.sd.23[i] = crq.fit[1,2]
    b1.crq.23[i] = crq.fit[2,1]
    b1.crq.sd.23[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.23[i] = NA
      b0.crq.sd.23[i] = NA
      b1.crq.23[i] = NA
      b1.crq.sd.23[i] = NA
    })
}
set.seed(21)
for (i in 1:2000){
  a<-data.gen(c.3,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.25, 100)
    b0.is.23[i] = ismb.fit[1,1]
    b0.is.sd.23[i] = ismb.fit[1,2]
    b1.is.23[i] = ismb.fit[2,1]
    b1.is.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.23[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.23[i,2]=ismb.fit[1,1]
    cover.23[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.23[i,4]=ind(1.22,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.23[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.23[i,6]=ismb.fit[2,1]
    cover.23[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.23[i,8]=ind(1.92,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.23[i] = NA
      b0.is.sd.23[i] = NA
      b1.is.23[i] = NA
      b1.is.sd.23[i] = NA
      # Coverage
      cover.23[i,1]=NA
      cover.23[i,2]=NA
      cover.23[i,3]=NA
      cover.23[i,4]=NA
      cover.23[i,5]=NA
      cover.23[i,6]=NA
      cover.23[i,7]=NA
      cover.23[i,8]=NA
    })
}
# Crq beta table
table2.crq[3,1]<-mean(b0.crq.23,na.rm=TRUE)
table2.crq[3,2]<-sd(b0.crq.23,na.rm=TRUE)
table2.crq[3,3]<-mean(b1.crq.23,na.rm=TRUE)
table2.crq[3,4]<-sd(b1.crq.23,na.rm=TRUE)

# IS beta table
table2.is[3,1]<-mean(b0.is.23,na.rm=TRUE)
table2.is[3,2]<-sd(b0.is.23,na.rm=TRUE)
table2.is[3,3]<-mean(b1.is.23,na.rm=TRUE)
table2.is[3,4]<-sd(b1.is.23,na.rm=TRUE)

# Variance table (Beta0)
b0var.table2[3,1]=sd(b0.is.23,na.rm=TRUE)
b0var.table2[3,2]=mean(b0.crq.sd.23,na.rm=TRUE)
b0var.table2[3,3]=mean(b0.is.sd.23,na.rm=TRUE)

# Variance table (Beta1)
b1var.table2[3,1]=sd(b1.is.23,na.rm=TRUE)
b1var.table2[3,2]=mean(b1.crq.sd.23,na.rm=TRUE)
b1var.table2[3,3]=mean(b1.is.sd.23,na.rm=TRUE)

# Coverage table
b0.coverage[3,3]=mean(cover.23[,4],na.rm=TRUE)
b1.coverage[3,3]=mean(cover.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.crq.25<-c()
b0.crq.sd.25<-c()
b1.crq.25<-c()
b1.crq.sd.25<-c()
b0.is.25<-c()
b0.is.sd.25<-c()
b1.is.25<-c()
b1.is.sd.25<-c()
cover.25=matrix(NA,2000,8)
colnames(cover.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(21)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.5,2)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.25)
    b0.crq.25[i] = crq.fit[1,1]
    b0.crq.sd.25[i] = crq.fit[1,2]
    b1.crq.25[i] = crq.fit[2,1]
    b1.crq.sd.25[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.25[i] = NA
      b0.crq.sd.25[i] = NA
      b1.crq.25[i] = NA
      b1.crq.sd.25[i] = NA
    })
}
set.seed(21)
for (i in 1:2000){
  a<-data.gen(c.5,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.25, 100)
    b0.is.25[i] = ismb.fit[1,1]
    b0.is.sd.25[i] = ismb.fit[1,2]
    b1.is.25[i] = ismb.fit[2,1]
    b1.is.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.25[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.25[i,2]=ismb.fit[1,1]
    cover.25[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.25[i,4]=ind(1.22,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.25[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.25[i,6]=ismb.fit[2,1]
    cover.25[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.25[i,8]=ind(1.92,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.25[i] = NA
      b0.is.sd.25[i] = NA
      b1.is.25[i] = NA
      b1.is.sd.25[i] = NA
      # Coverage
      cover.25[i,1]=NA
      cover.25[i,2]=NA
      cover.25[i,3]=NA
      cover.25[i,4]=NA
      cover.25[i,5]=NA
      cover.25[i,6]=NA
      cover.25[i,7]=NA
      cover.25[i,8]=NA
    })
}

# Crq beta table
table2.crq[4,1]<-mean(b0.crq.25,na.rm=TRUE)
table2.crq[4,2]<-sd(b0.crq.25,na.rm=TRUE)
table2.crq[4,3]<-mean(b1.crq.25,na.rm=TRUE)
table2.crq[4,4]<-sd(b1.crq.25,na.rm=TRUE)

# IS beta table
table2.is[4,1]<-mean(b0.is.25,na.rm=TRUE)
table2.is[4,2]<-sd(b0.is.25,na.rm=TRUE)
table2.is[4,3]<-mean(b1.is.25,na.rm=TRUE)
table2.is[4,4]<-sd(b1.is.25,na.rm=TRUE)

# Variance table (Beta0)
b0var.table2[4,1]=sd(b0.is.25,na.rm=TRUE)
b0var.table2[4,2]=mean(b0.crq.sd.25,na.rm=TRUE)
b0var.table2[4,3]=mean(b0.is.sd.25,na.rm=TRUE)

# Variance table (Beta1)
b1var.table2[4,1]=sd(b1.is.25,na.rm=TRUE)
b1var.table2[4,2]=mean(b1.crq.sd.25,na.rm=TRUE)
b1var.table2[4,3]=mean(b1.is.sd.25,na.rm=TRUE)

# Coverage table
b0.coverage[3,4]=mean(cover.25[,4],na.rm=TRUE)
b1.coverage[3,4]=mean(cover.25[,8],na.rm=TRUE)

b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
cover.27=matrix(NA,2000,8)
colnames(cover.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(21)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.7,2)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.25)
    b0.crq.27[i] = crq.fit[1,1]
    b0.crq.sd.27[i] = crq.fit[1,2]
    b1.crq.27[i] = crq.fit[2,1]
    b1.crq.sd.27[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.27[i] = NA
      b0.crq.sd.27[i] = NA
      b1.crq.27[i] = NA
      b1.crq.sd.27[i] = NA
    })
}
set.seed(21)
for (i in 1:2000){
  a<-data.gen(c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.25, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.27[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.27[i,2]=ismb.fit[1,1]
    cover.27[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.27[i,4]=ind(1.22,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.27[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.27[i,6]=ismb.fit[2,1]
    cover.27[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.27[i,8]=ind(1.92,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.27[i] = NA
      b0.is.sd.27[i] = NA
      b1.is.27[i] = NA
      b1.is.sd.27[i] = NA
      # Coverage
      cover.27[i,1]=NA
      cover.27[i,2]=NA
      cover.27[i,3]=NA
      cover.27[i,4]=NA
      cover.27[i,5]=NA
      cover.27[i,6]=NA
      cover.27[i,7]=NA
      cover.27[i,8]=NA
    })
}

# Crq beta table
table2.crq[5,1]<-mean(b0.crq.27,na.rm=TRUE)
table2.crq[5,2]<-sd(b0.crq.27,na.rm=TRUE)
table2.crq[5,3]<-mean(b1.crq.27,na.rm=TRUE)
table2.crq[5,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table2.is[5,1]<-mean(b0.is.27,na.rm=TRUE)
table2.is[5,2]<-sd(b0.is.27,na.rm=TRUE)
table2.is[5,3]<-mean(b1.is.27,na.rm=TRUE)
table2.is[5,4]<-sd(b1.is.27,na.rm=TRUE)

# Variance table (Beta0)
b0var.table2[5,1]=sd(b0.is.27,na.rm=TRUE)
b0var.table2[5,2]=mean(b0.crq.sd.27,na.rm=TRUE)
b0var.table2[5,3]=mean(b0.is.sd.27,na.rm=TRUE)

# Variance table (Beta1)
b1var.table2[5,1]=sd(b1.is.27,na.rm=TRUE)
b1var.table2[5,2]=mean(b1.crq.sd.27,na.rm=TRUE)
b1var.table2[5,3]=mean(b1.is.sd.27,na.rm=TRUE)

# Coverage table
b0.coverage[3,5]=mean(cover.27[,4],na.rm=TRUE)
b1.coverage[3,5]=mean(cover.27[,8],na.rm=TRUE)

#### t_0=3 & c=0% ####
b0.crq.30<-c()
b0.crq.sd.30<-c()
b1.crq.30<-c()
b1.crq.sd.30<-c()
b0.is.30<-c()
b0.is.sd.30<-c()
b1.is.30<-c()
b1.is.sd.30<-c()
cover.30=matrix(NA,2000,8)
colnames(cover.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(31)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,3)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.25)
    b0.crq.30[i] = crq.fit[1,1]
    b0.crq.sd.30[i] = crq.fit[1,2]
    b1.crq.30[i] = crq.fit[2,1]
    b1.crq.sd.30[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.30[i] = NA
      b0.crq.sd.30[i] = NA
      b1.crq.30[i] = NA
      b1.crq.sd.30[i] = NA
    })
}
set.seed(31)
for (i in 1:2000){
  a<-data.gen(c.0,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.25, 100)
    b0.is.30[i] = ismb.fit[1,1]
    b0.is.sd.30[i] = ismb.fit[1,2]
    b1.is.30[i] = ismb.fit[2,1]
    b1.is.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.30[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.30[i,2]=ismb.fit[1,1]
    cover.30[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.30[i,4]=ind(1.04,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.30[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.30[i,6]=ismb.fit[2,1]
    cover.30[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.30[i,8]=ind(2.06,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.30[i] = NA
      b0.is.sd.30[i] = NA
      b1.is.30[i] = NA
      b1.is.sd.30[i] = NA
      # Coverage
      cover.30[i,1]=NA
      cover.30[i,2]=NA
      cover.30[i,3]=NA
      cover.30[i,4]=NA
      cover.30[i,5]=NA
      cover.30[i,6]=NA
      cover.30[i,7]=NA
      cover.30[i,8]=NA
    })
}

# Crq beta table
table3.crq[1,1]<-mean(b0.crq.30,na.rm=TRUE)
table3.crq[1,2]<-sd(b0.crq.30,na.rm=TRUE)
table3.crq[1,3]<-mean(b1.crq.30,na.rm=TRUE)
table3.crq[1,4]<-sd(b1.crq.30,na.rm=TRUE)

# IS beta table
table3.is[1,1]<-mean(b0.is.30,na.rm=TRUE)
table3.is[1,2]<-sd(b0.is.30,na.rm=TRUE)
table3.is[1,3]<-mean(b1.is.30,na.rm=TRUE)
table3.is[1,4]<-sd(b1.is.30,na.rm=TRUE)

# Variance table (Beta0)
b0var.table3[1,1]=sd(b0.is.30,na.rm=TRUE)
b0var.table3[1,2]=mean(b0.crq.sd.30,na.rm=TRUE)
b0var.table3[1,3]=mean(b0.is.sd.30,na.rm=TRUE)

# Variance table (Beta1)
b1var.table3[1,1]=sd(b1.is.30,na.rm=TRUE)
b1var.table3[1,2]=mean(b1.crq.sd.30,na.rm=TRUE)
b1var.table3[1,3]=mean(b1.is.sd.30,na.rm=TRUE)

# Coverage table
b0.coverage[4,1]=mean(cover.30[,4],na.rm=TRUE)
b1.coverage[4,1]=mean(cover.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.crq.31<-c()
b0.crq.sd.31<-c()
b1.crq.31<-c()
b1.crq.sd.31<-c()
b0.is.31<-c()
b0.is.sd.31<-c()
b1.is.31<-c()
b1.is.sd.31<-c()
cover.31=matrix(NA,2000,8)
colnames(cover.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(31)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,3)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.25)
    b0.crq.31[i] = crq.fit[1,1]
    b0.crq.sd.31[i] = crq.fit[1,2]
    b1.crq.31[i] = crq.fit[2,1]
    b1.crq.sd.31[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.31[i] = NA
      b0.crq.sd.31[i] = NA
      b1.crq.31[i] = NA
      b1.crq.sd.31[i] = NA
    })
}
set.seed(31)
for (i in 1:2000){
  a<-data.gen(c.1,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.25, 100)
    b0.is.31[i] = ismb.fit[1,1]
    b0.is.sd.31[i] = ismb.fit[1,2]
    b1.is.31[i] = ismb.fit[2,1]
    b1.is.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.31[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.31[i,2]=ismb.fit[1,1]
    cover.31[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.31[i,4]=ind(1.04,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.31[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.31[i,6]=ismb.fit[2,1]
    cover.31[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.31[i,8]=ind(2.06,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.31[i] = NA
      b0.is.sd.31[i] = NA
      b1.is.31[i] = NA
      b1.is.sd.31[i] = NA
      # Coverage
      cover.31[i,1]=NA
      cover.31[i,2]=NA
      cover.31[i,3]=NA
      cover.31[i,4]=NA
      cover.31[i,5]=NA
      cover.31[i,6]=NA
      cover.31[i,7]=NA
      cover.31[i,8]=NA
    })
}

# Crq beta table
table3.crq[2,1]<-mean(b0.crq.31,na.rm=TRUE)
table3.crq[2,2]<-sd(b0.crq.31,na.rm=TRUE)
table3.crq[2,3]<-mean(b1.crq.31,na.rm=TRUE)
table3.crq[2,4]<-sd(b1.crq.31,na.rm=TRUE)

# IS beta table
table3.is[2,1]<-mean(b0.is.31,na.rm=TRUE)
table3.is[2,2]<-sd(b0.is.31,na.rm=TRUE)
table3.is[2,3]<-mean(b1.is.31,na.rm=TRUE)
table3.is[2,4]<-sd(b1.is.31,na.rm=TRUE)

# Variance table (Beta0)
b0var.table3[2,1]=sd(b0.is.31,na.rm=TRUE)
b0var.table3[2,2]=mean(b0.crq.sd.31,na.rm=TRUE)
b0var.table3[2,3]=mean(b0.is.sd.31,na.rm=TRUE)

# Variance table (Beta1)
b1var.table3[2,1]=sd(b1.is.31,na.rm=TRUE)
b1var.table3[2,2]=mean(b1.crq.sd.31,na.rm=TRUE)
b1var.table3[2,3]=mean(b1.is.sd.31,na.rm=TRUE)

# Coverage table
b0.coverage[4,2]=mean(cover.31[,4],na.rm=TRUE)
b1.coverage[4,2]=mean(cover.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
set.seed(31)
b0.crq.33<-c()
b0.crq.sd.33<-c()
b1.crq.33<-c()
b1.crq.sd.33<-c()
b0.is.33<-c()
b0.is.sd.33<-c()
b1.is.33<-c()
b1.is.sd.33<-c()
cover.33=matrix(NA,2000,8)
colnames(cover.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,3)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.25)
    b0.crq.33[i] = crq.fit[1,1]
    b0.crq.sd.33[i] = crq.fit[1,2]
    b1.crq.33[i] = crq.fit[2,1]
    b1.crq.sd.33[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.33[i] = NA
      b0.crq.sd.33[i] = NA
      b1.crq.33[i] = NA
      b1.crq.sd.33[i] = NA
    })
  ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.25, 100)
  b0.is.33[i] = ismb.fit[1,1]
  b0.is.sd.33[i] = ismb.fit[1,2]
  b1.is.33[i] = ismb.fit[2,1]
  b1.is.sd.33[i] = ismb.fit[2,2]
  # Coverage
  cover.33[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
  cover.33[i,2]=ismb.fit[1,1]
  cover.33[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
  cover.33[i,4]=ind(1.04,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
  cover.33[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
  cover.33[i,6]=ismb.fit[2,1]
  cover.33[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
  cover.33[i,8]=ind(2.06,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])
}

# Crq beta table
table3.crq[3,1]<-mean(b0.crq.33,na.rm=TRUE)
table3.crq[3,2]<-sd(b0.crq.33,na.rm=TRUE)
table3.crq[3,3]<-mean(b1.crq.33,na.rm=TRUE)
table3.crq[3,4]<-sd(b1.crq.33,na.rm=TRUE)

# IS beta table
table3.is[3,1]<-mean(b0.is.33,na.rm=TRUE)
table3.is[3,2]<-sd(b0.is.33,na.rm=TRUE)
table3.is[3,3]<-mean(b1.is.33,na.rm=TRUE)
table3.is[3,4]<-sd(b1.is.33,na.rm=TRUE)

# Variance table (Beta0)
b0var.table3[3,1]=sd(b0.is.33,na.rm=TRUE)
b0var.table3[3,2]=mean(b0.crq.sd.33,na.rm=TRUE)
b0var.table3[3,3]=mean(b0.is.sd.33,na.rm=TRUE)

# Variance table (Beta1)
b1var.table3[3,1]=sd(b1.is.33,na.rm=TRUE)
b1var.table3[3,2]=mean(b1.crq.sd.33,na.rm=TRUE)
b1var.table3[3,3]=mean(b1.is.sd.33,na.rm=TRUE)

# Coverage table
b0.coverage[4,3]=mean(cover.33[,4],na.rm=TRUE)
b1.coverage[4,3]=mean(cover.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.crq.35<-c()
b0.crq.sd.35<-c()
b1.crq.35<-c()
b1.crq.sd.35<-c()
b0.is.35<-c()
b0.is.sd.35<-c()
b1.is.35<-c()
b1.is.sd.35<-c()
cover.35=matrix(NA,2000,8)
colnames(cover.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(31)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.5,3)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.25)
    b0.crq.35[i] = crq.fit[1,1]
    b0.crq.sd.35[i] = crq.fit[1,2]
    b1.crq.35[i] = crq.fit[2,1]
    b1.crq.sd.35[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.35[i] = NA
      b0.crq.sd.35[i] = NA
      b1.crq.35[i] = NA
      b1.crq.sd.35[i] = NA
    })
}
set.seed(31)
for (i in 1:2000){
  a<-data.gen(c.5,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.25, 100)
    b0.is.35[i] = ismb.fit[1,1]
    b0.is.sd.35[i] = ismb.fit[1,2]
    b1.is.35[i] = ismb.fit[2,1]
    b1.is.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.35[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.35[i,2]=ismb.fit[1,1]
    cover.35[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.35[i,4]=ind(1.04,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.35[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.35[i,6]=ismb.fit[2,1]
    cover.35[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.35[i,8]=ind(2.06,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.35[i] = NA
      b0.is.sd.35[i] = NA
      b1.is.35[i] = NA
      b1.is.sd.35[i] = NA
      # Coverage
      cover.35[i,1]=NA
      cover.35[i,2]=NA
      cover.35[i,3]=NA
      cover.35[i,4]=NA
      cover.35[i,5]=NA
      cover.35[i,6]=NA
      cover.35[i,7]=NA
      cover.35[i,8]=NA
    })
}

# Crq beta table
table3.crq[4,1]<-mean(b0.crq.35,na.rm=TRUE)
table3.crq[4,2]<-sd(b0.crq.35,na.rm=TRUE)
table3.crq[4,3]<-mean(b1.crq.35,na.rm=TRUE)
table3.crq[4,4]<-sd(b1.crq.35,na.rm=TRUE)

# IS beta table
table3.is[4,1]<-mean(b0.is.35,na.rm=TRUE)
table3.is[4,2]<-sd(b0.is.35,na.rm=TRUE)
table3.is[4,3]<-mean(b1.is.35,na.rm=TRUE)
table3.is[4,4]<-sd(b1.is.35,na.rm=TRUE)

# Variance table (Beta0)
b0var.table3[4,1]=sd(b0.is.35,na.rm=TRUE)
b0var.table3[4,2]=mean(b0.crq.sd.35,na.rm=TRUE)
b0var.table3[4,3]=mean(b0.is.sd.35,na.rm=TRUE)

# Variance table (Beta1)
b1var.table3[4,1]=sd(b1.is.35,na.rm=TRUE)
b1var.table3[4,2]=mean(b1.crq.sd.35,na.rm=TRUE)
b1var.table3[4,3]=mean(b1.is.sd.35,na.rm=TRUE)

# Coverage table
b0.coverage[4,4]=mean(cover.35[,4],na.rm=TRUE)
b1.coverage[4,4]=mean(cover.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
cover.37=matrix(NA,2000,8)
colnames(cover.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(31)
for (i in 1:2000){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.7,3)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.25)
    b0.crq.37[i] = crq.fit[1,1]
    b0.crq.sd.37[i] = crq.fit[1,2]
    b1.crq.37[i] = crq.fit[2,1]
    b1.crq.sd.37[i] = crq.fit[2,2]}
    ,error=function(e){
      b0.crq.37[i] = NA
      b0.crq.sd.37[i] = NA
      b1.crq.37[i] = NA
      b1.crq.sd.37[i] = NA
    })
}
set.seed(31)
for (i in 1:2000){
  a<-data.gen(c.7,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.25, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.37[i,1]=ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.37[i,2]=ismb.fit[1,1]
    cover.37[i,3]=ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.37[i,4]=ind(1.04,ismb.fit[1,1]-1.96*ismb.fit[1,2],ismb.fit[1,1]+1.96*ismb.fit[1,2])
    cover.37[i,5]=ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.37[i,6]=ismb.fit[2,1]
    cover.37[i,7]=ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.37[i,8]=ind(2.06,ismb.fit[2,1]-1.96*ismb.fit[2,2],ismb.fit[2,1]+1.96*ismb.fit[2,2])}
    , error=function(e){
      b0.is.37[i] = NA
      b0.is.sd.37[i] = NA
      b1.is.37[i] = NA
      b1.is.sd.37[i] = NA
      # Coverage
      cover.37[i,1]=NA
      cover.37[i,2]=NA
      cover.37[i,3]=NA
      cover.37[i,4]=NA
      cover.37[i,5]=NA
      cover.37[i,6]=NA
      cover.37[i,7]=NA
      cover.37[i,8]=NA
    })
}

# Crq beta table
table3.crq[5,1]<-mean(b0.crq.37,na.rm=TRUE)
table3.crq[5,2]<-sd(b0.crq.37,na.rm=TRUE)
table3.crq[5,3]<-mean(b1.crq.37,na.rm=TRUE)
table3.crq[5,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table3.is[5,1]<-mean(b0.is.37,na.rm=TRUE)
table3.is[5,2]<-sd(b0.is.37,na.rm=TRUE)
table3.is[5,3]<-mean(b1.is.37,na.rm=TRUE)
table3.is[5,4]<-sd(b1.is.37,na.rm=TRUE)

# Variance table (Beta0)
b0var.table3[5,1]=sd(b0.is.37,na.rm=TRUE)
b0var.table3[5,2]=mean(b0.crq.sd.37,na.rm=TRUE)
b0var.table3[5,3]=mean(b0.is.sd.37,na.rm=TRUE)

# Variance table (Beta1)
b1var.table3[5,1]=sd(b1.is.37,na.rm=TRUE)
b1var.table3[5,2]=mean(b1.crq.sd.37,na.rm=TRUE)
b1var.table3[5,3]=mean(b1.is.sd.37,na.rm=TRUE)

# Coverage table
b0.coverage[4,5]=mean(cover.37[,4],na.rm=TRUE)
b1.coverage[4,5]=mean(cover.37[,8],na.rm=TRUE)

