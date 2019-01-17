library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)
library(xtable)
  
#### Given information ####
#Find C
#Generate T_i 1,000,000
#d.exp.beta.initial=5
#d.k=2
#d.r.initial=(log(2)^(1/d.k))/d.exp.beta.initial
#d.u=runif(n=1000000,min = 0,max = 1)
#d.T={-log(1-d.u)}^(1/d.k)/d.r.initial
#Find c which is derermined to achieve 0,10,20,30% censoring rate in case 1
#i=30
#while(1){
#  d.C<-runif(1000000,0,i)
#  i=i+0.01
#  if(sum(d.C<d.T)<=200000) break # 0% = 0, 10% = 1000000, 20% = 2000000, 30% = 3000000)
#}
#print(i)
#table(d.C<d.T)

c.0=5000000
c.1=53.03
c.2=26.58
c.3=17.73
exp.beta.initial=5
k=2
r.initial=(log(2))^(1/k)/exp.beta.initial

#### Data Generation function ####
data.gen<-function(censor, t_0){
  unif = runif(n=400,min = 0,max = 1)
  sim=matrix(NA,400,7)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5))
  sim[,1] = {{-log(1-unif)}^(1/k)}/r.initial
  # Generate C_i
  sim[,2] = runif(400,0,censor)
  # Generate Y_i (min(T,C))
  sim[,3] = apply(sim[,1:2], 1, FUN=min)
  sim[,4] = sim[,3]-(t_0)
  sim[,5] = log(sim[,4])
  # Covariates (Control=0, Treatment=1)
  sim[,6] = rbinom(400,size=1,p=0.5)
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,7]=I(sim[,1]<sim[,2])
  # Delete data that smaller than t_0
  sim = na.omit(sim)
  # Ordering
  sim = sim[order(sim[,4]),]
  n = nrow(sim)
  sim = as.data.frame(sim)
  ## Calcuation W_i
  for (i in 1:n){
    sim[i,8]<-sum(sim[,3]>=sim[i,3])
    sim[i,9]<-sum(sim[,3]==sim[i,3])
    if (sim[i,7]==0){
      sim[i,10]=(sim[i,8]-sim[i,9])/sim[i,8]
    } else {
      sim[i,10]=1
    }
    if (i==1){
      sim[i,11]<-sim[i,10]
    } else {
      sim[i,11]<-sim[i-1,11]*sim[i,10]
    }
    if (sim[i,11]==0){
      if (sim[i,7]==0){
        sim[i,12]=0
      } else {
        sim[i,12]<-sim[i-1,12]}
    }
    else {
      sim[i,12]<-sim[i,7]/(sim[i,11])}
  }
  # Column names
  colnames(sim) = c("T","C","Z","Z.diff","log(Z.diff)","X","censored","# at risk","# event","s/d","G_KM","Weight")
  return(sim)
}

#### Objective equation ####
objectF<-function(beta){
  beta<-as.matrix(beta)
  result=t(w*X)%*%(pnorm((a[,5]-X%*%beta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
  print(result)
}
#### revised object equation (with eta) ####
rev.objectF=function(beta){
  beta=as.matrix(beta)
  result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%beta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
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


#### Beta estimation and variance estimation and Coverage ####
# t_0=0
table0.crq<-matrix(NA,4,4)
rownames(table0.crq)<-c(0,10,20,30)
colnames(table0.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1.crq<-matrix(NA,4,4)
rownames(table1.crq)<-c(0,10,20,30)
colnames(table1.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table2.crq<-matrix(NA,4,4)
rownames(table2.crq)<-c(0,10,20,30)
colnames(table2.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3.crq<-matrix(NA,4,4)
rownames(table3.crq)<-c(0,10,20,30)
colnames(table3.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
b0.crq<-c()
b0.crq.sd<-c()
b1.crq<-c()
b1.crq.sd<-c()

table0.is<-matrix(NA,4,4)
rownames(table0.is)<-c(0,10,20,30)
colnames(table0.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1.is<-matrix(NA,4,4)
rownames(table1.is)<-c(0,10,20,30)
colnames(table1.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table2.is<-matrix(NA,4,4)
rownames(table2.is)<-c(0,10,20,30)
colnames(table2.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3.is<-matrix(NA,4,4)
rownames(table3.is)<-c(0,10,20,30)
colnames(table3.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
b0.is<-c()
b1.is<-c()

b0var.table0=matrix(NA,4,4)
colnames(b0var.table0)=c("true","MB","ISMB","Crq")
rownames(b0var.table0)=c(0,10,20,30)
b0var.table1=matrix(NA,4,4)
colnames(b0var.table1)=c("true","MB","ISMB","Crq")
rownames(b0var.table1)=c(0,10,20,30)
b0var.table2=matrix(NA,4,4)
colnames(b0var.table2)=c("true","MB","ISMB","Crq")
rownames(b0var.table2)=c(0,10,20,30)
b0var.table3=matrix(NA,4,4)
colnames(b0var.table3)=c("true","MB","ISMB","Crq")
rownames(b0var.table3)=c(0,10,20,30)

b1var.table0=matrix(NA,4,4)
colnames(b1var.table0)=c("true","MB","ISMB","Crq")
rownames(b1var.table0)=c(0,10,20,30)
b1var.table1=matrix(NA,4,4)
colnames(b1var.table1)=c("true","MB","ISMB","Crq")
rownames(b1var.table1)=c(0,10,20,30)
b1var.table2=matrix(NA,4,4)
colnames(b1var.table2)=c("true","MB","ISMB","Crq")
rownames(b1var.table2)=c(0,10,20,30)
b1var.table3=matrix(NA,4,4)
colnames(b1var.table3)=c("true","MB","ISMB","Crq")
rownames(b1var.table3)=c(0,10,20,30)

mb.varb0=c()
mb.varb1=c()
ismb.varb0=c()
ismb.varb1=c()

b0.coverage=matrix(NA,4,4)
colnames(b0.coverage)=c("0%","10%","20%","30%")
rownames(b0.coverage)=c("t0=0","t0=1","t0=2","t0=3")

b1.coverage=matrix(NA,4,4)
colnames(b1.coverage)=c("0%","10%","20%","30%")
rownames(b1.coverage)=c("t0=0","t0=1","t0=2","t0=3")

cover=matrix(NA,500,8)
colnames(cover)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

#### t_0=0 & c=0% ####
set.seed(1)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,0)
  a[,2]<-100
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.50))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.50))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.50))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.50))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.61,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.61,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.61,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table0.crq[1,1]<-mean(b0.crq)
table0.crq[1,2]<-sd(b0.crq)
table0.crq[1,3]<-mean(b1.crq)
table0.crq[1,4]<-sd(b1.crq)

# IS beta table
table0.is[1,1]<-mean(b0.is)
table0.is[1,2]<-sd(b0.is)
table0.is[1,3]<-mean(b1.is)
table0.is[1,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table0[1,1]=sd(b0.is)
b0var.table0[1,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table0[1,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table0[1,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table0[1,1]=sd(b1.is)
b1var.table0[1,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table0[1,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table0[1,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[1,1]=mean(cover[,4],na.rm=TRUE)
b1.coverage[1,1]=mean(cover[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
set.seed(2)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,0)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.61,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.61,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.61,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table0.crq[2,1]<-mean(b0.crq)
table0.crq[2,2]<-sd(b0.crq)
table0.crq[2,3]<-mean(b1.crq)
table0.crq[2,4]<-sd(b1.crq)

# IS beta table
table0.is[2,1]<-mean(b0.is)
table0.is[2,2]<-sd(b0.is)
table0.is[2,3]<-mean(b1.is)
table0.is[2,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table0[2,1]=sd(b0.is)
b0var.table0[2,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table0[2,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table0[2,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table0[2,1]=sd(b1.is)
b1var.table0[2,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table0[2,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table0[2,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[1,2]=mean(cover[,4],na.rm=TRUE)
b1.coverage[1,2]=mean(cover[,8],na.rm=TRUE)

#### t_0=0 & c=20% ####
set.seed(3)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.2,0)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.61,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.61,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.61,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table0.crq[3,1]<-mean(b0.crq)
table0.crq[3,2]<-sd(b0.crq)
table0.crq[3,3]<-mean(b1.crq)
table0.crq[3,4]<-sd(b1.crq)

# IS beta table
table0.is[3,1]<-mean(b0.is)
table0.is[3,2]<-sd(b0.is)
table0.is[3,3]<-mean(b1.is)
table0.is[3,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table0[3,1]=sd(b0.is)
b0var.table0[3,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table0[3,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table0[3,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table0[3,1]=sd(b1.is)
b1var.table0[3,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table0[3,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table0[3,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[1,3]=mean(cover[,4],na.rm=TRUE)
b1.coverage[1,3]=mean(cover[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
set.seed(4)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,0)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.61,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.61,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.61,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table0.crq[4,1]<-mean(b0.crq)
table0.crq[4,2]<-sd(b0.crq)
table0.crq[4,3]<-mean(b1.crq)
table0.crq[4,4]<-sd(b1.crq)

# IS beta table
table0.is[4,1]<-mean(b0.is)
table0.is[4,2]<-sd(b0.is)
table0.is[4,3]<-mean(b1.is)
table0.is[4,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table0[4,1]=sd(b0.is)
b0var.table0[4,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table0[4,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table0[4,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table0[4,1]=sd(b1.is)
b1var.table0[4,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table0[4,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table0[4,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[1,4]=mean(cover[,4],na.rm=TRUE)
b1.coverage[1,4]=mean(cover[,8],na.rm=TRUE)



#### t_0=1 & c=0% ####
set.seed(11)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,1)
  a[,2]<-100
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.41,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.41,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.41,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table1.crq[1,1]<-mean(b0.crq)
table1.crq[1,2]<-sd(b0.crq)
table1.crq[1,3]<-mean(b1.crq)
table1.crq[1,4]<-sd(b1.crq)

# IS beta table
table1.is[1,1]<-mean(b0.is)
table1.is[1,2]<-sd(b0.is)
table1.is[1,3]<-mean(b1.is)
table1.is[1,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table1[1,1]=sd(b0.is)
b0var.table1[1,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table1[1,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table1[1,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table1[1,1]=sd(b1.is)
b1var.table1[1,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table1[1,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table1[1,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[2,1]=mean(cover[,4],na.rm=TRUE)
b1.coverage[2,1]=mean(cover[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
set.seed(12)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,1)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.41,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.41,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.41,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table1.crq[2,1]<-mean(b0.crq)
table1.crq[2,2]<-sd(b0.crq)
table1.crq[2,3]<-mean(b1.crq)
table1.crq[2,4]<-sd(b1.crq)

# IS beta table
table1.is[2,1]<-mean(b0.is)
table1.is[2,2]<-sd(b0.is)
table1.is[2,3]<-mean(b1.is)
table1.is[2,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table1[2,1]=sd(b0.is)
b0var.table1[2,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table1[2,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table1[2,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table1[2,1]=sd(b1.is)
b1var.table1[2,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table1[2,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table1[2,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[2,2]=mean(cover[,4],na.rm=TRUE)
b1.coverage[2,2]=mean(cover[,8],na.rm=TRUE)

#### t_0=1 & c=20% ####
set.seed(13)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.2,1)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.41,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.41,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.41,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table1.crq[3,1]<-mean(b0.crq)
table1.crq[3,2]<-sd(b0.crq)
table1.crq[3,3]<-mean(b1.crq)
table1.crq[3,4]<-sd(b1.crq)

# IS beta table
table1.is[3,1]<-mean(b0.is)
table1.is[3,2]<-sd(b0.is)
table1.is[3,3]<-mean(b1.is)
table1.is[3,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table1[3,1]=sd(b0.is)
b0var.table1[3,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table1[3,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table1[3,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table1[3,1]=sd(b1.is)
b1var.table1[3,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table1[3,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table1[3,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[2,3]=mean(cover[,4],na.rm=TRUE)
b1.coverage[2,3]=mean(cover[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
set.seed(14)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,1)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.41,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.41,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.41,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table1.crq[4,1]<-mean(b0.crq)
table1.crq[4,2]<-sd(b0.crq)
table1.crq[4,3]<-mean(b1.crq)
table1.crq[4,4]<-sd(b1.crq)

# IS beta table
table1.is[4,1]<-mean(b0.is)
table1.is[4,2]<-sd(b0.is)
table1.is[4,3]<-mean(b1.is)
table1.is[4,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table1[4,1]=sd(b0.is)
b0var.table1[4,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table1[4,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table1[4,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table1[4,1]=sd(b1.is)
b1var.table1[4,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table1[4,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table1[4,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[2,4]=mean(cover[,4],na.rm=TRUE)
b1.coverage[2,4]=mean(cover[,8],na.rm=TRUE)



#### t_0=2 & c=0% ####
set.seed(21)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,2)
  a[,2]<-100
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.22,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.22,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.22,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table2.crq[1,1]<-mean(b0.crq)
table2.crq[1,2]<-sd(b0.crq)
table2.crq[1,3]<-mean(b1.crq)
table2.crq[1,4]<-sd(b1.crq)

# IS beta table
table2.is[1,1]<-mean(b0.is)
table2.is[1,2]<-sd(b0.is)
table2.is[1,3]<-mean(b1.is)
table2.is[1,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table2[1,1]=sd(b0.is)
b0var.table2[1,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table2[1,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table2[1,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table2[1,1]=sd(b1.is)
b1var.table2[1,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table2[1,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table2[1,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[3,1]=mean(cover[,4],na.rm=TRUE)
b1.coverage[3,1]=mean(cover[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
set.seed(22)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,2)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.22,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.22,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.22,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table2.crq[2,1]<-mean(b0.crq)
table2.crq[2,2]<-sd(b0.crq)
table2.crq[2,3]<-mean(b1.crq)
table2.crq[2,4]<-sd(b1.crq)

# IS beta table
table2.is[2,1]<-mean(b0.is)
table2.is[2,2]<-sd(b0.is)
table2.is[2,3]<-mean(b1.is)
table2.is[2,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table2[2,1]=sd(b0.is)
b0var.table2[2,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table2[2,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table2[2,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table2[2,1]=sd(b1.is)
b1var.table2[2,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table2[2,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table2[2,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[3,2]=mean(cover[,4],na.rm=TRUE)
b1.coverage[3,2]=mean(cover[,8],na.rm=TRUE)

#### t_0=2 & c=20% ####
set.seed(23)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.2,2)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.22,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.22,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.22,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table2.crq[3,1]<-mean(b0.crq)
table2.crq[3,2]<-sd(b0.crq)
table2.crq[3,3]<-mean(b1.crq)
table2.crq[3,4]<-sd(b1.crq)

# IS beta table
table2.is[3,1]<-mean(b0.is)
table2.is[3,2]<-sd(b0.is)
table2.is[3,3]<-mean(b1.is)
table2.is[3,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table2[3,1]=sd(b0.is)
b0var.table2[3,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table2[3,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table2[3,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table2[3,1]=sd(b1.is)
b1var.table2[3,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table2[3,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table2[3,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[3,3]=mean(cover[,4],na.rm=TRUE)
b1.coverage[3,3]=mean(cover[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
set.seed(24)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,2)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.22,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.22,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.22,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table2.crq[4,1]<-mean(b0.crq)
table2.crq[4,2]<-sd(b0.crq)
table2.crq[4,3]<-mean(b1.crq)
table2.crq[4,4]<-sd(b1.crq)

# IS beta table
table2.is[4,1]<-mean(b0.is)
table2.is[4,2]<-sd(b0.is)
table2.is[4,3]<-mean(b1.is)
table2.is[4,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table2[4,1]=sd(b0.is)
b0var.table2[4,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table2[4,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table2[4,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table2[4,1]=sd(b1.is)
b1var.table2[4,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table2[4,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table2[4,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[3,4]=mean(cover[,4],na.rm=TRUE)
b1.coverage[3,4]=mean(cover[,8],na.rm=TRUE)



#### t_0=3 & c=0% ####
set.seed(31)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.0,3)
  a[,2]<-100
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.04,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.04,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.04,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table3.crq[1,1]<-mean(b0.crq)
table3.crq[1,2]<-sd(b0.crq)
table3.crq[1,3]<-mean(b1.crq)
table3.crq[1,4]<-sd(b1.crq)

# IS beta table
table3.is[1,1]<-mean(b0.is)
table3.is[1,2]<-sd(b0.is)
table3.is[1,3]<-mean(b1.is)
table3.is[1,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table3[1,1]=sd(b0.is)
b0var.table3[1,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table3[1,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table3[1,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table3[1,1]=sd(b1.is)
b1var.table3[1,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table3[1,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table3[1,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[4,1]=mean(cover[,4],na.rm=TRUE)
b1.coverage[4,1]=mean(cover[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
set.seed(32)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.1,3)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.04,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.04,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.04,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table3.crq[2,1]<-mean(b0.crq)
table3.crq[2,2]<-sd(b0.crq)
table3.crq[2,3]<-mean(b1.crq)
table3.crq[2,4]<-sd(b1.crq)

# IS beta table
table3.is[2,1]<-mean(b0.is)
table3.is[2,2]<-sd(b0.is)
table3.is[2,3]<-mean(b1.is)
table3.is[2,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table3[2,1]=sd(b0.is)
b0var.table3[2,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table3[2,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table3[2,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table3[2,1]=sd(b1.is)
b1var.table3[2,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table3[2,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table3[2,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[4,2]=mean(cover[,4],na.rm=TRUE)
b1.coverage[4,2]=mean(cover[,8],na.rm=TRUE)

#### t_0=3 & c=20% ####
set.seed(33)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.2,3)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.04,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.04,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.04,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table3.crq[3,1]<-mean(b0.crq)
table3.crq[3,2]<-sd(b0.crq)
table3.crq[3,3]<-mean(b1.crq)
table3.crq[3,4]<-sd(b1.crq)

# IS beta table
table3.is[3,1]<-mean(b0.is)
table3.is[3,2]<-sd(b0.is)
table3.is[3,3]<-mean(b1.is)
table3.is[3,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table3[3,1]=sd(b0.is)
b0var.table3[3,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table3[3,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table3[3,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table3[3,1]=sd(b1.is)
b1var.table3[3,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table3[3,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table3[3,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[4,3]=mean(cover[,4],na.rm=TRUE)
b1.coverage[4,3]=mean(cover[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
set.seed(34)
for (i in 1:500){
  # Beta estimation method 1 : Crq package
  a<-data.gen(c.3,3)
  crq.fit<-crq(Surv(a[,5],a[,7])~a[,6],method='Portnoy')
  b0.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,1]
  b0.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[1,4]
  b1.crq[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,1]
  b1.crq.sd[i]<-summary(crq.fit,tau=c(0.25,0.5))[2][[1]]$coefficient[2,4]
  
  # Beta estimation method 2 : Our method
  n=nrow(a)
  w=a[,12]
  X=cbind(c(rep(1,nrow(a))),a[,6])
  G<-matrix(c(1/n,0,0,1/n),2,2)
  betastart<-c(1.04,0) # from table0
  is.fit<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  b0.is[i]<-is.fit$x[1]
  b1.is[i]<-is.fit$x[2]
  solbeta=is.fit$x
  
  # Variance estimation
  mb.b0=c()
  mb.b1=c()
  result.mb.b0=c()
  result.mb.b1=c()
  result.ismb=c()
  if (is.fit$fvec<1e-5&is.fit$fvec>-1e-5){
    for (j in 1:500){
      # Variance estimation method 1 : Multiplier bootstrap
      eta=rexp(n,1)
      betastart=c(1.04,0)
      mb.fit=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
      mb.b0=cbind(mb.b0,mb.fit$x[1])
      mb.b1=cbind(mb.b1,mb.fit$x[2])
      result.mb.b0=cbind(result.mb.b0,mb.fit$fvec[1])
      result.mb.b1=cbind(result.mb.b1,mb.fit$fvec[2])
      # Variance estimation method 2 : ISMB
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.25)
      result.ismb=cbind(result.ismb,result)
    }
    
    # Variance estimation method 1 : Multiplier bootstrap
    mb.b0=mb.b0[-1e-5<result.mb.b0&result.mb.b0<1e-5]
    mb.b1=mb.b1[-1e-5<result.mb.b1&result.mb.b1<1e-5]
    mb.varb0[i]=var(mb.b0)
    mb.varb1[i]=var(mb.b1)
    
    # Variance estimation method 2 : ISMB
    v=cov(t(result.ismb))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    ismb.varb0[i]=sigma[1,1]
    ismb.varb1[i]=sigma[2,2]
    
    # Coverage
    b0se=sqrt(sigma[1,1])
    b1se=sqrt(sigma[2,2])
    cover[i,1]=is.fit$x[1]-1.96*b0se
    cover[i,2]=is.fit$x[1]
    cover[i,3]=is.fit$x[1]+1.96*b0se
    cover[i,4]=ind(1.04,is.fit$x[1]-1.96*b0se,is.fit$x[1]+1.96*b0se)
    cover[i,5]=is.fit$x[2]-1.96*b1se
    cover[i,6]=is.fit$x[2]
    cover[i,7]=is.fit$x[2]+1.96*b1se
    cover[i,8]=ind(0,is.fit$x[2]-1.96*b1se,is.fit$x[2]+1.96*b1se)
  } else {
    mb.varb0[i]=NA
    mb.varb1[i]=NA
    ismb.varb0[i]=NA
    ismb.varb1[i]=NA
    cover[i,1:8]=NA
  }
}
# Crq beta table
table3.crq[4,1]<-mean(b0.crq)
table3.crq[4,2]<-sd(b0.crq)
table3.crq[4,3]<-mean(b1.crq)
table3.crq[4,4]<-sd(b1.crq)

# IS beta table
table3.is[4,1]<-mean(b0.is)
table3.is[4,2]<-sd(b0.is)
table3.is[4,3]<-mean(b1.is)
table3.is[4,4]<-sd(b1.is)

# Variance table (Beta0)
b0var.table3[4,1]=sd(b0.is)
b0var.table3[4,2]=sqrt(mean(mb.varb0,na.rm=TRUE))
b0var.table3[4,3]=sqrt(mean(ismb.varb0,na.rm=TRUE))
b0var.table3[4,4]=mean(b0.crq.sd)

# Variance table (Beta1)
b1var.table3[4,1]=sd(b1.is)
b1var.table3[4,2]=sqrt(mean(mb.varb1,na.rm=TRUE))
b1var.table3[4,3]=sqrt(mean(ismb.varb1,na.rm=TRUE))
b1var.table3[4,4]=mean(b1.crq.sd)

# Coverage table
b0.coverage[4,4]=mean(cover[,4],na.rm=TRUE)
b1.coverage[4,4]=mean(cover[,8],na.rm=TRUE)


