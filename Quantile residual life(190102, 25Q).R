library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)

#### Given information for quantile 25% ####
c.0=5000000
c.1=53.03
c.2=26.58
c.3=17.73
exp.beta.initial=5
k=2
r.initial=(log(4/3))^(1/k)/exp.beta.initial
true.b0t0=log(5)
true.b0t1=log((1/r.initial)*((log(4/3)+(r.initial*1)^k)^(1/k))-1)
true.b0t2=log((1/r.initial)*((log(4/3)+(r.initial*2)^k)^(1/k))-2)
true.b0t3=log((1/r.initial)*((log(4/3)+(r.initial*3)^k)^(1/k))-3)

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
  result=t(w*X)%*%(pnorm((a[,5]-X%*%beta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
  print(result)
}

#### revised object equation (with eta) ####
rev.objectF=function(beta){
  beta=as.matrix(beta)
  result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%beta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
}

#### Beta estimation ####
#### method 1 : Crq function (??�� package) ####
# t_0=0
table0.crq<-matrix(NA,4,6)

rownames(table0.crq)<-c(0,10,20,30)
colnames(table0.crq)<-c("beta_0","SE of beta_0(empirical)","SE of beta_0(avg)","beta_1","SE of beta_1(empirical)","SE of beta_1(avg)")
dummy0<-c()
dummy0.sd<-c()
dummy1<-c()
dummy1.sd<-c()

# t_0=0 & c=0%
for (i in 1:500){
  sim0<-data.gen(c.0,0)
  sim0[,2]<-100
  fit0<-crq(Surv(sim0[,5],sim0[,7])~sim0[,6],method='Portnoy')
  dummy0[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table0.crq[1,1]<-mean(dummy0)
table0.crq[1,2]<-sd(dummy0)
table0.crq[1,3]<-mean(dummy0.sd)
table0.crq[1,4]<-mean(dummy1)
table0.crq[1,5]<-sd(dummy1)
table0.crq[1,6]<-mean(dummy1.sd)

# t_0=0 & c=10%
for (i in 1:500){
  sim1<-data.gen(c.1,0)
  fit1<-crq(Surv(sim1[,5],sim1[,7])~sim1[,6],method='Portnoy')
  dummy0[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table0.crq[2,1]<-mean(dummy0)
table0.crq[2,2]<-sd(dummy0)
table0.crq[2,3]<-mean(dummy0.sd)
table0.crq[2,4]<-mean(dummy1)
table0.crq[2,5]<-sd(dummy1)
table0.crq[2,6]<-mean(dummy1.sd)

# t_0=0 & c=20%
for (i in 1:500){
  sim2<-data.gen(c.2,0)
  fit2<-crq(Surv(sim2[,5],sim2[,7])~sim2[,6],method='Portnoy')
  dummy0[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table0.crq[3,1]<-mean(dummy0)
table0.crq[3,2]<-sd(dummy0)
table0.crq[3,3]<-mean(dummy0.sd)
table0.crq[3,4]<-mean(dummy1)
table0.crq[3,5]<-sd(dummy1)
table0.crq[3,6]<-mean(dummy1.sd)

# t_0=0 & c=30%
for (i in 1:500){
  sim3<-data.gen(c.3,0)
  fit3<-crq(Surv(sim3[,5],sim3[,7])~sim3[,6],method='Portnoy')
  dummy0[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table0.crq[4,1]<-mean(dummy0)
table0.crq[4,2]<-sd(dummy0)
table0.crq[4,3]<-mean(dummy0.sd)
table0.crq[4,4]<-mean(dummy1)
table0.crq[4,5]<-sd(dummy1)
table0.crq[4,6]<-mean(dummy1.sd)

table0.crq

# t_0=1
table1.crq<-matrix(NA,4,6)
rownames(table1.crq)<-c(0,10,20,30)
colnames(table1.crq)<-c("beta_0","SE of beta_0(empirical)","SE of beta_0(avg)","beta_1","SE of beta_1(empirical)","SE of beta_1(avg)")
dummy0<-c()
dummy0.sd<-c()
dummy1<-c()
dummy1.sd<-c()

# t_0=1 & c=0%
for (i in 1:500){
  sim0<-data.gen(c.0,1)
  sim0[,2]<-100
  fit0<-crq(Surv(sim0[,5],sim0[,7])~sim0[,6],method='Portnoy')
  dummy0[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table1.crq[1,1]<-mean(dummy0)
table1.crq[1,2]<-sd(dummy0)
table1.crq[1,3]<-mean(dummy0.sd)
table1.crq[1,4]<-mean(dummy1)
table1.crq[1,5]<-sd(dummy1)
table1.crq[1,6]<-mean(dummy1.sd)

# t_0=1 & c=10%
for (i in 1:500){
  sim1<-data.gen(c.1,1)
  fit1<-crq(Surv(sim1[,5],sim1[,7])~sim1[,6],method='Portnoy')
  dummy0[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table1.crq[2,1]<-mean(dummy0)
table1.crq[2,2]<-sd(dummy0)
table1.crq[2,3]<-mean(dummy0.sd)
table1.crq[2,4]<-mean(dummy1)
table1.crq[2,5]<-sd(dummy1)
table1.crq[2,6]<-mean(dummy1.sd)

# t_0=1 & c=20%
for (i in 1:500){
  sim2<-data.gen(c.2,1)
  fit2<-crq(Surv(sim2[,5],sim2[,7])~sim2[,6],method='Portnoy')
  dummy0[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table1.crq[3,1]<-mean(dummy0)
table1.crq[3,2]<-sd(dummy0)
table1.crq[3,3]<-mean(dummy0.sd)
table1.crq[3,4]<-mean(dummy1)
table1.crq[3,5]<-sd(dummy1)
table1.crq[3,6]<-mean(dummy1.sd)

# t_0=1 & c=30%
for (i in 1:500){
  sim3<-data.gen(c.3,1)
  fit3<-crq(Surv(sim3[,5],sim3[,7])~sim3[,6],method='Portnoy')
  dummy0[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table1.crq[4,1]<-mean(dummy0)
table1.crq[4,2]<-sd(dummy0)
table1.crq[4,3]<-mean(dummy0.sd)
table1.crq[4,4]<-mean(dummy1)
table1.crq[4,5]<-sd(dummy1)
table1.crq[4,6]<-mean(dummy1.sd)

table1.crq

# t_0=2
table2.crq<-matrix(NA,4,6)
rownames(table2.crq)<-c(0,10,20,30)
colnames(table2.crq)<-c("beta_0","SE of beta_0(empirical)","SE of beta_0(avg)","beta_1","SE of beta_1(empirical)","SE of beta_1(avg)")
dummy0<-c()
dummy0.sd<-c()
dummy1<-c()
dummy1.sd<-c()

# t_0=2 & c=0%
for (i in 1:500){
  sim0<-data.gen(c.0,2)
  sim0[,2]<-100
  fit0<-crq(Surv(sim0[,5],sim0[,7])~sim0[,6],method='Portnoy')
  dummy0[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table2.crq[1,1]<-mean(dummy0)
table2.crq[1,2]<-sd(dummy0)
table2.crq[1,3]<-mean(dummy0.sd)
table2.crq[1,4]<-mean(dummy1)
table2.crq[1,5]<-sd(dummy1)
table2.crq[1,6]<-mean(dummy1.sd)

# t_0=2 & c=10%
for (i in 1:500){
  sim1<-data.gen(c.1,2)
  fit1<-crq(Surv(sim1[,5],sim1[,7])~sim1[,6],method='Portnoy')
  dummy0[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table2.crq[2,1]<-mean(dummy0)
table2.crq[2,2]<-sd(dummy0)
table2.crq[2,3]<-mean(dummy0.sd)
table2.crq[2,4]<-mean(dummy1)
table2.crq[2,5]<-sd(dummy1)
table2.crq[2,6]<-mean(dummy1.sd)

# t_0=2 & c=20%
for (i in 1:500){
  sim2<-data.gen(c.2,2)
  fit2<-crq(Surv(sim2[,5],sim2[,7])~sim2[,6],method='Portnoy')
  dummy0[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table2.crq[3,1]<-mean(dummy0)
table2.crq[3,2]<-sd(dummy0)
table2.crq[3,3]<-mean(dummy0.sd)
table2.crq[3,4]<-mean(dummy1)
table2.crq[3,5]<-sd(dummy1)
table2.crq[3,6]<-mean(dummy1.sd)

# t_0=2 & c=30%
for (i in 1:500){
  sim3<-data.gen(c.3,2)
  fit3<-crq(Surv(sim3[,5],sim3[,7])~sim3[,6],method='Portnoy')
  dummy0[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table2.crq[4,1]<-mean(dummy0)
table2.crq[4,2]<-sd(dummy0)
table2.crq[4,3]<-mean(dummy0.sd)
table2.crq[4,4]<-mean(dummy1)
table2.crq[4,5]<-sd(dummy1)
table2.crq[4,6]<-mean(dummy1.sd)

table2.crq

# t_0=3
table3.crq<-matrix(NA,4,6)
rownames(table3.crq)<-c(0,10,20,30)
colnames(table3.crq)<-c("beta_0","SE of beta_0(empirical)","SE of beta_0(avg)","beta_1","SE of beta_1(empirical)","SE of beta_1(avg)")
dummy0<-c()
dummy0.sd<-c()
dummy1<-c()
dummy1.sd<-c()

# t_0=3 & c=0%
for (i in 1:500){
  sim0<-data.gen(c.0,3)
  sim0[,2]<-100
  fit0<-crq(Surv(sim0[,5],sim0[,7])~sim0[,6],method='Portnoy')
  dummy0[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table3.crq[1,1]<-mean(dummy0)
table3.crq[1,2]<-sd(dummy0)
table3.crq[1,3]<-mean(dummy0.sd)
table3.crq[1,4]<-mean(dummy1)
table3.crq[1,5]<-sd(dummy1)
table3.crq[1,6]<-mean(dummy1.sd)

# t_0=3 & c=10%
for (i in 1:500){
  sim1<-data.gen(c.1,3)
  fit1<-crq(Surv(sim1[,5],sim1[,7])~sim1[,6],method='Portnoy')
  dummy0[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]
}
table3.crq[2,1]<-mean(dummy0)
table3.crq[2,2]<-sd(dummy0)
table3.crq[2,3]<-mean(dummy0.sd)
table3.crq[2,4]<-mean(dummy1)
table3.crq[2,5]<-sd(dummy1)
table3.crq[2,6]<-mean(dummy1.sd)

# t_0=3 & c=20%
for (i in 1:500){
  sim2<-data.gen(c.2,3)
  fit2<-crq(Surv(sim2[,5],sim2[,7])~sim2[,6],method='Portnoy')
  dummy0[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table3.crq[3,1]<-mean(dummy0)
table3.crq[3,2]<-sd(dummy0)
table3.crq[3,3]<-mean(dummy0.sd)
table3.crq[3,4]<-mean(dummy1)
table3.crq[3,5]<-sd(dummy1)
table3.crq[3,6]<-mean(dummy1.sd)

# t_0=3 & c=30%
for (i in 1:500){
  sim3<-data.gen(c.3,3)
  fit3<-crq(Surv(sim3[,5],sim3[,7])~sim3[,6],method='Portnoy')
  dummy0[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.10,0.25))[2][[1]]$coefficient[2,4]  
}
table3.crq[4,1]<-mean(dummy0)
table3.crq[4,2]<-sd(dummy0)
table3.crq[4,3]<-mean(dummy0.sd)
table3.crq[4,4]<-mean(dummy1)
table3.crq[4,5]<-sd(dummy1)
table3.crq[4,6]<-mean(dummy1.sd)

table3.crq

#### method 2 : Induced smooting + objective function (Our suggested method) ####
# t_0=0
table0<-matrix(NA,4,4)
rownames(table0)<-c(0,10,20,30)
colnames(table0)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
# true betas are (1.61, 0)

# t_0=0 & C=0%
b0.t0c0<-c()
b1.t0c0<-c()
for (m in 1:500){
  a<-data.gen(c.0,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c0<-nleqslv(betastart,objectF)
  b0.t0c0[m]<-z.t0c0$x[1]
  b1.t0c0[m]<-z.t0c0$x[2]
}

# find solution by nleqslv
table0[1,1]<-mean(b0.t0c0)
table0[1,2]<-sd(b0.t0c0)
table0[1,3]<-mean(b1.t0c0)
table0[1,4]<-sd(b1.t0c0)

# t_0=0 & C=10%
b0.t0c1<-c()
b1.t0c1<-c()
for (m in 1:500){
  a<-data.gen(c.1,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c1<-nleqslv(betastart,objectF)
  b0.t0c1[m]<-z.t0c1$x[1]
  b1.t0c1[m]<-z.t0c1$x[2]
}

# find solution by nleqslv
table0[2,1]<-mean(b0.t0c1)
table0[2,2]<-sd(b0.t0c1)
table0[2,3]<-mean(b1.t0c1)
table0[2,4]<-sd(b1.t0c1)

# t_0=0 & C=20
b0.t0c2<-c()
b1.t0c2<-c()
for (m in 1:500){
  a<-data.gen(c.2,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c2<-nleqslv(betastart,objectF)
  b0.t0c2[m]<-z.t0c2$x[1]
  b1.t0c2[m]<-z.t0c2$x[2]
}

# find solution by nleqslv
table0[3,1]<-mean(b0.t0c2)
table0[3,2]<-sd(b0.t0c2)
table0[3,3]<-mean(b1.t0c2)
table0[3,4]<-sd(b1.t0c2)

# t_0=0 & C=30
b0.t0c3<-c()
b1.t0c3<-c()
for (m in 1:500){
  a<-data.gen(c.3,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c3<-nleqslv(betastart,objectF)
  b0.t0c3[m]<-z.t0c3$x[1]
  b1.t0c3[m]<-z.t0c3$x[2]
}

# find solution by nleqslv
table0[4,1]<-mean(b0.t0c3)
table0[4,2]<-sd(b0.t0c3)
table0[4,3]<-mean(b1.t0c3)
table0[4,4]<-sd(b1.t0c3)

# t_0=1
table1<-matrix(NA,4,4)
rownames(table1)<-c(0,10,20,30)
colnames(table1)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

# true betas are (1.41, 0)
b0.t1c0<-c()
b1.t1c0<-c()
for (m in 1:500){
  a<-data.gen(c.0,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c0<-nleqslv(betastart,objectF)
  b0.t1c0[m]<-z.t1c0$x[1]
  b1.t1c0[m]<-z.t1c0$x[2]
}

# find solution by nleqslv
table1[1,1]<-mean(b0.t1c0)
table1[1,2]<-sd(b0.t1c0)
table1[1,3]<-mean(b1.t1c0)
table1[1,4]<-sd(b1.t1c0)

# t_0=1 & C=10%
b0.t1c1<-c()
b1.t1c1<-c()
for (m in 1:500){
  a<-data.gen(c.1,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c1<-nleqslv(betastart,objectF)
  b0.t1c1[m]<-z.t1c1$x[1]
  b1.t1c1[m]<-z.t1c1$x[2]
}

# find solution by nleqslv
table1[2,1]<-mean(b0.t1c1)
table1[2,2]<-sd(b0.t1c1)
table1[2,3]<-mean(b1.t1c1)
table1[2,4]<-sd(b1.t1c1)

# t_0=1 & C=20
b0.t1c2<-c()
b1.t1c2<-c()
for (m in 1:500){
  a<-data.gen(c.2,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c2<-nleqslv(betastart,objectF)
  b0.t1c2[m]<-z.t1c2$x[1]
  b1.t1c2[m]<-z.t1c2$x[2]
}

# find solution by nleqslv
table1[3,1]<-mean(b0.t1c2)
table1[3,2]<-sd(b0.t1c2)
table1[3,3]<-mean(b1.t1c2)
table1[3,4]<-sd(b1.t1c2)

# t_0=1 & C=30
b0.t1c3<-c()
b1.t1c3<-c()
for (m in 1:500){
  a<-data.gen(c.3,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c3<-nleqslv(betastart,objectF)
  b0.t1c3[m]<-z.t1c3$x[1]
  b1.t1c3[m]<-z.t1c3$x[2]
}

# find solution by nleqslv
table1[4,1]<-mean(b0.t1c3)
table1[4,2]<-sd(b0.t1c3)
table1[4,3]<-mean(b1.t1c3)
table1[4,4]<-sd(b1.t1c3)

# t_0=2
table2<-matrix(NA,4,4)
rownames(table2)<-c(0,10,20,30)
colnames(table2)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

# true betas are (1.22, 0)
# t_0=2 & C=0%
b0.t2c0<-c()
b1.t2c0<-c()
for (m in 1:500){
  a<-data.gen(c.0,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c0<-nleqslv(betastart,objectF)
  b0.t2c0[m]<-z.t2c0$x[1]
  b1.t2c0[m]<-z.t2c0$x[2]
}

# find solution by nleqslv
table2[1,1]<-mean(b0.t2c0)
table2[1,2]<-sd(b0.t2c0)
table2[1,3]<-mean(b1.t2c0)
table2[1,4]<-sd(b1.t2c0)

# t_0=2 & C=10%
b0.t2c1<-c()
b1.t2c1<-c()
for (m in 1:500){
  a<-data.gen(c.1,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c1<-nleqslv(betastart,objectF)
  b0.t2c1[m]<-z.t2c1$x[1]
  b1.t2c1[m]<-z.t2c1$x[2]
}

# find solution by nleqslv
table2[2,1]<-mean(b0.t2c1)
table2[2,2]<-sd(b0.t2c1)
table2[2,3]<-mean(b1.t2c1)
table2[2,4]<-sd(b1.t2c1)

# t_0=2 & C=20
b0.t2c2<-c()
b1.t2c2<-c()
for (m in 1:500){
  a<-data.gen(c.2,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c2<-nleqslv(betastart,objectF)
  b0.t2c2[m]<-z.t2c2$x[1]
  b1.t2c2[m]<-z.t2c2$x[2]
}

# find solution by nleqslv
table2[3,1]<-mean(b0.t2c2)
table2[3,2]<-sd(b0.t2c2)
table2[3,3]<-mean(b1.t2c2)
table2[3,4]<-sd(b1.t2c2)

# t_0=2 & C=30
b0.t2c3<-c()
b1.t2c3<-c()
for (m in 1:500){
  a<-data.gen(c.3,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c3<-nleqslv(betastart,objectF)
  b0.t2c3[m]<-z.t2c3$x[1]
  b1.t2c3[m]<-z.t2c3$x[2]
}

# find solution by nleqslv
table2[4,1]<-mean(b0.t2c3)
table2[4,2]<-sd(b0.t2c3)
table2[4,3]<-mean(b1.t2c3)
table2[4,4]<-sd(b1.t2c3)

# t_0=3
table3<-matrix(NA,4,4)
rownames(table3)<-c(0,10,20,30)
colnames(table3)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

# true betas are (1.04, 0)
# t_0=3 & C=0%
b0.t3c0<-c()
b1.t3c0<-c()
for (m in 1:500){
  a<-data.gen(c.0,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c0<-nleqslv(betastart,objectF)
  b0.t3c0[m]<-z.t3c0$x[1]
  b1.t3c0[m]<-z.t3c0$x[2]
}

# find solution by nleqslv
table3[1,1]<-mean(b0.t3c0)
table3[1,2]<-sd(b0.t3c0)
table3[1,3]<-mean(b1.t3c0)
table3[1,4]<-sd(b1.t3c0)

# t_0=3 & C=10%
b0.t3c1<-c()
b1.t3c1<-c()
for (m in 1:500){
  a<-data.gen(c.1,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c1<-nleqslv(betastart,objectF)
  b0.t3c1[m]<-z.t3c1$x[1]
  b1.t3c1[m]<-z.t3c1$x[2]
}

# find solution by nleqslv
table3[2,1]<-mean(b0.t3c1)
table3[2,2]<-sd(b0.t3c1)
table3[2,3]<-mean(b1.t3c1)
table3[2,4]<-sd(b1.t3c1)

# t_0=3 & C=20
b0.t3c2<-c()
b1.t3c2<-c()
for (m in 1:500){
  a<-data.gen(c.2,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c2<-nleqslv(betastart,objectF)
  b0.t3c2[m]<-z.t3c2$x[1]
  b1.t3c2[m]<-z.t3c2$x[2]
}

# find solution by nleqslv
table3[3,1]<-mean(b0.t3c2)
table3[3,2]<-sd(b0.t3c2)
table3[3,3]<-mean(b1.t3c2)
table3[3,4]<-sd(b1.t3c2)

# t_0=3 & C=30
b0.t3c3<-c()
b1.t3c3<-c()
for (m in 1:500){
  a<-data.gen(c.3,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c3<-nleqslv(betastart,objectF)
  b0.t3c3[m]<-z.t3c3$x[1]
  b1.t3c3[m]<-z.t3c3$x[2]
}

# find solution by nleqslv
table3[4,1]<-mean(b0.t3c3)
table3[4,2]<-sd(b0.t3c3)
table3[4,3]<-mean(b1.t3c3)
table3[4,4]<-sd(b1.t3c3)

#### Variance estimation ####
b0var.table0=matrix(NA,4,3)
colnames(b0var.table0)=c("true","MB","ISMB(simple)")
rownames(b0var.table0)=c(0,10,20,30)
b0var.table1=matrix(NA,4,3)
colnames(b0var.table1)=c("true","MB","ISMB(simple)")
rownames(b0var.table1)=c(0,10,20,30)
b0var.table2=matrix(NA,4,3)
colnames(b0var.table2)=c("true","MB","ISMB(simple)")
rownames(b0var.table2)=c(0,10,20,30)
b0var.table3=matrix(NA,4,3)
colnames(b0var.table3)=c("true","MB","ISMB(simple)")
rownames(b0var.table3)=c(0,10,20,30)

b1var.table0=matrix(NA,4,3)
colnames(b1var.table0)=c("true","MB","ISMB(simple)")
rownames(b1var.table0)=c(0,10,20,30)
b1var.table1=matrix(NA,4,3)
colnames(b1var.table1)=c("true","MB","ISMB(simple)")
rownames(b1var.table1)=c(0,10,20,30)
b1var.table2=matrix(NA,4,3)
colnames(b1var.table2)=c("true","MB","ISMB(simple)")
rownames(b1var.table2)=c(0,10,20,30)
b1var.table3=matrix(NA,4,3)
colnames(b1var.table3)=c("true","MB","ISMB(simple)")
rownames(b1var.table3)=c(0,10,20,30)

#### Method 1 : True variance by beta estimation ####
# Beta 0
b0var.table0[1,1]<-sd(b0.t0c0)
b0var.table0[2,1]<-sd(b0.t0c1)
b0var.table0[3,1]<-sd(b0.t0c2)
b0var.table0[4,1]<-sd(b0.t0c3)
b0var.table1[1,1]<-sd(b0.t1c0)
b0var.table1[2,1]<-sd(b0.t1c1)
b0var.table1[3,1]<-sd(b0.t1c2)
b0var.table1[4,1]<-sd(b0.t1c3)
b0var.table2[1,1]<-sd(b0.t2c0)
b0var.table2[2,1]<-sd(b0.t2c1)
b0var.table2[3,1]<-sd(b0.t2c2)
b0var.table2[4,1]<-sd(b0.t2c3)
b0var.table3[1,1]<-sd(b0.t3c0)
b0var.table3[2,1]<-sd(b0.t3c1)
b0var.table3[3,1]<-sd(b0.t3c2)
b0var.table3[4,1]<-sd(b0.t3c3)

# Beta 1
b1var.table0[1,1]<-sd(b1.t0c0)
b1var.table0[2,1]<-sd(b1.t0c1)
b1var.table0[3,1]<-sd(b1.t0c2)
b1var.table0[4,1]<-sd(b1.t0c3)
b1var.table1[1,1]<-sd(b1.t1c0)
b1var.table1[2,1]<-sd(b1.t1c1)
b1var.table1[3,1]<-sd(b1.t1c2)
b1var.table1[4,1]<-sd(b1.t1c3)
b1var.table2[1,1]<-sd(b1.t2c0)
b1var.table2[2,1]<-sd(b1.t2c1)
b1var.table2[3,1]<-sd(b1.t2c2)
b1var.table2[4,1]<-sd(b1.t2c3)
b1var.table3[1,1]<-sd(b1.t3c0)
b1var.table3[2,1]<-sd(b1.t3c1)
b1var.table3[3,1]<-sd(b1.t3c2)
b1var.table3[4,1]<-sd(b1.t3c3)


#### method 2 : Multiplier Bootstrap ####
# t_0=0
# t_0=0 & C=0%
tic()
vart0c0.b0=c()
vart0c0.b1=c()
for (i in 1:500){
  a=data.gen(c.0,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.00.b0=c()
  mb.00.b1=c()
  result.mb.00.b0=c()
  result.mb.00.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c0=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.00.b0[j]=z.c0$x[1]
    mb.00.b1[j]=z.c0$x[2]
    result.mb.00.b0[j]=z.c0$fvec[1]
    result.mb.00.b1[j]=z.c0$fvec[2]
  }
  mb.00.b0=mb.00.b0[-1e-5<result.mb.00.b0&result.mb.00.b0<1e-5]
  mb.00.b1=mb.00.b1[-1e-5<result.mb.00.b1&result.mb.00.b1<1e-5]
  vart0c0.b0[i]=var(mb.00.b0)
  vart0c0.b1[i]=var(mb.00.b1)
}
b0var.table0[1,2]=sqrt(mean(vart0c0.b0))
b1var.table0[1,2]=sqrt(mean(vart0c0.b1))
toc()

# t_0=0 & C=10%
vart0c1.b0=c()
vart0c1.b1=c()
for (i in 1:500){
  a=data.gen(c.1,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.01.b0=c()
  mb.01.b1=c()
  result.mb.01.b0=c()
  result.mb.01.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c1=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.01.b0[j]=z.c1$x[1]
    mb.01.b1[j]=z.c1$x[2]
    result.mb.01.b0[j]=z.c1$fvec[1]
    result.mb.01.b1[j]=z.c1$fvec[2]
  }
  mb.01.b0=mb.01.b0[-1e-5<result.mb.01.b0&result.mb.01.b0<1e-5]
  mb.01.b1=mb.01.b1[-1e-5<result.mb.01.b1&result.mb.01.b1<1e-5]
  vart0c1.b0[i]=var(mb.01.b0)
  vart0c1.b1[i]=var(mb.01.b1)
}
b0var.table0[2,2]=sqrt(mean(vart0c1.b0))
b1var.table0[2,2]=sqrt(mean(vart0c1.b1))

# t_0=0 & C=20%
vart0c2.b0=c()
vart0c2.b1=c()
for (i in 1:500){
  a=data.gen(c.2,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.02.b0=c()
  mb.02.b1=c()
  result.mb.02.b0=c()
  result.mb.02.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c2=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.02.b0[j]=z.c2$x[1]
    mb.02.b1[j]=z.c2$x[2]
    result.mb.02.b0[j]=z.c2$fvec[1]
    result.mb.02.b1[j]=z.c2$fvec[2]
  }
  mb.02.b0=mb.02.b0[-1e-5<result.mb.02.b0&result.mb.02.b0<1e-5]
  mb.02.b1=mb.02.b1[-1e-5<result.mb.02.b1&result.mb.02.b1<1e-5]
  vart0c2.b0[i]=var(mb.02.b0)
  vart0c2.b1[i]=var(mb.02.b1)
}
b0var.table0[3,2]=sqrt(mean(vart0c2.b0))
b1var.table0[3,2]=sqrt(mean(vart0c2.b1))

# t_0=0 & C=30%
vart0c3.b0=c()
vart0c3.b1=c()
for (i in 1:500){
  a=data.gen(c.3,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.03.b0=c()
  mb.03.b1=c()
  result.mb.03.b0=c()
  result.mb.03.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c3=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.03.b0[j]=z.c3$x[1]
    mb.03.b1[j]=z.c3$x[2]
    result.mb.03.b0[j]=z.c3$fvec[1]
    result.mb.03.b1[j]=z.c3$fvec[2]
  }
  mb.03.b0=mb.03.b0[-1e-5<result.mb.03.b0&result.mb.03.b0<1e-5]
  mb.03.b1=mb.03.b1[-1e-5<result.mb.03.b1&result.mb.03.b1<1e-5]
  vart0c3.b0[i]=var(mb.03.b0)
  vart0c3.b1[i]=var(mb.03.b1)
}
b0var.table0[4,2]=sqrt(mean(vart0c3.b0))
b1var.table0[4,2]=sqrt(mean(vart0c3.b1))

# t_0=1
# t_0=1 & C=0%
vart1c0.b0=c()
vart1c0.b1=c()
for (i in 1:500){
  a=data.gen(c.0,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.10.b0=c()
  mb.10.b1=c()
  result.mb.10.b0=c()
  result.mb.10.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.44,0)
    z.c0=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.10.b0[j]=z.c0$x[1]
    mb.10.b1[j]=z.c0$x[2]
    result.mb.10.b0[j]=z.c0$fvec[1]
    result.mb.10.b1[j]=z.c0$fvec[2]
  }
  mb.10.b0=mb.10.b0[-1e-5<result.mb.10.b0&result.mb.10.b0<1e-5]
  mb.10.b1=mb.10.b1[-1e-5<result.mb.10.b1&result.mb.10.b1<1e-5]
  vart1c0.b0[i]=var(mb.10.b0)
  vart1c0.b1[i]=var(mb.10.b1)
}
b0var.table1[1,2]=sqrt(mean(vart1c0.b0))
b1var.table1[1,2]=sqrt(mean(vart1c0.b1))

# t_0=1 & C=10%
vart1c1.b0=c()
vart1c1.b1=c()
for (i in 1:500){
  a=data.gen(c.1,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.11.b0=c()
  mb.11.b1=c()
  result.mb.11.b0=c()
  result.mb.11.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.44,0)
    z.c1=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.11.b0[j]=z.c1$x[1]
    mb.11.b1[j]=z.c1$x[2]
    result.mb.11.b0[j]=z.c1$fvec[1]
    result.mb.11.b1[j]=z.c1$fvec[2]
  }
  mb.11.b0=mb.11.b0[-1e-5<result.mb.11.b0&result.mb.11.b0<1e-5]
  mb.11.b1=mb.11.b1[-1e-5<result.mb.11.b1&result.mb.11.b1<1e-5]
  vart1c1.b0[i]=var(mb.11.b0)
  vart1c1.b1[i]=var(mb.11.b1)
}
b0var.table1[2,2]=sqrt(mean(vart1c1.b0))
b1var.table1[2,2]=sqrt(mean(vart1c1.b1))

# t_0=1 & C=20%
vart1c2.b0=c()
vart1c2.b1=c()
for (i in 1:500){
  a=data.gen(c.2,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.12.b0=c()
  mb.12.b1=c()
  result.mb.12.b0=c()
  result.mb.12.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.44,0)
    z.c2=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.12.b0[j]=z.c2$x[1]
    mb.12.b1[j]=z.c2$x[2]
    result.mb.12.b0[j]=z.c2$fvec[1]
    result.mb.12.b1[j]=z.c2$fvec[2]
  }
  mb.12.b0=mb.12.b0[-1e-5<result.mb.12.b0&result.mb.12.b0<1e-5]
  mb.12.b1=mb.12.b1[-1e-5<result.mb.12.b1&result.mb.12.b1<1e-5]
  vart1c2.b0[i]=var(mb.12.b0)
  vart1c2.b1[i]=var(mb.12.b1)
}
b0var.table1[3,2]=sqrt(mean(vart1c2.b0))
b1var.table1[3,2]=sqrt(mean(vart1c2.b1))

# t_0=1 & C=30%
vart1c3.b0=c()
vart1c3.b1=c()
for (i in 1:500){
  a=data.gen(c.3,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.13.b0=c()
  mb.13.b1=c()
  result.mb.13.b0=c()
  result.mb.13.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.44,0)
    z.c3=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.13.b0[j]=z.c3$x[1]
    mb.13.b1[j]=z.c3$x[2]
    result.mb.13.b0[j]=z.c3$fvec[1]
    result.mb.13.b1[j]=z.c3$fvec[2]
  }
  mb.13.b0=mb.13.b0[-1e-5<result.mb.13.b0&result.mb.13.b0<1e-5]
  mb.13.b1=mb.13.b1[-1e-5<result.mb.13.b1&result.mb.13.b1<1e-5]
  vart1c3.b0[i]=var(mb.13.b0)
  vart1c3.b1[i]=var(mb.13.b1)
}
b0var.table1[4,2]=sqrt(mean(vart1c3.b0))
b1var.table1[4,2]=sqrt(mean(vart1c3.b1))

# t_0=2
# t_0=2 & C=0%
vart2c0.b0=c()
vart2c0.b1=c()
for (i in 1:500){
  a=data.gen(c.0,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.20.b0=c()
  mb.20.b1=c()
  result.mb.20.b0=c()
  result.mb.20.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c0=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.20.b0[j]=z.c0$x[1]
    mb.20.b1[j]=z.c0$x[2]
    result.mb.20.b0[j]=z.c0$fvec[1]
    result.mb.20.b1[j]=z.c0$fvec[2]
  }
  mb.20.b0=mb.20.b0[-1e-5<result.mb.20.b0&result.mb.20.b0<1e-5]
  mb.20.b1=mb.20.b1[-1e-5<result.mb.20.b1&result.mb.20.b1<1e-5]
  vart2c0.b0[i]=var(mb.20.b0)
  vart2c0.b1[i]=var(mb.20.b1)
}
b0var.table2[1,2]=sqrt(mean(vart2c0.b0))
b1var.table2[1,2]=sqrt(mean(vart2c0.b1))

# t_0=2 & C=10%
vart2c1.b0=c()
vart2c1.b1=c()
for (i in 1:500){
  a=data.gen(c.1,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.21.b0=c()
  mb.21.b1=c()
  result.mb.21.b0=c()
  result.mb.21.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c1=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.21.b0[j]=z.c1$x[1]
    mb.21.b1[j]=z.c1$x[2]
    result.mb.21.b0[j]=z.c1$fvec[1]
    result.mb.21.b1[j]=z.c1$fvec[2]
  }
  mb.21.b0=mb.21.b0[-1e-5<result.mb.21.b0&result.mb.21.b0<1e-5]
  mb.21.b1=mb.21.b1[-1e-5<result.mb.21.b1&result.mb.21.b1<1e-5]
  vart2c1.b0[i]=var(mb.21.b0)
  vart2c1.b1[i]=var(mb.21.b1)
}
b0var.table2[2,2]=sqrt(mean(vart2c1.b0))
b1var.table2[2,2]=sqrt(mean(vart2c1.b1))

# t_0=2 & C=20%
vart2c2.b0=c()
vart2c2.b1=c()
for (i in 1:500){
  a=data.gen(c.2,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.22.b0=c()
  mb.22.b1=c()
  result.mb.22.b0=c()
  result.mb.22.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c2=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.22.b0[j]=z.c2$x[1]
    mb.22.b1[j]=z.c2$x[2]
    result.mb.22.b0[j]=z.c2$fvec[1]
    result.mb.22.b1[j]=z.c2$fvec[2]
  }
  mb.22.b0=mb.22.b0[-1e-5<result.mb.22.b0&result.mb.22.b0<1e-5]
  mb.22.b1=mb.22.b1[-1e-5<result.mb.22.b1&result.mb.22.b1<1e-5]
  vart2c2.b0[i]=var(mb.22.b0)
  vart2c2.b1[i]=var(mb.22.b1)
}
b0var.table2[3,2]=sqrt(mean(vart2c2.b0))
b1var.table2[3,2]=sqrt(mean(vart2c2.b1))

# t_0=2 & C=30%
vart2c3.b0=c()
vart2c3.b1=c()
for (i in 1:500){
  a=data.gen(c.3,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.23.b0=c()
  mb.23.b1=c()
  result.mb.23.b0=c()
  result.mb.23.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c3=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.23.b0[j]=z.c3$x[1]
    mb.23.b1[j]=z.c3$x[2]
    result.mb.23.b0[j]=z.c3$fvec[1]
    result.mb.23.b1[j]=z.c3$fvec[2]
  }
  mb.23.b0=mb.23.b0[-1e-5<result.mb.23.b0&result.mb.23.b0<1e-5]
  mb.23.b1=mb.23.b1[-1e-5<result.mb.23.b1&result.mb.23.b1<1e-5]
  vart2c3.b0[i]=var(mb.23.b0)
  vart2c3.b1[i]=var(mb.23.b1)
}
b0var.table2[4,2]=sqrt(mean(vart2c3.b0))
b1var.table2[4,2]=sqrt(mean(vart2c3.b1))

# t_0=3
# t_0=3 & C=0%
vart3c0.b0=c()
vart3c0.b1=c()
for (i in 1:500){
  a=data.gen(c.0,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.30.b0=c()
  mb.30.b1=c()
  result.mb.30.b0=c()
  result.mb.30.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c0=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.30.b0[j]=z.c0$x[1]
    mb.30.b1[j]=z.c0$x[2]
    result.mb.30.b0[j]=z.c0$fvec[1]
    result.mb.30.b1[j]=z.c0$fvec[2]
  }
  mb.30.b0=mb.30.b0[-1e-5<result.mb.30.b0&result.mb.30.b0<1e-5]
  mb.30.b1=mb.30.b1[-1e-5<result.mb.30.b1&result.mb.30.b1<1e-5]
  vart3c0.b0[i]=var(mb.30.b0)
  vart3c0.b1[i]=var(mb.30.b1)
}
b0var.table3[1,2]=sqrt(mean(vart3c0.b0))
b1var.table3[1,2]=sqrt(mean(vart3c0.b1))

# t_0=3 & C=10%
vart3c1.b0=c()
vart3c1.b1=c()
for (i in 1:500){
  a=data.gen(c.1,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.31.b0=c()
  mb.31.b1=c()
  result.mb.31.b0=c()
  result.mb.31.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c1=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.31.b0[j]=z.c1$x[1]
    mb.31.b1[j]=z.c1$x[2]
    result.mb.31.b0[j]=z.c1$fvec[1]
    result.mb.31.b1[j]=z.c1$fvec[2]
  }
  mb.31.b0=mb.31.b0[-1e-5<result.mb.31.b0&result.mb.31.b0<1e-5]
  mb.31.b1=mb.31.b1[-1e-5<result.mb.31.b1&result.mb.31.b1<1e-5]
  vart3c1.b0[i]=var(mb.31.b0)
  vart3c1.b1[i]=var(mb.31.b1)
}
b0var.table3[2,2]=sqrt(mean(vart3c1.b0))
b1var.table3[2,2]=sqrt(mean(vart3c1.b1))

# t_0=3 & C=20%
vart3c2.b0=c()
vart3c2.b1=c()
for (i in 1:500){
  a=data.gen(c.2,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.32.b0=c()
  mb.32.b1=c()
  result.mb.32.b0=c()
  result.mb.32.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c2=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.32.b0[j]=z.c2$x[1]
    mb.32.b1[j]=z.c2$x[2]
    result.mb.32.b0[j]=z.c2$fvec[1]
    result.mb.32.b1[j]=z.c2$fvec[2]
  }
  mb.32.b0=mb.32.b0[-1e-5<result.mb.32.b0&result.mb.32.b0<1e-5]
  mb.32.b1=mb.32.b1[-1e-5<result.mb.32.b1&result.mb.32.b1<1e-5]
  vart3c2.b0[i]=var(mb.32.b0)
  vart3c2.b1[i]=var(mb.32.b1)
}
b0var.table3[3,2]=sqrt(mean(vart3c2.b0))
b1var.table3[3,2]=sqrt(mean(vart3c2.b1))

# t_0=3 & C=30%
vart3c3.b0=c()
vart3c3.b1=c()
for (i in 1:500){
  a=data.gen(c.3,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.33.b0=c()
  mb.33.b1=c()
  result.mb.33.b0=c()
  result.mb.33.b1=c()
  w=a[,12]
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c3=nleqslv(betastart,rev.objectF,control=list(ftol=1e-5))
    mb.33.b0[j]=z.c3$x[1]
    mb.33.b1[j]=z.c3$x[2]
    result.mb.33.b0[j]=z.c3$fvec[1]
    result.mb.33.b1[j]=z.c3$fvec[2]
  }
  mb.33.b0=mb.33.b0[-1e-5<result.mb.33.b0&result.mb.33.b0<1e-5]
  mb.33.b1=mb.33.b1[-1e-5<result.mb.33.b1&result.mb.33.b1<1e-5]
  vart3c3.b0[i]=var(mb.33.b0)
  vart3c3.b1[i]=var(mb.33.b1)
}
b0var.table3[4,2]=sqrt(mean(vart3c3.b0))
b1var.table3[4,2]=sqrt(mean(vart3c3.b1))

#### method 3 : Induced Smoothing Multiplier Bootstrap ####
# t_0=0
# t_0=0 & C=0%

svart0c0.b0=c()
svart0c0.b1=c()
for (m in 1:500){
  a<-data.gen(c.0,0)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c0$x
  u=c()
  if (z.t0c0$fvec<1e-5&z.t0c0$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart0c0.b0[m]=sigma[1,1]
    svart0c0.b1[m]=sigma[2,2]
  } else {
    svart0c0.b0[m]=NA
    svart0c0.b1[m]=NA
  }
}
b0var.table0[1,3]=sqrt(mean(svart0c0.b0,na.rm=TRUE))
b1var.table0[1,3]=sqrt(mean(svart0c0.b1,na.rm=TRUE))

# t_0=0 & C=10%
svart0c1.b0=c()
svart0c1.b1=c()
for (i in 1:500){
  a<-data.gen(c.1,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c1$x
  u=c()
  for (j in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
    u=cbind(u,result)
  }
  v=cov(t(u))
  a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
  inva.beta=qr.solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart0c1.b0[i]=sigma[1,1]
  svart0c1.b1[i]=sigma[2,2]
}

b0var.table0[2,3]=sqrt(mean(svart0c1.b0,na.rm=TRUE))
b1var.table0[2,3]=sqrt(mean(svart0c1.b1,na.rm=TRUE))

# t_0=0 & C=20%
svart0c2.b0=c()
svart0c2.b1=c()
for (i in 1:500){
  a<-data.gen(c.2,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c2$x
  if (z.t0c2$fvec<1e-5&z.t0c2$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart0c2.b0[i]=sigma[1,1]
    svart0c2.b1[i]=sigma[2,2]
  } else {
    svart0c2.b0[i]=NA
    svart0c2.b1[i]=NA
  }
}
b0var.table0[3,3]=sqrt(mean(svart0c2.b0,na.rm=TRUE))
b1var.table0[3,3]=sqrt(mean(svart0c2.b1,na.rm=TRUE))

# t_0=0 & C=30%
svart0c3.b0=c()
svart0c3.b1=c()
for (i in 1:500){
  a<-data.gen(c.3,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c3$x
  if (z.t0c3$fvec<1e-5&z.t0c3$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart0c3.b0[i]=sigma[1,1]
    svart0c3.b1[i]=sigma[2,2]
  } else {
    svart0c3.b0[i]=NA
    svart0c3.b1[i]=NA
  }
}
b0var.table0[4,3]=sqrt(mean(svart0c3.b0,na.rm=TRUE))
b1var.table0[4,3]=sqrt(mean(svart0c3.b1,na.rm=TRUE))

# t_0=1
# t_0=1 & C=0%
svart1c0.b0=c()
svart1c0.b1=c()
for (m in 1:500){
  a<-data.gen(c.0,1)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c0$x
  u=c()
  if (z.t1c0$fvec<1e-5&z.t1c0$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart1c0.b0[m]=sigma[1,1]
    svart1c0.b1[m]=sigma[2,2]
  } else {
    svart1c0.b0[m]=NA
    svart1c0.b1[m]=NA
  }
}
b0var.table1[1,3]=sqrt(mean(svart1c0.b0,na.rm=TRUE))
b1var.table1[1,3]=sqrt(mean(svart1c0.b1,na.rm=TRUE))

# t_0=1 & C=10%
svart1c1.b0=c()
svart1c1.b1=c()
for (m in 1:500){
  a<-data.gen(c.1,1)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c1$x
  u=c()
  if (z.t1c1$fvec<1e-5&z.t1c1$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart1c1.b0[m]=sigma[1,1]
    svart1c1.b1[m]=sigma[2,2]
  } else {
    svart1c1.b0[m]=NA
    svart1c1.b1[m]=NA
  }
}
b0var.table1[2,3]=sqrt(mean(svart1c1.b0,na.rm=TRUE))
b1var.table1[2,3]=sqrt(mean(svart1c1.b1,na.rm=TRUE))

# t_0=1 & C=20%
svart1c2.b0=c()
svart1c2.b1=c()
for (m in 1:500){
  a<-data.gen(c.2,1)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c2$x
  u=c()
  if (z.t1c2$fvec<1e-5&z.t1c2$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart1c2.b0[m]=sigma[1,1]
    svart1c2.b1[m]=sigma[2,2]
  } else {
    svart1c2.b0[m]=NA
    svart1c2.b1[m]=NA
  }
}
b0var.table1[3,3]=sqrt(mean(svart1c2.b0,na.rm=TRUE))
b1var.table1[3,3]=sqrt(mean(svart1c2.b1,na.rm=TRUE))

# t_0=1 & C=30%
svart1c3.b0=c()
svart1c3.b1=c()
for (m in 1:500){
  a<-data.gen(c.3,1)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c3$x
  u=c()
  if (z.t1c3$fvec<1e-5&z.t1c3$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart1c3.b0[m]=sigma[1,1]
    svart1c3.b1[m]=sigma[2,2]
  } else {
    svart1c3.b0[m]=NA
    svart1c3.b1[m]=NA
  }
}
b0var.table1[4,3]=sqrt(mean(svart1c3.b0,na.rm=TRUE))
b1var.table1[4,3]=sqrt(mean(svart1c3.b1,na.rm=TRUE))

# t_0=2
# t_0=2 & C=0%
svart2c0.b0=c()
svart2c0.b1=c()
for (m in 1:500){
  a<-data.gen(c.0,2)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c0$x
  u=c()
  if (z.t2c0$fvec<1e-5&z.t2c0$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart2c0.b0[m]=sigma[1,1]
    svart2c0.b1[m]=sigma[2,2]
  } else {
    svart2c0.b0[m]=NA
    svart2c0.b1[m]=NA
  }
}
b0var.table2[1,3]=sqrt(mean(svart2c0.b0,na.rm=TRUE))
b1var.table2[1,3]=sqrt(mean(svart2c0.b1,na.rm=TRUE))

# t_0=2 & C=10%
svart2c1.b0=c()
svart2c1.b1=c()
for (m in 1:500){
  a<-data.gen(c.1,2)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c1$x
  u=c()
  if (z.t2c1$fvec<1e-5&z.t2c1$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart2c1.b0[m]=sigma[1,1]
    svart2c1.b1[m]=sigma[2,2]
  } else {
    svart2c1.b0[m]=NA
    svart2c1.b1[m]=NA
  }
}
b0var.table2[2,3]=sqrt(mean(svart2c1.b0,na.rm=TRUE))
b1var.table2[2,3]=sqrt(mean(svart2c1.b1,na.rm=TRUE))

# t_0=2 & C=20%
svart2c2.b0=c()
svart2c2.b1=c()
for (m in 1:500){
  a<-data.gen(c.2,2)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c2$x
  u=c()
  if (z.t2c2$fvec<1e-5&z.t2c2$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart2c2.b0[m]=sigma[1,1]
    svart2c2.b1[m]=sigma[2,2]
  } else {
    svart2c2.b0[m]=NA
    svart2c2.b1[m]=NA
  }
}
b0var.table2[3,3]=sqrt(mean(svart2c2.b0,na.rm=TRUE))
b1var.table2[3,3]=sqrt(mean(svart2c2.b1,na.rm=TRUE))

# t_0=2 & C=30%
svart2c3.b0=c()
svart2c3.b1=c()
for (m in 1:500){
  a<-data.gen(c.3,2)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c3$x
  u=c()
  if (z.t2c3$fvec<1e-5&z.t2c3$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart2c3.b0[m]=sigma[1,1]
    svart2c3.b1[m]=sigma[2,2]
  } else {
    svart2c3.b0[m]=NA
    svart2c3.b1[m]=NA
  }
}
b0var.table2[4,3]=sqrt(mean(svart2c3.b0,na.rm=TRUE))
b1var.table2[4,3]=sqrt(mean(svart2c3.b1,na.rm=TRUE))

# t_0=3
# t_0=3 & C=0%
svart3c0.b0=c()
svart3c0.b1=c()
for (m in 1:500){
  a<-data.gen(c.0,3)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c0$x
  u=c()
  if (z.t3c0$fvec<1e-5&z.t3c0$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart3c0.b0[m]=sigma[1,1]
    svart3c0.b1[m]=sigma[2,2]
  } else {
    svart3c0.b0[m]=NA
    svart3c0.b1[m]=NA
  }
}
b0var.table3[1,3]=sqrt(mean(svart3c0.b0,na.rm=TRUE))
b1var.table3[1,3]=sqrt(mean(svart3c0.b1,na.rm=TRUE))

# t_0=3 & C=10%
svart3c1.b0=c()
svart3c1.b1=c()
for (m in 1:500){
  a<-data.gen(c.1,3)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c1$x
  u=c()
  if (z.t3c1$fvec<1e-5&z.t3c1$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart3c1.b0[m]=sigma[1,1]
    svart3c1.b1[m]=sigma[2,2]
  } else {
    svart3c1.b0[m]=NA
    svart3c1.b1[m]=NA
  }
}
b0var.table3[2,3]=sqrt(mean(svart3c1.b0,na.rm=TRUE))
b1var.table3[2,3]=sqrt(mean(svart3c1.b1,na.rm=TRUE))

# t_0=3 & C=20%
svart3c2.b0=c()
svart3c2.b1=c()
for (m in 1:500){
  a<-data.gen(c.2,3)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c2$x
  u=c()
  if (z.t3c2$fvec<1e-5&z.t3c2$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart3c2.b0[m]=sigma[1,1]
    svart3c2.b1[m]=sigma[2,2]
  } else {
    svart3c2.b0[m]=NA
    svart3c2.b1[m]=NA
  }
}
b0var.table3[3,3]=sqrt(mean(svart3c2.b0,na.rm=TRUE))
b1var.table3[3,3]=sqrt(mean(svart3c2.b1,na.rm=TRUE))

# t_0=3 & C=30%
svart3c3.b0=c()
svart3c3.b1=c()
for (m in 1:500){
  a<-data.gen(c.3,3)
  ndata=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,ndata)),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/ndata,0,0,1/ndata),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c3$x
  u=c()
  if (z.t3c3$fvec<1e-5&z.t3c3$fvec>-1e-5){
    for (j in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart3c3.b0[m]=sigma[1,1]
    svart3c3.b1[m]=sigma[2,2]
  } else {
    svart3c3.b0[m]=NA
    svart3c3.b1[m]=NA
  }
}
b0var.table3[4,3]=sqrt(mean(svart3c3.b0,na.rm=TRUE))
b1var.table3[4,3]=sqrt(mean(svart3c3.b1,na.rm=TRUE))

#### Coverage check (95% Confidence Interval) ####
b0.coverage=matrix(NA,4,4)
colnames(b0.coverage)=c("0%","10%","20%","30%")
rownames(b0.coverage)=c("t0=0","t0=1","t0=2","t0=3")

b1.coverage=matrix(NA,4,4)
colnames(b1.coverage)=c("0%","10%","20%","30%")
rownames(b1.coverage)=c("t0=0","t0=1","t0=2","t0=3")

#### Indicator function ####
ind=function(a,b,c){
  if (a>=b&a<=c) {
    result=1
  } else {
    result=0
  }
  print(result)
}

## t_0=0 & c=0%
covert0c0=matrix(NA,500,8)
colnames(covert0c0)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.0,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c0$x
  u=c()
  if (z.t0c0$fvec[1]<1e-5&z.t0c0$fvec[1]>-1e-5&z.t0c0$fvec[2]<1e-5&z.t0c0$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t0c0=sqrt(sigma[1,1])
    b1se.t0c0=sqrt(sigma[2,2])
    covert0c0[i,1]=z.t0c0$x[1]-1.96*b0se.t0c0
    covert0c0[i,2]=z.t0c0$x[1]
    covert0c0[i,3]=z.t0c0$x[1]+1.96*b0se.t0c0
    covert0c0[i,4]=ind(1.61,z.t0c0$x[1]-1.96*b0se.t0c0,z.t0c0$x[1]+1.96*b0se.t0c0)
    covert0c0[i,5]=z.t0c0$x[2]-1.96*b1se.t0c0
    covert0c0[i,6]=z.t0c0$x[2]
    covert0c0[i,7]=z.t0c0$x[2]+1.96*b1se.t0c0
    covert0c0[i,8]=ind(0,z.t0c0$x[2]-1.96*b1se.t0c0,z.t0c0$x[2]+1.96*b1se.t0c0)
  } else {
    covert0c0[i,1:8]=NA
  }
}

b0.coverage[1,1]=mean(covert0c0[,4],na.rm=TRUE)
b1.coverage[1,1]=mean(covert0c0[,8],na.rm=TRUE)

## t_0=0 & c=10%
covert0c1=matrix(NA,500,8)
colnames(covert0c1)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.1,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c1$x
  u=c()
  if (z.t0c1$fvec[1]<1e-5&z.t0c1$fvec[1]>-1e-5&z.t0c1$fvec[2]<1e-5&z.t0c1$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t0c1=sqrt(sigma[1,1])
    b1se.t0c1=sqrt(sigma[2,2])
    covert0c1[i,1]=z.t0c1$x[1]-1.96*b0se.t0c1
    covert0c1[i,2]=z.t0c1$x[1]
    covert0c1[i,3]=z.t0c1$x[1]+1.96*b0se.t0c1
    covert0c1[i,4]=ind(1.61,z.t0c1$x[1]-1.96*b0se.t0c1,z.t0c1$x[1]+1.96*b0se.t0c1)
    covert0c1[i,5]=z.t0c1$x[2]-1.96*b1se.t0c1
    covert0c1[i,6]=z.t0c1$x[2]
    covert0c1[i,7]=z.t0c1$x[2]+1.96*b1se.t0c1
    covert0c1[i,8]=ind(0,z.t0c1$x[2]-1.96*b1se.t0c1,z.t0c1$x[2]+1.96*b1se.t0c1)
  } else {
    covert0c1[i,1:8]=NA
  }
}

b0.coverage[1,2]=mean(covert0c1[,4],na.rm=TRUE)
b1.coverage[1,2]=mean(covert0c1[,8],na.rm=TRUE)

## t_0=0 & c=20%
covert0c2=matrix(NA,500,8)
colnames(covert0c2)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.2,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c2$x
  u=c()
  if (z.t0c2$fvec[1]<1e-5&z.t0c2$fvec[1]>-1e-5&z.t0c2$fvec[2]<1e-5&z.t0c2$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t0c2=sqrt(sigma[1,1])
    b1se.t0c2=sqrt(sigma[2,2])
    covert0c2[i,1]=z.t0c2$x[1]-1.96*b0se.t0c2
    covert0c2[i,2]=z.t0c2$x[1]
    covert0c2[i,3]=z.t0c2$x[1]+1.96*b0se.t0c2
    covert0c2[i,4]=ind(1.61,z.t0c2$x[1]-1.96*b0se.t0c2,z.t0c2$x[1]+1.96*b0se.t0c2)
    covert0c2[i,5]=z.t0c2$x[2]-1.96*b1se.t0c2
    covert0c2[i,6]=z.t0c2$x[2]
    covert0c2[i,7]=z.t0c2$x[2]+1.96*b1se.t0c2
    covert0c2[i,8]=ind(0,z.t0c2$x[2]-1.96*b1se.t0c2,z.t0c2$x[2]+1.96*b1se.t0c2)
  } else {
    covert0c2[i,1:8]=NA
  }
}

b0.coverage[1,3]=mean(covert0c2[,4],na.rm=TRUE)
b1.coverage[1,3]=mean(covert0c2[,8],na.rm=TRUE)

## t_0=0 & c=30%
covert0c3=matrix(NA,500,8)
colnames(covert0c3)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.3,0)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.t0c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t0c3$x
  u=c()
  if (z.t0c3$fvec[1]<1e-5&z.t0c3$fvec[1]>-1e-5&z.t0c3$fvec[2]<1e-5&z.t0c3$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t0c3=sqrt(sigma[1,1])
    b1se.t0c3=sqrt(sigma[2,2])
    covert0c3[i,1]=z.t0c3$x[1]-1.96*b0se.t0c3
    covert0c3[i,2]=z.t0c3$x[1]
    covert0c3[i,3]=z.t0c3$x[1]+1.96*b0se.t0c3
    covert0c3[i,4]=ind(1.61,z.t0c3$x[1]-1.96*b0se.t0c3,z.t0c3$x[1]+1.96*b0se.t0c3)
    covert0c3[i,5]=z.t0c3$x[2]-1.96*b1se.t0c3
    covert0c3[i,6]=z.t0c3$x[2]
    covert0c3[i,7]=z.t0c3$x[2]+1.96*b1se.t0c3
    covert0c3[i,8]=ind(0,z.t0c3$x[2]-1.96*b1se.t0c3,z.t0c3$x[2]+1.96*b1se.t0c3)
  } else {
    covert0c3[i,1:8]=NA
  }
}

b0.coverage[1,4]=mean(covert0c3[,4],na.rm=TRUE)
b1.coverage[1,4]=mean(covert0c3[,8],na.rm=TRUE)

## t_0=1 & c=0%
covert1c0=matrix(NA,500,8)
colnames(covert1c0)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.0,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c0$x
  u=c()
  if (z.t1c0$fvec[1]<1e-5&z.t1c0$fvec[1]>-1e-5&z.t1c0$fvec[2]<1e-5&z.t1c0$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t1c0=sqrt(sigma[1,1])
    b1se.t1c0=sqrt(sigma[2,2])
    covert1c0[i,1]=z.t1c0$x[1]-1.96*b0se.t1c0
    covert1c0[i,2]=z.t1c0$x[1]
    covert1c0[i,3]=z.t1c0$x[1]+1.96*b0se.t1c0
    covert1c0[i,4]=ind(1.41,z.t1c0$x[1]-1.96*b0se.t1c0,z.t1c0$x[1]+1.96*b0se.t1c0)
    covert1c0[i,5]=z.t1c0$x[2]-1.96*b1se.t1c0
    covert1c0[i,6]=z.t1c0$x[2]
    covert1c0[i,7]=z.t1c0$x[2]+1.96*b1se.t1c0
    covert1c0[i,8]=ind(0,z.t1c0$x[2]-1.96*b1se.t1c0,z.t1c0$x[2]+1.96*b1se.t1c0)
  } else {
    covert1c0[i,1:8]=NA
  }
}

b0.coverage[2,1]=mean(covert1c0[,4],na.rm=TRUE)
b1.coverage[2,1]=mean(covert1c0[,8],na.rm=TRUE)

## t_0=1 & c=10%
covert1c1=matrix(NA,500,8)
colnames(covert1c1)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.1,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c1$x
  u=c()
  if (z.t1c1$fvec[1]<1e-5&z.t1c1$fvec[1]>-1e-5&z.t1c1$fvec[2]<1e-5&z.t1c1$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t1c1=sqrt(sigma[1,1])
    b1se.t1c1=sqrt(sigma[2,2])
    covert1c1[i,1]=z.t1c1$x[1]-1.96*b0se.t1c1
    covert1c1[i,2]=z.t1c1$x[1]
    covert1c1[i,3]=z.t1c1$x[1]+1.96*b0se.t1c1
    covert1c1[i,4]=ind(1.41,z.t1c1$x[1]-1.96*b0se.t1c1,z.t1c1$x[1]+1.96*b0se.t1c1)
    covert1c1[i,5]=z.t1c1$x[2]-1.96*b1se.t1c1
    covert1c1[i,6]=z.t1c1$x[2]
    covert1c1[i,7]=z.t1c1$x[2]+1.96*b1se.t1c1
    covert1c1[i,8]=ind(0,z.t1c1$x[2]-1.96*b1se.t1c1,z.t1c1$x[2]+1.96*b1se.t1c1)
  } else {
    covert1c1[i,1:8]=NA
  }
}

b0.coverage[2,2]=mean(covert1c1[,4],na.rm=TRUE)
b1.coverage[2,2]=mean(covert1c1[,8],na.rm=TRUE)

## t_0=1 & c=20%
covert1c2=matrix(NA,500,8)
colnames(covert1c2)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.2,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c2$x
  u=c()
  if (z.t1c2$fvec[1]<1e-5&z.t1c2$fvec[1]>-1e-5&z.t1c2$fvec[2]<1e-5&z.t1c2$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t1c2=sqrt(sigma[1,1])
    b1se.t1c2=sqrt(sigma[2,2])
    covert1c2[i,1]=z.t1c2$x[1]-1.96*b0se.t1c2
    covert1c2[i,2]=z.t1c2$x[1]
    covert1c2[i,3]=z.t1c2$x[1]+1.96*b0se.t1c2
    covert1c2[i,4]=ind(1.41,z.t1c2$x[1]-1.96*b0se.t1c2,z.t1c2$x[1]+1.96*b0se.t1c2)
    covert1c2[i,5]=z.t1c2$x[2]-1.96*b1se.t1c2
    covert1c2[i,6]=z.t1c2$x[2]
    covert1c2[i,7]=z.t1c2$x[2]+1.96*b1se.t1c2
    covert1c2[i,8]=ind(0,z.t1c2$x[2]-1.96*b1se.t1c2,z.t1c2$x[2]+1.96*b1se.t1c2)
  } else {
    covert1c2[i,1:8]=NA
  }
}

b0.coverage[2,3]=mean(covert1c2[,4],na.rm=TRUE)
b1.coverage[2,3]=mean(covert1c2[,8],na.rm=TRUE)

## t_0=1 & c=30%
covert1c3=matrix(NA,500,8)
colnames(covert1c3)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.3,1)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.41,0) # from table0
  z.t1c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t1c3$x
  u=c()
  if (z.t1c3$fvec[1]<1e-5&z.t1c3$fvec[1]>-1e-5&z.t1c3$fvec[2]<1e-5&z.t1c3$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t1c3=sqrt(sigma[1,1])
    b1se.t1c3=sqrt(sigma[2,2])
    covert1c3[i,1]=z.t1c3$x[1]-1.96*b0se.t1c3
    covert1c3[i,2]=z.t1c3$x[1]
    covert1c3[i,3]=z.t1c3$x[1]+1.96*b0se.t1c3
    covert1c3[i,4]=ind(1.41,z.t1c3$x[1]-1.96*b0se.t1c3,z.t1c3$x[1]+1.96*b0se.t1c3)
    covert1c3[i,5]=z.t1c3$x[2]-1.96*b1se.t1c3
    covert1c3[i,6]=z.t1c3$x[2]
    covert1c3[i,7]=z.t1c3$x[2]+1.96*b1se.t1c3
    covert1c3[i,8]=ind(0,z.t1c3$x[2]-1.96*b1se.t1c3,z.t1c3$x[2]+1.96*b1se.t1c3)
  } else {
    covert1c3[i,1:8]=NA
  }
}

b0.coverage[2,4]=mean(covert1c3[,4],na.rm=TRUE)
b1.coverage[2,4]=mean(covert1c3[,8],na.rm=TRUE)

## t_0=2 & c=0%
covert2c0=matrix(NA,500,8)
colnames(covert2c0)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.0,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c0$x
  u=c()
  if (z.t2c0$fvec[1]<1e-5&z.t2c0$fvec[1]>-1e-5&z.t2c0$fvec[2]<1e-5&z.t2c0$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t2c0=sqrt(sigma[1,1])
    b1se.t2c0=sqrt(sigma[2,2])
    covert2c0[i,1]=z.t2c0$x[1]-1.96*b0se.t2c0
    covert2c0[i,2]=z.t2c0$x[1]
    covert2c0[i,3]=z.t2c0$x[1]+1.96*b0se.t2c0
    covert2c0[i,4]=ind(1.22,z.t2c0$x[1]-1.96*b0se.t2c0,z.t2c0$x[1]+1.96*b0se.t2c0)
    covert2c0[i,5]=z.t2c0$x[2]-1.96*b1se.t2c0
    covert2c0[i,6]=z.t2c0$x[2]
    covert2c0[i,7]=z.t2c0$x[2]+1.96*b1se.t2c0
    covert2c0[i,8]=ind(0,z.t2c0$x[2]-1.96*b1se.t2c0,z.t2c0$x[2]+1.96*b1se.t2c0)
  } else {
    covert2c0[i,1:8]=NA
  }
}

b0.coverage[3,1]=mean(covert2c0[,4])
b1.coverage[3,1]=mean(covert2c0[,8])

## t_0=2 & c=10%
covert2c1=matrix(NA,500,8)
colnames(covert2c1)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.1,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c1$x
  u=c()
  if (z.t2c1$fvec[1]<1e-5&z.t2c1$fvec[1]>-1e-5&z.t2c1$fvec[2]<1e-5&z.t2c1$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t2c1=sqrt(sigma[1,1])
    b1se.t2c1=sqrt(sigma[2,2])
    covert2c1[i,1]=z.t2c1$x[1]-1.96*b0se.t2c1
    covert2c1[i,2]=z.t2c1$x[1]
    covert2c1[i,3]=z.t2c1$x[1]+1.96*b0se.t2c1
    covert2c1[i,4]=ind(1.22,z.t2c1$x[1]-1.96*b0se.t2c1,z.t2c1$x[1]+1.96*b0se.t2c1)
    covert2c1[i,5]=z.t2c1$x[2]-1.96*b1se.t2c1
    covert2c1[i,6]=z.t2c1$x[2]
    covert2c1[i,7]=z.t2c1$x[2]+1.96*b1se.t2c1
    covert2c1[i,8]=ind(0,z.t2c1$x[2]-1.96*b1se.t2c1,z.t2c1$x[2]+1.96*b1se.t2c1)
  } else {
    covert2c1[i,1:8]=NA
  }
}

b0.coverage[3,2]=mean(covert2c1[,4])
b1.coverage[3,2]=mean(covert2c1[,8])

## t_0=2 & c=20%
covert2c2=matrix(NA,500,8)
colnames(covert2c2)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.2,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c2$x
  u=c()
  if (z.t2c2$fvec[1]<1e-5&z.t2c2$fvec[1]>-1e-5&z.t2c2$fvec[2]<1e-5&z.t2c2$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t2c2=sqrt(sigma[1,1])
    b1se.t2c2=sqrt(sigma[2,2])
    covert2c2[i,1]=z.t2c2$x[1]-1.96*b0se.t2c2
    covert2c2[i,2]=z.t2c2$x[1]
    covert2c2[i,3]=z.t2c2$x[1]+1.96*b0se.t2c2
    covert2c2[i,4]=ind(1.22,z.t2c2$x[1]-1.96*b0se.t2c2,z.t2c2$x[1]+1.96*b0se.t2c2)
    covert2c2[i,5]=z.t2c2$x[2]-1.96*b1se.t2c2
    covert2c2[i,6]=z.t2c2$x[2]
    covert2c2[i,7]=z.t2c2$x[2]+1.96*b1se.t2c2
    covert2c2[i,8]=ind(0,z.t2c2$x[2]-1.96*b1se.t2c2,z.t2c2$x[2]+1.96*b1se.t2c2)
  } else {
    covert2c2[i,1:8]=NA
  }
}

b0.coverage[3,3]=mean(covert2c2[,4])
b1.coverage[3,3]=mean(covert2c2[,8])

## t_0=2 & c=30%
covert2c3=matrix(NA,500,8)
colnames(covert2c3)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.3,2)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.22,0) # from table0
  z.t2c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t2c3$x
  u=c()
  if (z.t2c3$fvec[1]<1e-5&z.t2c3$fvec[1]>-1e-5&z.t2c3$fvec[2]<1e-5&z.t2c3$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t2c3=sqrt(sigma[1,1])
    b1se.t2c3=sqrt(sigma[2,2])
    covert2c3[i,1]=z.t2c3$x[1]-1.96*b0se.t2c3
    covert2c3[i,2]=z.t2c3$x[1]
    covert2c3[i,3]=z.t2c3$x[1]+1.96*b0se.t2c3
    covert2c3[i,4]=ind(1.22,z.t2c3$x[1]-1.96*b0se.t2c3,z.t2c3$x[1]+1.96*b0se.t2c3)
    covert2c3[i,5]=z.t2c3$x[2]-1.96*b1se.t2c3
    covert2c3[i,6]=z.t2c3$x[2]
    covert2c3[i,7]=z.t2c3$x[2]+1.96*b1se.t2c3
    covert2c3[i,8]=ind(0,z.t2c3$x[2]-1.96*b1se.t2c3,z.t2c3$x[2]+1.96*b1se.t2c3)
  } else {
    covert2c3[i,1:8]=NA
  }
}

b0.coverage[3,4]=mean(covert2c3[,4])
b1.coverage[3,4]=mean(covert2c3[,8])

## t_0=3 & c=0%
covert3c0=matrix(NA,500,8)
colnames(covert3c0)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.0,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c0<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c0$x
  u=c()
  if (z.t3c0$fvec[1]<1e-5&z.t3c0$fvec[1]>-1e-5&z.t3c0$fvec[2]<1e-5&z.t3c0$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t3c0=sqrt(sigma[1,1])
    b1se.t3c0=sqrt(sigma[2,2])
    covert3c0[i,1]=z.t3c0$x[1]-1.96*b0se.t3c0
    covert3c0[i,2]=z.t3c0$x[1]
    covert3c0[i,3]=z.t3c0$x[1]+1.96*b0se.t3c0
    covert3c0[i,4]=ind(1.04,z.t3c0$x[1]-1.96*b0se.t3c0,z.t3c0$x[1]+1.96*b0se.t3c0)
    covert3c0[i,5]=z.t3c0$x[2]-1.96*b1se.t3c0
    covert3c0[i,6]=z.t3c0$x[2]
    covert3c0[i,7]=z.t3c0$x[2]+1.96*b1se.t3c0
    covert3c0[i,8]=ind(0,z.t3c0$x[2]-1.96*b1se.t3c0,z.t3c0$x[2]+1.96*b1se.t3c0)
  } else {
    covert3c0[i,1:8]=NA
  }
}

b0.coverage[4,1]=mean(covert3c0[,4])
b1.coverage[4,1]=mean(covert3c0[,8])

## t_0=3 & c=10%
covert3c1=matrix(NA,500,8)
colnames(covert3c1)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.1,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c1<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c1$x
  u=c()
  if (z.t3c1$fvec[1]<1e-5&z.t3c1$fvec[1]>-1e-5&z.t3c1$fvec[2]<1e-5&z.t3c1$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t3c1=sqrt(sigma[1,1])
    b1se.t3c1=sqrt(sigma[2,2])
    covert3c1[i,1]=z.t3c1$x[1]-1.96*b0se.t3c1
    covert3c1[i,2]=z.t3c1$x[1]
    covert3c1[i,3]=z.t3c1$x[1]+1.96*b0se.t3c1
    covert3c1[i,4]=ind(1.04,z.t3c1$x[1]-1.96*b0se.t3c1,z.t3c1$x[1]+1.96*b0se.t3c1)
    covert3c1[i,5]=z.t3c1$x[2]-1.96*b1se.t3c1
    covert3c1[i,6]=z.t3c1$x[2]
    covert3c1[i,7]=z.t3c1$x[2]+1.96*b1se.t3c1
    covert3c1[i,8]=ind(0,z.t3c1$x[2]-1.96*b1se.t3c1,z.t3c1$x[2]+1.96*b1se.t3c1)
  } else {
    covert3c1[i,1:8]=NA
  }
}

b0.coverage[4,2]=mean(covert3c1[,4])
b1.coverage[4,2]=mean(covert3c1[,8])

## t_0=3 & c=20%
covert3c2=matrix(NA,500,8)
colnames(covert3c2)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.2,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c2<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c2$x
  u=c()
  if (z.t3c2$fvec[1]<1e-5&z.t3c2$fvec[1]>-1e-5&z.t3c2$fvec[2]<1e-5&z.t3c2$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t3c2=sqrt(sigma[1,1])
    b1se.t3c2=sqrt(sigma[2,2])
    covert3c2[i,1]=z.t3c2$x[1]-1.96*b0se.t3c2
    covert3c2[i,2]=z.t3c2$x[1]
    covert3c2[i,3]=z.t3c2$x[1]+1.96*b0se.t3c2
    covert3c2[i,4]=ind(1.04,z.t3c2$x[1]-1.96*b0se.t3c2,z.t3c2$x[1]+1.96*b0se.t3c2)
    covert3c2[i,5]=z.t3c2$x[2]-1.96*b1se.t3c2
    covert3c2[i,6]=z.t3c2$x[2]
    covert3c2[i,7]=z.t3c2$x[2]+1.96*b1se.t3c2
    covert3c2[i,8]=ind(0,z.t3c2$x[2]-1.96*b1se.t3c2,z.t3c2$x[2]+1.96*b1se.t3c2)
  } else {
    covert3c2[i,1:8]=NA
  }
}

b0.coverage[4,3]=mean(covert3c2[,4])
b1.coverage[4,3]=mean(covert3c2[,8])

## t_0=3 & c=30%
covert3c3=matrix(NA,500,8)
colnames(covert3c3)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
for (i in 1:500) {
  a<-data.gen(c.3,3)
  n=nrow(a)
  w=a[,12]
  # Matrix X_i (x[,1]=intercept, x[,2]=covariate)  
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G<-matrix(c(1/n,0,0,1/n),2,2)
  # Find estimated value of beta
  betastart<-c(1.04,0) # from table0
  z.t3c3<-nleqslv(betastart,objectF,control=list(ftol=1e-5))
  solbeta=z.t3c3$x
  u=c()
  if (z.t3c3$fvec[1]<1e-5&z.t3c3$fvec[1]>-1e-5&z.t3c3$fvec[2]<1e-5&z.t3c3$fvec[2]>-1e-5){
    for (j in 1:500){
      eta=rexp(n,1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.75)
      u=cbind(u,result)
    }
    v=cov(t(u))
    a.beta=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X)))))*X)%*%(X/sqrt(diag(X%*%G%*%t(X))))
    inva.beta=qr.solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    b0se.t3c3=sqrt(sigma[1,1])
    b1se.t3c3=sqrt(sigma[2,2])
    covert3c3[i,1]=z.t3c3$x[1]-1.96*b0se.t3c3
    covert3c3[i,2]=z.t3c3$x[1]
    covert3c3[i,3]=z.t3c3$x[1]+1.96*b0se.t3c3
    covert3c3[i,4]=ind(1.04,z.t3c3$x[1]-1.96*b0se.t3c3,z.t3c3$x[1]+1.96*b0se.t3c3)
    covert3c3[i,5]=z.t3c3$x[2]-1.96*b1se.t3c3
    covert3c3[i,6]=z.t3c3$x[2]
    covert3c3[i,7]=z.t3c3$x[2]+1.96*b1se.t3c3
    covert3c3[i,8]=ind(0,z.t3c3$x[2]-1.96*b1se.t3c3,z.t3c3$x[2]+1.96*b1se.t3c3)
  } else {
    covert3c3[i,1:8]=NA
  }
}

b0.coverage[4,4]=mean(covert3c3[,4])
b1.coverage[4,4]=mean(covert3c3[,8])