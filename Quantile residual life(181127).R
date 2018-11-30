library(quantreg)
library(survival)
library(BB)
library(nleqslv)
library(tictoc)
  
#### Given information ####
c.0=5000000
c.1=53.03
c.2=26.58
c.3=17.73
exp.beta.initial=5
k=2
r.initial=(log(2))^(1/k)/exp.beta.initial

#### Data Generation function ####
data.gen<-function(censor, t_0){
  u = runif(n=400,min = 0,max = 1)
  sim=matrix(NA,400,7)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5))
  sim[,1] = {{-log(1-u)}^(1/k)}/r.initial
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

#### Beta estimation ####
#### method 1 : Crq function (±âÁ¸ package) ####
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
  dummy0[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]
  dummy1[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit0,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]  
  dummy1[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit1,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]
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
  dummy0[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit2,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
  dummy0[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,1]
  dummy0.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[1,4]   
  dummy1[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,1]
  dummy1.sd[i]<-summary(fit3,tau=c(0.4,0.5))[2][[1]]$coefficient[2,4]  
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
var.table0=matrix(NA,4,4)
colnames(var.table0)=c("true","MB","ISMB(simple)","ISMB(Iter)")
rownames(var.table0)=c(0,10,20,30)
var.table1=matrix(NA,4,4)
colnames(var.table1)=c("true","MB","ISMB(simple)","ISMB(Iter)")
rownames(var.table1)=c(0,10,20,30)
var.table2=matrix(NA,4,4)
colnames(var.table2)=c("true","MB","ISMB(simple)","ISMB(Iter)")
rownames(var.table2)=c(0,10,20,30)
var.table3=matrix(NA,4,4)
colnames(var.table3)=c("true","MB","ISMB(simple)","ISMB(Iter)")
rownames(var.table3)=c(0,10,20,30)
#### Method 1 : True variance by beta estimation ####
var.table0[1,1]<-sd(b0.t0c0)
var.table0[2,1]<-sd(b0.t0c1)
var.table0[3,1]<-sd(b0.t0c2)
var.table0[4,1]<-sd(b0.t0c3)
var.table1[1,1]<-sd(b0.t1c0)
var.table1[2,1]<-sd(b0.t1c1)
var.table1[3,1]<-sd(b0.t1c2)
var.table1[4,1]<-sd(b0.t1c3)
var.table2[1,1]<-sd(b0.t2c0)
var.table2[2,1]<-sd(b0.t2c1)
var.table2[3,1]<-sd(b0.t2c2)
var.table2[4,1]<-sd(b0.t2c3)
var.table3[1,1]<-sd(b0.t3c0)
var.table3[2,1]<-sd(b0.t3c1)
var.table3[3,1]<-sd(b0.t3c2)
var.table3[4,1]<-sd(b0.t3c3)

#### method 2 : Multiplier Bootstrap ####
# t_0=0
# t_0=0 & C=0
vart0c0=c()
for (k in 1:500){
  a=data.gen(c.0,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.00.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c0=nleqslv(betastart,rev.objectF)
    mb.00.b0[m]=z.c0$x[1]
  }
  mb.00.b0=mb.00.b0[0<mb.00.b0&mb.00.b0<2]
  vart0c0[k]=sd(mb.00.b0)
}
vart0c0=na.omit(vart0c0)
var.table0[1,2]=mean(vart0c0)

# t_0=0 & C=10%
vart0c1=c()
for (k in 1:500){
  a=data.gen(c.1,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.01.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c1=nleqslv(betastart,rev.objectF)
    mb.01.b0[m]=z.c1$x[1]
  }
  mb.01.b0=mb.01.b0[0<mb.01.b0&mb.01.b0<2]
  vart0c1[k]=sd(mb.01.b0)
}
vart0c1=na.omit(vart0c1)
var.table0[2,2]=mean(vart0c1)

# t_0=0 & C=20%
vart0c2=c()
for (k in 1:500){
  a=data.gen(c.2,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.02.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c2=nleqslv(betastart,rev.objectF)
    mb.02.b0[m]=z.c2$x[1]
  }
  mb.02.b0=mb.02.b0[0<mb.02.b0&mb.02.b0<2]
  vart0c2[k]=sd(mb.02.b0)
}
vart0c2=na.omit(vart0c2)
var.table0[3,2]=mean(vart0c2)

# t_0=0 & C=30%
vart0c3=c()
for (k in 1:500){
  a=data.gen(c.3,0)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.03.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.61,0)
    z.c3=nleqslv(betastart,rev.objectF)
    mb.03.b0[m]=z.c3$x[1]
  }
  mb.03.b0=mb.03.b0[0<mb.03.b0&mb.03.b0<2]
  vart0c3[k]=sd(mb.03.b0)
}
vart0c3=na.omit(vart0c3)
var.table0[4,2]=mean(vart0c3)

# t_0=1
# t_0=1 & C=0
vart1c0=c()
for (k in 1:500){
  a=data.gen(c.0,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.10.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.41,0)
    z.c0=nleqslv(betastart,rev.objectF)
    mb.10.b0[m]=z.c0$x[1]
  }
  mb.10.b0=mb.10.b0[0<mb.10.b0&mb.10.b0<2]
  vart1c0[k]=sd(mb.10.b0)
}
vart1c0=na.omit(vart1c0)
var.table1[1,2]=mean(vart1c0)

# t_0=1 & C=10%
vart1c1=c()
for (k in 1:500){
  a=data.gen(c.1,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.11.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.41,0)
    z.c1=nleqslv(betastart,rev.objectF)
    mb.11.b0[m]=z.c1$x[1]
  }
  mb.11.b0=mb.11.b0[0<mb.11.b0&mb.11.b0<2]
  vart1c1[k]=sd(mb.11.b0)
}
vart1c1=na.omit(vart1c1)
var.table1[2,2]=mean(vart1c1)

# t_0=1 & C=20%
vart1c2=c()
for (k in 1:500){
  a=data.gen(c.2,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.12.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.41,0)
    z.c2=nleqslv(betastart,rev.objectF)
    mb.12.b0[m]=z.c2$x[1]
  }
  mb.12.b0=mb.12.b0[0<mb.12.b0&mb.12.b0<2]
  vart1c2[k]=sd(mb.12.b0)
}
vart1c2=na.omit(vart1c2)
var.table1[3,2]=mean(vart1c2)

# t_0=1 & C=30%
vart1c3=c()
for (k in 1:500){
  a=data.gen(c.3,1)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.13.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.41,0)
    z.c3=nleqslv(betastart,rev.objectF)
    mb.13.b0[m]=z.c3$x[1]
  }
  mb.13.b0=mb.13.b0[0<mb.13.b0&mb.13.b0<2]
  vart1c3[k]=sd(mb.13.b0)
}
vart1c3=na.omit(vart1c3)
var.table1[4,2]=mean(vart1c3)

# t_0=2
# t_0=2 & C=0
vart2c0=c()
for (k in 1:500){
  a=data.gen(c.0,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.20.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c0=nleqslv(betastart,rev.objectF)
    mb.20.b0[m]=z.c0$x[1]
  }
  mb.20.b0=mb.20.b0[0<mb.20.b0&mb.20.b0<2]
  vart2c0[k]=sd(mb.20.b0)
}
vart2c0=na.omit(vart2c0)
var.table2[1,2]=mean(vart2c0)

# t_0=2 & C=10%
vart2c1=c()
for (k in 1:500){
  a=data.gen(c.1,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.21.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c1=nleqslv(betastart,rev.objectF)
    mb.21.b0[m]=z.c1$x[1]
  }
  mb.21.b0=mb.21.b0[0<mb.21.b0&mb.21.b0<2]
  vart2c1[k]=sd(mb.21.b0)
}
vart2c1=na.omit(vart2c1)
var.table2[2,2]=mean(vart2c1)

# t_0=2 & C=20%
vart2c2=c()
for (k in 1:500){
  a=data.gen(c.2,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.22.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c2=nleqslv(betastart,rev.objectF)
    mb.22.b0[m]=z.c2$x[1]
  }
  mb.22.b0=mb.22.b0[0<mb.22.b0&mb.22.b0<2]
  vart2c2[k]=sd(mb.22.b0)
}
vart2c2=na.omit(vart2c2)
var.table2[3,2]=mean(vart2c2)

# t_0=2 & C=30%
vart2c3=c()
for (k in 1:500){
  a=data.gen(c.3,2)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.23.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.22,0)
    z.c3=nleqslv(betastart,rev.objectF)
    mb.23.b0[m]=z.c3$x[1]
  }
  mb.23.b0=mb.23.b0[0<mb.23.b0&mb.23.b0<2]
  vart2c3[k]=sd(mb.23.b0)
}
vart2c3=na.omit(vart2c3)
var.table2[4,2]=mean(vart2c3)

# t_0=3
# t_0=3 & C=0
vart3c0=c()
for (k in 1:500){
  a=data.gen(c.0,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.30.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c0=nleqslv(betastart,rev.objectF)
    mb.30.b0[m]=z.c0$x[1]
  }
  mb.30.b0=mb.30.b0[0<mb.30.b0&mb.30.b0<2]
  vart3c0[k]=sd(mb.30.b0)
}
vart3c0=na.omit(vart3c0)
var.table3[1,2]=mean(vart3c0)

# t_0=3 & C=10%
vart3c1=c()
for (k in 1:500){
  a=data.gen(c.1,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.31.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c1=nleqslv(betastart,rev.objectF)
    mb.31.b0[m]=z.c1$x[1]
  }
  mb.31.b0=mb.31.b0[0<mb.31.b0&mb.31.b0<2]
  vart3c1[k]=sd(mb.31.b0)
}
vart3c1=na.omit(vart3c1)
var.table3[2,2]=mean(vart3c1)

# t_0=3 & C=20%
vart3c2=c()
for (k in 1:500){
  a=data.gen(c.2,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.32.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c2=nleqslv(betastart,rev.objectF)
    mb.32.b0[m]=z.c2$x[1]
  }
  mb.32.b0=mb.32.b0[0<mb.32.b0&mb.32.b0<2]
  vart3c2[k]=sd(mb.32.b0)
}
vart3c2=na.omit(vart3c2)
var.table3[3,2]=mean(vart3c2)

# t_0=3 & C=30%
vart3c3=c()
for (k in 1:500){
  a=data.gen(c.3,3)
  n=nrow(a)
  X=cbind(c(rep(1,nrow(a))),a[,6])
  # Defining Gamma (1/n)I_p
  G=matrix(c(1/n,0,0,1/n),2,2)
  mb.33.b0=c()
  w=a[,12]
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    betastart=c(1.04,0)
    z.c3=nleqslv(betastart,rev.objectF)
    mb.33.b0[m]=z.c3$x[1]
  }
  mb.33.b0=mb.33.b0[0<mb.33.b0&mb.33.b0<2]
  vart3c3[k]=sd(mb.33.b0)
}
vart3c3=na.omit(vart3c3)
var.table3[4,2]=mean(vart3c)

#### method 3 : Induced Smoothing Multiplier Bootstrap ####
# t_0=0
# t_0=0 & C=0
svart0c0=c()
for (k in 1:500){
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
  solbeta=z.t0c0$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart0c0[k]=sqrt(sigma[1,1])
}
var.table0[1,3]=mean(svart0c0,na.rm=TRUE)

# t_0=0 & C=10%
svart0c1=c()
  for (k in 1:500){
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
    solbeta=z.t0c1$x
    u=c()
    for (m in 1:500){
      eta=rexp(nrow(a),1)
      result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
      u=cbind(u,result)
    }
    v=apply(u,1,var)
    v=diag(v,2,2)
    a.beta=matrix(c(0,0,0,0),2,2)
    a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
    a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
    a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
    a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
    inva.beta=solve(a.beta)
    sigma=inva.beta%*%v%*%t(inva.beta)
    svart0c1[k]=sqrt(sigma[1,1])
  }
var.table0[2,3]=mean(svart0c1,na.rm=TRUE)

# t_0=0 & C=20%
svart0c2=c()
for (k in 1:500){
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
  solbeta=z.t0c2$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart0c2[k]=sqrt(sigma[1,1])
}
var.table0[3,3]=mean(svart0c2,na.rm=TRUE)

# t_0=0 & C=30%
svart0c3=c()
for (k in 1:500){
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
  solbeta=z.t0c3$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart0c3[k]=sqrt(sigma[1,1])
}
var.table0[4,3]=mean(svart0c3,na.rm=TRUE)

# t_0=1
# t_0=1 & C=0%
svart1c0=c()
for (k in 1:500){
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
  solbeta=z.t1c0$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart1c0[k]=sqrt(sigma[1,1])
}
var.table1[1,3]=mean(svart1c0,na.rm=TRUE)

# t_0=1 & C=10%
svart1c1=c()
for (k in 1:500){
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
  solbeta=z.t1c1$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart1c1[k]=sqrt(sigma[1,1])
}
var.table1[2,3]=mean(svart1c1,na.rm=TRUE)

# t_0=1 & C=20%
svart1c2=c()
for (k in 1:500){
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
  solbeta=z.t1c2$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart1c2[k]=sqrt(sigma[1,1])
}
var.table1[3,3]=mean(svart1c2,na.rm=TRUE)

# t_0=1 & C=30%
svart1c3=c()
for (k in 1:500){
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
  solbeta=z.t1c3$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart1c3[k]=sqrt(sigma[1,1])
}
var.table1[4,3]=mean(svart1c3,na.rm=TRUE)

# t_0=2
# t_0=2 & C=0%
svart2c0=c()
for (k in 1:500){
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
  solbeta=z.t2c0$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart2c0[k]=sqrt(sigma[1,1])
}
var.table2[1,3]=mean(svart2c0,na.rm=TRUE)

# t_0=2 & C=10%
svart2c1=c()
for (k in 1:500){
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
  solbeta=z.t2c1$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart2c1[k]=sqrt(sigma[1,1])
}
var.table2[2,3]=mean(svart2c1,na.rm=TRUE)

# t_0=2 & C=20%
svart2c2=c()
for (k in 1:500){
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
  solbeta=z.t2c2$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart2c2[k]=sqrt(sigma[1,1])
}
var.table2[3,3]=mean(svart2c2,na.rm=TRUE)

# t_0=2 & C=30%
svart2c3=c()
for (k in 1:500){
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
  solbeta=z.t2c3$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart2c3[k]=sqrt(sigma[1,1])
}
var.table2[4,3]=mean(svart2c3,na.rm=TRUE)

# t_0=3
# t_0=3 & C=0%
svart3c0=c()
for (k in 1:500){
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
  solbeta=z.t3c0$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart3c0[k]=sqrt(sigma[1,1])
}
var.table3[1,3]=mean(svart3c0,na.rm=TRUE)

# t_0=3 & C=10%
svart3c1=c()
for (k in 1:500){
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
  solbeta=z.t3c1$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart3c1[k]=sqrt(sigma[1,1])
}
var.table3[2,3]=mean(svart3c1,na.rm=TRUE)

# t_0=03& C=20%
svart3c2=c()
for (k in 1:500){
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
  solbeta=z.t3c2$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart3c2[k]=sqrt(sigma[1,1])
}
var.table3[3,3]=mean(svart3c2,na.rm=TRUE)

# t_0=3 & C=30%
svart3c3=c()
for (k in 1:500){
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
  solbeta=z.t3c3$x
  u=c()
  for (m in 1:500){
    eta=rexp(nrow(a),1)
    result=t(eta*w*X)%*%(pnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))-0.5)
    u=cbind(u,result)
  }
  v=apply(u,1,var)
  v=diag(v,2,2)
  a.beta=matrix(c(0,0,0,0),2,2)
  a.beta[1,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,1])
  a.beta[1,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,1]*X[,2])
  a.beta[2,1]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,1])
  a.beta[2,2]=t(w*as.vector(dnorm((a[,5]-X%*%solbeta)/sqrt(diag(X%*%G%*%t(X))))))%*%(X[,2]*X[,2])
  inva.beta=solve(a.beta)
  sigma=inva.beta%*%v%*%t(inva.beta)
  svart3c3[k]=sqrt(sigma[1,1])
}
var.table3[4,3]=mean(svart3c3,na.rm=TRUE)
