library(quantreg)
library(survival)

#### Data generation of JUNG and KIM #####
#Find C
#Generate T_i 1,000,000
d.exp.beta.initial=5
d.k=2
d.r.initial=(log(2)^(1/d.k))/d.exp.beta.initial
d.u=runif(n=1000000,min = 0,max = 1)
d.T={-log(1-d.u)}^(1/d.k)/d.r.initial
# Find c which is derermined to achieve 0,10,20,30% censoring rate in case 1
i=70
while(1){
  d.C<-runif(1000000,0,i)
  i=i+0.01
  if(sum(d.C<d.T)<=0) break # 0% = 0, 10% = 100000, 20% = 200000, 30% = 300000)
}
print(i)
table(d.C<d.T)


# result of finding C
c.0=5000000
c.1=53.03
c.2=26.58
c.3=17.73
exp.beta.initial=5
k=2
r.initial=(log(2))^(1/k)/exp.beta.initial

#### Data Generation function ####
data.gen<-function(censor, t_0){
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5))
  u=runif(n=200,min = 0,max = 1)
  T={{-log(1-u)}^(1/k)}/r.initial
  # Generate C_i
  c<-runif(200,0,censor)
  #Covariates (Control=0, Treatment=1)
  X=rbinom(200,size=1,p=0.5)
  #Combine data
  sim<-matrix(NA,200,7)
  sim[,1]=T
  sim[,2]=c
  for (i in 1:200)
    {
      sim[i,3]<-min(sim[i,1],sim[i,2])
      sim[i,4]<-sim[i,3]-(t_0)
      sim[i,5]<-log(sim[i,4])
  sim[,6]<-X
  if(sim[i,1]<=sim[i,2])
    {
    sim[i,7]=1
    }
  else{
    sim[i,7]=0
  }
  }
  colnames(sim)<-c("T","C","Z","Z.diff","log(Z.diff)","X","censored")
  sim<-as.data.frame(sim)
  return(sim)
}
std<-function(x){
  sd(x)/sqrt(length(x))
}

# t_0=0
table0<-matrix(NA,4,4)
rownames(table0)<-c(0,10,20,30)
colnames(table0)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
dummy0<-c()
dummy1<-c()

# t_0=0 & c=0%
for (i in 1:100){
sim0<-data.gen(c.0,0)
sim0[,2]<-100
fit0<-rq(sim0[,5]~sim0[,6],data=sim0,tau=0.5,method='fn')
dummy0[i]<-fit0$coefficient[1]
dummy1[i]<-fit0$coefficient[2]
}
table0[1,1]<-mean(dummy0)
table0[1,2]<-std(dummy0)
table0[1,3]<-mean(dummy1)
table0[1,4]<-std(dummy1)

# t_0=0 & c=10%
for (i in 1:100){
sim1<-data.gen(c.1,0)
sim1<-sim1[sim1$censored==1,]
fit1<-rq(sim1[,5]~sim1[,6],data=sim1,tau=0.5,method='fn')
dummy0[i]<-fit1$coefficient[1]
dummy1[i]<-fit1$coefficient[2]
}
table0[2,1]<-mean(dummy0)
table0[2,2]<-std(dummy0)
table0[2,3]<-mean(dummy1)
table0[2,4]<-std(dummy1)

# t_0=0 & c=20%
for (i in 1:100){
sim2<-data.gen(c.2,0)
sim2<-sim2[sim2$censored==1,]
fit2<-rq(sim2[,5]~sim2[,6],data=sim2,tau=0.5,method='fn')
dummy0[i]<-fit2$coefficient[1]
dummy1[i]<-fit2$coefficient[2]
}
table0[3,1]<-mean(dummy0)
table0[3,2]<-std(dummy0)
table0[3,3]<-mean(dummy1)
table0[3,4]<-std(dummy1)

# t_0=0 & c=30%
for (i in 1:100){
sim3<-data.gen(c.3,0)
sim3<-sim3[sim3$censored==1,]
fit3<-rq(sim3[,5]~sim3[,6],data=sim3,tau=0.5,method='fn')
dummy0[i]<-fit3$coefficient[1]
dummy1[i]<-fit3$coefficient[2]
}
table0[4,1]<-mean(dummy0)
table0[4,2]<-std(dummy0)
table0[4,3]<-mean(dummy1)
table0[4,4]<-std(dummy1)

table0


#### Comparison to linear regression ####
# t_0=0
table1<-matrix(NA,4,4)
rownames(table1)<-c(0,10,20,30)
colnames(table1)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
dummy0<-c()
dummy1<-c()

# t_0=0 & c=0%
for (i in 1:100){
  sim0<-data.gen(c.0,0)
  sim0[,2]<-100
  fit0<-lm(sim0[,5]~sim0[,6],data=sim0)
  dummy0[i]<-fit0$coefficient[1]
  dummy1[i]<-fit0$coefficient[2]
}
table1[1,1]<-mean(dummy0)
table1[1,2]<-std(dummy0)
table1[1,3]<-mean(dummy1)
table1[1,4]<-std(dummy1)

# t_0=0 & c=10%
for (i in 1:100){
  sim1<-data.gen(c.1,0)
  sim1<-sim1[sim1$censored==1,]
  fit1<-lm(sim1[,5]~sim1[,6],data=sim1)
  dummy0[i]<-fit1$coefficient[1]
  dummy1[i]<-fit1$coefficient[2]
}
table1[2,1]<-mean(dummy0)
table1[2,2]<-std(dummy0)
table1[2,3]<-mean(dummy1)
table1[2,4]<-std(dummy1)

# t_0=0 & c=20%
for (i in 1:100){
  sim2<-data.gen(c.2,0)
  sim2<-sim2[sim2$censored==1,]
  fit2<-lm(sim2[,5]~sim2[,6],data=sim2)
  dummy0[i]<-fit2$coefficient[1]
  dummy1[i]<-fit2$coefficient[2]
}
table1[3,1]<-mean(dummy0)
table1[3,2]<-std(dummy0)
table1[3,3]<-mean(dummy1)
table1[3,4]<-std(dummy1)

# t_0=0 & c=30%
for (i in 1:100){
  sim3<-data.gen(c.3,0)
  sim3<-sim3[sim3$censored==1,]
  fit3<-lm(sim3[,5]~sim3[,6],data=sim3)
  dummy0[i]<-fit3$coefficient[1]
  dummy1[i]<-fit3$coefficient[2]
}
table1[4,1]<-mean(dummy0)
table1[4,2]<-std(dummy0)
table1[4,3]<-mean(dummy1)
table1[4,4]<-std(dummy1)

table1

#### Inducded smoothing + Kim's article ####
#### object equation ####
objectF<-function(b){
  logz<-as.matrix(a.order[,5])
  b<-as.matrix(b)
  x0<-as.matrix(x[,1])
  x1<-as.matrix(x[,2])
  normal<-matrix(0,n,1)
  for (i in 1:n){
    normal[i,1]<-(0.5-1+pnorm((logz-x%*%b)[i,1]/(sqrt(t(x[i,])%*%G%*%(as.matrix(x[i,]))))))
  }
  result0<-t(w)%*%(I*x0*normal)
  result1<-t(w)%*%(I*x1*normal)
  print(c(result0,result1))
}

library(BB)
library(nleqslv)

#### Method 1 : BB Solver (200 data, 100 simulation) ####
table0<-matrix(NA,4,4)
rownames(table0)<-c(0,10,20,30)
colnames(table0)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
## Censoring proportion = 0%
## true betas are (1.61, 0)
b0.c0<-c()
b1.c0<-c()

for (m in 1:100){
  a<-data.gen(c.0,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calculation w_i with survfit
  a.order$censored.rev<-1+a.order$censored*-1
  a.km.fit<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order)
  a.order[,9]<-a.km.fit$surv
  for (i in 1:n)
  {
    if (a.order[i,9]==0){
      a.order[i,10]<-a.order[i-1,10]}
    else{
      a.order[i,10]<-a.order[i,7]/(a.order[i,9])}
  }
  w<-a.order[,10]
  ## Indicator Z_i>T_0
  I<-matrix(1,n,1)
  for (l in 1:n){
    if (a.order[l,4]>0)
    {
      I[l,1]=1
    }
    else {I[l,1]=0}
  }
  ## Matrix X_i (x[,1]=intercept, x[,2]=covariate)
  x<-matrix(1,n,2)
  x[,1]<-1
  x[,2]<-a.order[,6]
  
  # Defining Gamma (1/sqrt(n))I_p
  G<-matrix(c(1/sqrt(n),0,0,1/sqrt(n)),2,2)
  
  # Find estimated value of beta
  betastart<-c(1.61,0) # from table0
  z.c0<-BBsolve(betastart,objectF)
  b0.c0[m]<-z.c0$par[1]
  b1.c0[m]<-z.c0$par[2]
}

# find solution by BBsolve
table0[1,1]<-mean(b0.c0)
table0[1,2]<-std(b0.c0)
table0[1,3]<-mean(b1.c0)
table0[1,4]<-std(b1.c0)

## Censoring proportion = 10%
## true betas are (1.57, -0.02)

b0.c1<-c()
b1.c1<-c()

for (m in 1:100){
  a<-data.gen(c.1,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calculation w_i with survfit
  a.order$censored.rev<-1+a.order$censored*-1
  a.km.fit<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order)
  a.order[,9]<-a.km.fit$surv
  for (i in 1:n)
  {
    if (a.order[i,9]==0){
    a.order[i,10]<-a.order[i-1,10]}
    else{
    a.order[i,10]<-a.order[i,7]/(a.order[i,9])}
  }
  w<-a.order[,10]

  ## Indicator Z_i>T_0
  I<-matrix(1,n,1)
  for (l in 1:n){
    if (a.order[l,4]>0)
    {
      I[l,1]=1
    }
    else {I[l,1]=0}
  }
  ## Matrix X_i (x[,1]=intercept, x[,2]=covariate)
  x<-matrix(1,n,2)
  x[,1]<-1
  x[,2]<-a.order[,6]
  
  # Defining Gamma (1/sqrt(n))I_p
  G<-matrix(c(1/sqrt(n),0,0,1/sqrt(n)),2,2)
  
  # Find estimated value of beta
  betastart<-c(1.57,-.015)
  z.c1<-BBsolve(betastart,objectF)
  b0.c1[m]<-z.c1$par[1]  
  b1.c1[m]<-z.c1$par[2]
}
# find solution by BBsolve
table0[2,1]<-mean(b0.c2)
table0[2,2]<-std(b0.c2)
table0[2,3]<-mean(b1.c2)
table0[2,4]<-std(b1.c2)


## Censoring proportion = 20%
## true betas are (1.52, 0.006)
b0.c2<-c()
b1.c2<-c()

for (m in 1:100){
  a<-data.gen(c.2,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calculation w_i with survfit
  a.order$censored.rev<-1+a.order$censored*-1
  a.km.fit<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order)
  a.order[,9]<-a.km.fit$surv
  for (i in 1:n)
  {
    if (a.order[i,9]==0){
      a.order[i,10]<-a.order[i-1,10]}
    else{
      a.order[i,10]<-a.order[i,7]/(a.order[i,9])}
  }
  w<-a.order[,10]
  
  ## Indicator Z_i>T_0
  I<-matrix(1,n,1)
  for (l in 1:n){
    if (a.order[l,4]>0)
    {
      I[l,1]=1
    }
    else {I[l,1]=0}
  }
  ## Matrix X_i (x[,1]=intercept, x[,2]=covariate)
  x<-matrix(1,n,2)
  x[,1]<-1
  x[,2]<-a.order[,6]
  
  # Defining Gamma (1/sqrt(n))I_p
  G<-matrix(c(1/sqrt(n),0,0,1/sqrt(n)),2,2)
  
  # Find estimated value of beta
  betastart<-c(1.52,.006)
  z3<-BBsolve(betastart,objectF)
  b0.c2[m]<-z3$par[1]  
  b1.c2[m]<-z3$par[2]
}

# find solution by BBsolve
table0[3,1]<-mean(b0.c2)
table0[3,2]<-std(b0.c2)
table0[3,3]<-mean(b1.c2)
table0[3,4]<-std(b1.c2)

## Censoring proportion = 30%
## true betas are (1.47, -0.011)
b0.c3<-c()
b1.c3<-c()

for (m in 1:100){
  a<-data.gen(c.3,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calculation w_i with survfit
  a.order$censored.rev<-1+a.order$censored*-1
  a.km.fit<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order)
  a.order[,9]<-a.km.fit$surv
  for (i in 1:n)
  {
    if (a.order[i,9]==0){
      a.order[i,10]<-a.order[i-1,10]}
    else{
      a.order[i,10]<-a.order[i,7]/(a.order[i,9])}
  }
  w<-a.order[,10]
  
  ## Indicator Z_i>T_0
  I<-matrix(1,n,1)
  for (l in 1:n){
    if (a.order[l,4]>0)
    {
      I[l,1]=1
    }
    else {I[l,1]=0}
  }
  ## Matrix X_i (x[,1]=intercept, x[,2]=covariate)
  x<-matrix(1,n,2)
  x[,1]<-1
  x[,2]<-a.order[,6]
  
  # Defining Gamma (1/sqrt(n))I_p
  G<-matrix(c(1/sqrt(n),0,0,1/sqrt(n)),2,2)
  
  # Find estimated value of beta
  betastart<-c(1.48,0.016)
  z4<-BBsolve(betastart,objectF)
  b0.c3[m]<-z4$par[1]  
  b1.c3[m]<-z4$par[2]
}

# find solution by BBsolve
table0[4,1]<-mean(b0.c3)
table0[4,2]<-std(b0.c3)
table0[4,3]<-mean(b1.c3)
table0[4,4]<-std(b1.c3)

## solve without induced smoothing
# t_0=0
table0.xis<-matrix(NA,4,4)
rownames(table0.xis)<-c(0,10,20,30)
colnames(table0.xis)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
dummy0<-c()
dummy1<-c()

# t_0=0 & c=0%
for (i in 1:100){
  a0<-data.gen(c.0,0)
  a.order0<-a0[order(a0$Z.diff),]
  n<-nrow(a.order0)
  ## Calculation w_i with survfit
  a.order0$censored.rev<-1+a.order0$censored*-1
  a.km.fit0<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order0)
  a.order0$G_KM<-a.km.fit0$surv
  for (j in 1:n)
  {
    if (a.order0[j,9]==0)
    {
      a.order0[j,10]<-a.order0[j-1,10]
    }
    else{a.order0[j,10]<-a.order0[j,7]/(a.order0[j,9])}
  }
  fit0<-rq(a.order0[,5]~a.order0[,6],data=a.order0,tau=0.5, weight=a.order0[,10],method='fn')
  dummy0[i]<-fit0$coefficient[1]
  dummy1[i]<-fit0$coefficient[2]
}
table0.xis[1,1]<-mean(dummy0)
table0.xis[1,2]<-std(dummy0)
table0.xis[1,3]<-mean(dummy1)
table0.xis[1,4]<-std(dummy1)

# t_0=0 & c=10%
for (i in 1:100){
  a1<-data.gen(c.1,0)
  a.order1<-a1[order(a1$Z.diff),]
  n<-nrow(a.order1)
  ## Calculation w_i with survfit
  a.order1$censored.rev<-1+a.order1$censored*-1
  a.km.fit1<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order1)
  a.order1$G_KM<-a.km.fit1$surv
  for (j in 1:n)
  {
    if (a.order1[j,9]==0)
    {
      a.order1[j,10]<-a.order1[j-1,10]
    }
    else{a.order1[j,10]<-a.order1[j,7]/(a.order1[j,9])}
  }
  fit1<-rq(a.order1[,5]~a.order1[,6],data=a.order1,tau=0.5, weight=a.order1[,10],method='fn')
  dummy0[i]<-fit1$coefficient[1]
  dummy1[i]<-fit1$coefficient[2]
}
table0.xis[2,1]<-mean(dummy0)
table0.xis[2,2]<-std(dummy0)
table0.xis[2,3]<-mean(dummy1)
table0.xis[2,4]<-std(dummy1)


# t_0=0 & c=20%
for (i in 1:100){
  a2<-data.gen(c.2,0)
  a.order2<-a2[order(a2$Z.diff),]
  n<-nrow(a.order2)
  ## Calculation w_i with survfit
  a.order2$censored.rev<-1+a.order2$censored*-1
  a.km.fit2<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order2)
  a.order2$G_KM<-a.km.fit2$surv
  for (j in 1:n)
  {
    if (a.order2[j,9]==0)
    {
      a.order2[j,10]<-a.order2[j-1,10]
    }
    else{a.order2[j,10]<-a.order2[j,7]/(a.order2[j,9])}
  }
  fit2<-rq(a.order2[,5]~a.order2[,6],data=a.order2,tau=0.5, weight=a.order2[,10],method='fn')
  dummy0[i]<-fit2$coefficient[1]
  dummy1[i]<-fit2$coefficient[2]
}
table0.xis[3,1]<-mean(dummy0)
table0.xis[3,2]<-std(dummy0)
table0.xis[3,3]<-mean(dummy1)
table0.xis[3,4]<-std(dummy1)

# t_0=0 & c=30%
for (i in 1:100){
  a3<-data.gen(c.3,0)
  a.order3<-a3[order(a3$Z.diff),]
  n<-nrow(a.order3)
  ## Calculation w_i with survfit
  a.order3$censored.rev<-1+a.order3$censored*-1
  a.km.fit3<-survfit(formula = Surv(Z.diff,censored.rev) ~ 1, data = a.order3)
  a.order3$G_KM<-a.km.fit3$surv
  for (j in 1:n)
  {
    if (a.order3[j,9]==0)
    {
      a.order3[j,10]<-a.order3[j-1,10]
    }
    else{a.order3[j,10]<-a.order3[j,7]/(a.order3[j,9])}
  }
  fit3<-rq(a.order3[,5]~a.order3[,6],data=a.order3,tau=0.5, weight=a.order3[,10],method='fn')
  dummy0[i]<-fit3$coefficient[1]
  dummy1[i]<-fit3$coefficient[2]
}
table0.xis[4,1]<-mean(dummy0)
table0.xis[4,2]<-std(dummy0)
table0.xis[4,3]<-mean(dummy1)
table0.xis[4,4]<-std(dummy1)

table0
table0.xis
