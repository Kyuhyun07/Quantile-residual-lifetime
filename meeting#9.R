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

#Data Generation function

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

## Inducded smoothing + Kim's article

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

library(nleqslv)
library(BB)

#### 100 data w/ no censoring####
## BBsolve
b.20<-c()
b.210<-c()

for (m in 1:100){
  a<-data.gen(c.0,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calculation w_i with survfit
  w<-matrix(1,200,1)
  
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
  z2<-BBsolve(betastart,objectF)
  b.20[m]<-z2$par[1]
  b.210[m]<-z2$par[2]
}

# find solution by BBsolve
est.20<-c(sum(b.20)/length(b.20),sum(b.210)/length(b.210))

## previous beta are 1.61, 0

#### 100 data w/ 10% censoring####
## BBsolve
b.2<-c()
b.21<-c()

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
    if (a.order[i,9]==1)
    {
      a.order[i,10]<-0
    }
    else{a.order[i,10]<-a.order[i,7]/(1-a.order[i,9])}
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
  betastart<-c(1.58,-.012)
  z2<-BBsolve(betastart,objectF)
  b.2[m]<-z2$par[1]  
  b.21[m]<-z2$par[2]
}

# find solution by BBsolve
est.21<-c(sum(b.2)/length(b.2),sum(b.21)/length(b.21))
est.21
b.2


## previous beta are 1.564346 0.005893996
