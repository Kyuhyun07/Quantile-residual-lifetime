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
r.initial=log(2)^(1/k)/exp.beta.initial

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
  sim<-matrix(NA,200,11)
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
  colnames(sim)<-c("T","C","Z","Z.diff","log(Z.diff)","X","censored","n","d","prod","G_KM")
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
fit0<-rq(log(sim0[,1])~sim0[,6],data=sim0,tau=0.5)
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
fit1<-rq(sim1[,5]~sim1[,6],data=sim1,tau=0.5)
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
fit2<-rq(sim2[,5]~sim2[,6],data=sim2,tau=0.5)
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
fit3<-rq(sim3[,5]~sim3[,6],data=sim3,tau=0.5)
dummy0[i]<-fit3$coefficient[1]
dummy1[i]<-fit3$coefficient[2]
}
table0[4,1]<-mean(dummy0)
table0[4,2]<-std(dummy0)
table0[4,3]<-mean(dummy1)
table0[4,4]<-std(dummy1)

table0

## Inducded smoothing + Kim's article

#### one data ####
## Making data censoring 10%
a<-data.gen(c.1,0)
a.order<-a[order(a$Z.diff),]
n<-nrow(a.order)
## Calcuation W_i
for (i in 1:n){
  a.order[i,8]<-sum(a.order[,3]>=a.order[i,3])
  a.order[i,9]<-sum(a.order[,3]==a.order[i,3])
  if (a.order[i,7]==0){
    a.order[i,10]=(a.order[i,8]-a.order[i,9])/a.order[i,8]
  }
  else {
    a.order[i,10]=1
  }
  if (i==1){
    a.order[i,11]<-a.order[i,10]
  }
  else{
  a.order[i,11]<-a.order[i-1,11]*a.order[i,10]
  }
}

w<-matrix(NA,n,1)
for (i in 1:n){
  if (a.order[i,11]==1){
    w[i,1]=1
  }
  else{
    w[i,1]<-a.order[i,7]/((1-a.order[i,11])*n)
}
}
## make weight
w<-w/sum(w)

## Indicator Z_i>T_0
I<-matrix(1,n,1)
for (i in 1:200){
  if (a.order[i,4]>0)
    {
    I[i,1]=1
  }
  else {I[i,1]=0}
}

## Matrix X_i (x[,1]=intercept, x[,2]=covariate)
x<-matrix(1,200,2)
x[,1]<-1
x[,2]<-a.order[,6]

# Defining Gamma (1/sqrt(n))I_p
G<-matrix(c(1/sqrt(n),0,0,1/sqrt(n)),2,2)

#### object equation ####
objectF<-function(b){
logz<-as.matrix(a.order[,5])
b<-as.matrix(b)
x0<-as.matrix(x[,1])
x1<-as.matrix(x[,2])
normal<-matrix(0,n,1)
for (i in 1:n){
  normal[i,1]<-(0.5-1+pnorm((logz-x%*%b)[i,1]/(sqrt(t(x[i,])%*%G%*%x[i,]))))
}
result0<-t(w)%*%(I*x0*normal)
result1<-t(w)%*%(I*x1*normal)
print(c(result0,result1))
}

library(nleqslv)
library(BB)

betastart<-c(1.525536,0.006749514)

z<-nleqslv(betastart,objectF,control=list(btol=.01))

#### 100 data w/ 10% censoring####
## nleqslv
b.1<-c()
b.11<-c()

## BBsolve
b.2<-c()
b.21<-c()

## dfsane
b.3<-c()
b.31<-c()

for (m in 1:100){
  a<-data.gen(c.1,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calcuation W_i
  for (i in 1:n){
    a.order[i,8]<-sum(a.order[,3]>=a.order[i,3])
    a.order[i,9]<-sum(a.order[,3]==a.order[i,3])
    if (a.order[i,7]==0){
      a.order[i,10]=(a.order[i,8]-a.order[i,9])/a.order[i,8]
    }
    else {
      a.order[i,10]=1
    }
    if (i==1){
      a.order[i,11]<-a.order[i,10]
    }
    else{
      a.order[i,11]<-a.order[i-1,11]*a.order[i,10]
    }
  }
  w<-matrix(NA,n,1)
  for (i in 1:n){
    if (a.order[i,11]==1){
      w[i,1]=1
    }
    else{
      w[i,1]<-a.order[i,7]/((1-a.order[i,11])*n)
    }
  }
  # make weight
  w<-w/sum(w)
  
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
  betastart<-c(1.537047,-0.0006810171)
  z1<-nleqslv(betastart,objectF,control=list(btol=.01))
  z2<-BBsolve(betastart,objectF)
  z3<-dfsane(betastart,objectF)
  
  b.1[m]<-z1$x[1]
  b.2[m]<-z2$par[1]  
  b.3[m]<-z3$par[1]   
  b.11[m]<-z1$x[2]
  b.21[m]<-z2$par[2]
  b.31[m]<-z3$par[2] 
}

# find solution by nleqslv 
c(sum(b.1)/length(b.1),sum(b.11)/length(b.11))
# find solution by BBsolve
c(sum(b.2)/length(b.2),sum(b.21)/length(b.21))
# find solution by dfsane
c(sum(b.3)/length(b.3),sum(b.31)/length(b.31))

## previous beta are 1.537047 -0.0006810171

#### 100 data w/ no censoring####
## nleqslv
b.1<-c()
b.11<-c()

## BBsolve
b.2<-c()
b.21<-c()

## dfsane
b.3<-c()
b.31<-c()

for (m in 1:100){
  a<-data.gen(c.0,0)
  a.order<-a[order(a$Z.diff),]
  n<-nrow(a.order)
  ## Calcuation W_i
  for (i in 1:n){
    a.order[i,8]<-sum(a.order[,3]>=a.order[i,3])
    a.order[i,9]<-sum(a.order[,3]==a.order[i,3])
    if (a.order[i,7]==0){
      a.order[i,10]=(a.order[i,8]-a.order[i,9])/a.order[i,8]
    }
    else {
      a.order[i,10]=1
    }
    if (i==1){
      a.order[i,11]<-a.order[i,10]
    }
    else{
      a.order[i,11]<-a.order[i-1,11]*a.order[i,10]
    }
  }
  w<-matrix(NA,n,1)
  for (i in 1:n){
    if (a.order[i,11]==1){
      w[i,1]=1
    }
    else{
      w[i,1]<-a.order[i,7]/((1-a.order[i,11])*n)
    }
  }
  # make weight
  w<-w/sum(w)
  
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
  betastart<-c(1.611372,-0.0031985300) # from table0
  z1<-nleqslv(betastart,objectF,control=list(btol=.01))
  z2<-BBsolve(betastart,objectF)
  z3<-dfsane(betastart,objectF)
  
  b.1[m]<-z1$x[1]
  b.2[m]<-z2$par[1]
  b.3[m]<-z3$par[1]   
  b.11[m]<-z1$x[2]
  b.21[m]<-z2$par[2]
  b.31[m]<-z3$par[2] 
}

# find solution by nleqslv 
c(sum(b.1)/length(b.1),sum(b.11)/length(b.11))
# find solution by BBsolve
c(sum(b.2)/length(b.2),sum(b.21)/length(b.21))
# find solution by dfsane
c(sum(b.3)/length(b.3),sum(b.31)/length(b.31))

## previous beta are 1.611372 -0.0031985300
