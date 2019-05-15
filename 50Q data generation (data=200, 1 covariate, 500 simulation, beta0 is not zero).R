library(quantreg)
library(survival)
library(nleqslv)
library(tictoc)
library(xtable)
  
### Given information ####
##Find C
##Generate T_i 1,000,000
#d.exp.beta.initial.0=5
#d.exp.beta.initial.1=25
#d.k=2
#d.r.initial.0=(log(2)^(1/d.k))/d.exp.beta.initial.0
#d.r.initial.1=(log(2)^(1/d.k))/d.exp.beta.initial.1
#d.u=runif(n=10000000,min = 0,max = 1)
#d.x=rbinom(10000000,size=1,p=0.5)
#d.T<-c()
#for (q in 1:10000000){
#  if (d.x[q]==0){
#    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.0
#  } else {
#    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.1
#  }
#}
##Find c which is derermined to achieve 0,10,20,30% censoring rate in case 1
#i=52.68
#while(1){
#  d.C<-runif(10000000,0,i)
#  i=i+0.01
#  if(sum(d.C<d.T)<=2000000) break # 0% = 0, 10% = 100000, 20% = 200000, 30% = 300000)
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

## True Beta
# When T_0=0, beta_0=1.61, beta_1=0
# When T_0=1, beta_0=1.41, beta_1=0
# When T_0=2, beta_0=1.22, beta_1=0
# When T_0=3, beta_0=1.04, beta_1=0


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

#### t_0=0 & c=0% ####
set.seed(1)
for (i in 1:500){
  a<-data.gen(c.0,0)
  a[,2]<-100
}
#### t_0=0 & c=10% ####
set.seed(1)
for (i in 1:500){
  a<-data.gen(c.1,0)
}

#### t_0=0 & c=20% ####
set.seed(1)
for (i in 1:500){
  a<-data.gen(c.2,0)
}

#### t_0=0 & c=30% ####
set.seed(1)
for (i in 1:500){
  a<-data.gen(c.3,0)
}

#### t_0=1 & c=0% ####
set.seed(11)
for (i in 1:500){
  a<-data.gen(c.0,1)
  a[,2]<-100
}

#### t_0=1 & c=10% ####
set.seed(11)
for (i in 1:500){
  a<-data.gen(c.1,1)
}

#### t_0=1 & c=20% ####
set.seed(11)
for (i in 1:500){
  a<-data.gen(c.2,1)
}

#### t_0=1 & c=30% ####
set.seed(11)
for (i in 1:500){
  a<-data.gen(c.3,1)
}

#### t_0=2 & c=0% ####
set.seed(21)
for (i in 1:500){
  a<-data.gen(c.0,2)
  a[,2]<-100
}

#### t_0=2 & c=10% ####
set.seed(21)
for (i in 1:500){
  a<-data.gen(c.1,2)
}

#### t_0=2 & c=20% ####
set.seed(21)
for (i in 1:500){
  a<-data.gen(c.2,2)
}

#### t_0=2 & c=30% ####
set.seed(21)
for (i in 1:500){
  a<-data.gen(c.3,2)
}

#### t_0=3 & c=0% ####
set.seed(31)
for (i in 1:500){
  a<-data.gen(c.0,3)
  a[,2]<-100
}

#### t_0=3 & c=10% ####
set.seed(31)
for (i in 1:500){
  a<-data.gen(c.1,3)
}

#### t_0=3 & c=20% ####
set.seed(31)
for (i in 1:500){
  a<-data.gen(c.2,3)
}

#### t_0=3 & c=30% ####
set.seed(31)
for (i in 1:500){
  a<-data.gen(c.3,3)
}