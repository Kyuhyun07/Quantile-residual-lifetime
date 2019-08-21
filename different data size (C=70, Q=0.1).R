#### Condition ####
# data size = 200 / 400 / 800 / 1600 / 3200
# beta0, beta1 effective
# Quantile 10%
# Censoring = 70%
# simulation 200
# eta = 100

library(quantreg)
library(survival)
library(nleqslv)
library(xtable)

#### True Beta ####
#beta_0    beta_1
#t_0=0 1.609438 0.6931472
#t_0=1 1.410748 0.7974189
#t_0=2 1.219403 0.9070615
#t_0=3 1.040613 1.0174711

#### Find censoring point ####
d.exp.beta.initial.0=5
d.exp.beta.initial.1=10
d.k=2
d.r.initial.0=(log(10/9)^(1/d.k))/d.exp.beta.initial.0
d.r.initial.1=(log(10/9)^(1/d.k))/d.exp.beta.initial.1
d.u=runif(n=1000000,min = 0,max = 1)
d.x=rbinom(1000000,size=1,p=0.5)
d.T<-c()
for (q in 1:1000000){
  if (d.x[q]==0){
    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.0
  } else {
    d.T[q]={{-log(1-d.u[q])}^(1/d.k)}/d.r.initial.1
  }
}
##Find c which is derermined to achieve 0,10,20,30% censoring rate in case 1
i=20
censor.t3=c()
for (m in 1:2){
  while(1){
    d.C<-runif(1000000,0,i)
    dummy = matrix(NA, 1000000, 3)
    dummy[,1] = d.C
    dummy[,2] = d.T        
    dummy[,3] = apply(dummy[,1:2], 1, FUN=min)
    newdummy = dummy[dummy[,3]>=0,]
    i=i+0.01
    if(sum(newdummy[,1]<newdummy[,2])/nrow(newdummy)<0.1*(9-2*m)) break
  }
  print(i)
  table(newdummy[,1]<newdummy[,2])[2]/nrow(newdummy)
  censor.t3[m]=i
}

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
  
  colnames(data)[1:2]=c("Y-t_0", "log(y-t_0)")
  covar = paste("covariate",1:nc,sep="")
  colnames(data)[3:(nc+2)]=covar
  #for (j in 3:(nc+2)){
  #  colnames(data)[j]="covariate"
  #}
  colnames(data)[(nc+3)] = c("delta")
  
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
  G = diag(1/n, nc+1, nc+1)
  
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

#### tables ####
table200.crq<-matrix(NA,4,4)
rownames(table200.crq)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table200.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table400.crq<-matrix(NA,4,4)
rownames(table400.crq)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table400.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table800.crq<-matrix(NA,4,4)
rownames(table800.crq)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table800.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1600.crq<-matrix(NA,4,4)
rownames(table1600.crq)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table1600.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3200.crq<-matrix(NA,4,4)
rownames(table3200.crq)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table3200.crq)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

table200.is<-matrix(NA,4,4)
rownames(table200.is)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table200.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table400.is<-matrix(NA,4,4)
rownames(table400.is)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table400.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table800.is<-matrix(NA,4,4)
rownames(table800.is)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table800.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table1600.is<-matrix(NA,4,4)
rownames(table1600.is)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table1600.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")
table3200.is<-matrix(NA,4,4)
rownames(table3200.is)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table3200.is)<-c("beta_0","SE of beta_0","beta_1","SE of beta_1")

table.event<-matrix(NA,4,5)
rownames(table.event)<-c("t_0=0","t_0=1","t_0=2","t_0=3")
colnames(table.event)<-c("data=200","data=400","data=800","data=1600","data=3200")


#### Censoring Q=0.1, t_0=0####
c.0=5000000
#c.1=
#c.3=
c.5=38.94
c.7=23.48

#### My assumption ####
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/9))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/9))^(1/k)/exp.beta.initial.1

#### t_0=0 & c=70% & data = 200####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
avg.event.0.200<-c()
set.seed(1)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(200,c.7,0)
  avg.event.0.200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.10)
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
for (i in 1:200){
  a<-data.gen(200,c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.10, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
    })
  }
# Crq beta table
table200.crq[1,1]<-mean(b0.crq.07,na.rm=TRUE)
table200.crq[1,2]<-sd(b0.crq.07,na.rm=TRUE)
table200.crq[1,3]<-mean(b1.crq.07,na.rm=TRUE)
table200.crq[1,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table200.is[1,1]<-mean(b0.is.07,na.rm=TRUE)
table200.is[1,2]<-sd(b0.is.07,na.rm=TRUE)
table200.is[1,3]<-mean(b1.is.07,na.rm=TRUE)
table200.is[1,4]<-sd(b1.is.07,na.rm=TRUE)

table.event[1,1] = mean(avg.event.0.200)

#### t_0=0 & c=70% & data = 400####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
avg.event.0.400<-c()
set.seed(1)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(400,c.7,0)
  avg.event.0.400[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.10)
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
for (i in 1:200){
  a<-data.gen(400,c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.10, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
    })
}
# Crq beta table
table400.crq[1,1]<-mean(b0.crq.07,na.rm=TRUE)
table400.crq[1,2]<-sd(b0.crq.07,na.rm=TRUE)
table400.crq[1,3]<-mean(b1.crq.07,na.rm=TRUE)
table400.crq[1,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table400.is[1,1]<-mean(b0.is.07,na.rm=TRUE)
table400.is[1,2]<-sd(b0.is.07,na.rm=TRUE)
table400.is[1,3]<-mean(b1.is.07,na.rm=TRUE)
table400.is[1,4]<-sd(b1.is.07,na.rm=TRUE)

table.event[1,2] = mean(avg.event.0.400)

#### t_0=0 & c=70% & data = 800####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
avg.event.0.800<-c()
set.seed(1)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(800,c.7,0)
  avg.event.0.800[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.10)
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
for (i in 1:200){
  a<-data.gen(800,c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.10, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
    })
}
# Crq beta table
table800.crq[1,1]<-mean(b0.crq.07,na.rm=TRUE)
table800.crq[1,2]<-sd(b0.crq.07,na.rm=TRUE)
table800.crq[1,3]<-mean(b1.crq.07,na.rm=TRUE)
table800.crq[1,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table800.is[1,1]<-mean(b0.is.07,na.rm=TRUE)
table800.is[1,2]<-sd(b0.is.07,na.rm=TRUE)
table800.is[1,3]<-mean(b1.is.07,na.rm=TRUE)
table800.is[1,4]<-sd(b1.is.07,na.rm=TRUE)

table.event[1,3] = mean(avg.event.0.800)

#### t_0=0 & c=70% & data = 1600####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
avg.event.0.1600<-c()
set.seed(1)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(1600,c.7,0)
  avg.event.0.1600[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.10)
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
for (i in 1:200){
  a<-data.gen(1600,c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.10, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
    })
}
# Crq beta table
table1600.crq[1,1]<-mean(b0.crq.07,na.rm=TRUE)
table1600.crq[1,2]<-sd(b0.crq.07,na.rm=TRUE)
table1600.crq[1,3]<-mean(b1.crq.07,na.rm=TRUE)
table1600.crq[1,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table1600.is[1,1]<-mean(b0.is.07,na.rm=TRUE)
table1600.is[1,2]<-sd(b0.is.07,na.rm=TRUE)
table1600.is[1,3]<-mean(b1.is.07,na.rm=TRUE)
table1600.is[1,4]<-sd(b1.is.07,na.rm=TRUE)

table.event[1,4] = mean(avg.event.0.1600)

#### t_0=0 & c=70% & data = 3200####
b0.crq.07<-c()
b0.crq.sd.07<-c()
b1.crq.07<-c()
b1.crq.sd.07<-c()
b0.is.07<-c()
b0.is.sd.07<-c()
b1.is.07<-c()
b1.is.sd.07<-c()
avg.event.0.3200<-c()
set.seed(1)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(3200,c.7,0)
  avg.event.0.3200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 0, 0.10)
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
for (i in 1:200){
  a<-data.gen(3200,c.7,0)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 0, 0.10, 100)
    b0.is.07[i] = ismb.fit[1,1]
    b0.is.sd.07[i] = ismb.fit[1,2]
    b1.is.07[i] = ismb.fit[2,1]
    b1.is.sd.07[i] = ismb.fit[2,2]}
    , error=function(e){
      b0.is.07[i] = NA
      b0.is.sd.07[i] = NA
      b1.is.07[i] = NA
      b1.is.sd.07[i] = NA
    })
}
# Crq beta table
table3200.crq[1,1]<-mean(b0.crq.07,na.rm=TRUE)
table3200.crq[1,2]<-sd(b0.crq.07,na.rm=TRUE)
table3200.crq[1,3]<-mean(b1.crq.07,na.rm=TRUE)
table3200.crq[1,4]<-sd(b1.crq.07,na.rm=TRUE)

# IS beta table
table3200.is[1,1]<-mean(b0.is.07,na.rm=TRUE)
table3200.is[1,2]<-sd(b0.is.07,na.rm=TRUE)
table3200.is[1,3]<-mean(b1.is.07,na.rm=TRUE)
table3200.is[1,4]<-sd(b1.is.07,na.rm=TRUE)

table.event[1,5] = mean(avg.event.0.3200)


#### Given Information Q=0.1, t_0=1####
c.0=5000000
c.1=195.41
c.3=65.89
c.5=37.71
c.7=22.75

#### t_0=1 & c=70% & data = 200####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
avg.event.1.200<-c()
set.seed(11)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(200,c.7,1)
  avg.event.1.200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.10)
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
for (i in 1:200){
  a<-data.gen(200,c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.10, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
    }
    , error=function(e){
      b0.is.17[i] = NA
      b0.is.sd.17[i] = NA
      b1.is.17[i] = NA
      b1.is.sd.17[i] = NA
    })
  }

# Crq beta table
table200.crq[2,1]<-mean(b0.crq.17,na.rm=TRUE)
table200.crq[2,2]<-sd(b0.crq.17,na.rm=TRUE)
table200.crq[2,3]<-mean(b1.crq.17,na.rm=TRUE)
table200.crq[2,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table200.is[2,1]<-mean(b0.is.17,na.rm=TRUE)
table200.is[2,2]<-sd(b0.is.17,na.rm=TRUE)
table200.is[2,3]<-mean(b1.is.17,na.rm=TRUE)
table200.is[2,4]<-sd(b1.is.17,na.rm=TRUE)

table.event[2,1] = mean(avg.event.1.200)

#### t_0=1 & c=70% & data = 400####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
avg.event.1.400<-c()
set.seed(11)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(400,c.7,1)
  avg.event.1.400[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.10)
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
for (i in 1:200){
  a<-data.gen(400,c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.10, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
    }
    , error=function(e){
      b0.is.17[i] = NA
      b0.is.sd.17[i] = NA
      b1.is.17[i] = NA
      b1.is.sd.17[i] = NA
    })
  }

# Crq beta table
table400.crq[2,1]<-mean(b0.crq.17,na.rm=TRUE)
table400.crq[2,2]<-sd(b0.crq.17,na.rm=TRUE)
table400.crq[2,3]<-mean(b1.crq.17,na.rm=TRUE)
table400.crq[2,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table400.is[2,1]<-mean(b0.is.17,na.rm=TRUE)
table400.is[2,2]<-sd(b0.is.17,na.rm=TRUE)
table400.is[2,3]<-mean(b1.is.17,na.rm=TRUE)
table400.is[2,4]<-sd(b1.is.17,na.rm=TRUE)

table.event[2,2] = mean(avg.event.1.400)

#### t_0=1 & c=70% & data = 800####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
avg.event.1.800<-c()
set.seed(11)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(800,c.7,1)
  avg.event.1.800[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.10)
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
for (i in 1:200){
  a<-data.gen(800,c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.10, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.17[i] = NA
    b0.is.sd.17[i] = NA
    b1.is.17[i] = NA
    b1.is.sd.17[i] = NA
  })
}

# Crq beta table
table800.crq[2,1]<-mean(b0.crq.17,na.rm=TRUE)
table800.crq[2,2]<-sd(b0.crq.17,na.rm=TRUE)
table800.crq[2,3]<-mean(b1.crq.17,na.rm=TRUE)
table800.crq[2,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table800.is[2,1]<-mean(b0.is.17,na.rm=TRUE)
table800.is[2,2]<-sd(b0.is.17,na.rm=TRUE)
table800.is[2,3]<-mean(b1.is.17,na.rm=TRUE)
table800.is[2,4]<-sd(b1.is.17,na.rm=TRUE)

table.event[2,3] = mean(avg.event.1.800)

#### t_0=1 & c=70% & data = 1600####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
avg.event.1.1600<-c()
set.seed(11)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(1600,c.7,1)
  avg.event.1.1600[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.10)
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
for (i in 1:200){
  a<-data.gen(1600,c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.10, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.17[i] = NA
    b0.is.sd.17[i] = NA
    b1.is.17[i] = NA
    b1.is.sd.17[i] = NA
  })
}

# Crq beta table
table1600.crq[2,1]<-mean(b0.crq.17,na.rm=TRUE)
table1600.crq[2,2]<-sd(b0.crq.17,na.rm=TRUE)
table1600.crq[2,3]<-mean(b1.crq.17,na.rm=TRUE)
table1600.crq[2,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table1600.is[2,1]<-mean(b0.is.17,na.rm=TRUE)
table1600.is[2,2]<-sd(b0.is.17,na.rm=TRUE)
table1600.is[2,3]<-mean(b1.is.17,na.rm=TRUE)
table1600.is[2,4]<-sd(b1.is.17,na.rm=TRUE)

table.event[2,4] = mean(avg.event.1.1600)

#### t_0=1 & c=70% & data = 3200####
b0.crq.17<-c()
b0.crq.sd.17<-c()
b1.crq.17<-c()
b1.crq.sd.17<-c()
b0.is.17<-c()
b0.is.sd.17<-c()
b1.is.17<-c()
b1.is.sd.17<-c()
avg.event.1.3200<-c()
set.seed(11)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(3200,c.7,1)
  avg.event.1.3200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 1, 0.10)
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
for (i in 1:200){
  a<-data.gen(3200,c.7,1)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 1, 0.10, 100)
    b0.is.17[i] = ismb.fit[1,1]
    b0.is.sd.17[i] = ismb.fit[1,2]
    b1.is.17[i] = ismb.fit[2,1]
    b1.is.sd.17[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.17[i] = NA
    b0.is.sd.17[i] = NA
    b1.is.17[i] = NA
    b1.is.sd.17[i] = NA
  })
}

# Crq beta table
table3200.crq[2,1]<-mean(b0.crq.17,na.rm=TRUE)
table3200.crq[2,2]<-sd(b0.crq.17,na.rm=TRUE)
table3200.crq[2,3]<-mean(b1.crq.17,na.rm=TRUE)
table3200.crq[2,4]<-sd(b1.crq.17,na.rm=TRUE)

# IS beta table
table3200.is[2,1]<-mean(b0.is.17,na.rm=TRUE)
table3200.is[2,2]<-sd(b0.is.17,na.rm=TRUE)
table3200.is[2,3]<-mean(b1.is.17,na.rm=TRUE)
table3200.is[2,4]<-sd(b1.is.17,na.rm=TRUE)

table.event[2,5] = mean(avg.event.1.3200)



#### Given Information Q=0.1, t_0=2####
c.0=5000000
c.1=188.21
c.3=64.11
c.5=36.89
c.7=22.24

#### t_0=2 & c=70% & data = 200####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
avg.event.2.200<-c()
set.seed(21)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(200,c.7,2)
  avg.event.2.200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.10)
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
for (i in 1:200){
  a<-data.gen(200,c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.10, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.27[i] = NA
    b0.is.sd.27[i] = NA
    b1.is.27[i] = NA
    b1.is.sd.27[i] = NA
  })
}

# Crq beta table
table200.crq[3,1]<-mean(b0.crq.27,na.rm=TRUE)
table200.crq[3,2]<-sd(b0.crq.27,na.rm=TRUE)
table200.crq[3,3]<-mean(b1.crq.27,na.rm=TRUE)
table200.crq[3,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table200.is[3,1]<-mean(b0.is.27,na.rm=TRUE)
table200.is[3,2]<-sd(b0.is.27,na.rm=TRUE)
table200.is[3,3]<-mean(b1.is.27,na.rm=TRUE)
table200.is[3,4]<-sd(b1.is.27,na.rm=TRUE)

table.event[3,1] = mean(avg.event.2.200)

#### t_0=2 & c=70% & data = 400####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
avg.event.2.400<-c()
set.seed(21)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(400,c.7,2)
  avg.event.2.400[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.10)
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
for (i in 1:200){
  a<-data.gen(400,c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.10, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.27[i] = NA
    b0.is.sd.27[i] = NA
    b1.is.27[i] = NA
    b1.is.sd.27[i] = NA
  })
}

# Crq beta table
table400.crq[3,1]<-mean(b0.crq.27,na.rm=TRUE)
table400.crq[3,2]<-sd(b0.crq.27,na.rm=TRUE)
table400.crq[3,3]<-mean(b1.crq.27,na.rm=TRUE)
table400.crq[3,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table400.is[3,1]<-mean(b0.is.27,na.rm=TRUE)
table400.is[3,2]<-sd(b0.is.27,na.rm=TRUE)
table400.is[3,3]<-mean(b1.is.27,na.rm=TRUE)
table400.is[3,4]<-sd(b1.is.27,na.rm=TRUE)

table.event[3,2] = mean(avg.event.2.400)

#### t_0=2 & c=70% & data = 800####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
avg.event.2.800<-c()
set.seed(21)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(800,c.7,2)
  avg.event.2.800[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.10)
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
for (i in 1:200){
  a<-data.gen(800,c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.10, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.27[i] = NA
    b0.is.sd.27[i] = NA
    b1.is.27[i] = NA
    b1.is.sd.27[i] = NA
  })
}

# Crq beta table
table800.crq[3,1]<-mean(b0.crq.27,na.rm=TRUE)
table800.crq[3,2]<-sd(b0.crq.27,na.rm=TRUE)
table800.crq[3,3]<-mean(b1.crq.27,na.rm=TRUE)
table800.crq[3,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table800.is[3,1]<-mean(b0.is.27,na.rm=TRUE)
table800.is[3,2]<-sd(b0.is.27,na.rm=TRUE)
table800.is[3,3]<-mean(b1.is.27,na.rm=TRUE)
table800.is[3,4]<-sd(b1.is.27,na.rm=TRUE)

table.event[3,3] = mean(avg.event.2.800)

#### t_0=2 & c=70% & data = 1600####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
avg.event.2.1600<-c()
set.seed(21)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(1600,c.7,2)
  avg.event.2.1600[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.10)
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
for (i in 1:200){
  a<-data.gen(1600,c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.10, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.27[i] = NA
    b0.is.sd.27[i] = NA
    b1.is.27[i] = NA
    b1.is.sd.27[i] = NA
  })
}

# Crq beta table
table1600.crq[3,1]<-mean(b0.crq.27,na.rm=TRUE)
table1600.crq[3,2]<-sd(b0.crq.27,na.rm=TRUE)
table1600.crq[3,3]<-mean(b1.crq.27,na.rm=TRUE)
table1600.crq[3,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table1600.is[3,1]<-mean(b0.is.27,na.rm=TRUE)
table1600.is[3,2]<-sd(b0.is.27,na.rm=TRUE)
table1600.is[3,3]<-mean(b1.is.27,na.rm=TRUE)
table1600.is[3,4]<-sd(b1.is.27,na.rm=TRUE)

table.event[3,4] = mean(avg.event.2.1600)

#### t_0=2 & c=70% & data = 3200####
b0.crq.27<-c()
b0.crq.sd.27<-c()
b1.crq.27<-c()
b1.crq.sd.27<-c()
b0.is.27<-c()
b0.is.sd.27<-c()
b1.is.27<-c()
b1.is.sd.27<-c()
avg.event.2.3200<-c()
set.seed(21)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(3200,c.7,2)
  avg.event.2.3200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 2, 0.10)
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
for (i in 1:200){
  a<-data.gen(3200,c.7,2)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 2, 0.10, 100)
    b0.is.27[i] = ismb.fit[1,1]
    b0.is.sd.27[i] = ismb.fit[1,2]
    b1.is.27[i] = ismb.fit[2,1]
    b1.is.sd.27[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.27[i] = NA
    b0.is.sd.27[i] = NA
    b1.is.27[i] = NA
    b1.is.sd.27[i] = NA
  })
}

# Crq beta table
table3200.crq[3,1]<-mean(b0.crq.27,na.rm=TRUE)
table3200.crq[3,2]<-sd(b0.crq.27,na.rm=TRUE)
table3200.crq[3,3]<-mean(b1.crq.27,na.rm=TRUE)
table3200.crq[3,4]<-sd(b1.crq.27,na.rm=TRUE)

# IS beta table
table3200.is[3,1]<-mean(b0.is.27,na.rm=TRUE)
table3200.is[3,2]<-sd(b0.is.27,na.rm=TRUE)
table3200.is[3,3]<-mean(b1.is.27,na.rm=TRUE)
table3200.is[3,4]<-sd(b1.is.27,na.rm=TRUE)

table.event[3,5] = mean(avg.event.2.3200)

#### Given Information Q=0.1, t_0=3####
c.0=5000000
c.1=181.05
c.3=62.4
c.5=36.06
c.7=21.84

#### t_0=3 & c=70% & data = 200####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
avg.event.3.200<-c()
set.seed(31)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(200,c.7,3)
  avg.event.3.200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.10)
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
for (i in 1:200){
  a<-data.gen(200,c.7,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.10, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.37[i] = NA
    b0.is.sd.37[i] = NA
    b1.is.37[i] = NA
    b1.is.sd.37[i] = NA
  })
}

# Crq beta table
table200.crq[4,1]<-mean(b0.crq.37,na.rm=TRUE)
table200.crq[4,2]<-sd(b0.crq.37,na.rm=TRUE)
table200.crq[4,3]<-mean(b1.crq.37,na.rm=TRUE)
table200.crq[4,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table200.is[4,1]<-mean(b0.is.37,na.rm=TRUE)
table200.is[4,2]<-sd(b0.is.37,na.rm=TRUE)
table200.is[4,3]<-mean(b1.is.37,na.rm=TRUE)
table200.is[4,4]<-sd(b1.is.37,na.rm=TRUE)

table.event[4,1] = mean(avg.event.3.200)

#### t_0=3 & c=70% & data = 400####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
avg.event.3.400<-c()
set.seed(31)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(400,c.7,3)
  avg.event.3.400[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.10)
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
for (i in 1:200){
  a<-data.gen(400,c.7,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.10, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.37[i] = NA
    b0.is.sd.37[i] = NA
    b1.is.37[i] = NA
    b1.is.sd.37[i] = NA
  })
}

# Crq beta table
table400.crq[4,1]<-mean(b0.crq.37,na.rm=TRUE)
table400.crq[4,2]<-sd(b0.crq.37,na.rm=TRUE)
table400.crq[4,3]<-mean(b1.crq.37,na.rm=TRUE)
table400.crq[4,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table400.is[4,1]<-mean(b0.is.37,na.rm=TRUE)
table400.is[4,2]<-sd(b0.is.37,na.rm=TRUE)
table400.is[4,3]<-mean(b1.is.37,na.rm=TRUE)
table400.is[4,4]<-sd(b1.is.37,na.rm=TRUE)

table.event[4,2] = mean(avg.event.3.400)

#### t_0=3 & c=70% & data = 800####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
avg.event.3.800<-c()
set.seed(31)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(800,c.7,3)
  avg.event.3.800[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.10)
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
for (i in 1:200){
  a<-data.gen(800,c.7,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.10, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.37[i] = NA
    b0.is.sd.37[i] = NA
    b1.is.37[i] = NA
    b1.is.sd.37[i] = NA
  })
}

# Crq beta table
table800.crq[4,1]<-mean(b0.crq.37,na.rm=TRUE)
table800.crq[4,2]<-sd(b0.crq.37,na.rm=TRUE)
table800.crq[4,3]<-mean(b1.crq.37,na.rm=TRUE)
table800.crq[4,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table800.is[4,1]<-mean(b0.is.37,na.rm=TRUE)
table800.is[4,2]<-sd(b0.is.37,na.rm=TRUE)
table800.is[4,3]<-mean(b1.is.37,na.rm=TRUE)
table800.is[4,4]<-sd(b1.is.37,na.rm=TRUE)

table.event[4,3] = mean(avg.event.3.800)

#### t_0=3 & c=70% & data = 1600####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
avg.event.3.1600<-c()
set.seed(31)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(1600,c.7,3)
  avg.event.3.1600[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.10)
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
for (i in 1:200){
  a<-data.gen(1600,c.7,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.10, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.37[i] = NA
    b0.is.sd.37[i] = NA
    b1.is.37[i] = NA
    b1.is.sd.37[i] = NA
  })
}

# Crq beta table
table1600.crq[4,1]<-mean(b0.crq.37,na.rm=TRUE)
table1600.crq[4,2]<-sd(b0.crq.37,na.rm=TRUE)
table1600.crq[4,3]<-mean(b1.crq.37,na.rm=TRUE)
table1600.crq[4,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table1600.is[4,1]<-mean(b0.is.37,na.rm=TRUE)
table1600.is[4,2]<-sd(b0.is.37,na.rm=TRUE)
table1600.is[4,3]<-mean(b1.is.37,na.rm=TRUE)
table1600.is[4,4]<-sd(b1.is.37,na.rm=TRUE)

table.event[4,4] = mean(avg.event.3.1600)

#### t_0=3 & c=70% & data = 3200####
b0.crq.37<-c()
b0.crq.sd.37<-c()
b1.crq.37<-c()
b1.crq.sd.37<-c()
b0.is.37<-c()
b0.is.sd.37<-c()
b1.is.37<-c()
b1.is.sd.37<-c()
avg.event.3.3200<-c()
set.seed(31)
for (i in 1:200){
  # Beta estimation method 1 : Crq package
  a<-data.gen(3200,c.7,3)
  avg.event.3.3200[i]=sum(a[,7]==1)
  tryCatch({
    crq.fit = crq.est(a[,3], 1, a[,6], a[,7], 3, 0.10)
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
for (i in 1:200){
  a<-data.gen(3200,c.0,3)
  tryCatch({
    ismb.fit = ismb.est(a[,3], 1, a[,6], a[,7], 3, 0.10, 100)
    b0.is.37[i] = ismb.fit[1,1]
    b0.is.sd.37[i] = ismb.fit[1,2]
    b1.is.37[i] = ismb.fit[2,1]
    b1.is.sd.37[i] = ismb.fit[2,2]
  }
  , error=function(e){
    b0.is.37[i] = NA
    b0.is.sd.37[i] = NA
    b1.is.37[i] = NA
    b1.is.sd.37[i] = NA
  })
}

# Crq beta table
table3200.crq[4,1]<-mean(b0.crq.37,na.rm=TRUE)
table3200.crq[4,2]<-sd(b0.crq.37,na.rm=TRUE)
table3200.crq[4,3]<-mean(b1.crq.37,na.rm=TRUE)
table3200.crq[4,4]<-sd(b1.crq.37,na.rm=TRUE)

# IS beta table
table3200.is[4,1]<-mean(b0.is.37,na.rm=TRUE)
table3200.is[4,2]<-sd(b0.is.37,na.rm=TRUE)
table3200.is[4,3]<-mean(b1.is.37,na.rm=TRUE)
table3200.is[4,4]<-sd(b1.is.37,na.rm=TRUE)

table.event[4,5] = mean(avg.event.3.3200)