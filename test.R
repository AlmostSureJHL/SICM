rm(list=ls(all=TRUE))
setwd('E:/ECNU/Reaserch/SSDR/Code_SSDR')
source('functions_set.R')
source('functions_set_forTarma2023.R')

library(MASS)
library(ggplot2)
library(ggpubr)

gap<-c(4,20,20,0.8,0.7,0.6)
Ra<-c(10,10,1,1,1)

Mpower=matrix(0,5,2)

for (k in 1:2){
  
  delta=gap[k]
  n<-60
  p<-60#floor(exp(n^{0.4}))+230
  
  rhop<-8
  rep_time<-300
  
  
  eta<-0.1
  alpha=0.05
  
  
  q=gap[k]#floor(0.5*p^{delta})
  non0<-sqrt(eta)
  beta<-c(rep(non0,q),rep(0,(p-q)))
  
  
  set.seed(1000)
  rho1<-runif(rhop,0,1)
  #rho<-c(rho1,rep(0,p-rhop))
  #Gamma<-cbind(diag(rep(rho[1],p)),matrix(0,p,m-p))
  #for (i in 2:(p-1))
  #{Gamma<-Gamma+cbind(matrix(0,p,i-1),diag(rep(rho[i],p)),matrix(0,p,m-p-i+1))}
  #Gamma<-Gamma+cbind(matrix(0,p,m-p),diag(rep(rho[p],p)))
  
  SCMIstat=seq(0,0,length=rep_time) 
  siindSCMI=seq(0,0,length=rep_time) 
  siindKLR2023=seq(0,0,length=rep_time) 
  
  for (j in 1:rep_time){
    print(j)
    z=matrix(rnorm((p+rhop-1)*n,0,1),n,(p+rhop-1))
    X=matrix(NA,n,p)#matrix(runif(n*p,0,0.5),n,p)#
    for (i in 1:p){
      X[,i]=z[,i:(i+rhop-1)]%*%rho1
    }
    Hazard<-X^4%*%beta
    epsion<-runif(n,0,1)
    T<-exp(-log(1-epsion)/Hazard)
    C<-rexp(n,1/Ra[k])
    D<-(T<=C)*1  #observe_rat:sum(D)/n
    Y<-pmin(T,C)
    sum(D)/n
    
    SCMIstat[j]=T_n(X,Y,D,'Gauss')$Z_alpha
    siindSCMI[j]=(SCMIstat[j]>qnorm(1-alpha,0,1))
    
    siindKLR2023[j]=(wild_boot_logrank(X, Y, D, 'Gauss', 'Gauss', num_bootstrap=999)$pvalue<0.05)
    # print(c(siindCCov[j],siindChen2010[j],siindMDD2018[j]))
    print(c(mean(siindSCMI[1:j]),mean(siindKLR2023[1:j])))
  }
  Mpower[k,]=c(mean(siindSCMI),mean(siindKLR2023))
}

