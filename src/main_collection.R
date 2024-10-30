library(quantreg) # Quantile Regression
library(lpSolve) # Linear and Integer Programming
library(Hmisc)
library(parallel) # for detectCores() 
library(snowfall) # parallel programming 
library(purrr) # for cross2()
library(LowRankQP) # for quadratic programming
library(quadprog)
source('functions.R')

# read data
data_sampled<-read.csv('data_sampled.csv',header = TRUE, fileEncoding = "GBK")
data_sampled<-data_sampled[,-1]
sales<-data_sampled[,1]
covariate<-data_sampled[,-1]

# intialize
n.all<-dim(data_sampled)[1]
R.r<-500 
M<-ncol(covariate)+1 
settings<-list(n=c(50,100,150),tau=c(0.95,0.9,0.5,0.1,0.05))
settings<-cross2(settings$n,settings$tau,.filter = NULL)
settings <- lapply(settings, unlist, 1)

#################### main p=1 ####################
sfInit(parallel = TRUE, cpus = detectCores()-10) 
sfLibrary(LowRankQP)
sfLibrary(quantreg)
sfLibrary(lpSolve)
sfExport('Indic','turnbits_cir','FPE_qr','JCVMA_qr')
sfExport('M','n.all','R.r','settings','sales','covariate')

main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  n.s<-n.all-n 
  FPEr.JCVMA5<-numeric(R.r)
  FPEr.JCVMA10<-numeric(R.r)
  FPEr.JMA<-numeric(R.r)
  FPEr.SAIC<-numeric(R.r)
  FPEr.SBIC<-numeric(R.r)
  FPEr.AIC<-numeric(R.r)
  FPEr.BIC<-numeric(R.r)
  FPEr.MAP5<-numeric(R.r)
  FPEr.MAP10<-numeric(R.r)
  FPEr.EWA<-numeric(R.r)
  FPEr.LM<-numeric(R.r)
  
  #######################################################
  for (r in 1:R.r){
  ################ divide training and testing set ####################
  train<-sample(n.all,n)
  X<-covariate[train,]
  X<-as.matrix(cbind(1,X))
  X.s<-covariate[-train,]
  X.s<-as.matrix(cbind(1,X.s))
  y<-sales[train]
  y.s<-sales[-train]
  data<-data.frame(y,X)
  # theta.hat
  theta.hat<-matrix(0,M,M)
  for (m in 1:M) {
    fit<-rq(data=data.frame(X[,1:m]), y~0+.,tau=tau)
    theta.hat[1:m,m]<-coef(fit)
  }
  
  ####################  JCVMA ######################
  JCVMA5<-JCVMA_qr(X,y,J=5,tau=tau)
  w.JCVMA5<-JCVMA5$JCVMA
  w.MAP5<-JCVMA5$MAP
  JCVMA10<-JCVMA_qr(X,y,J=10,tau=tau)
  w.JCVMA10<-JCVMA10$JCVMA
  w.MAP10<-JCVMA10$MAP
  w.JMA<-JCVMA_qr(X,y,J=n,tau=tau)$JCVMA
  
  ######################### (S)AIC & (S)BIC    ###################
  AIC<-rep(0,len=M)
  BIC<-rep(0,len=M)
  for (m in 1:M) {
    AIC[m]<-2*n*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+2*m
    BIC[m]<-2*n*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+m*log(n)
  }
  # AIC & BIC weights
  AIC<-AIC-min(AIC)
  BIC<-BIC-min(BIC)
  w.AIC<-rep(0,len=M)
  w.BIC<-rep(0,len=M)
  w.AIC[which.min(AIC)]<-1
  w.BIC[which.min(BIC)]<-1
  # SAIC & SBIC weights
  w.SAIC<-rep(0,len=M)
  w.SBIC<-rep(0,len=M)
  for (m in 1:M) {
    w.SAIC[m]<-exp(-0.5*AIC[m])/sum(exp(-0.5*AIC))
    w.SBIC[m]<-exp(-0.5*BIC[m])/sum(exp(-0.5*BIC))
  }
  
  ###################### EWA & Largest Model  ###################
  w.EWA<-rep(1/M,M)
  w.LM<-c(rep(0,M-1),1)
  
  #################### FPE ############################
  muhat.s<-X.s %*% theta.hat
  FPEr.JCVMA5[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JCVMA5)
  FPEr.JCVMA10[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JCVMA10)
  FPEr.JMA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JMA)
  FPEr.SAIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SAIC)
  FPEr.SBIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SBIC)
  FPEr.AIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.AIC)
  FPEr.BIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.BIC)
  FPEr.MAP5[r]<-FPE_qr(tau,y.s,muhat.s%*%w.MAP5)
  FPEr.MAP10[r]<-FPE_qr(tau,y.s,muhat.s%*%w.MAP10)
  FPEr.EWA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.EWA)
  FPEr.LM[r]<-FPE_qr(tau,y.s,muhat.s%*%w.LM)
  write(r, file = "test.txt",append = TRUE, sep =" ") 
}
FPE.JCVMA5<-mean(FPEr.JCVMA5)
FPE.JCVMA10<-mean(FPEr.JCVMA10)
FPE.JMA<-mean(FPEr.JMA)
FPE.SAIC<-mean(FPEr.SAIC)
FPE.SBIC<-mean(FPEr.SBIC)
FPE.AIC<-mean(FPEr.AIC)
FPE.BIC<-mean(FPEr.BIC)
FPE.MAP5<-mean(FPEr.MAP5)
FPE.MAP10<-mean(FPEr.MAP10)
FPE.EWA<-mean(FPEr.EWA)
FPE.LM<-mean(FPEr.LM)
return(c(n,tau,FPE.JCVMA5,FPE.JCVMA10,FPE.JMA,
         FPE.SAIC,FPE.SBIC,FPE.MAP5,FPE.MAP10,
         FPE.EWA,FPE.AIC,FPE.BIC,FPE.LM))
}

#R.r=5
#FPE_parallel<-lapply(1:length(settings), main_parallel) 
FPE_parallel<-sfLapply(1:length(settings), main_parallel) 
FPE<-matrix(unlist(FPE_parallel),nrow=length(settings),byrow = TRUE)
FPE<-data.frame(FPE)
colnames(FPE)<-c('n','tau','JCVMA5','JCVMA10','JMA','SAIC','SBIC',
                  'MAP5','MAP10','EWA','AIC','BIC','LM')
unlink("test.txt")
sfStop()

FPE_result<-data.frame(p=1,FPE)

#################### main p=2 ####################
sfInit(parallel = TRUE, cpus = detectCores()-10) 
sfLibrary(LowRankQP)
sfExport('Indic','turnbits_cir','er','FPE_er','JCVMA_er')
sfExport('M','n.all','R.r','settings','sales','covariate')

main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  n.s<-n.all-n
  FPEr.JCVMA5<-numeric(R.r)
  FPEr.JCVMA10<-numeric(R.r)
  FPEr.JMA<-numeric(R.r)
  FPEr.SAIC<-numeric(R.r)
  FPEr.SBIC<-numeric(R.r)
  FPEr.AIC<-numeric(R.r)
  FPEr.BIC<-numeric(R.r)
  FPEr.MAP5<-numeric(R.r)
  FPEr.MAP10<-numeric(R.r)
  FPEr.EWA<-numeric(R.r)
  FPEr.LM<-numeric(R.r)
  
  #######################################################
  for (r in 1:R.r){
    ################ divide training and testing set ####################
    train<-sample(n.all,n)
    X<-covariate[train,]
    X<-as.matrix(cbind(1,X))
    X.s<-covariate[-train,]
    X.s<-as.matrix(cbind(1,X.s))
    y<-sales[train]
    y.s<-sales[-train]
    data<-data.frame(y,X)
    theta.hat<-matrix(0,M,M)
    for (m in 1:M) {
      theta.hat[1:m,m]<-er(X=as.matrix(X[,1:m]),y=data[,1],
                           tau=tau,max.iter = 100,tol = 1)
    }
    
    ####################  JCVMA ######################
    JCVMA5<-JCVMA_er(X,y,J=5,tau=tau)
    w.JCVMA5<-JCVMA5$JCVMA
    w.MAP5<-JCVMA5$MAP
    JCVMA10<-JCVMA_er(X,y,J=10,tau=tau)
    w.JCVMA10<-JCVMA10$JCVMA
    w.MAP10<-JCVMA10$MAP
    w.JMA<-JCVMA_er(X,y,J=n,tau=tau)$JCVMA
    
    ######################### (S)AIC & (S)BIC    ###################
    AIC<-rep(0,len=M)
    BIC<-rep(0,len=M)
    for (m in 1:M) {
      AIC[m]<-n*log(FPE_er(tau,y,X%*%theta.hat[,m]))+2*m
      BIC[m]<-n*log(FPE_er(tau,y,X%*%theta.hat[,m]))+m*log(n)
    }
    # AIC & BIC weights
    AIC<-AIC-min(AIC)
    BIC<-BIC-min(BIC)
    w.AIC<-rep(0,len=M)
    w.BIC<-rep(0,len=M)
    w.AIC[which.min(AIC)]<-1
    w.BIC[which.min(BIC)]<-1
    # SAIC & SBIC weights
    w.SAIC<-rep(0,len=M)
    w.SBIC<-rep(0,len=M)
    for (m in 1:M) {
      w.SAIC[m]<-exp(-0.5*AIC[m])/sum(exp(-0.5*AIC))
      w.SBIC[m]<-exp(-0.5*BIC[m])/sum(exp(-0.5*BIC))
    }
    
    ###################### EWA & Largest Model  ###################
    w.EWA<-rep(1/M,M)
    w.LM<-c(rep(0,M-1),1)
    
    #################### FPE ############################
    muhat.s<-matrix(0,n.s,M)
    for (m in 1:M) {
      muhat.s[,m]<-X.s %*% theta.hat[,m]
    }
    FPEr.JCVMA5[r]<-FPE_er(tau,y.s,muhat.s%*%w.JCVMA5)
    FPEr.JCVMA10[r]<-FPE_er(tau,y.s,muhat.s%*%w.JCVMA10)
    FPEr.JMA[r]<-FPE_er(tau,y.s,muhat.s%*%w.JMA)
    FPEr.SAIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SAIC)
    FPEr.SBIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SBIC)
    FPEr.AIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.AIC)
    FPEr.BIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.BIC)
    FPEr.MAP5[r]<-FPE_er(tau,y.s,muhat.s%*%w.MAP5)
    FPEr.MAP10[r]<-FPE_er(tau,y.s,muhat.s%*%w.MAP10)
    FPEr.EWA[r]<-FPE_er(tau,y.s,muhat.s%*%w.EWA)
    FPEr.LM[r]<-FPE_er(tau,y.s,muhat.s%*%w.LM)
    write(r, file = "test.txt",append = TRUE, sep =" ") 
  }
  FPE.JCVMA5<-mean(FPEr.JCVMA5, na.rm = TRUE) #mean
  FPE.JCVMA10<-mean(FPEr.JCVMA10, na.rm = TRUE)
  FPE.JMA<-mean(FPEr.JMA, na.rm = TRUE)
  FPE.SAIC<-mean(FPEr.SAIC)
  FPE.SBIC<-mean(FPEr.SBIC)
  FPE.AIC<-mean(FPEr.AIC)
  FPE.BIC<-mean(FPEr.BIC)
  FPE.MAP5<-mean(FPEr.MAP5)
  FPE.MAP10<-mean(FPEr.MAP10)
  FPE.EWA<-mean(FPEr.EWA)
  FPE.LM<-mean(FPEr.LM)
  return(c(n,tau,FPE.JCVMA5,FPE.JCVMA10,FPE.JMA,
           FPE.SAIC,FPE.SBIC,FPE.MAP5,FPE.MAP10,
           FPE.EWA,FPE.AIC,FPE.BIC,FPE.LM))
}

#R.r=5
#FPE_parallel<-lapply(1:length(settings), main_parallel) 
FPE_parallel<-sfLapply(1:length(settings), main_parallel) 
FPE<-matrix(unlist(FPE_parallel),nrow=length(settings),byrow = TRUE)
FPE<-data.frame(FPE)
colnames(FPE)<-c('n','tau','JCVMA5','JCVMA10','JMA','SAIC','SBIC',
                 'MAP5','MAP10','EWA','AIC','BIC','LM')
unlink("test.txt")
sfStop()

FPE_result<-rbind(FPE_result,data.frame(p=2,FPE))
write.csv(FPE_result, 'FPE_result.csv',fileEncoding = "GBK")

FPE_result<-read.csv('FPE_result.csv')
FPE_result<-FPE_result[,-1]
FPE_normalized<-cbind(FPE_result[,c(1,2,3)],FPE_result[,-c(1,2,3,12,13,14)]/FPE_result$LM)
FPE_final<-FPE_normalized[,-which(colnames(FPE_normalized)%in%c('JCVMA10','MAP10') )]
