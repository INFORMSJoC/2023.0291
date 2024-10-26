library(quantreg) # Quantile Regression
library(lpSolve) # Linear and Integer Programming
setwd('/Users/gudieqi/Desktop/MA under AL/simulation/design2')
library(parallel) # for detectCores() 
library(snowfall) # parallel programming 
library(purrr) # for cross2()
source('functions.R')

################################################
# initialization
k<-5
R.r<-500 
k.unsure<-4 # the first 1 variables are fixed, and the later 4 variables are unsure
# Z is the selection matrix
Z<-matrix(unlist(lapply(0:(2^k.unsure-1),turnbits_cir,k.unsure)),
          ncol = k.unsure ,byrow = T)
Z<-cbind(matrix(TRUE,2^k.unsure,1),Z) # 3 fixed variable, 8 & 16 are correct models
M<-dim(Z)[1] 
beta<-c(1,1,1,1,0)
settings<-list(n=c(100,200,400,800,1600),tau=c(0.95,0.9,0.5,0.1,0.05))
settings<-cross2(settings$n,settings$tau,.filter = NULL)
settings <- lapply(settings, unlist, 1)

#################### main p=1 ####################
sfInit(parallel = TRUE, cpus = detectCores()-30) 
sfLibrary(quantreg)
sfLibrary(lpSolve)
sfExport('Indic','turnbits_cir','FPE_qr','JCVMA_qr')
sfExport('settings','R.r','k.unsure','Z','M','k','beta')

main_parallel <-function(ii){
  ###################### initialization ###########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  theta<-beta+c(1,0,0,0,0)*qnorm(tau)
  MSEr.JCVMA5<-numeric(R.r)
  MSEr.JCVMA10<-numeric(R.r)
  MSEr.true<-numeric(R.r)
  wr.correct.JCVMA5<-numeric(R.r)
  wr.correct.JCVMA10<-numeric(R.r)
  ######################  main ###########################
  for (r in 1:R.r){
  X<-matrix(0,n,k) 
  X[,1] <- 1 
  X[,-1]<-rnorm(n*(k-1),0,1)
  eplison<-X%*%c(1,0,0,0,0)
  eplison<-eplison*rnorm(n,0,1)
  y<-X %*% beta+eplison
  data<-data.frame(y,X)
  theta.hat<-matrix(0,k,M)
  for (m in 1:M) {
    fit<-rq(data=data.frame(X[,Z[m,]]), y~0+.,tau=tau, method = 'fn')
    theta.hat[Z[m,],m]<-coef(fit)
  }
  
  ####################  JCVMA ######################
  w.JCVMA5<-JCVMA_qr(X,y,J=5,tau=tau,nest = FALSE,fixed = 1)$JCVMA
  w.JCVMA10<-JCVMA_qr(X,y,J=10,tau=tau,nest = FALSE,fixed = 1)$JCVMA
  # w.JMA<-JCVMA_qr(X,y,J=n,tau=tau,nest = FALSE,fixed = 1)$JCVMA

  ####################  MSE and weights ######################
  w.true<-rep(0,M)
  w.true[8]<-1
  theta.hat.JCVMA5<-theta.hat %*% w.JCVMA5
  theta.hat.JCVMA10<-theta.hat %*% w.JCVMA10
  theta.hat.true<-theta.hat %*% w.true
  MSEr.JCVMA5[r]<- t(theta.hat.JCVMA5-theta) %*% (theta.hat.JCVMA5-theta)
  MSEr.JCVMA10[r]<- t(theta.hat.JCVMA10-theta) %*% (theta.hat.JCVMA10-theta)
  MSEr.true[r]<-t(theta.hat.true-theta) %*% (theta.hat.true-theta)
  wr.correct.JCVMA5<-w.JCVMA5[8]+w.JCVMA5[16]
  wr.correct.JCVMA10<-w.JCVMA10[8]+w.JCVMA10[16]
  write(r, file = "test.txt",append = TRUE, sep =" ") 
  }
  MSE.JCVMA5<-mean(MSEr.JCVMA5)
  MSE.JCVMA10<-mean(MSEr.JCVMA10)
  MSE.true<-mean(MSEr.true)
  w.correct.JCVMA5<-mean(wr.correct.JCVMA5)
  w.correct.JCVMA10<-mean(wr.correct.JCVMA10)
  return(c(n,tau,MSE.JCVMA5, MSE.JCVMA10, MSE.true, 
           w.correct.JCVMA5, w.correct.JCVMA10))
}
# R.r<-5
# MSE_parallel<-lapply(1:length(settings), main_parallel) 
MSE_parallel<-sfLapply(1:length(settings), main_parallel) 
MSE<-matrix(unlist(MSE_parallel),nrow=length(settings),byrow = TRUE)
MSE<-data.frame(MSE)
colnames(MSE)<-c('n','tau','MSE.JCVMA5','MSE.JCVMA10','MSE.true',
                 'w.correct.JCVMA5','w.correct.JCVMA10')
unlink("test.txt")
sfStop()
MSE_result<-data.frame(p=1,MSE)
View(FPE_result)
write.csv(MSE_result, 'MSE_result.csv',fileEncoding = "GBK")

#################### main p=2 ####################
sfInit(parallel = TRUE, cpus = detectCores()-30) 
sfLibrary(LowRankQP)
sfLibrary(expectreg)
sfExport('Indic','turnbits_cir','er','FPE_er','JCVMA_er')
sfExport('settings','R.r','k.unsure','Z','M','k','beta')

main_parallel <-function(ii){
  ###################### initialization ###########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  theta<-beta+c(1,0,0,0,0)*enorm(tau)
  MSEr.JCVMA5<-numeric(R.r)
  MSEr.JCVMA10<-numeric(R.r)
  MSEr.true<-numeric(R.r)
  wr.correct.JCVMA5<-numeric(R.r)
  wr.correct.JCVMA10<-numeric(R.r)
  ######################  main ###########################
  for (r in 1:R.r){
    X<-matrix(0,n,k) 
    X[,1] <- 1 
    X[,-1]<-rnorm(n*(k-1),0,1)
    eplison<-X%*%c(1,0,0,0,0)
    eplison<-eplison*rnorm(n,0,1)
    y<-X %*% beta+eplison
    data<-data.frame(y,X)
    theta.hat<-matrix(0,k,M)
    for (m in 1:M) {
      theta.hat[Z[m,],m]<-er(X=as.matrix(X[,Z[m,]]),y=data[,1],
                             tau=tau,max.iter = 100,tol = 1)
    }  
    
    ####################  JCVMA ######################
    w.JCVMA5<-JCVMA_er(X,y,J=5,tau=tau,nest = FALSE,fixed = 1)$JCVMA
    w.JCVMA10<-JCVMA_er(X,y,J=10,tau=tau,nest = FALSE,fixed = 1)$JCVMA
    # w.JMA<-JCVMA_er(X,y,J=n,tau=tau,nest = FALSE,fixed = 1)$JCVMA
    
    ####################  MSE and weights ######################
    w.true<-rep(0,M)
    w.true[8]<-1
    theta.hat.JCVMA5<-theta.hat %*% w.JCVMA5
    theta.hat.JCVMA10<-theta.hat %*% w.JCVMA10
    theta.hat.true<-theta.hat %*% w.true
    MSEr.JCVMA5[r]<- t(theta.hat.JCVMA5-theta) %*% (theta.hat.JCVMA5-theta)
    MSEr.JCVMA10[r]<- t(theta.hat.JCVMA10-theta) %*% (theta.hat.JCVMA10-theta)
    MSEr.true[r]<-t(theta.hat.true-theta) %*% (theta.hat.true-theta)
    wr.correct.JCVMA5<-w.JCVMA5[8]+w.JCVMA5[16]
    wr.correct.JCVMA10<-w.JCVMA10[8]+w.JCVMA10[16]
    write(r, file = "test.txt",append = TRUE, sep =" ") # 测试是否正常运行
  }
  MSE.JCVMA5<-mean(MSEr.JCVMA5)
  MSE.JCVMA10<-mean(MSEr.JCVMA10)
  MSE.true<-mean(MSEr.true)
  w.correct.JCVMA5<-mean(wr.correct.JCVMA5)
  w.correct.JCVMA10<-mean(wr.correct.JCVMA10)
  return(c(n,tau,MSE.JCVMA5, MSE.JCVMA10, MSE.true, 
           w.correct.JCVMA5, w.correct.JCVMA10))
}

# R.r<-5
# MSE_parallel<-lapply(1:length(settings), main_parallel) 
MSE_parallel<-sfLapply(1:length(settings), main_parallel) 
MSE<-matrix(unlist(MSE_parallel),nrow=length(settings),byrow = TRUE)
MSE<-data.frame(MSE)
colnames(MSE)<-c('n','tau','MSE.JCVMA5','MSE.JCVMA10','MSE.true',
                 'w.correct.JCVMA5','w.correct.JCVMA10')
unlink("test.txt")
sfStop()
MSE_result<-data.frame(p=2,MSE)
View(MSE_result)
write.csv(MSE_result, 'MSE_result.csv',fileEncoding = "GBK")




