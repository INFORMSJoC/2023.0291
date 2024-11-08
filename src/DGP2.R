#####################  library the packages  ################
library(quantreg) # rq, Quantile Regression
library(lpSolve) # lp, Linear and Integer Programming
library(LowRankQP) # for quadratic programming
library(expectreg)
library(parallel) # for detectCores() 
library(snowfall) # parallel programming 
library(purrr) # for cross3()
library(quadprog)
setwd('/Users/gudieqi/Desktop/MA under AL/simulation')
source('functions.R')

################################################
# initialization
k.all<-30
R.r<-500 # replications
n.s<-100
k.unsure<-5 # the first 3 variables are fixed, and the later 5 variables are unsure
# Z is the selection matrix
Z<-matrix(unlist(lapply(0:(2^k.unsure-1),turnbits_cir,k.unsure)),ncol = k.unsure ,byrow = T)
Z<-cbind(matrix(TRUE,2^k.unsure,3),Z) # 3 fixed variable
M<-dim(Z)[1] # number of models 
k<-3+k.unsure # covarible dimension

# simulation settings
settings<-list(R.square=seq(0.1,0.9,0.1),n=c(50,100,200,400),
               tau=c(0.95,0.9,0.5,0.1,0.05))
settings<-cross3(settings$R.square,settings$n, settings$tau)
settings<-lapply(settings, unlist,1)

#################### DGP1 for p=1 ####################
## initialize snowfall
sfInit(parallel = TRUE, cpus = detectCores()-10) 
# library packages
sfLibrary(LowRankQP)
sfLibrary(quantreg)
sfLibrary(lpSolve)
# library functions
sfExport('Indic','turnbits_cir','FPE_qr','JCVMA_qr')
# library variables
sfExport('settings','k.all','R.r','n.s','k.unsure','Z','M','k')

###################### main p=1 ##########################
main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  R.square<-setting[1]
  n<-setting[2]
  tau<-setting[3]
  q<-qnorm(tau) # tau quantile
  FPE<-1/sqrt(2*pi)*exp(-q^2/2)*7 # the item to be subtracted
  # construct the coefficients
  theta<- sqrt((63/(1-R.square)-63)/38)
  beta<-theta*c(1,1,1,0,0,1,2,3,rep(1,k.all-k))
  EFPEr.JCVMA5<-numeric(R.r)
  EFPEr.JCVMA10<-numeric(R.r)
  EFPEr.JMA<-numeric(R.r)
  EFPEr.SAIC<-numeric(R.r)
  EFPEr.SBIC<-numeric(R.r)
  EFPEr.AIC<-numeric(R.r)
  EFPEr.BIC<-numeric(R.r)
  EFPEr.MAP5<-numeric(R.r)
  EFPEr.MAP10<-numeric(R.r)
  EFPEr.EWA<-numeric(R.r)
  EFPE0<-matrix(0,R.r,M)
  
  #######################################################
  for (r in 1:R.r) {
    ###################### DGP ############################
    # construct samples yi,xi and out-of-sample data ys,xs
    X<-matrix(0,n,k.all) 
    X.s<-matrix(0,n.s,k.all) 
    X[,1] <- 1 # the first column of the design matrix is 1
    X.s[,1] <- 1 
    X[,-1]<-rnorm(n*(k.all-1),0,1) # the 2 to k columns are N(0,1)
    X.s[,-1]<-rnorm(n.s*(k.all-1),0,1)
    eplison<-matrix(0,n,1)
    eplison.s<-matrix(0,n.s,1)
    for (j in 2:8) {
      eplison<-eplison+X[,j]^2
      eplison.s<-eplison.s+X.s[,j]^2
    }
    eplison<-eplison*rnorm(n,0,1)
    eplison.s<-eplison.s*rnorm(n.s,0,1)
    y<-X %*% beta+eplison
    y.s<-X.s %*% beta+eplison.s
    X<-X[,1:k]
    X.s<-X.s[,1:k]
    data<-data.frame(y,X)
    # coefficient estimate theta.hat
    theta.hat<-matrix(0,k,M)
    for (m in 1:M) {
      fit<-rq(data=data.frame(X[,Z[m,]]), y~0+.,tau=tau)
      theta.hat[Z[m,],m]<-coef(fit)
    }
    
    ####################  JCVMA ######################
    JCVMA5<-JCVMA_qr(X,y,J=5,tau=tau,nest = FALSE,fixed = 3)
    w.JCVMA5<-JCVMA5$JCVMA
    w.MAP5<-JCVMA5$MAP
    JCVMA10<-JCVMA_qr(X,y,J=10,tau=tau,nest = FALSE,fixed = 3)
    w.JCVMA10<-JCVMA10$JCVMA
    w.MAP10<-JCVMA10$MAP
    w.JMA<-JCVMA_qr(X,y,J=n,tau=tau,nest = FALSE,fixed = 3)$JCVMA

    ######################### (S)AIC & (S)BIC    ###################
    AIC<-rep(0,len=M)
    BIC<-rep(0,len=M)
    for (m in 1:M) {
      AIC[m]<-2*n*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+2*sum(Z[m,])
      BIC[m]<-2*n*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+sum(Z[m,])*log(n)
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
    
    ###################### Equal weights averaging  ###################
    w.EWA<-rep(1/M,M)
    
    #################### EFPE ############################
    muhat.s<-matrix(0,n.s,M)
    for (m in 1:M) {
      muhat.s[,m]<-X.s %*% theta.hat[,m]
      EFPE0[r,m]<-FPE_qr(tau,y.s,muhat.s[,m])-FPE
    }
    EFPEr.JCVMA5[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JCVMA5)-FPE
    EFPEr.JCVMA10[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JCVMA10)-FPE
    EFPEr.JMA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JMA)-FPE
    EFPEr.SAIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SAIC)-FPE
    EFPEr.SBIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SBIC)-FPE
    EFPEr.AIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.AIC)-FPE
    EFPEr.BIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.BIC)-FPE
    EFPEr.MAP5[r]<-FPE_qr(tau,y.s,muhat.s%*%w.MAP5)-FPE
    EFPEr.MAP10[r]<-FPE_qr(tau,y.s,muhat.s%*%w.MAP10)-FPE
    EFPEr.EWA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.EWA)-FPE
    # test to see if it is working properly
    write(r, file = "test.txt",append = TRUE, sep =" ")
  }
  EFPE<-min(colMeans(EFPE0))
  EFPE.JCVMA5<-mean(EFPEr.JCVMA5)/EFPE
  EFPE.JCVMA10<-mean(EFPEr.JCVMA10)/EFPE
  EFPE.JMA<-mean(EFPEr.JMA)/EFPE
  EFPE.SAIC<-mean(EFPEr.SAIC)/EFPE
  EFPE.SBIC<-mean(EFPEr.SBIC)/EFPE
  EFPE.AIC<-mean(EFPEr.AIC)/EFPE
  EFPE.BIC<-mean(EFPEr.BIC)/EFPE
  EFPE.MAP5<-mean(EFPEr.MAP5)/EFPE
  EFPE.MAP10<-mean(EFPEr.MAP10)/EFPE
  EFPE.EWA<-mean(EFPEr.EWA)/EFPE
  return(c(n,M,tau,R.square,EFPE.JCVMA5,EFPE.JCVMA10,EFPE.JMA,
           EFPE.SAIC,EFPE.SBIC,EFPE.MAP5,EFPE.MAP10,
           EFPE.EWA,EFPE.AIC,EFPE.BIC))
}
# R.r<-5
# EFPE_parallel<-lapply(1:length(settings), main_parallel) 
EFPE_parallel<-sfLapply(1:length(settings), main_parallel) 
EFPE<-matrix(unlist(EFPE_parallel),nrow=length(settings),byrow = TRUE)
EFPE<-data.frame(EFPE)
colnames(EFPE)<-c('n','M','tau','R2','JCVMA5','JCVMA10','JMA','SAIC','SBIC',
                  'MAP5','MAP10','EWA','AIC','BIC')
unlink("test.txt")
sfStop()

#################### DGP1 for p=2 ####################
## initialize snowfall
sfInit(parallel = TRUE, cpus = detectCores()-20) 
sfLibrary(LowRankQP)
sfLibrary(expectreg)
sfExport('Indic','turnbits_cir','FPE_er','JCVMA_er','er')
sfExport('settings','k.all','R.r','n.s','k.unsure','Z','M','k')

###################### main p=2 ##########################
main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  R.square<-setting[1]
  n<-setting[2]
  tau<-setting[3]
  e<-enorm(tau) # tau expectile
  FPE<-((1-2*tau)*e*1/sqrt(2*pi)*exp(-e^2/2)+(1-2*tau)*pnorm(e)+
          tau+(1-2*tau)*e^2*pnorm(e)+tau*e^2)*63 
  # construct the coefficients
  theta<- sqrt((63/(1-R.square)-63)/38)
  beta<-theta*c(1,1,1,0,0,1,2,3,rep(1,k.all-k))
  EFPEr.JCVMA5<-numeric(R.r)
  EFPEr.JCVMA10<-numeric(R.r)
  EFPEr.JMA<-numeric(R.r)
  EFPEr.SAIC<-numeric(R.r)
  EFPEr.SBIC<-numeric(R.r)
  EFPEr.AIC<-numeric(R.r)
  EFPEr.BIC<-numeric(R.r)
  EFPEr.MAP5<-numeric(R.r)
  EFPEr.MAP10<-numeric(R.r)
  EFPEr.EWA<-numeric(R.r)
  EFPE0<-matrix(0,R.r,M)
  
  #######################################################
  for (r in 1:R.r) {
    ###################### DGP ############################
    # construct samples yi,xi and out-of-sample data ys,xs
    X<-matrix(0,n,k.all) 
    X.s<-matrix(0,n.s,k.all) 
    X[,1] <- 1 # the first column of the design matrix is 1
    X.s[,1] <- 1 
    X[,-1]<-rnorm(n*(k.all-1),0,1) # the 2 to k columns are N(0,1)
    X.s[,-1]<-rnorm(n.s*(k.all-1),0,1)
    eplison<-matrix(0,n,1)
    eplison.s<-matrix(0,n.s,1)
    for (j in 2:8) {
      eplison<-eplison+X[,j]^2
      eplison.s<-eplison.s+X.s[,j]^2
    }
    eplison<-eplison*rnorm(n,0,1)
    eplison.s<-eplison.s*rnorm(n.s,0,1)
    y<-X %*% beta+eplison
    y.s<-X.s %*% beta+eplison.s
    X<-X[,1:k]
    X.s<-X.s[,1:k]
    data<-data.frame(y,X)
    # coefficient estimate theta.hat
    theta.hat<-matrix(0,k,M)
    for (m in 1:M) {
      theta.hat[Z[m,],m]<-er(X=as.matrix(X[,Z[m,]]),y=data[,1],
                           tau=tau,max.iter = 100,tol = 1)
    }  
    
    ####################  JCVMA ######################
    JCVMA5<-JCVMA_er(X,y,J=5,tau=tau,nest = FALSE,fixed = 3)
    w.JCVMA5<-JCVMA5$JCVMA
    w.MAP5<-JCVMA5$MAP
    JCVMA10<-JCVMA_er(X,y,J=10,tau=tau,nest = FALSE,fixed = 3)
    w.JCVMA10<-JCVMA10$JCVMA
    w.MAP10<-JCVMA10$MAP
    w.JMA<-JCVMA_er(X,y,J=n,tau=tau,nest = FALSE,fixed = 3)$JCVMA
    
    ######################### (S)AIC & (S)BIC ###################
    AIC<-rep(0,len=M)
    BIC<-rep(0,len=M)
    for (m in 1:M) {
      AIC[m]<-n*log(FPE_er(tau,y,X%*%theta.hat[,m]))+2*sum(Z[m,])
      BIC[m]<-n*log(FPE_er(tau,y,X%*%theta.hat[,m]))+sum(Z[m,])*log(n)
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
    
    ###################### Equal weights averaging  ###################
    w.EWA<-rep(1/M,M)
    
    #################### EFPE ############################
    muhat.s<-matrix(0,n.s,M)
    for (m in 1:M) {
      muhat.s[,m]<-X.s %*% theta.hat[,m]
      EFPE0[r,m]<-FPE_er(tau,y.s,muhat.s[,m])-FPE
    }
    EFPEr.JCVMA5[r]<-FPE_er(tau,y.s,muhat.s%*%w.JCVMA5)-FPE
    EFPEr.JCVMA10[r]<-FPE_er(tau,y.s,muhat.s%*%w.JCVMA10)-FPE
    EFPEr.JMA[r]<-FPE_er(tau,y.s,muhat.s%*%w.JMA)-FPE
    EFPEr.SAIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SAIC)-FPE
    EFPEr.SBIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SBIC)-FPE
    EFPEr.AIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.AIC)-FPE
    EFPEr.BIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.BIC)-FPE
    EFPEr.MAP5[r]<-FPE_er(tau,y.s,muhat.s%*%w.MAP5)-FPE
    EFPEr.MAP10[r]<-FPE_er(tau,y.s,muhat.s%*%w.MAP10)-FPE
    EFPEr.EWA[r]<-FPE_er(tau,y.s,muhat.s%*%w.EWA)-FPE
    # test to see if it is working properly
    write(r, file = "test.txt",append = TRUE, sep =" ")
  }
  EFPE<-min(colMeans(EFPE0))
  EFPE.JCVMA5<-mean(EFPEr.JCVMA5)/EFPE
  EFPE.JCVMA10<-mean(EFPEr.JCVMA10)/EFPE
  EFPE.JMA<-mean(EFPEr.JMA)/EFPE
  EFPE.SAIC<-mean(EFPEr.SAIC)/EFPE
  EFPE.SBIC<-mean(EFPEr.SBIC)/EFPE
  EFPE.AIC<-mean(EFPEr.AIC)/EFPE
  EFPE.BIC<-mean(EFPEr.BIC)/EFPE
  EFPE.MAP5<-mean(EFPEr.MAP5)/EFPE
  EFPE.MAP10<-mean(EFPEr.MAP10)/EFPE
  EFPE.EWA<-mean(EFPEr.EWA)/EFPE
  return(c(n,M,tau,R.square,EFPE.JCVMA5,EFPE.JCVMA10,EFPE.JMA,
           EFPE.SAIC,EFPE.SBIC,EFPE.MAP5,EFPE.MAP10,
           EFPE.EWA,EFPE.AIC,EFPE.BIC))
}
#R.r<-2
#EFPE_parallel<-lapply(1:3, main_parallel) 
#EFPE<-matrix(unlist(EFPE_parallel),nrow=3,byrow = TRUE)
EFPE_parallel<-sfLapply(1:length(settings), main_parallel) 
EFPE<-matrix(unlist(EFPE_parallel),nrow=length(settings),byrow = TRUE)
EFPE<-data.frame(EFPE)
colnames(EFPE)<-c('n','M','tau','R2','JCVMA5','JCVMA10','JMA','SAIC','SBIC',
                  'MAP5','MAP10','EWA','AIC','BIC')
unlink("test.txt")
sfStop()

###################### Plot ##########################
###################### p=1 ##########################
for (tau in unique(EFPE$tau)) {
  pdf(paste('DGP2_QR_',tau,'.pdf',sep=""),width = 10,height = 8)
  par(mfrow = c(2, 2),oma = c(3.2, 0.5, 0.1, 0.6),xpd = NA,
      cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.4)  
  for (n in unique(EFPE$n)) {
    EFPEn<-EFPE[EFPE$n==n & EFPE$tau==tau,]
    # png(paste('DGP2_QR',n,tau,'.png'))
    plot(x = EFPEn$R2,y = EFPEn$JCVMA5,
         xlab =expression(R^2), ylab="normalized EFPE",col="red",
         lwd=2.5,type='o',lty=1,ylim = c(0.7,1.4),pch=4, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$JCVMA10,col="orange",lwd=2.5,type='o',lty=1,
          pch=2, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$JMA,col="purple",lwd=2.5,type='o',lty=1,
          pch=3, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$SAIC,col="blue",lwd=2.5,type='l',lty=2, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$SBIC,col="green",lwd=2.5,type='l',lty=3, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$MAP5,col="pink",lwd=2.5,type='l',lty=4, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$MAP10,col="lightblue",lwd=2.5,type='l',lty=5, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$EWA,col="black",lwd=2.5,type='l',lty=6, xpd=FALSE)
    title(main = substitute(paste('p=1, n=',n,', M=',M,', ',tau,'=',t), 
                            list(n=EFPEn$n,t=EFPEn$tau,M=EFPEn$M)))
    # dev.off()
    # legend('topright',c("JCVMA5", "JCVMA10",'JMA',"SAIC","SBIC",'MAP5','MAP10',"EWA"),
    #        lwd = 2,lty = c(1,1,1,2,3,4,5,6), pch = c(4,2,3,rep(NA,5)),
    #        col = c("red","orange",'purple','blue','green','pink','lightblue','black'))
  }
  legend(-0.65,0.44,c("JCVMA5", "JCVMA10",'JMA',"SAIC","SBIC",'MAP5','MAP10',"EWA"),
         lwd = 2.5,lty = c(1,1,1,2,3,4,5,6), pch = c(4,2,3,rep(NA,5)), ncol = 4,
         col = c("red","orange",'purple','blue','green','pink','lightblue','black'),
         cex=1.3)
  dev.off()
}

###################### p=2 ##########################
for (tau in unique(EFPE$tau)) {
  pdf(paste('DGP2_ER_',tau,'.pdf',sep=""),width = 10,height = 8)
  par(mfrow = c(2, 2),oma = c(3.2, 0.5, 0.1, 0.6),xpd = NA,
      cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.4)  
  for (n in unique(EFPE$n)) {
    EFPEn<-EFPE[EFPE$n==n & EFPE$tau==tau,]
    # png(paste('DGP1_QR',n,tau,'.png'))
    plot(x = EFPEn$R2,y = EFPEn$JCVMA5,
         xlab =expression(R^2), ylab="normalized EFPE",col="red",
         lwd=2.5,type='o',lty=1,ylim = c(0.7,1.4),pch=4, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$JCVMA10,col="orange",lwd=2.5,type='o',lty=1,
          pch=2, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$JMA,col="purple",lwd=2.5,type='o',lty=1,
          pch=3, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$SAIC,col="blue",lwd=2.5,type='l',lty=2, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$SBIC,col="green",lwd=2.5,type='l',lty=3, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$MAP5,col="pink",lwd=2.5,type='l',lty=4, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$MAP10,col="lightblue",lwd=2.5,type='l',lty=5, xpd=FALSE)
    lines(x = EFPEn$R2, y = EFPEn$EWA,col="black",lwd=2.5,type='l',lty=6, xpd=FALSE)
    title(main = substitute(paste('p=2, n=',n,', M=',M,', ',tau,'=',t), 
                            list(n=EFPEn$n,t=EFPEn$tau,M=EFPEn$M)))
    # dev.off()
    # legend('topright',c("JCVMA5", "JCVMA10",'JMA',"SAIC","SBIC",'MAP5','MAP10',"EWA"),
    #        lwd = 2,lty = c(1,1,1,2,3,4,5,6), pch = c(4,2,3,rep(NA,5)),
    #        col = c("red","orange",'purple','blue','green','pink','lightblue','black'))
  }
  legend(-0.65,0.44,c("JCVMA5", "JCVMA10",'JMA',"SAIC","SBIC",'MAP5','MAP10',"EWA"),
         lwd = 2.5,lty = c(1,1,1,2,3,4,5,6), pch = c(4,2,3,rep(NA,5)), ncol = 4,
         col = c("red","orange",'purple','blue','green','pink','lightblue','black'),
         cex=1.3)
  dev.off()
}




