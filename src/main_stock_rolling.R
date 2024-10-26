#####################  导入需要的包  ################
library(quantreg) 
library(lpSolve) 
library(Hmisc)
library(parallel) # for detectCores() 
library(snowfall) # parallel programming 
library(purrr) # for cross2()
library(LowRankQP) # for quadratic programming
library(expectreg)
library(quadprog)
source('functions.R')

# 读取数据
data_sampled<-read.csv('data_sampled.csv',header = TRUE, fileEncoding = "GBK")
data_sampled<-data_sampled[,-1]
sales<-data_sampled[,1]
covariate<-data_sampled[,-1]

# 赋初值
n.all<-dim(data_sampled)[1]
M<-ncol(covariate)+1 # 模型个数
settings<-list(n=c(60,240),
               tau=c(0.95,0.9,0.5,0.1,0.05))
settings<-cross2(settings$n,settings$tau,.filter = NULL)
settings<-lapply(settings, unlist, 1)

#################### main p=1 ####################
sfInit(parallel = TRUE, cpus = detectCores()-40) 
sfLibrary(LowRankQP)
sfLibrary(quantreg)
sfLibrary(lpSolve)
sfExport('Indic','turnbits_cir','FPE_qr','JCVMA_qr')
sfExport('M','n.all','settings','sales','covariate')

main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  R.r<-n.all-n
  FPEr.JCVMA5<-numeric(R.r)
  FPEr.JMA<-numeric(R.r)
  FPEr.SAIC<-numeric(R.r)
  FPEr.SBIC<-numeric(R.r)
  FPEr.AIC<-numeric(R.r)
  FPEr.BIC<-numeric(R.r)
  FPEr.MAP5<-numeric(R.r)
  FPEr.EWA<-numeric(R.r)
  FPEr.LM<-numeric(R.r)
  FPEr.AVE<-numeric(R.r)
  #######################################################
  for (r in 1:R.r){
    ################ Recursive window  ####################
    X<-covariate[r:(n-1+r),]
    # X<-jitter(X)
    X<-as.matrix(cbind(1,X))
    X.s<-covariate[n+r,]
    X.s<-as.matrix(cbind(1,X.s))
    y<-sales[r:(n-1+r)]
    y.s<-sales[n+r]
    data<-data.frame(y,X)
    # 参数估计theta.hat
    theta.hat<-matrix(0,M,M)
    for (m in 1:M) {
      fit<-rq(data=data.frame(X[,1:m]), y~0+.,tau=tau)
      theta.hat[1:m,m]<-coef(fit)
    }
    
    ####################  JCVMA ######################
    JCVMA5<-JCVMA_qr(X,y,J=5,tau=tau)
    w.JCVMA5<-JCVMA5$JCVMA
    w.MAP5<-JCVMA5$MAP
    w.JMA<-JCVMA_qr(X,y,J=n,tau=tau)$JCVMA
    
    ######################### (S)AIC & (S)BIC    ###################
    AIC<-rep(0,len=M)
    BIC<-rep(0,len=M)
    for (m in 1:M) {
      AIC[m]<-2*(n)*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+2*m
      BIC[m]<-2*(n)*log(FPE_qr(tau,y,X%*%theta.hat[,m]))+m*log(n)
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
    FPEr.JMA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.JMA)
    FPEr.SAIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SAIC)
    FPEr.SBIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.SBIC)
    FPEr.AIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.AIC)
    FPEr.BIC[r]<-FPE_qr(tau,y.s,muhat.s%*%w.BIC)
    FPEr.MAP5[r]<-FPE_qr(tau,y.s,muhat.s%*%w.MAP5)
    FPEr.EWA[r]<-FPE_qr(tau,y.s,muhat.s%*%w.EWA)
    FPEr.LM[r]<-FPE_qr(tau,y.s,muhat.s%*%w.LM)
    FPEr.AVE[r]<-FPE_qr(tau,y.s,quantile(y,tau)) #historical quantile/expectile returns
    write(r, file = "test.txt",append = TRUE, sep =" ") # 测试是否正常运行
  }
  FPE.JCVMA5<-mean(FPEr.JCVMA5)
  FPE.JMA<-mean(FPEr.JMA)
  FPE.SAIC<-mean(FPEr.SAIC)
  FPE.SBIC<-mean(FPEr.SBIC)
  FPE.AIC<-mean(FPEr.AIC)
  FPE.BIC<-mean(FPEr.BIC)
  FPE.MAP5<-mean(FPEr.MAP5)
  FPE.EWA<-mean(FPEr.EWA)
  FPE.LM<-mean(FPEr.LM)
  FPE.AVE<-mean(FPEr.AVE)
  return(c(n,tau,FPE.JCVMA5,FPE.JMA,
           FPE.SAIC,FPE.SBIC,FPE.MAP5,
           FPE.EWA,FPE.AIC,FPE.BIC,FPE.LM,FPE.AVE))
}

#FPE_parallel<-lapply(1:length(settings), main_parallel) 
FPE_parallel<-sfLapply(1:length(settings), main_parallel) 
FPE<-matrix(unlist(FPE_parallel),nrow=length(settings),byrow = TRUE)
FPE<-data.frame(FPE)
colnames(FPE)<-c('n','tau','JCVMA5','JMA','SAIC','SBIC',
                 'MAP5','EWA','AIC','BIC','LM','AVE')
unlink("test.txt")
sfStop()
FPE_result<-data.frame(p=1,FPE)
View(FPE_result)
write.csv(FPE_result, 'FPE_result.csv',fileEncoding = "GBK")
FPE_normalized<-cbind(FPE_result[,c(1,2,3)],FPE_result[,-c(1,2,3,10,11,12,13)]/FPE_result$LM)

#################### main p=2 ####################
sfInit(parallel = TRUE, cpus = detectCores()-20) 
sfLibrary(LowRankQP)
sfLibrary(expectreg)
sfExport('Indic','turnbits_cir','er','FPE_er','JCVMA_er')
sfExport('M','n.all','settings','sales','covariate')

main_parallel <-function(ii){
  ################### initialization ########################
  setting<-settings[[ii]]
  n<-setting[1]
  tau<-setting[2]
  R.r<-n.all-n
  FPEr.JCVMA5<-numeric(R.r)
  FPEr.JMA<-numeric(R.r)
  FPEr.SAIC<-numeric(R.r)
  FPEr.SBIC<-numeric(R.r)
  FPEr.AIC<-numeric(R.r)
  FPEr.BIC<-numeric(R.r)
  FPEr.MAP5<-numeric(R.r)
  FPEr.EWA<-numeric(R.r)
  FPEr.LM<-numeric(R.r)
  FPEr.AVE<-numeric(R.r)
  #######################################################
  for (r in 1:R.r){
    ################ Recursive window  ####################
    X<-covariate[r:(n-1+r),]
    X<-as.matrix(cbind(1,X))
    # X[,9]<-jitter(X[,9])
    X.s<-covariate[n+r,]
    X.s<-as.matrix(cbind(1,X.s))
    y<-sales[r:(n-1+r)]
    y.s<-sales[n+r]
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
    w.JMA<-JCVMA_er(X,y,J=n,tau=tau)$JCVMA
    
    ######################### (S)AIC & (S)BIC    ###################
    AIC<-rep(0,len=M)
    BIC<-rep(0,len=M)
    for (m in 1:M) {
      AIC[m]<-(n)*log(FPE_er(tau,y,X%*%theta.hat[,m]))+2*m
      BIC[m]<-(n)*log(FPE_er(tau,y,X%*%theta.hat[,m]))+m*log(n)
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
    FPEr.JCVMA5[r]<-FPE_er(tau,y.s,muhat.s%*%w.JCVMA5)
    FPEr.JMA[r]<-FPE_er(tau,y.s,muhat.s%*%w.JMA)
    FPEr.SAIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SAIC)
    FPEr.SBIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.SBIC)
    FPEr.AIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.AIC)
    FPEr.BIC[r]<-FPE_er(tau,y.s,muhat.s%*%w.BIC)
    FPEr.MAP5[r]<-FPE_er(tau,y.s,muhat.s%*%w.MAP5)
    FPEr.EWA[r]<-FPE_er(tau,y.s,muhat.s%*%w.EWA)
    FPEr.LM[r]<-FPE_er(tau,y.s,muhat.s%*%w.LM)
    FPEr.AVE[r]<-FPE_er(tau,y.s,expectile(y,tau)) #historical quantile/expectile returns
    write(r, file = "test.txt",append = TRUE, sep =" ") # 测试是否正常运行
  }
  FPE.JCVMA5<-mean(FPEr.JCVMA5, na.rm = TRUE)
  FPE.JMA<-mean(FPEr.JMA, na.rm = TRUE)
  FPE.SAIC<-mean(FPEr.SAIC)
  FPE.SBIC<-mean(FPEr.SBIC)
  FPE.AIC<-mean(FPEr.AIC)
  FPE.BIC<-mean(FPEr.BIC)
  FPE.MAP5<-mean(FPEr.MAP5)
  FPE.EWA<-mean(FPEr.EWA)
  FPE.LM<-mean(FPEr.LM)
  FPE.AVE<-mean(FPEr.AVE)
  return(c(n,tau,FPE.JCVMA5,FPE.JMA,
           FPE.SAIC,FPE.SBIC,FPE.MAP5,
           FPE.EWA,FPE.AIC,FPE.BIC,FPE.LM,FPE.AVE))
}

#FPE_parallel<-lapply(7, main_parallel) 
FPE_parallel<-sfLapply(1:length(settings), main_parallel) 
FPE<-matrix(unlist(FPE_parallel),nrow=length(settings),byrow = TRUE)
FPE<-data.frame(FPE)
colnames(FPE)<-c('n','tau','JCVMA5','JMA','SAIC','SBIC',
                 'MAP5','EWA','AIC','BIC','LM','AVE')
unlink("test.txt")
sfStop()

FPE_result<-rbind(FPE_result,data.frame(p=2,FPE))
View(FPE_result)
write.csv(FPE_result, 'FPE_result.csv',fileEncoding = "GBK")
FPE_normalized<-cbind(FPE_result[,c(1,2,3)],FPE_result[,-c(1,2,3,10,11,12,13)]/FPE_result$LM)

# FPE_final<-FPE_normalized[,-which(colnames(FPE_normalized)%in%c('JCVMA10','MAP10') )]


######################### latex table  #######################
latex_table<-function(data){
  cat("\\begin{table}[htbp]\n")
  cat("\\centering\n")
  cat("\\caption{relative FPE of forecast of excess stock returns 
      in Design I}\n")
  cat("\\label{table:Rsquare} \n")
  cat("\\begin{tabular}{lllllllll}\n")
  cat("\\hline\n")
  cat("$p$ & $\\tau$  & $n_1$ & JCVMA5  & JMA & SAIC 
		& SBIC & MAP5  & EWA   \\\\\n")
  cat("\\hline\n")
  for (i in 1:nrow(data)) {
    cat(paste(data[i,1],"& ", data[i,2],"& ", 
              data[i,3],"& "))
    first<-order(as.numeric(data[i,-c(1,2,3)]), decreasing = FALSE)[1]
    second<-order(as.numeric(data[i,-c(1,2,3)]), decreasing = FALSE)[2]
    third<-order(as.numeric(data[i,-c(1,2,3)]), decreasing = FALSE)[3]
    for (j in 4:8) {
      if(j==first+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[1]}$ & ", sep="")
      else if(j==second+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[2]}$ & ", sep="")
      else if(j==third+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[3]}$ & ", sep="")
      else
        cat(sprintf("%.3f", round(data[i,j],3))," & ",sep="")
    }
    for (j in 9) {
      if(j==first+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[1]}$", sep="")
      else if(j==second+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[2]}$", sep="")
      else if(j==third+3)
        cat(sprintf("%.3f", round(data[i,j],3)),"$^{[3]}$", sep="")
      else
        cat(sprintf("%.3f", round(data[i,j],3))," ",sep="")
    }
    cat("\\\\\n")
  }
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}
# latex_table(FPE_final)
latex_table(FPE_normalized)




