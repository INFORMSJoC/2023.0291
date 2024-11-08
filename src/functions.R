# define the indicator function
Indic<-function(e){
  g<-rep(0,length(e))
  g[which(e<=0)]=1
  g
}

##turnbits_cir function is used to convert decimal to p-dimensional binary vector
##n is a decimal number, p is the binary vector dimension
##Output p-dimensional logical row vector
turnbits_cir <- function(n,p){
  z<-rep(0,p)           
  tn<-n                 
  z[1]<-tn%%2          
  for(j in 2:p){         
    tn<- (tn-z[j-1])/2  
    if(tn == 0) break    
    z[j]<-tn%%2         
  }
  as.logical(z)         
}

# FPE for p=1
FPE_qr<-function(tau,y,yhat){
  (tau-Indic(y-yhat))%*%(y-yhat)/length(y)
}

# conduct JCVMA & CV & MAP for p=1
JCVMA_qr<-function(X, y, J=5, tau, nest = TRUE, fixed = 1){
  ## nest: whether use nested candidate model
  ## fixed: the number of fixed covariates
  ## the intercept term should be included
  ######################  Data Preparation  #########################
  n<-nrow(X)
  k<-ncol(X)
  Q<-floor(n/J)
  if(nest){
    M<-k
    Z<-matrix(TRUE,k,k)
    Z[upper.tri(Z)] <-FALSE
  } else{
    k.unsure<-k-fixed
    Z<-matrix(unlist(lapply(0:(2^k.unsure-1),turnbits_cir,k.unsure)),ncol = k.unsure ,byrow = T)
    Z<-cbind(matrix(TRUE,2^k.unsure,fixed),Z) 
    M<-nrow(Z)
  }
  data<-data.frame(y,X)
  
  ######################  Divide Data  ###########################
  # divide the data into J groups, index represents the group index
  index <- rep(1:J,Q) 
  index <- sample(index,n)
  # initialize the matrix composed of mu_i(m)
  mu<-matrix(0,n,M)
  # initialize y used to conduct CV after rearranging
  y.CV<-rep(0,n) 
  # prepare data for JCV
  for (j in 1:J) {
    data.train<-data[index!=j,]
    data.test<-data[index==j,]
    for (m in 1:M) {
      X.train.m<-data.frame(data.train[,-1][,Z[m,]])
      X.test.m<-data.test[,-1][,Z[m,]]
      rq.m<-rq(data.train[,1]~0+., tau=tau, data=X.train.m)
      mu[((j-1)*Q+1):(j*Q),m]<-as.matrix(X.test.m)%*% coef(rq.m)
    }
    y.CV[((j-1)*Q+1):(j*Q)]<-data.test[,1]
  } 
  
  ####################   JCV weights    ###########################
  CV<-rep(0,len=M)
  w.CV<-rep(0,len=M)
  for (m in 1:M) {
    CV[m]<-FPE_qr(tau,y.CV,mu[,m])
  }
  w.CV[which.min(CV)]<-1
  
  ####################  JCVMA weights use LP ######################
  obj1<-rep(tau,len=n);obj2<-rep(1-tau,len=n);obj3<-rep(0,len=M)
  f.obj<-c(obj1,obj2,obj3)
  # ui,vi,wi>=0
  con1<-diag(2*n+M)
  # wi<=1
  con21<-matrix(0,nrow = M,ncol = 2*n);
  con22<-diag(M);
  con2<-cbind(con21,con22) 
  # sum(w)=1
  con3<-c(rep(0,len=2*n),rep(1,len=M))
  # sum(w*mu_i(m)) + u(j-1)Q+q - v(j-1)Q+q=y(j-1)Q+q
  con4<-cbind(diag(n),-diag(n),mu)
  f.con<-rbind(con1,con2,con3,con4) # constraint matrix
  f.dir<-c(rep(">=",len=2*n+M),rep("<=",len=M),rep("=",len=n+1) ) # direction
  f.rhs<-c(rep(0,len=2*n+M),rep(1,len=M+1), y.CV) # the right hand side of the constraints
  # solve LP
  lp.result<-lp("min", f.obj, f.con, f.dir, f.rhs)
  w.JCVMA<-lp.result$solution[(2*n+1):(2*n+M)]
  
  ####################  MAP weights use QP ######################
  a1<-t(y.CV-mu)%*%(y.CV-mu)
  MAP<-LowRankQP(Vmat=a1,dvec=rep(0,M),Amat=matrix(1,nrow=1,ncol=M),bvec=1,
                 uvec=rep(1,M),method="LU")
  w.MAP<-MAP$alpha
  result<-list(CV=w.CV,JCVMA=w.JCVMA,MAP=w.MAP)
  return(result)
}

# define expectile regression function
er <- function(X, y, tau, max.iter, tol){
  # X is the model matrix
  # y is the response vector of observed proportion
  # maxIter is the maximum number of iterations
  # tol is a convergence criterion
  b <- bLast <- rep(0, ncol(X)) # initialize
  it <- 1 # iteration index
  while (it <= max.iter){
    ypred <- c(X %*% b)
    w <- as.vector(tau *(y>= ypred) + (1-tau)* (y<ypred))
    b <- lsfit(X, y, w, intercept=FALSE)$coef
    if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) break
    bLast <- b
    it <- it + 1 # increment index
  }
  if (it > max.iter) warning('maximum iterations exceeded')
  # ## the loss function
  # loss <- sum(w*(y-c(X%*%b))^2)
  # 
  # ## the variance
  # Am <- t(X) %*% diag(w) %*% X
  # if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
  # 
  # A <- solve(Am)
  # H <- w * diag(X %*% A %*% t(X))
  # B <- t(X) %*% diag(w^2*(y-ypred)^2/(1-H)) %*% X
  # Vb <- A %*% B %*% A  # variance function
  # list(coefficients = b, variance = Vb, loss = loss, it = it)
  return(b)
}

# FPE for p=2
FPE_er<-function(tau,y,yhat){
  abs(tau-Indic(y-yhat))%*%(y-yhat)^2/length(y)
}

# conduct JCVMA & CV & MAP for p=2
JCVMA_er<-function(X, y, J=5, tau, nest = TRUE, fixed = 1){
  ## nest: whether use nested candidate model
  ## fixed: the number of fixed covariates
  ## the intercept term should be included
  ######################  Data Preparation  #########################
  n<-nrow(X)
  k<-ncol(X)
  Q<-floor(n/J)
  if(nest){
    M<-k
    Z<-matrix(TRUE,k,k)
    Z[upper.tri(Z)] <-FALSE
  } else{
    k.unsure<-k-fixed
    Z<-matrix(unlist(lapply(0:(2^k.unsure-1),turnbits_cir,k.unsure)),
              ncol = k.unsure ,byrow = T)
    Z<-cbind(matrix(TRUE,2^k.unsure,fixed),Z) 
    M<-nrow(Z)
  }
  data<-data.frame(y,X)
  
  ######################  Divide Data  ###########################
  # divide the data into J groups, index represents the group index
  index <- rep(1:J,Q) 
  index <- sample(index,n)
  # initialize the matrix composed of mu_i(m)
  mu<-matrix(0,n,M)
  # initialize y used to conduct CV after rearranging
  y.CV<-rep(0,n) 
  # prepare data for JCV
  for (j in 1:J) {
    data.train<-data[index!=j,]
    data.test<-data[index==j,]
    for (m in 1:M) {
      X.train.m<-as.matrix(data.train[,-1][,Z[m,]])
      X.test.m<-data.test[,-1][,Z[m,]]
      eq.m<-er(X=X.train.m,y=data.train[,1],tau=tau,max.iter = 100,tol = 1)
      mu[((j-1)*Q+1):(j*Q),m]<-as.matrix(X.test.m)%*%as.matrix(eq.m)
    }
    y.CV[((j-1)*Q+1):(j*Q)]<-data.test[,1]
  } 
  
  ####################   JCV weights    ###########################
  CV<-rep(0,len=M)
  w.CV<-rep(0,len=M)
  for (m in 1:M) {
    CV[m]<-FPE_er(tau,y.CV,mu[,m])
  }
  w.CV[which.min(CV)]<-1
  
  ####################  JCVMA weights use QP ######################
  # coefficient matrix
  H<-diag(c(rep(2*tau,n),rep(2*(1-tau),n),rep(0,M)))
  d<-rep(0,2*n+M)
  # linear constraint matrix
  # sum(w)=1
  con1<-c(rep(0,len=2*n),rep(1,len=M))
  # sum(w*mu_i(m)) + u(j-1)Q+q - v(j-1)Q+q=y(j-1)Q+q
  con2<-cbind(diag(n),-diag(n),mu)
  coneq<-rbind(con1,con2)
  # ui,vi,wi>=0,wi<=1
  beq<-c(1,y.CV)
  uvec<-c(rep(100,2*n),rep(1,M))
  # solve QP
  qp.result<-LowRankQP(Vmat=H,dvec=d,Amat=coneq,bvec=beq,uvec=uvec,method="LU")
  w.JCVMA<-qp.result$alpha[(2*n+1):(2*n+M)]
  
  ####################  MAP weights use QP ######################
  a1<-t(y.CV-mu)%*%(y.CV-mu)
  MAP<-LowRankQP(Vmat=a1,dvec=rep(0,M),Amat=matrix(1,nrow=1,ncol=M),bvec=1,
            uvec=rep(1,M),method="LU")
  w.MAP <- MAP$alpha
  result<-list(CV=w.CV,JCVMA=w.JCVMA,MAP=w.MAP)
  return(result)
}



