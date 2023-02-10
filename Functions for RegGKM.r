

###################################################
#  inner  functions for RegGKM
##################################################

### Data generating ###
# N for training
# Nvalid for validation
# Ntest for prediction
# X<-matrix(runif((N+Ntest+Nvalid)*P,-0.01,0.01),(N+Ntest+Nvalid),P);
datf <- function(datid,P,Q,PT,QT,N,Nvalid,Ntest){
	X<-matrix(rnorm((N+Ntest+Nvalid)*P,0,1),(N+Ntest+Nvalid),P);
	Z<-matrix(rnorm((N+Ntest+Nvalid)*Q,0,1),(N+Ntest+Nvalid),Q);
	true_beta<-c(rep(1,PT),rep(0,P-PT))
	XL<-cbind(X,Z)
	if(datid == 1){
		str <- 'H = 0.01 * (cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] )'
	}else if(datid == 2){
		str <- 'H = cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] + 0.9 * Z[,6] * sin(Z[,7]) - 0.8 * cos(Z[,6]) * Z[,7] + 2 * Z[,8] * sin(Z[,9]) * sin(Z[,10]) - 1.5 * Z[,8]^3 - Z[,8] * Z[,9] - 0.1 * exp(Z[,10]) * cos(Z[,10])'
	}else if(datid == 3){
		str <- 'H = 10*cos(Z[,1]) - 15*(Z[,2])^2 + 10*exp(-Z[,3])*Z[,4] - 8 * sin(Z[,5]) * cos(Z[,3]) + 20 * Z[,1] * Z[,5]'
	}else if(datid == 4){
		str <- 'H = - 0.01 * (cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] ) '   
	}else if(datid == 599){
		str <- 'H = 0.06*(10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2)' 
	}else if(datid == 699){
		str <- 'H = 5*cos(Z[,1])* Z[,2]^2 - 10* (Z[,1])^3 + 8*exp(Z[,1]) * cos(Z[,2]) - 8* sin(Z[,2]) * cos(Z[,3]) + 7* exp(Z[,3]) * sin(Z[,4])+6* Z[,2] ^2* sin(Z[,4])- 10* Z[,4]^3 -6* exp(Z[,4]) * cos(Z[,5])-6*sin(Z[,4])* Z[,5]^2' 
	}else if(datid == 899){
		str <- 'H = 0.06*(10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2)' 
	}else if(datid==799){
		str <- 'H = 2*cos(Z[,1])* Z[,2] - 2* (Z[,1])^2 + 3.5*exp(Z[,1]) * Z[,2] - 3* sin(Z[,2]) * cos(Z[,3]) + 2* exp(Z[,3]) * sin(Z[,4])+5* Z[,2] * sin(Z[,4])- 1.5* cos(Z[,3])*Z[,4]^2 -3* exp(Z[,4]) * cos(Z[,5])-4*sin(Z[,4])* Z[,5]^2' 
	}else if(datid == 100){
		str <- 'H = 10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2' 
	}else if(datid == 101){
		str <- 'H = 10 * Z[,1] + Z[,2]^2 + 20 * Z[,3] * Z[,4] + Z[,5]' 
	}
	eval(parse(text=str))
	fx <- as.vector(X%*%true_beta) # X * beta
	y <- fx + H                    # X * beta + h(Z)
	g <- exp(y)                    # hazard function
	D <- (unique(rexp((N+Ntest+Nvalid)*2,g)))[1:(N+Ntest+Nvalid)]
	U <- runif(N+Ntest+Nvalid,2000000000,3000000000)# censoring rate is about 0%
	C <- rexp(N+Ntest+Nvalid, 1/(U*g)) 
	TIM <- ifelse(D <= C, D, C)      # survival time
	TAU <- ifelse(D <= C, 1, 0)    # death indicator 
	train <- 1:N                                                # taining data id
	valid <- (N+1):(N+Nvalid)                                   # validation data id
	test <- (N + Nvalid+1):(N + Nvalid + Ntest)                 # testing data id
	ZT <- as.matrix(Z[test,]);Zvalid <- as.matrix(Z[valid,]);Z <- as.matrix(Z[train,]) # train data must be listed at last
	XT <- as.matrix(X[test,]);Xvalid <- as.matrix(X[valid,]);X <- as.matrix(X[train,])
	XLT <- as.matrix(XL[test,]);XLvalid <- as.matrix(XL[valid,]);XL <- as.matrix(XL[train,])
	TIMvalid <- TIM[valid]; TIMT <- TIM[test]; TIM <- TIM[train]
	TAUvalid <- TAU[valid]; TAUT <- TAU[test]; TAU <- TAU[train]
	fxT <- as.matrix(fx[test]);fxvalid <- as.matrix(fx[valid]);fx <- as.matrix(fx[train]) 
	HT <- as.matrix(H[test]);Hvalid <- as.matrix(H[valid]);H <- as.matrix(H[train]) 
	orde_valid <- order(TIMvalid); orde_t <- order(TIMT); orde <- order(TIM) # order for Tvalid, TT, T
	TIMvalid <- sort(TIMvalid); TIMT <- sort(TIMT); TIM<- sort(TIM)                # sort Tvalid, TT, T in increasing order
	ZT <- as.matrix(ZT[orde_t,]); Zvalid <- as.matrix(Zvalid[orde_valid,]); Z <- as.matrix(Z[orde,])    # sort Z based on T increasing order
	XT <- as.matrix(XT[orde_t,]); Xvalid <- as.matrix(Xvalid[orde_valid,]); X <- as.matrix(X[orde,])
	XLT <- as.matrix(XLT[orde_t,]); XLvalid <- as.matrix(XLvalid[orde_valid,]); XL <- as.matrix(XL[orde,])
	TAUT <- TAUT[orde_t]; TAUvalid <- TAUvalid[orde_valid]; TAU <- TAU[orde]
	HT <- as.matrix(HT[orde_t]); Hvalid <- as.matrix(Hvalid[orde_valid]); H <- as.matrix(H[orde])
	fxT <- as.matrix(fxT[orde_t]); fxvalid <- as.matrix(fxvalid[orde_valid]); fx <- as.matrix(fx[orde])  
	res <- list(X=X,Z=Z,XL=XL,TIM=TIM,TAU=TAU,XT=XT,ZT=ZT,XLT=XLT,TIMT=TIMT,TAUT=TAUT,Xvalid=Xvalid,Zvalid=Zvalid,XLvalid=XLvalid,TIMvalid=TIMvalid,TAUvalid=TAUvalid,HT=HT,Hvalid=Hvalid,H=H,fxT=fxT,fxvalid=fxvalid,fx=fx)
	res
}

# uf, used in gradient of f function in updating delta, should be listed in data.R
# u is a Q*N*N tensor
uf <- function(Z){  
	N <- nrow(Z)
	Q <- ncol(Z)
	u <- array(0,c(Q,N,N))
	tZ <-t(Z)
	for(i in 1:N){
		u[,,i]<-(Z[i,]-tZ)^2
	}
	u <- u/Q     
	return(u)
}

# Kernel matrix: N*N 
Kf <- function(Z,delta){
	N <- nrow(Z)
	Q <- ncol(Z)
	Zsqu <- Z^2 %*% delta %*% t(rep(1,N))
	res <- Zsqu -2 * (Z*matrix(delta,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
	K <- exp(-res/Q) # modi Q = 1000
	return(K)  # N*N
}

# h(\cdot): a N*1 vector
hf <- function(Z,alpha,delta){
	h <- Kf(Z,delta) %*% alpha 
	return(h) # N *1
}

# output:eta=X %*% betaa + hf(Z,alpha,delta)
# eta: a N*1 vector
etaf <- function(X,Z,alpha,betaa,delta){  # eta n*1 matrix!
	eta <- X %*% betaa + Kf(Z,delta) %*% alpha # modi 20220808
	return(eta)
}

#output:R=R, len=len, E=E, len_E=len_E
Rf <- function(TIM){   # Rf should be listed in data.R
	N <- length(TIM)
	R <- matrix(0,N,N)  # risk set
	len <- c()          # length of R[i,]
	for (i in 1:N){
		Ri <- which(TIM >= TIM[i]) 
		R[i,1:length(Ri)] <- Ri; 
	} 
	len <- apply(R,1,function(x) sum(x!=0))
	E <- matrix(0,N,N) # E_i: the risk sets containing individual i
	len_E <- c()
	for (i in 1:N){
		Ei <- c()
		for (m in 1:N){
			if( is.element(i,R[m,])){
				Ei=c(Ei,m)
			}else{
				Ei=Ei
			}	# needed to modified: ifelse
		}
		len.e <- length(Ei)
		E[i,1:len.e] <- Ei; 
		E[i,-(1:len.e)] <- 0
		len_E <- c(len_E,len.e)
	}
	return(list( R=R, len=len, E=E, len_E=len_E ))
}

# output: sum_exp is a N*1 vector, where i-th elements is sum_{l in R_i} exp(X_l beta+K_l alpha) for j in R_i
# similar to Simon et al. (2011)
sum.exp <- function(X,Z,R,len,alpha,betaa,delta){
	exp.eta <- exp(etaf(X,Z,alpha,betaa,delta)) + 1e-6
	Ra <- ifelse(R>0,1,0)
	sum_exp <- as.vector(Ra %*% exp.eta)
	return(sum_exp)  # N*1 vector 
}

#output: updated beta for given lambda_1
up.beta<-function(Y,X,Z,alpha,delta,lambda_1,weig){
	yy <- Y - hf(Z,alpha,delta)
	if(any(is.na(yy))|(sd(yy)==0)){ # 
		betaa <- rep(0,ncol(X))
	}else if(ncol(X)==1){
		fit <- lm(yy~X,weights=weig)  # 
		betaa <- as.vector(fit$coef[2])
		inter <- fit$coef[1]
	}else{
		glmfit <- glmnet(x=X,y=yy,family="gaussian",type.gaussian="naive",weights=weig,alpha=1,lambda=lambda_1) # ,standardize=FALSE)
		betaa <- as.vector(glmfit$beta)
		inter <- glmfit$a0
	}
	return(betaa)
}

#output: updated alpha for given lambda_3
up.alpha <- function(Y,X,Z,alpha,betaa,delta,lambda_3,weig){
	N <- length(Y)
	K = Kf(Z,delta)
	W <- diag(weig,nrow=N,ncol=N) #N*N weight matrix
	A = (1/N) * K %*% W %*% K + lambda_3 * K + diag(1e-6,N)
	B = (1/N)* K %*% W %*% (Y - X %*% betaa)
	alpha_n <- try(solve(A,B),silent=TRUE)
	if("try-error" %in% class(alpha_n)){
		alpha_n <- alpha
	}
	if(sum(alpha_n)!=0){
		alpha_n <- alpha_n/sum(alpha_n)
	}
	return(alpha_n)  # N*1 matrix
}

##### output: updated delta for given lambda_2,lambda_3 # spg
# spg function from BB
up.delta.spg <- function(X,Z,u,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	Q <- ncol(Z)
	u <- u; TAU <- TAU; R <- R; len <- len; alpha <- alpha; betaa <- betaa; delta <- delta; lambda_2 <- lambda_2;lambda_3<- lambda_3;X <- X; Z <- Z
	
	f <- function(delta_n,X,Z,u,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){
		sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
		A <- log(sum_exp)
		res <- (1/N) * (t(TAU) %*% (X %*% betaa + hf(Z,alpha,delta_n) - A)) - lambda_2 * sum(delta_n) - (1/2) *lambda_3 * (t(alpha) %*% Kf(Z,delta_n) %*% alpha)
		res
	}

	fg <- function(delta_n,X,Z,u,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){      		## gradient of function 'f'
		eta <- etaf(X,Z,alpha,betaa,delta_n)
		sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
		I = rep(1,Q)
		V <- matrix(1,N,N)
		V[lower.tri(V)]=0
		G <- t(TAU) * (1/t(sum_exp)) %*% V * exp(t(eta)) * t(alpha)
		IKF <-  I %*% t(as.vector(Kf(Z,delta_n))) * array(u,c(Q,N*N))
		res <-  as.vector( -(1/N)* IKF %*% kronecker(TAU,alpha) + (1/N)*IKF %*% (rep(G,each=N) * rep(alpha,N)) - lambda_2 * I  + (1/2) * lambda_3 *IKF %*% kronecker(alpha,alpha) )	  
		res
	}
	fg_res <- fg(delta,X,Z,u,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q)
	if ((max(fg_res) < 0) & (max(delta) ==0)){
		delta_n <- delta
	} else if((min(fg_res) > 0) & (min(delta)==1)) {
		delta_n <- delta
	} else {
		spg_res <- spg(par=delta,fn=f,gr=fg,X=X,Z=Z,u=u,TAU=TAU,R=R,len=len,alpha=alpha,betaa=betaa,lambda_2=lambda_2,lambda_3=lambda_3,N=N,Q=Q,lower=rep(0,Q),upper=rep(1,Q),control=list(maximize=TRUE,checkGrad=FALSE,eps=1e-5,ftol=1e-9,trace=FALSE),quiet=TRUE) #
		delta_n <- spg_res$par
	}
	delta_n <- ifelse(delta_n < 1e-4,0,delta_n)  
	return(delta_n)
}

# optimize for one q step
up.delta.optim <- function(X,Z,u,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	Q <- ncol(Z)
	delta_n <- delta
	for(q in 1:Q){ 
		f <- function(delta_hat){
			delta_n[q] <- delta_hat
			sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
			A <- log(sum_exp)
			res <- (1/N) * (t(TAU) %*% (X %*% betaa + hf(Z,alpha,delta_n) - A)) - lambda_2 * sum(delta_n) - (1/2) *lambda_3 * (t(alpha) %*% Kf(Z,delta_n) %*% alpha)
			-res  # optim default minimized
		}
		fg <- function(delta_hat){      		## gradient of function 'f'
			delta_n[q] <- delta_hat
			eta <- etaf(X,Z,alpha,betaa,delta_n)
			sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
			V<-matrix(1,N,N)
			V[lower.tri(V)]=0
			G <- V%*%(exp(eta) * ( ((Kf(Z,delta_n)) * ((-1) * (u[q,,])) ) %*% alpha ))
			res <- (1/N)*(t(TAU) %*% (Kf(Z,delta_n) * (-u[q,,])) %*% alpha) - (1/N)*(t(TAU) %*% (G/sum_exp))- lambda_2 - (1/2) * lambda_3 *(t(alpha) %*% (Kf(Z,delta_n) * (-u[q,,])) %*% alpha) 
			-res 
		}
		delta_n[q] <- optimize(f,c(0,1))$minimum
	}
	delta <- ifelse(delta_n<1e-4,0,delta_n)
	return(delta)
}

# optimParallel function from 
up.delta.optimParallel <- function(X,Z,u,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	Q <- ncol(Z)
	f <- function(delta_n,X,Z,u,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){
		Zsqu <- Z^2 %*% delta_n %*% t(rep(1,N))
		res <- Zsqu -2 * (Z*matrix(delta_n,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
		K <- exp(-res/Q) # modi Q = 1000
		h <- K %*% alpha
		eta <- X %*% betaa + h
		exp.eta <- exp(eta) + 1e-5
		sum_exp <- as.vector(ifelse(R>0,1,0) %*% exp.eta)
		A <- log(sum_exp)
		res <- (1/N) * (t(TAU) %*% (X %*% betaa + h - A)) - lambda_2 * sum(delta_n)
		- (1/2) *lambda_3 * (t(alpha) %*% K %*% alpha)
		-res
	}
	fg <- function(delta_n,X,Z,u,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){      		## gradient of function 'f'
		Zsqu <- Z^2 %*% delta_n %*% t(rep(1,N))
		res <- Zsqu -2 * (Z*matrix(delta_n,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
		K <- exp(-res/Q) # modi Q = 1000
		h <- K %*% alpha
		eta <- X %*% betaa + h
		exp.eta <- exp(eta) + 1e-5
		sum_exp <- as.vector(ifelse(R>0,1,0) %*% exp.eta)
		I = rep(1,Q)
		V <- matrix(1,N,N)
		V[lower.tri(V)]=0
		G <- t(TAU) * (1/t(sum_exp)) %*% V * exp(t(eta)) * t(alpha)
		IKF <-  I %*% t(as.vector(K)) * array(u,c(Q,N*N))
		res <-  as.vector( -(1/N)* IKF %*% kronecker(TAU,alpha) + (1/N)*IKF %*% (rep(G,each=N) * rep(alpha,N)) - lambda_2 * I  + (1/2) * lambda_3 *IKF %*% kronecker(alpha,alpha) )	  
		-res
	}
	cl <- makeCluster(1)
	setDefaultCluster(cl=cl)
	opt.par <-optimParallel(par=delta,fn=f,gr=fg,X=X,Z=Z,u=u,TAU=TAU,R=R,len=len,alpha=alpha,betaa=betaa,lambda_2=lambda_2,lambda_3=lambda_3,N=N,Q=Q,lower=rep(0,Q),upper=rep(1,Q),control = list(factr=1e-5))	
	delta_n <- opt.par$par
	print(paste(opt.par$value,opt.par$message,opt.par$convergence,sep=','))
	setDefaultCluster(cl=NULL)
	stopCluster(cl)
	delta <- ifelse(delta_n<1e-4,0,delta_n)
	return(delta)
}

#input: X,Z,TAU,R,len,E,len_E,alpha,betaa,delta
#output: Y, weigths, A
#first and second derivatives of log-partial likelihood, working response, weight matrix
#first derivative: grad, n*1 vector; second derivative: hess
upda.esti <- function(X,Z,TAU,R,len,E,len_E,alpha,betaa,delta){
	N <- dim(X)[1]
	eta <- as.vector(etaf(X,Z,alpha,betaa,delta))
	expeta <- exp(eta) + 1e-5 # modi 20220809 
	Ea <- ifelse(E!=0,1,0)
	Ra <- ifelse(R!=0,1,0)
	sum_E <- (Ea %*% TAU) / (Ra %*% expeta)
	sum_Esqua <- (Ea %*% TAU) / (Ra %*% expeta)^2
	grad <- TAU - expeta * sum_E
	hess <- expeta * (expeta * sum_Esqua - sum_E)   
	hess <- ifelse(abs(hess) < 1e-5,(-1)*1e-5, hess)
	Y <- eta - as.vector(grad / hess)
	weig <- as.vector((-1) * hess)   # N*1 vector: whose elements could not be zero! 
	return(list(Y=Y, weig=weig))#, A=A))   
}

### update \beta, \alpha and \delta
Qf <- function(X,Z,u,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,lambda_1,lambda_2,lambda_3){
	N <- length(TAU)
	P <- ncol(X)
	Q <- ncol(Z)
	beta_n <- beta_nn
	delta_n <- delta_nn
	alpha_n <- alpha_nn
	A_n  <- log(sum.exp(X,Z,R,len,alpha_n,beta_n,delta_n))# old
	like_new <- (1/N) * (t(TAU)%*% (X %*% beta_n + hf(Z,alpha_n,delta_n) - A_n)) -  lambda_1  * sum(abs(beta_n)) - lambda_2 * sum(delta_n) - 1/2 *lambda_3 * t(alpha_n) %*% Kf(Z,delta_n) %*% alpha_n
	iter = 0   ### iteration
	if(is.na(like_new)){
		mess = 1
		break
	}
	repeat{
		alpha <- alpha_n
		betaa <- beta_n
		delta <- delta_n
		delta <- delta_n
		A <- A_n
		like_old <- like_new
		iter = iter + 1
		upda_esti_r <- upda.esti(X,Z,TAU,R,len,E,len_E,alpha,betaa,delta) 		# update Y,weights
		Y <- upda_esti_r$Y; 
		weig <- upda_esti_r$weig
		beta_n <- up.beta(Y,X,Z,alpha,delta,lambda_1,weig)
		alpha_n <- up.alpha(Y,X,Z,alpha,betaa,delta,lambda_3,weig)
		delta_n <- up.delta(X,Z,u,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3)
		A_n <- log(sum.exp(X,Z,R,len,alpha_n,beta_n,delta_n))# new
		like_new <- (1/N) * (t(TAU)%*% (X %*% beta_n + hf(Z,alpha_n,delta_n) - A_n)) -  lambda_1  * sum(abs(beta_n)) - lambda_2 * sum(delta_n) - 1/2 *lambda_3 * t(alpha_n) %*% Kf(Z,delta_n) %*% alpha_n
		if(is.na(like_new)){
			mess = 1
			break
		}
		like_diff <- like_new - like_old  
		if(abs(like_diff) <= 1e-04){# | (iter>100)){# | (par_diff < 1e-4)){   modi 20230207 不考虑参数大小。
			mess = 0
			break
		}else if(iter > 100){
			mess = 1
			break
		}
	}
	result<-list(mess=mess,like_new=like_new,iter=iter,betaa=beta_n,alpha=alpha_n,delta=delta_n)  
	return(result)
}

##############
###caculate the estimate of h(Z) and A on testing ZT or validation Zvalid data for gaussian kernel
##############
# R,len for validation or testing time (Tvalid,TT)
hathaf <- function(Xvalid,Z,Zvalid,Rvalid,lenvalid,alpha,betaa,delta){ 
	Nvalid <- nrow(Zvalid)
	N <- nrow(Z)
	zf <- t(Z^2 %*% delta %*% t(rep(1,Nvalid)))
	ztf <- Zvalid^2 %*% delta %*% t(rep(1,N))
	zz <- 2*Zvalid %*% diag(delta,length(delta)) %*% t(Z)
	hat_h <- exp(-(zf + ztf - zz)) %*% alpha #Nvalid*1 matrix= Nvalid*N matrix %*% N*1 vector
	sum_exp <- c()
	for (i in 1:Nvalid) {
		l <- Rvalid[i, (1:lenvalid[i])]   # lenT[i] vector
		sum_expi <- sum(exp( as.matrix(Xvalid[l,]) %*% betaa + hat_h[l]))
		sum_exp <-c(sum_exp,sum_expi)
	}
	hat_A <- as.matrix(log(sum_exp)) # Nvalid*1 matrix
	return(list( hat_h=hat_h, hat_A=hat_A ) ) 
}

# validated partial log-likelihood vpll=(1/Nvalid)*(t(TAU) %*% (Xvalid %*% betaa + hathf(Z,Zvalid,alpha,delta) - A_valid)) (in real data analysis, use cross-validated partial log-likelihood)
# Q.func: generalized error function to be minimized, -1* validated partial log-likelihood
#                  traning: X,Z,u,TAU,R,len,E,len_E,
#                  initial value: alpha_nn,beta_nn,delta_nn,
#                  validation: Xvalid,Zvalid,TAUvalid,Rvalid,lenvalid
# arguments for Q.func are lambda_1,lambda_2,lambda_3
# Q.func must satisfy:(1) containing arguments(bounds=bounds,parms.coding=parms.coding,seed=seed,maxevals=maxevals,verbose=verbose)
#                     (2) must use q.val  the generalized error and must be returned as "list" format
# Q.func <- function(lambda_1,lambda_2,lambda_3,
#                    X,Z,u,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,
#                    Xvalid,Zvalid,TAUvalid,Rvalid,lenvalid,
#                    bounds=bounds,parms.coding=parms.coding,seed=seed,maxevals=maxevals,verbose=verbose)
lamf<- function(X,Z,TAU,TIM,Xvalid,Zvalid,TAUvalid,TIMvalid,lambda_1_seq,lambda_2_seq,lambda_3_seq){
	N<-nrow(Z)# number of training observations
	Nvalid <- nrow(Zvalid)# number of testing observations
	if(!is.matrix(Z)){Z <- matrix(Z,nrow=N)}
	if(!is.matrix(X)){X <- matrix(X,nrow=N)}
	P <- ncol(X)
	Q <- ncol(Z)
	u <- uf(Z)
	LH <- c()
	ITER<-c()
	BETA <- matrix(ncol=P)
	DELTA <- matrix(ncol=Q)
	ALPHA <- matrix(ncol=N)
	LAMBDA <- matrix(ncol=3)
	R<-Rf(TIM)$R
	len<-Rf(TIM)$len
	E<-Rf(TIM)$E
	len_E<-Rf(TIM)$len_E
	beta_nn <- runif(P,0,1) # initial
	alpha_nn <- runif(N,0,1)
	delta_nn <- runif(Q,0,1)
	for(lambda_1 in lambda_1_seq){
		for(lambda_2 in lambda_2_seq){
			for(lambda_3 in lambda_3_seq){
				a_b_d <- Qf(X,Z,u,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,lambda_1,lambda_2,lambda_3)
				alpha <- as.vector(a_b_d$alpha) 
				betaa <- as.vector(a_b_d$betaa)
				delta <- as.vector(a_b_d$delta)
				iter  <- a_b_d$iter
				ITER<-c(ITER,iter)
				Rfvalid <- Rf(TIMvalid)
				Rvalid<-Rfvalid$R
				lenvalid<-Rfvalid$len
				hath <- hathaf(Xvalid,Z,Zvalid,Rvalid,lenvalid,alpha,betaa,delta)	
				hat_h <- hath$hat_h # hat_h is a estimate of h(Zvalid)
				hat_A <- hath$hat_A # hat_A is a estimate of A_valid
				LH<-c(LH,c((-1/Nvalid) * (t(TAUvalid) %*% ( Xvalid%*%betaa+hat_h-hat_A ))))
				BETA <- rbind(BETA,betaa)
				DELTA <- rbind(DELTA,delta)
				ALPHA <- rbind(ALPHA,alpha)
				LAMBDA <- rbind(LAMBDA,c(lambda_1,lambda_2,lambda_3))
				alpha_nn <- alpha  # warmstart
				beta_nn <- betaa
				delta_nn <- delta
			}
		}
	}
	id <- which.min(LH) 
	res <- list(LH = LH[id],BETA = BETA[id+1,],DELTA = DELTA[id+1,],ALPHA=ALPHA[id+1,],LAMBDA = LAMBDA[id+1,],ITER = ITER[id]) 
	res
}

cou <- function(x,PT){
	a <- all(x[1:PT]!=0)
	b <- all(x[-(1:PT)]==0)
	c <- any(x[1:PT]==0)
	c(a*b,a-a*b,c)
}

#lasso-COX
hata_lasso <- function(XLvalid,Rvalid,lenvalid,betaa){ 
	Nvalid <- nrow(XLvalid)
	sum_exp <- c()
	for (i in 1:Nvalid) {
		l <- Rvalid[i, (1:lenvalid[i])]   # lenT[i] vector
		sum_expi <- sum(exp(as.matrix(XLvalid [l,]) %*% betaa))
		sum_exp <-c(sum_exp,sum_expi)
	}
	hat_A <- as.matrix(log(sum_exp)) # Nvalid*1 matrix
	return(list(hat_A=hat_A)) 
}

lamf_lasso<- function(XL,TIM,TAU,XLvalid,TIMvalid,TAUvalid,lambda_lasso){
	Nvalid<-length(TIMvalid) # number of validation observations
	D <- ncol(XL)
	LH<-c()
	BETA <- matrix(ncol=D)
	LAMBDA <- matrix(ncol=1)
	S<-cbind(time = TIM, status = TAU)
	for(lambda_1 in lambda_lasso){
		a_b_c<-glmnet(x=XL,y=S,family = "cox",alpha=1,lambda = lambda_1,standardize=FALSE)
		betaa <- as.vector(a_b_c$beta)
		Rfvalid <- Rf(TIMvalid)
		Rvalid<-Rfvalid$R
		lenvalid<-Rfvalid$len
		hath <- hata_lasso(XLvalid,Rvalid,lenvalid,betaa)	
		hat_A <- hath$hat_A 
		LH<-c(LH,c((-1/Nvalid)* (t(TAUvalid) %*% (XLvalid%*%betaa-hat_A))))
		BETA <- rbind(BETA,betaa)
		LAMBDA <- rbind(LAMBDA,c(lambda_1))
	}
	id <- which.min(LH)
	res <- list(LH = LH[id],BETA = BETA[id+1,],LAMBDA = LAMBDA[id+1,]) 
	res
}

# cou for lasso
coul <- function(x,PT,P,QT,Q){
	a <- all(x[1:PT]!=0 & x[(P+1):(P+QT)]!=0)
	b <- all(x[(PT+1):P]==0 & x[(P+QT+1):(P+Q)]==0)
	c <- any(x[1:PT]==0 | x[(P+1):(P+QT)]==0)
	c(a*b,a-a*b,c)
}

# standardize funciton
ST<-function(dat){ 
	dat <- apply(dat,2,function(x) (x-min(x))/(max(x)-min(x)))
	dat
}

hata_cosso <- function(ceta,RT,lenT){ 
	NT <- nrow(ceta)
	sum_exp <- c()
	for (i in 1:NT) {
		l <- RT[i, (1:lenT[i])]   # lenT[i] vector
		sum_expi <- sum(exp(as.matrix(ceta[l,])))
		sum_exp <-c(sum_exp,sum_expi)
	}
	hat_A <- as.matrix(log(sum_exp)) # Nvalid*1 matrix
	return(list(hat_A=hat_A)) 
}