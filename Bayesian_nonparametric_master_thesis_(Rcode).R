rm(list=ls())
library(MCMCpack) #to generate sig.j rinvgamma
library(geoR) #to generate sig rinvchisq

#--------------------------------------------------------------------------
#basic function
#--------------------------------------------------------------------------

#log(PI(V_j))
PI <- function(V,j,vec=TRUE){
  out <- cumsum(log(c(1,1-V[1:(j-1)])))+log(V[1:j])
  if(vec==FALSE) out <- out[j]
  return(out)
}

#--------------------------------------------------------------------------
#simulation data
#jump distribution
samp.I2 <- function(t,m,V.j){ #(V.j,m)
  pr    <- PI(V.j,m) 
  pr    <- exp(pr)   
  return(sample(1:m,t,prob=pr,replace=TRUE))
}

#--------------------------------------------------------------------------
#simulation data
samp.tr.J.temp<-function(x){ #c(tr.j)
  j<-x
  return(rnorm(1,tr.mu.j[j],sqrt(tr.sig.j[j])))  
}

#--------------------------------------------------------------------------
samp.I <- function(x){ #(J,Z)
  J     <- x[1]
  Z     <- x[2]
  pr    <- PI(V.j,m,vec=TRUE)+Z*(dnorm(J,mu.j,sqrt(sig.j),log=T)+1/2*log(sig*sig.j/(sig.j+sig)))
  pr    <- exp(pr-max(pr))
  return(sample(1:m,1,prob=pr))
}
#--------------------------------------------------------------------------
#samp Z
#> vtom(z)
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    0    0    1    0    0


vtom <- function(z){ 
  out <- matrix(0,t,m)
  out[cbind(1:t,z)] <- 1
  return(out)
}

#--------------------------------------------------------------------------
# sample V
samp.V<-function(x){#c(J-1)
  j<-x
  return(rbeta(1,B[j]+1,t-sum(B[1:j])+gam))
}


#--------------------------------------------------------------------------
# sample J
samp.J<-function(x){ #c(Yobs,II)
  Yobs  <-x[1]
  j     <-x[2]
  return(rnorm(1,(sig.j[j]*(Yobs-mu)+sig*mu.j[j])/(sig.j[j]+sig),sqrt(sig*sig.j[j]/(sig.j[j]+sig))))
}

#--------------------------------------------------------------------------
samp.mu.j <-function(x){ #c(J)
  j     <-x
  denom <-tau2.1*length(J[(II==j)&(Z==1)])+sig.j[j]
  return(rnorm(1,(tau2.1*sum(J[(II==j)&(Z==1)])+sig.j[j]*xi)/denom,sqrt(sig.j[j]*tau2.1/denom)))
}

#--------------------------------------------------------------------------
samp.sig.j<-function(x){ #c(II,Z,J,mu.j,m)
  j     <-x
  return(rinvgamma(1,length(J[(II==j)&(Z==1)])/2+1,sum((J[(II==j)&(Z==1)]-mu.j[j])^2)/2+beta)) #두번째가 되는지 확인할 것
}

#after check simulation data and save it
Yobs<-tr.Yobs

#--------------------------------------------------------------------------
#prior distributions for Hyper-Parameters
#--------------------------------------------------------------------------
                    #whose / Prior dist / hyper parameter
a0        <- 0      #mu         normal           mean
b0        <- 10     #mu         normal           scale

a.sig2    <- 0.001  #varinace  Inv-gamma         1st
b.sig2    <- 0.001  #variance  Inv-gamma         2nd

#scaled inv-chi 
#nu	       <- 10  #sigma      scaled inv-chi    d.f. 
#tau2.0    <- 0.1   #sigma      scaled inv-chi    scale

a.gamma   <- 1      #gamma      gamma             1st
b.gamma   <- 1      #gamma      gamma             2nd

a.lambda  <- 1      #lambda     beta              1st
b.lambda  <- 1      #lambda     beta              2nd

a.xi      <- 0      #m of mu.j  normal            mean
b.xi      <- 10     #m of mu.j  normal            scale
a.tau2.1  <- 0.001  #v of mu.j  Inv-gamma         1st
b.tau2.1  <- 0.001  #v of mu.j  Inv-gamma         2nd

a.beta    <- 0      #sigma.j    gamma             1st
b.beta    <- 0      #sigma.j    gamma             2nd

#--------------------------------------------------------------------------
#Settings
#--------------------------------------------------------------------------
n.iters       <- 32000
m             <- 20#simulation:10 realdata:20                            
k             <- 1                            
 
#Save the result
burn           <- n.iters/2                       
gibbs.theta1   <- matrix(0,nrow=burn, ncol=9+m*3) 
gibbs.theta2  <- matrix(0,nrow=burn, ncol=9+m*3)

set.seed(2009121282)
#====================================================
# Gibbs Sampler
#====================================================
for(start in 1:2){
start<-1
#----------------------------------------------------
# Setting Starting Value
#----------------------------------------------------
  if(start==1){ #Starting Values 1
    set.seed(2009121282)
    mu      <- 0.001
    sig     <- 2
    lamb    <- 0.025
    Z       <- rbinom(t,1,lamb)
    
    ##dirichlet
    gam    <- 2
    xi     <- 0.1
    tau2.1 <- 1.5
    beta   <- 1
    
    J      <- rep(0.1,t)
    mu.j   <- rnorm(m,xi,sqrt(tau2.1))
    sig.j  <- rinvgamma(m,1,beta)
    V.j    <- c(rbeta(m-1,1,gam),1)
    II     <- sample(1:m,t,replace=T)   
    I      <- matrix(0,t,m)             
    for(i in 1:t) I[i,II[i]]<-1         
  }#start==1

  if(start==2){ #Starting Values 2
    mu     <- 0.0011#0.002
    sig    <- 2.0001
    lamb   <- 0.0251
    Z      <- rbinom(t,1,lamb)
    
    ##dirichlet
    gam    <- 2.0001
    xi     <- 0.1001
    tau2.1 <- 1.5001
    beta   <- 1.0001
    
    J      <- rep(0.1,t)
    mu.j   <- rnorm(m,xi,sqrt(tau2.1))
    sig.j  <- rinvgamma(m,1,beta)
    V.j    <- c(rbeta(m-1,1,gam),1)
    II     <- sample(1:m,t,replace=T)   # jump group 
    I      <- matrix(0,t,m)             # group matrix
    for(i in 1:t) I[i,II[i]]<-1  
    k             <- 1        
  }# start==2

#----------------------------------------------------
# Iteration
#----------------------------------------------------
  for(iter in 1:n.iters){ 
    
    if(iter%%100==0){cat("Iteration",iter,"\n")}
    J_hat   <-(sig.j[II]*(Yobs-mu)+sig*mu.j[II])/(sig.j[II]+sig)

    # Sample Z
    tmp1	<-dnorm(Yobs,mu+J_hat,sqrt(sig),log=T)+log(lamb)+dnorm(J_hat,mu.j[II],sqrt(sig.j[II]),log=T)+1/2*log(sig*sig.j[II]/(sig.j[II]+sig))
    tmp2	<-dnorm(Yobs,mu,sqrt(sig),log=T)+log(1-lamb)
    
    tmp1	<-exp(tmp1-max(tmp1,tmp2))
    tmp2	<-exp(tmp2-max(tmp1,tmp2))
    
    p1    <-tmp1/(tmp1+tmp2)
    Z	    <-rbinom(t,1,p1)
    
    # Sample I
    A     <-cbind(J_hat,Z)
    II    <-apply(A,1,samp.I)
    I     <-vtom(II)
    
    III<-Z*II
    IIII<-vtom(III)
    G<-apply(IIII,2,sum)
    
    if(iter%%100==1) cat("- # of clusters = [",length(G[G!=0]),"]:",c(1:m)[G!=0],"\n")
    
    M<-length(G[G!=0])
    
    # Sample J
    E     <-cbind(Yobs,II)
    J     <-apply(E,1,samp.J)

    # Sample V
    B     <-apply(I,2,sum)
    H     <-cbind(1:(m-1))
    V.j   <-c(apply(H,1,samp.V),1)
    
    # V_J 
    while(length(V.j[V.j==1])>1){
      idx   <-cbind((1:m)[V.j==1][-length(V.j[V.j==1])])
      V.j[idx]<-apply(idx,1,samp.V)
    }
    
    # Sample gamma
    gam   <-rgamma(1,a.gamma+m-1,rate=(b.gamma-sum(log(1-V.j[1:(m-1)]))))
    
    # Sample mu.j
    C     <-as.matrix(c(1:m))
    mu.j  <-apply(C,1,samp.mu.j)
    
    # Sample xi
    xi    <-rnorm(1,(a.xi+b.xi*sum(mu.j))/(b.xi*m+1),sqrt((b.xi*tau2.1/(b.xi*m+1))))
    
    # Sample tau2.1
    tau2.1<-rinvgamma(1,a.tau2.1+(m+1)/2,b.tau2.1+sum((mu.j-xi)^2)/2+(xi-a.xi)^2/(2*b.xi))
    
    # Sample sig.j
    D     <-as.matrix(c(1:m))
    sig.j <-apply(D,1,samp.sig.j)
    
    # Sample beta
    beta  <-rgamma(1,a.beta+m,rate=(sum(1/sig.j)+b.beta))
    
    # Sample mu
    S1	  <-sum(Yobs-(Z*J))
    mu	  <-rnorm(1,(b0*S1+a0)/(b0*t+1),sqrt(b0*sig/(b0*t+1)))
    
    # Sample sig
    S2	  <-sum((Yobs-mu-(Z*J))^2) #vector
    sig   <-rinvgamma(1,a.sig2+(t+1)/2,b.sig2+S2/2+(mu-a0)^2/(2*b0))
    
    # Sample lambda
    lamb	<-rbeta(1,sum(Z)+a.lambda,t-sum(Z)+b.lambda)
    
    # B: jump 횟수 확인
    B	<-sum(Z)
    
    # Save
    if(iter>burn){
      if(start==1){gibbs.theta1[k,]<-c(mu,sig,lamb,gam,xi,tau2.1,beta,M,B,exp(PI(V.j,m)),mu.j,sig.j)}
      if(start==2){gibbs.theta2[k,]<-c(mu,sig,lamb,gam,xi,tau2.1,beta,M,B,exp(PI(V.j,m)),mu.j,sig.j)}
      k <- k+1
    }# iter>burn
  }# iteration in 1:n.iters
}# start in 1:2

#---------------------------------------------------------
# Gibbs sumup function
#---------------------------------------------------------
Gibbs.Sumup <- function(theta,half=F){
  if(half) theta <- theta[,(length(theta[1,,1])/2+1):length(theta[1,,1]),]
  J <- length(theta[,1,1])  # Number of Chains
  n <- length(theta[1,,1])  # Number of Iterations
  p <- length(theta[1,1,])  # Number of Paramters
  output <- matrix(0,p,6)
  stats <- array(0,c(J,p,2))
  for(i in 1:J){
    stats[i,,1] <- apply(theta[i,,],2,mean)
    stats[i,,2] <- apply(theta[i,,],2,var)
  }
  Between <- apply(stats[,,1],2,var)*n
  Within <- apply(stats[,,2],2,mean)
  for(i in 1:p){
    output[i,1] <- mean(theta[,,i])
    output[i,2] <- sqrt(var(as.vector(theta[,,i])))
    output[i,3] <- sort(as.vector(theta[,,i]))[as.integer(n*J*.025)]
    output[i,4] <- median(theta[,,i])
    output[i,5] <- sort(as.vector(theta[,,i]))[as.integer(n*J*.975)+1]
    output[i,6] <- sqrt(((n-1)*Within/n + Between/n)/Within)[i]
  }
  dimnames(output) <- list(NULL,c("mean","sd","2.5%","median","97.5%","Rhat"))
  rownames(output) <- colnames(theta[1,,])
  round(output,6)
}
dim(gibbs.theta1)
head(gibbs.theta1)

#=========================================================
gibbs1 <- gibbs.theta1[,1:3]  # gibbs1.fst <- gibbs1 # gibbs1 <- gibbs1.fst
gibbs2 <- gibbs.theta2[,1:3] # gibbs2.fst <- gibbs2 # gibbs2 <- gibbs2.fst
#---------------------------------------------------------

pppp<-3
theta.gib <- array(0, c(2,n.iters/2,pppp))
theta.gib[1,,] <- gibbs1; theta.gib[2,,] <- gibbs2

dim(theta.gib)
dimnames(theta.gib)=list(1:2,1:n.iters,c("mu","sig2","lamb"))
Gibbs.Sumup(theta.gib,half=FALSE)


