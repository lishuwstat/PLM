## code for one covariate with non-linear effect


#fun==1 Simulation setup I (\Lambda(t)=t/10)
#fun==2 Simulation setup II (\Lambda(t)=log(1+3t/2))

## n denotes sample size

## r is the parameter in the transformation function G

## gt: the length of two consecutive observations (\mu) follows the exponential distribution with mean 1/gt

## tau is the length of study

gen.data<-function(n,fun,r,gt,tau){
  
  if(fun==1){
    base<-function(t){
      res<-t/10
      return(res)
    }}


  if(fun==2){
    base<-function(t){
      res<-log(1+3/2*t)
      return(res)
    }}
  
  
  if(fun==1){
    phi<-function(w){
      res<-(w^3+1)
      return(res)
    }}
  
  if(fun==2){
    phi<-function(w){
      res<-w^2-1
      return(res)
    }}

  
  # gen.t<-function(tt,lamb){
  #   dif2<-(base(tt)-lamb)^2
  #   return(dif2)
  # }
  
  b<-matrix(c(0.5,-0.5))
  X2<-rbinom(n,1,0.5)
  X1<-rnorm(n,0,1)
  X<-cbind(X1,X2)
  W<-runif(n,-1,1)
  xb<-X%*%b+phi(W)   
  exb<-exp(xb)
  U<-runif(n,0,1)
  Lambt<-G.inv(-log(U),r)*exb^(-1)  
 
  T<-numeric()
  # for(i in 1:n){
  #   T[i]<-optimize(gen.t,interval=c(0,1e+8),lamb=Lambt[i])$minimum
  # }
  if(fun==1){T<-10*Lambt}
  if(fun==2){T<-2/3*(exp(Lambt)-1)}
  
 
  L = R = c()
  
  d1 = d2 = d3 = rep(0,n)
  
  i = 1

  for(i in 1:n){
    
    Tn = rexp(1000, gt)      # gap times
    
    Sn = cumsum(Tn)          # observed time
    
    N = min(which(Sn>tau))-1 # number of observations
    
    Cn =c(0,  Sn[1:N]) 
    
    L.ind =  max(which(T[i] >= Cn))
    
    if(T[i] <= Cn[2] ){  #left censor
      d1[i] = 1
      L[i] = 0
      R[i] = Cn[2]
    }
    
    if (Cn[2] < T[i] &&  T[i] < max(Cn) ){  #interval censor
      
      L[i] = Cn[L.ind ]
      R[i] = Cn[L.ind +1]
      d2[i] = 1
    } 
    
    if(T[i] >= max(Cn)){ #right censor
      
      L[i] = max(Cn) 
      R[i] = Inf
      d3[i] = 1
    }
  }
  
  L.prob = mean(d1)
  I.prob = mean(d2)
  R.prob = mean(d3)
    
  L[d1==1]<-NA
  R[d3==1]<-NA
  sdata<-data.frame(Li=L,Ri=R,d1=d1,d2=d2,d3=d3,X1=X1,X2=X2,Ti=T,W=W)
  return(sdata)
}


## Li and Ri correspond to the left and right endpoints of the interval that brackets the failure time, respectively.

##d1, d2 and d3 correspond to the left, interval and right indicators, respectively.






