################# necessary functions ################

G.fun <- function(u,r){
  if(r > 0){
    results = log(((1+r*u)<=0)*10^(-3)+(1+r*u>0)*(1+r*u))/r
  } else results = u
  return(results)
}

G.inv <- function(u,r){
  if (r > 0) {
    results = (exp(u*r)-1)/r
  } else results = u
  return(results)
}

##############################################################

fphi <- function(mu,r){
  
  results = (1/r)^(1/r)*mu^(1/r-1)*exp(-1/r*mu)/gamma(1/r)
  
  return(results)
}

#####################################################

Partial.Linear.PG<-function(sdata, n.phi1,n.phi2, n.lam, order, degree, t.seq, w1.seq,w2.seq,max.iter,cov.rate){
 
  X<-cbind(sdata$X1,sdata$X2)
  d<-dim(X)[2]
  
  W1=sdata$W1
  W2=sdata$W2
  g0<-rep(0.01,order+n.lam)
  b0<-rep(0,d)
  a10<-rep(0,degree+n.phi1)
  a20<-rep(0,degree+n.phi2)
  
  #B-spline
  id1 <- seq(0, 1, length.out = (n.phi1 + 2))
  id1 <- id1[-c(1, (n.phi1 + 2))]
  id2 <- seq(0, 1, length.out = (n.phi2 + 2))
  id2 <- id2[-c(1, (n.phi2 + 2))]
  w1.max <- max(W1) + 1e-05
  w1.min <- min(W1) - 1e-05
  knots1 <- quantile(W1, id1)
  w2.max <- max(W2) + 1e-05
  w2.min <- min(W2) - 1e-05
  knots2 <- quantile(W2, id2)
  Bw1<-bs(x=W1, knots = knots1, degree = degree, Boundary.knots=c(w1.min,w1.max))
  Bw2<-bs(x=W2, knots = knots2, degree = degree, Boundary.knots=c(w2.min,w2.max))
  Bws1<-bs(x=w1.seq, knots = knots1, degree = degree, Boundary.knots=c(w1.min,w1.max))
  Bws2<-bs(x=w2.seq, knots = knots2, degree = degree, Boundary.knots=c(w2.min,w2.max))
  
  fit<-partial.linear.PG.EM(g0,b0,a10,a20,sdata=sdata,X,Bw1,Bw2,r,d,n.phi1,n.phi2,n.lam,order=order, t.seq=t.seq,max.iter,cov.rate,equal = FALSE)
  
  fit$psi1<-Bws1%*%fit$a1
  fit$psi2<-Bws2%*%fit$a2
  
  return(fit)
  
}

##########################################################

Q1<-function(bb,aa1,aa2,X,Bw1,Bw2,Ezil,Eyil,Efi,Ml.Li,Ml.Ri,sdata){ 
  N<-nrow(sdata)
  L<-ncol(Ml.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  exb<-exp(X%*%bb+Bw1%*%aa1+Bw2%*%aa2)
  numer<-Ezil+r.mat(1-sdata$d1)*Eyil 
  denom<-(r.mat(1-sdata$d3)*t(Ml.Ri)+r.mat(sdata$d3)*t(Ml.Li))*r.mat(exb*Efi)   
  g1<- apply(numer,1,sum)/apply(denom,1,sum)
  return(g1)
}

##########################################################

Q2<-function(bb,aa1,aa2,X,Bw1,Bw2,Ezil,Eyil,Efi,Ml.Li,Ml.Ri,sdata){  
  N<-nrow(sdata)
  L<-ncol(Ml.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  exb<-exp(X%*%bb+Bw1%*%aa1+Bw2%*%aa2)
  p1<-sum((Ezil+r.mat(sdata$d2+sdata$d3)*Eyil)*r.mat(X%*%bb+Bw1%*%aa1+Bw2%*%aa2))
  g1<- Q1(bb,aa1,aa2,X,Bw1,Bw2,Ezil,Eyil,Efi,Ml.Li,Ml.Ri,sdata)
  g1[g1==0]=0.01
  p2<-sum(Ezil*c.mat(log(g1))+r.mat(sdata$d2+sdata$d3)*Eyil*c.mat(log(g1)))
  p3<-sum(c.mat(g1)*r.mat(exb*Efi)*(r.mat(1-sdata$d3)*t(Ml.Ri)+r.mat(sdata$d3)*t(Ml.Li)))
  Q<-(p1+p2-p3)
  return(-Q)
}

##########################################################

loglik <-
  function(b0,a10,a20,g0,sdata,X,Bw1,Bw2,r,Ml.Li,Ml.Ri){
    GRi<-Ml.Ri%*%g0
    GLi<-Ml.Li%*%g0
    xb<-X%*%b0+Bw1%*%a10+Bw2%*%a20
    l1<-(1-exp(-G.fun(GRi*exp(xb),r)))^sdata$d1
    l1[which(l1==0)] <-10^-3
    l2<-((1-exp(-G.fun(GRi*exp(xb),r)))-(1-exp(-G.fun(GLi*exp(xb),r))))^sdata$d2
    l2[which(l2<=0)] <-10^-3
    l3<-(1-(1-exp(-G.fun(GLi*exp(xb),r))))^sdata$d3
    lpart1<-log(l1)
    lpart2<-log(l2)
    lpart3<-log(l3)
    
    return(sum(lpart1+lpart2+lpart3))
  }

###################################################
# The EM algorithm

partial.linear.PG.EM <-
  function(g0,b0,a10,a20,sdata,X,Bw1,Bw2,r,d,n.phi1,n.phi2,n.lam,order,t.seq,max.iter,cov.rate,equal = FALSE){

    B <- length(b0)
    A1 <- length(a10)
    A2 <- length(a20)
    G <- length(g0)
    N <- length(sdata$d1)
    
    r.mat<-function(x) matrix(x,ncol=N,nrow=G,byrow=TRUE)
    c.mat<-function(x) matrix(x,ncol=N,nrow=G)
    
    sdata$Li[sdata$d1 == 1] <- sdata$Ri[sdata$d1 == 1]  
    sdata$Ri[sdata$d3 == 1] <- sdata$Li[sdata$d3 == 1]  
    ti <- c(sdata$Li[sdata$d1 == 0], sdata$Ri[sdata$d3 == 0])
    
    if (equal == TRUE) {         
      ti.max <- max(ti) + .00001
      ti.min <- min(ti) - .00001
      knots <- seq(ti.min, ti.max, length.out = (n.lam + 2))
    }
    if (equal == FALSE) {       #
      id <- seq(0, 1, length.out = (n.lam + 2))
      id <- id[-c(1, (n.lam + 2))]
      ti.max <- max(ti) + .00001
      ti.min <- min(ti) - .00001
      knots <- c(ti.min, quantile(ti, id), ti.max)
    }
    
    
    Ml.Li<-Ispline(sdata$Li,order=order,knots=knots) 
    Ml.Ri<-Ispline(sdata$Ri,order=order,knots=knots) 
    Mt <- t(Ispline(t.seq,order=order,knots=knots) )
    Ml.Li=t(Ml.Li)
    Ml.Ri=t(Ml.Ri)
    
    ###############################################
    # Lets start iterating
    ll<-numeric()
    dd<-1
    ii<-0
   
    while(dd>cov.rate & ii<max.iter){
    ii<-ii+1
    
    # E-step
    exb0<-exp(X%*%b0+Bw1%*%a10+Bw2%*%a20)
    Lamd.Li<-Ml.Li%*%g0
    Lamd.Ri<-Ml.Ri%*%g0
    dz<- 1-exp(-G.fun(Lamd.Ri*exb0,r))
    dz[which(dz==0)] <- 1
    
    
    EZi<-sdata$d1*exb0*Lamd.Ri/dz  
    Ezil<-(r.mat(sdata$d1)*t(Ml.Ri))*c.mat(g0)*r.mat(exb0)/r.mat(dz)

    dy2<- exp(-G.fun(Lamd.Li*exb0,r))-exp(-G.fun(Lamd.Ri*exb0,r))
    dy2[which(dy2==0)] <- 1
    dy1<-(Lamd.Ri-Lamd.Li)*exb0
    dy1[which(dy1==0)] <- 1
    
    
    if(r>0){
    EYi = matrix(1,nrow = N,ncol=1)
    cc  <-  gaussLaguerre(17) 
    for(i in 1:N)
    {
      EYi[i] <- sum((sdata$d2[i]*cc$x*fphi(cc$x,r)*(Lamd.Ri[i]-Lamd.Li[i])*exb0[i]*(1-exp(-dy1[i]*cc$x))^(-1)*(exp(-Lamd.Li[i]*exb0[i]*cc$x)-exp(-Lamd.Ri[i]*exb0[i]*cc$x))/dy2[i]*exp(cc$x))*cc$w)
    }
    EYi2 = matrix(1,nrow = N,ncol=1)
    for(i in 1:N)
    {
      EYi2[i] <- sum((sdata$d2[i]*cc$x*fphi(cc$x,r)*exb0[i]*(1-exp(-dy1[i]*cc$x))^(-1)*(exp(-Lamd.Li[i]*exb0[i]*cc$x)-exp(-Lamd.Ri[i]*exb0[i]*cc$x))/dy2[i]*exp(cc$x))*cc$w)
    }
    Eyil <-t(Ml.Ri-Ml.Li)*c.mat(g0)*r.mat(EYi2)
    }
   

    if(r==0){
      dy<-1-exp(-(Lamd.Ri-Lamd.Li)*exb0)
      dy[which(dy==0)] <- 1
      EYi<- sdata$d2*exb0*(Lamd.Ri-Lamd.Li)/dy
      Eyil <-(r.mat(sdata$d2)*t(Ml.Ri-Ml.Li)*c.mat(g0)*r.mat(exb0))/r.mat(dy)
    }

    
    pc1<-sdata$d1*(1-(1+r*Lamd.Ri*exb0)^(-1/r-1))/dz    
    pc2<-sdata$d2*((1+r*Lamd.Li*exb0)^(-1/r-1)-(1+r*Lamd.Ri*exb0)^(-1/r-1))/dy2
    pc3<-sdata$d3*(1+r*Lamd.Li*exb0)^(-1) 
    Efi<-(r>0)*(pc1+pc2+pc3)+(r==0)*rep(1,n)
    
    
    fit_b<-nlm(Q2,p=b0,aa1=a10,aa2=a20,X=X,Bw1=Bw1,Bw2=Bw2,Ezil=Ezil,Eyil=Eyil,Efi=Efi,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,sdata=sdata,iterlim = 1,print.level=0)
    fit_a1<-nlm(Q2,p=a10,bb=b0,aa2=a20,X=X,Bw1=Bw1,Bw2=Bw2,Ezil=Ezil,Eyil=Eyil,Efi=Efi,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,sdata=sdata,iterlim = 1,print.level=0)
    fit_a2<-nlm(Q2,p=a20,bb=b0,aa1=a10,X=X,Bw1=Bw1,Bw2=Bw2,Ezil=Ezil,Eyil=Eyil,Efi=Efi,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,sdata=sdata,iterlim = 1,print.level=0)
    
    
    if((!is.character(fit_b))&&(!is.character(fit_a1))){
      b0<-fit_b$estimate
      a10<-fit_a1$estimate
      a20<-fit_a2$estimate
      g0<-Q1(bb=b0,aa1=a10,aa2=a20,X=X,Bw1=Bw1,Bw2=Bw2,Ezil=Ezil,Eyil=Eyil,Efi=Efi,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,sdata=sdata)
      ll<-c(ll,loglik(b0,a10,a20,g0,sdata=sdata,X=X,Bw1=Bw1,Bw2=Bw2,r=r,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri))
    }
    if(ii>2) dd<-abs(ll[ii]-ll[ii-1])
    
    
    }
    loglik<-ll[ii]
    AIC<-2*(A1+A2+B+G)-2*loglik
    BIC<-(A1+A2+B+G)*log(N)-2*loglik
    
    Lam<-Mt%*%g0
    
    return(list("b"=b0, "g"=g0,"a1"=a10,"a2"=a20, "Lam"=Lam,"ll"=loglik, "AIC"=AIC, "BIC"=BIC, "loops"=ii))
    
  }


