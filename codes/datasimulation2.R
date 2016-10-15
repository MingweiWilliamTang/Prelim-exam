sampleCTMC_BM=function(pi1,k,t,w,s){
  
  # t total length of the time
  # return of dataframe has three colomns:
  # the 1st col is the jumping time, the 2nd col is the states 
  
  # construct the transition rate 
  k12 = k*pi1
  k21 = k*(1-pi1)
  ks = c(k12,k21)
  
  # generate the initial states
  state = numeric(1) 
  state[1] = sample(x = c(1,2),size = 1,prob = c(pi1,1-pi1))
  ts=numeric(1)
  i=1
  x=1
  y=1
  z=1
  
  alpha = exp(-1/2*((x^2+y^2)/w[1]^2+(z^2)/w[2]^2))
  
  while( ts[i]<t ){
    k=ks[state[i]]
    t.temp=rexp(n = 1, rate = k)
    x.temp=rnorm(1,x[i],s*sqrt(t.temp))
    y.temp=rnorm(1,y[i],s*sqrt(t.temp))
    z.temp=rnorm(1,z[i],s*sqrt(t.temp))  
   ts=c(ts,ts[i]+t.temp)
    x=c(x,x.temp)
    y=c(y,y.temp)
    z=c(z,z.temp)
    alpha.temp=exp(-1/2*((x.temp^2+y.temp^2)/w[1]^2+z.temp^2/w[2]^2))
   alpha=c(alpha,alpha.temp)
    state=c(state,ifelse(state[i]==1,2,1))
    i=i+1
  }
  
  return(data.frame(ts,state,alpha))
}











SampleM2=function(a,b,pi1,k,A0,t,w,s){
  
  # return of data frame for ti and tay
  
  # sample the state space space
  twostates = sampleCTMC(pi1,k,t)
  n=dim(twostates)[1]
  
  # brownian motion
  x=0
  y=0
  z=0
  xs=x;ys=y;zs=z;
  alpha = exp(-1/2*((x^2+y^2)/w[1]^2+(z^2)/w[2]^2))
  
  
  
  gamm = (twostates$state==1)*(a)+b *(twostates$state==2)
  rate = A0/gamm
  # sample from homogenous possion process first
  L=2*A0*(1/a+1/b)
  ts=numeric(1)
  tt=ts
  taus= rexp(1,rate=gamm[1])
  i=1
  print("done")
  while( tt<t ){
    t.temp=rexp(n = 1, rate = L)
    tt=tt+t.temp
    
    x.temp=rnorm(1,x,s*sqrt(t.temp))
    y.temp=rnorm(1,y,s*sqrt(t.temp))
    z.temp=rnorm(1,z,s*sqrt(t.temp))
      
    alpha.temp=exp(-1/2*((x.temp^2+y.temp^2)/w[1]^2+z.temp^2/w[2]^2))
      
    id=max(which(tt>twostates$ts))
    lt=rate[id]*alpha.temp
    
    x=x.temp
    y=y.temp
    z=z.temp
    
    if(runif(1)<lt/L){ 
      ts=c(ts,tt)
      tau.temp=rexp(1,rate=gamm[id])
      taus=c(taus,tau.temp)
      alpha=c(alpha,alpha.temp)
      xs=c(xs,x.temp)
      ys=c(ys,y.temp)
      zs=c(zs,z.temp)
      i=i+1
    }
  }
  return(data.frame(ts,taus,alpha,xs,ys,zs))
}

#SampleM2(1.5,0.25,0.2,3000,37000,0.01,c(1000,500),10000)
