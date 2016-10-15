###########
# sample from CTMC
sampleCTMC=function(pi1,k,t){
  
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
  
  while( ts[i]<t ){
  k=ks[state[i]]
  t.temp=rexp(n = 1, rate = k)
  ts=c(ts,ts[i]+t.temp)
  state=c(state,ifelse(state[i]==1,2,1))
  i=i+1
  }
  
  return(data.frame(ts,state))
}


############
# sample observation time ti and delay time tau

SampleM1=function(a,b,pi1,k,A0,t){
  
  # return of data frame for ti and tay
  
  # sample the state space space
  twostates = sampleCTMC(pi1,k,t)
  n=dim(twostates)[1]
  
  gamm = (twostates$state==1)*(a)+b *(twostates$state==2)
  rate1 = A0/gamm
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
    id=max(which(tt>twostates$ts))
    lt=rate1[id]
    if(runif(1)<lt/L){ 
    ts=c(ts,tt)
    tau.temp=rexp(1,rate=gamm[id])
    taus=c(taus,tau.temp)
    i=i+1
    }
  }
  return(data.frame(ts,taus))
}





