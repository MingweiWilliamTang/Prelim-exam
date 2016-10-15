ag_likelihood0=function(datasets,alphas,a,b,pi1,k,A0)
{
  N=length(datasets)
  L=0
  for(i in 1:N){
    dt=datasets[[i]]$dt
    taus=datasets[[i]]$taus
    alpha=alphas[[i]]$alpha
#    cat(likelihood_2(dt,taus,alpha,a,b,pi1,k,A0),'\n')
    L=L+likelihood_2(dt,taus,alpha,a,b,pi1,k,A0)
  }
  return(L)
}

ag_likelihood=function(datasets,alphas,xts,a,b,pi1,k,A0)
{
 N=length(datasets)
 L=0
 for(i in 1:N){
   dt=datasets[[i]]$dt
   taus=datasets[[i]]$taus
   alpha=alphas[[i]]$alpha
   xt=xts[[i]]$xt
   L=L+likelihood_3(dt,taus,alpha,xt,a,b,pi1,k,A0)
 }
  return(L)
}

agg_path=function(datasets,xts,lambda,kxi)
{
  N=length(datasets)
  L=0
  for(i in 1:N){
    dt=datasets[[i]]$dt
    xt=xts[[i]]$xt
    L=L+logPath(xt,dt,lambda,kxi)
  }
  return(L)
}

updateTheta=function(datasets,alphas,xts,cs,a0,b0,pi10,k0,A00,iter){
  L0=ag_likelihood(datasets,alphas,xts,a0,b0,pi10,k0,A00)
  for(i in 1:iter){
    
    ### update a,b
    
    temp1=rgamma(1,1/cs[1],scale=a0*cs[1]) 
    temp2=rgamma(1,1/cs[2],scale=b0*cs[2])
    a=max(temp1,temp2)
    b=min(temp1,temp2)
    
    L1=ag_likelihood(datasets,alphas,xts,a,b,pi10,k0,A00)
    
    r = L1+dgamma(a0,1/cs[1],scale=a*cs[1],log=T)+dgamma(b0,1/cs[2],scale=(b*cs[2]),log=T)+ 
      dgamma(a,1,scale=2,log=T)+dgamma(b,1.5625,scale =1/1.5625,log=T)-
      (L0+dgamma(a,1/cs[1],scale=a0*cs[1],log=T)+dgamma(b,1/cs[2],scale=(b0*cs[2]),log=T)+ 
         dgamma(a0,1,scale=2,log=T)+dgamma(b0,1.5625,scale =1/1.5625,log=T))
    r=exp(r)
    if(is.na(r)){cat(c(a,b,pi10,k0,A00,L1))}
    if(runif(1,0,1)<r) { 
      a0=a
      b0=b  
      L0 = L1
    }
    ### update pi1
    
    pi1 = rbeta(1,cs[3]*pi10,cs[3]*(1-pi10))
    if(pi1==0 || pi1==1){pi1=0.5}
    L1 = ag_likelihood(datasets,alphas,xts,a0,b0,pi1,k0,A00)
    
    r = (L1 + dbeta(pi10,cs[3]*pi1,cs[3]*(1-pi1),log =T) + dbeta(pi1,0.89,0.89,log=T))-
    (L0 + dbeta(pi1,cs[3]*pi10,cs[3]*(1-pi10),log=T) + dbeta(pi10,0.89,0.89,log=T))

    r=exp(r)
    if(is.na(r)){cat(c(a0,b0,pi1,k0,A00,L1))}
    if(runif(1,0,1)<r) { 
      pi10 = pi1
       L0 = L1
    }
    
    # update k
    k = rgamma(1,1/cs[4],scale=(k0*cs[4]))
    
    L1 = ag_likelihood(datasets,alphas,xts,a0,b0,pi10,k,A00)
    
    r = (L1 + dgamma(k0,1/cs[4],scale=(k*cs[4]),log=T) + dexp(k,1/40000.0,log=T) )-
      (L0 + dgamma(k,1/cs[4],scale=(k0*cs[4]),log=T) + dexp(k0,1/40000.0,log=T) )
    
    r=exp(r)
    if(is.na(r)){cat(c(a0,b0,pi10,k,A00,L1))}
    if(runif(1,0,1)<r){ 
      k0 = k
      L0 = L1
    }
    
    A0 = rgamma(1,1/cs[5],scale =(A00*cs[5]))
    L1 = ag_likelihood(datasets,alphas,xts,a0,b0,pi10,k0,A0)
    
    r =( L1 + dgamma(A00,1/cs[5],scale=(A0*cs[5]),log=T) + 
      dgamma(A0,1.96,scale =100000/5.6,log=T) ) -
      (L0 + dgamma(A0,1/cs[5],scale = (A00*cs[5]),log=T) + 
         dgamma(A00,1.96,scale =100000/5.6,log=T) )
    
    r=exp(r)
    if(is.na(r)){cat(c(a0,b0,pi10,k0,A0,L1))}
    if(runif(1,0,1)<r){ 
      A00 = A0
      L0 = L1
    }
  }
  return(c(a0,b0,pi10,k0,A00,L0))
}


agg_log_path=function(datasets,xts,lambda,kxi)
{
  N=length(datasets)
  p=0
  for(i in 1:N)
  {
    dt=datasets[[i]]$dt
    xt=xts[[i]]$xt
    p=p+logPath(xt,dt,lambda,kxi)
  }
  return(p)
}

agg_OU_para=function(datasets,xts,lambda0,kxi0,c1,c2,iter)
{
  L0 = agg_log_path(datasets,xts ,lambda0, kxi0)
  for(i in 1:iter){
  lambda = rgamma(1,1/c1,scale=c1* lambda0)
  
  L1 = agg_log_path(datasets,xts ,lambda, kxi0)
  
  r = L1 + dgamma(lambda0,1/c1,scale=c1*lambda,log=T) + dgamma(lambda,40,scale=5,log=T)-
  L0- dgamma(lambda, 1/c1,scale=c1*lambda0,log=T)-dgamma(lambda0,40,scale=5,log=T)
  if(is.na(r)){
    cat(c(lambda,kxi0))
               r=-100
               lambda0=1
  }
  if(runif(1,0,1)< exp(r))
  {
    lambda0=lambda
    L0=L1
  }
  kxi = rgamma(1,1/c2,scale =c2*kxi0);
  L1 = agg_log_path(datasets,xts,lambda0,kxi)
  r = L1 + dgamma(kxi0,1/c2,scale= c2*kxi,log=T) + dgamma(kxi,1,scale=2,log=T)-
  L0- dgamma(kxi, 1/c2,scale= c2*kxi0,log=T)-dgamma(kxi0,1,scale=2,log=T)
  if(is.na(r)){cat(c(lambda0,kxi))
               r=-100
               kxi0=1}
  if(runif(1,0,1)<exp(r))
  {
    kxi0=kxi
    L0=L1
  }
  }
  return(c(lambda0,kxi0))
}


agg_Br=function(datasets,alphas,Bs,xts,a,b,pi1,k,A0,c1,iter=1,maxlik=0)
{
  N=length(datasets)
  for(i in 1:N){
    dt=datasets[[i]]$dt
    taus=datasets[[i]]$taus
    alpha0=alphas[[i]]$alpha
    xt=xts[[i]]$xt
    Bx=Bs[[i]]$Bx
    By=Bs[[i]]$By
    Bz=Bs[[i]]$Bz
    res=updateBrowian(dt,taus,alpha0,Bx,By,Bz,xt,a,b,pi1,k,A0,c1,iter,maxlik)
    Bs[[i]]$Bx=res[1,]
    Bs[[i]]$By=res[2,]
    Bs[[i]]$Bz=res[3,]
    alphas[[i]]$alpha=res[4,]
  }
  return(list(r1=Bs,r2=alphas))
}

agg_OU_xt=function(datasets,alphas,xts,lambda,kxi,a,b,pi1,k,A0,c1,iter=1)
{
  N=length(datasets)
  xt2=xts
  for(i in 1:N){
    dt=datasets[[i]]$dt
    taus=datasets[[i]]$taus
    alpha=alphas[[i]]$alpha
    xt0=xts[[i]]$xt
    xt2[[i]]$xt=updatext(dt,taus,alpha,xt0,lambda,kxi,a,b,pi1,k,A0,c1,iter)
  }
return(xt2)
}


###
agg_cont_diff=function(datasets,cs,iter,burnin,thin)
{
 ### initialize
  N=length(datasets)
  cat('number of datasets N =',N,'\n')
  w=c(310,1760)
  sigma=1000;
  ag_res=matrix(rep(0,iter*8),ncol=8)
  
 # a0 = rgamma(1,1.0,scale=2.0);
#  b0 = rgamma(1,1.5625,scale=1.0/1.5625);
#  k0 = rexp(1,1/40000.0);
 # pi10 = rbeta(1,0.89,0.89);
  #A00 = rgamma(1,1.96,scale=100000/5.6);
a0 = 1.5#rgamma(1,1.0,scale=2.0);
b0 = 0.25#rgamma(1,1.5625,scale=1.0/1.5625);
k0 = 3000#rexp(1,1/40000.0);
pi10 = 0.3#rbeta(1,0.89,0.89);
A00 = 25000#rgamma(1,1.96,scale=100000/5.6);  
a=max(a0,b0)
  b=min(a0,b0)
  a0=a
  b0=b
  pi1=pi10
  k=k0
  A0=A00
  lambda0 = rgamma(1,100 , scale=5.0);
  kxi0 = rgamma(1,1.0 , scale=2.0);
  kxi=0
  lambda=lambda0
  
  Bs=list()
  alphas=list()
  xts=list()
  for(j in 1:N){
    dt=datasets[[j]]$dt
    taus=datasets[[j]]$taus
     n = length(taus)
    {  
    #Bx=numeric(n)
    #xt=alpha=By=Bz=Bx
  #  Bx[1] = By[1] = Bz[1] = 0
    
   # alpha[1] = exp(-( (Bx[1]*Bx[1] +By[1]*By[1] )/ (2*w[1]*w[1])) - 
    #               ((Bz[1]*Bz[1])/(2*w[2]*w[2])))
    #for(i in 2:n)
    #{
    #  Bx[i]=rnorm(1,Bx[i-1],sigma*sqrt(dt[i-1]));
    #  By[i]=rnorm(1,By[i-1],sigma*sqrt(dt[i-1]));
    #  Bz[i]=rnorm(1,Bz[i-1],sigma*sqrt(dt[i-1]));
    #  alpha[i]= exp(-( (Bx[i]*Bx[i] +By[i]*By[i] )/ (2*w[1]*w[1])) - 
    #                ((Bz[i]*Bz[i])/(2*w[2]*w[2])));
    #}
     }
    alpha=datasets[[j]]$alpha
    Bx=datasets[[j]]$Bx
    By=datasets[[j]]$By
    Bz=datasets[[j]]$Bz
      # xt
    # xt[1] = rnorm(1,0 , sqrt(kxi0));
  xt=numeric(n)  
#  xt[1] = 0
 #   for (i in 2:n)
#    {
 #     xt[i] = rnorm(1,xt[i-1] * exp(-lambda * dt[i-1]),sqrt(kxi *
  #                                                           (1 - exp(-2*lambda * dt[i-1] ))))
  #  }
    xts[[j]]=list(xt=xt)
    Bs[[j]]=list(Bx=Bx,By=By,Bz=Bz)
    alphas[[j]]=list(alpha=alpha)
  }
    for(i in 1:(burnin+iter))
    {
      if(i%%100==0){cat("finish ",i," simulations","\n")}
      th=ifelse(i>burnin,thin,1)
      for(j in 1:th){
        res1=updateTheta(datasets,alphas,xts,cs[1:5],a0,b0,pi10,k0,A00,iter=10)
        a=a0=res1[1]
        b=b0=res1[2]
        pi1=pi10=res1[3]
        k=k0=res1[4]
        A0=A00=res1[5]
        
        res2=agg_OU_para(datasets,xts,lambda0,kxi0,cs[6],cs[7],iter=10)
        lambda=lambda0=res2[1]
        kxi=kxi0=res2[2]
        
        res3=agg_Br(datasets,alphas,Bs,xts,a,b,pi1,k,A0,cs[8],iter=1)
        Bs=res3$r1
        alphas=res3$r2
        
        xts=agg_OU_xt(datasets,alphas,xts,lambda,kxi,a,b,pi1,k,A0,cs[9],iter=1)
      }
      if(i>burnin)
      {
      BF=exp(ag_likelihood0(datasets,alphas,a,b,pi1,k,A0)-
               ag_likelihood(datasets,alphas,xts,a,b,pi1,k,A0))
      ag_res[i-burnin,]=c(res1[1:5],res2,BF)
      }
    }
  return(list(r1=ag_res,r2=xts,r3=alphas))
}


###############
agg_two_states=function(datasets,cs,iter,burnin,thin)
{
  ### initialize
  N=length(datasets)
  cat('number of datasets N =',N,'\n')
  w=c(310,1760)
  sigma=1000;
  ag_res=matrix(rep(0,iter*5),ncol=5)
  
  a0 = 1.5#rgamma(1,1.0,scale=2.0);
  b0 = 0.25#rgamma(1,1.5625,scale=1.0/1.5625);
  k0 = 10000#rexp(1,1/40000.0);
  pi10 = 0.4#rbeta(1,0.89,0.89);
  A00 = 35000#rgamma(1,1.96,scale=100000/5.6);
  a=max(a0,b0)
  b=min(a0,b0)
  a0=a
  b0=b
  pi1=pi10
  k=k0
  A0=A00
  
  Bs=list()
  alphas=list()
  xts=list()
  for(j in 1:N){
    dt=datasets[[j]]$dt
    taus=datasets[[j]]$taus
    n = length(taus)
{  
      Bx=numeric(n)
      xt=alpha=By=Bz=Bx
        Bx[1] = By[1] = Bz[1] = 0
      
       alpha[1] = exp(-( (Bx[1]*Bx[1] +By[1]*By[1] )/ (2*w[1]*w[1])) - 
                     ((Bz[1]*Bz[1])/(2*w[2]*w[2])))
      for(i in 2:n)
      {
        Bx[i]=rnorm(1,Bx[i-1],sigma*sqrt(dt[i-1]));
        By[i]=rnorm(1,By[i-1],sigma*sqrt(dt[i-1]));
        Bz[i]=rnorm(1,Bz[i-1],sigma*sqrt(dt[i-1]));
        alpha[i]= exp(-( (Bx[i]*Bx[i] +By[i]*By[i] )/ (2*w[1]*w[1])) - 
                      ((Bz[i]*Bz[i])/(2*w[2]*w[2])));
      }
    }
#alpha=datasets[[j]]$alpha
#Bx=datasets[[j]]$Bx
#By=datasets[[j]]$By
#Bz=datasets[[j]]$Bz
# xt
# xt[1] = rnorm(1,0 , sqrt(kxi0));
xt=numeric(n)
#  xt[1] = 0
#   for (i in 2:n)
#    {
#     xt[i] = rnorm(1,xt[i-1] * exp(-lambda * dt[i-1]),sqrt(kxi *
#                                                           (1 - exp(-2*lambda * dt[i-1] ))))
#  }
xts[[j]]=list(xt=xt)
Bs[[j]]=list(Bx=Bx,By=By,Bz=Bz)
alphas[[j]]=list(alpha=alpha)
  }
for(i in 1:(burnin+iter))
{
  if(i%%100==0){cat("finish ",i," simulations","\n")}
  th=ifelse(i>burnin,thin,1)
  for(j in 1:th){
    res1=updateTheta(datasets,alphas,xts,cs[1:5],a0,b0,pi10,k0,A00,iter=10)
    a=a0=res1[1]
    b=b0=res1[2]
    pi1=pi10=res1[3]
    k=k0=res1[4]
    A0=A00=res1[5]
    
    res3=agg_Br(datasets,alphas,Bs,xts,a,b,pi1,k,A0,cs[6],iter=1,1)
    Bs=res3$r1
    alphas=res3$r2
  }
  if(i>burnin)
  {
    ag_res[i-burnin,]=c(res1[1:5])
  }
}
return(list(r1=ag_res,r3=alphas))
}
##################
agg_cont_diff1=function(datasets,cs,iter,burnin,thin)
{
  ### initialize
  N=length(datasets)
  cat('number of datasets N =',N,'\n')
  w=c(310,1760)
  sigma=1000;
  ag_res=matrix(rep(0,iter*10),ncol=10)
  
  # a0 = rgamma(1,1.0,scale=2.0);
  #  b0 = rgamma(1,1.5625,scale=1.0/1.5625);
  #  k0 = rexp(1,1/40000.0);
  # pi10 = rbeta(1,0.89,0.89);
  #A00 = rgamma(1,1.96,scale=100000/5.6);
  a0 = 1.7#rgamma(1,1.0,scale=2.0);
  b0 = 0.28#rgamma(1,1.5625,scale=1.0/1.5625);
  k0 = 10000#rexp(1,1/40000.0);
  pi10 = 0.4#rbeta(1,0.89,0.89);
  A00 =25000 #rgamma(1,1.96,scale=100000/5.6);  
  a=max(a0,b0)
  b=min(a0,b0)
  a0=a
  b0=b
  pi1=pi10
  k=k0
  A0=A00
  lambda0 = 200;
  #rgamma(1,100 , scale=5.0);
  kxi0 = rgamma(1,1.0 , scale=2.0);
  kxi=5
  lambda=lambda0
  
  Bs=list()
  alphas=list()
  xts=list()
  for(j in 1:N){
    dt=datasets[[j]]$dt
    taus=datasets[[j]]$taus
    n = length(taus)
  
      Bx=numeric(n)
      xt=alpha=By=Bz=Bx
      alpha=rep(1,n)
        Bx[1] = By[1] = Bz[1] = 0
      
       alpha[1] = exp(-( (Bx[1]*Bx[1] +By[1]*By[1] )/ (2*w[1]*w[1])) - 
                     ((Bz[1]*Bz[1])/(2*w[2]*w[2])))
     for(i in 2:n)
      {
        Bx[i]=rnorm(1,Bx[i-1],sigma*sqrt(dt[i-1]));
        By[i]=rnorm(1,By[i-1],sigma*sqrt(dt[i-1]));
        Bz[i]=rnorm(1,Bz[i-1],sigma*sqrt(dt[i-1]));
        alpha[i]= exp(-( (Bx[i]*Bx[i] +By[i]*By[i] )/ (2*w[1]*w[1])) - 
                      ((Bz[i]*Bz[i])/(2*w[2]*w[2])));
      }
#alpha=datasets[[j]]$alpha
#Bx=datasets[[j]]$Bx
#By=datasets[[j]]$By
#Bz=datasets[[j]]$Bz
# xt
# xt[1] = rnorm(1,0 , sqrt(kxi0));
xt=numeric(n)  
#  xt[1] = rnorm(1,0,sqrt(kxi))

for (i in 2:n)
{
     xt[i] = rnorm(1,xt[i-1] * exp(-lambda * dt[i-1]),sqrt(kxi *
                                                           (1 - exp(-2*lambda * dt[i-1] ))))
  }
xts[[j]]=list(xt=xt/10)
Bs[[j]]=list(Bx=Bx,By=By,Bz=Bz)
alphas[[j]]=list(alpha=alpha)
  }
for(i in 1:(burnin+iter))
{
  if(i%%100==0){cat("finish ",i," simulations","\n")}
  th=ifelse(i>burnin,thin,1)
 # if(i==1500){
  #  a0=1.7
  #  b0=0.29
  #  pi10=0.38
  #  k0=10000
  #  A00=32000
  #  lambda0=100
  #}
  for(j in 1:th){
    res1=updateTheta(datasets,alphas,xts,cs[1:5],a0,b0,pi10,k0,A00,iter=10)
    a=a0=res1[1]
    b=b0=res1[2]
    pi1=pi10=res1[3]
    k=k0=res1[4]
    A0=A00=res1[5]
    
    res2=agg_OU_para(datasets,xts,lambda0,kxi0,cs[6],cs[7],iter=10)
    lambda=lambda0=res2[1]
    kxi=kxi0=res2[2]
    
    res3=agg_Br(datasets,alphas,Bs,xts,a,b,pi1,k,A0,cs[8],iter=1,1)
    Bs=res3$r1
   alphas=res3$r2
    
    xts=agg_OU_xt(datasets,alphas,xts,lambda,kxi,a,b,pi1,k,A0,cs[9],iter=1)
  }
  if(i>burnin)
  {
    l1=ag_likelihood0(datasets,alphas,a,b,pi1,k,A0)
    l2=ag_likelihood(datasets,alphas,xts,a,b,pi1,k,A0)
    BF=exp(ag_likelihood0(datasets,alphas,a,b,pi1,k,A0)-
             ag_likelihood(datasets,alphas,xts,a,b,pi1,k,A0))
    ag_res[i-burnin,]=c(res1[1:5],res2,l1,l2,BF)
  }
}
return(list(r1=ag_res,r2=xts,r3=alphas))
}

