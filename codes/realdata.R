# data cleaning
setwd("DNAhairpinData/")
d=read.table("trajectory11.dat")


datasets2=list()
for(i in 1:50){
  filename=paste0("trajectory",as.character(i),".dat")
  data=read.table(filename)
  dt=diff(data[,1])
  dt[dt==0]=10^(-8)
  taus=data[,2]
  taus[taus<=0]=13.2
  datasets2[[i]] <- list(dt=dt,taus=taus)
}

xts2=list()
xts=list()
alphas=list()
datasets=list()
for(i in 1:1){
  a=SampleM2(1.5,0.25,0.3,3000,25000,0.1,c(310,1760),1000)
  dt=diff(a[,1])
  alphas[[i]]=list(alpha=a[,3])
  datasets[[i]]=list(dt=dt,taus=a[,2],alpha=a[,3],Bx=a[,4],By=a[,5],Bz=a[,6]) 
  xts[[i]]=list(xt=rep(0,length(a[,3])))  
  xts2[[i]]=list(xt=rnorm(length(a[,3]),0,0.1))
}


