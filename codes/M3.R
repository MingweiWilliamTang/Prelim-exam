M2a=a
res20=ContdiffMCMC(dt,taus,c(0.01,0.1,10,0.01,0.001,1,0.1,0.5),c(310,1760),1000,3000,1,1,a[,4],a[,5],a[,6])

mcmc=res2$MCMC[2500:3000,]
par(mfrow=c(2,3))
hist(mcmc[,1],main="posterior histogram of a",xlab="a");abline(v=1.5,col="red",lwd=2)
hist(mcmc[,2],main="posterior histogram of b",xlab="b");abline(v=0.25,col="red",lwd=2)
hist(mcmc[,3],breaks=15,main="posterior histogram of pi1",xlab="pi1");abline(v=0.3,col="red",lwd=2)
hist(mcmc[,4],breaks=30,main="posterior histogram of k",xlab="k");abline(v=3000,col="red",lwd=2)
hist(mcmc[,5],breaks=20,main="posterior histogram of A0",xlab="A0");abline(v=25000,col="red",lwd=2)
#hist(mcmc[,6])
hist(mcmc[,7],breaks=5,xlim=c(0,0.7*10^(-4)),main="posterior histogram of xi",xlab="xi")
plot(res2$MCMC[,8],type="l")

## path

plot(a[,1],a[,3],type="l",ylim=c(0,1))


plot(a[,1],a[,3],type="l",ylim=c(0,1))

for(i in 2500:3000){
  lines(a[,1],res2$path[i,],col="red",lwd=0.1)
}

likelihood_3(dt,taus,res$path[2000,],res$energy[2000,],1.5,0.25,0.2,1000,2500)


likelihood_3(dt,taus,a[,3],res$energy[2000,],1.5,0.25,0.2,1000,25000)

likelihood_3(dt,taus,res2$path[2500,],rep(0,length(taus)),1.5,0.25,0.3,3000,25000)
likelihood_3(dt,taus,res2$path[2500,],res2$energy[2500,],1.5,0.25,0.3,3000,25000)


likelihood_3(dt,taus,a[,3],res2$energy[2000,],1.5,0.25,0.3,3000,25000)
likelihood_3(dt,taus,a[,3],rep(0,length(taus)),1.5,0.25,0.3,3000,25000)
likelihood_3(dt,taus,a[,3],rep(0,length(taus)),1.5,0.25,0.2,1000,25000)

plot(a[,1],rep(0,length(a[,1])),type="l",ylim=c(-1,1))
for(i in 1500:2000)
{lines(a[,1],res20$energy[i,],col="red",lwd=0.5)}
res$energy[1,]-res$energy[500,]
plot(res$energy[,])
logPath(res$energy[1,],dt,0.1,0.1)
logPath(rep(0,length(taus)),dt,0.1,0.1)



logPath(p,rep(1,10),0.1,1)-logPath(q,rep(1,10),0.1,1)
OUratio(4,q,rep(1,10),10,0.1,1)

q=rnorm(11)
p=q;
p[5]=10


likelihood_3(dt,taus,a[,3],res$energy[2000,],1.5,0.25,0.2,3000,25000)
