simu.pair1=SampleM1(1.5,0.25,0.3,4200,35000,0.3)
dt=diff(simu.pair1$ts)
taus=simu.pair1$taus
a=SimpleMCMC(dt,taus,c(0.001,0.01,10,0.01,0.001),4000,1,1)


likelihood_1(dt,taus,1.5,0.25,0.3,4000,23000)

plot(a[1:4000,6],type="l",main="Trace plot for loglikelihood")

draw=seq(1500,4000,by=5)
par(mfrow=c(2,3))
hist(a[draw,1],breaks=8,main="posterior a",xlab="a");abline(v=1.5,col="red",lwd=2,lty=2)
hist(a[draw,2]-0.003,main="posterior b",xlab="b");abline(v=0.25,col="red",lwd=2,lty=2)
hist(a[draw,3],main="posterior pi1",xlab="pi1");abline(v=0.3,col="red",lwd=2,lty=2)

hist(a[draw,4],ylim=c(0,120),main="posterior of k",xlab="k");abline(v=4200,col="red",lwd=2,lty=2)
hist(a[draw,5],ylim=c(0,120),main="posterior of A0",xlab="A0");abline(v=35000,col="red",lwd=3,lty=2)
mean
maxlog1=which.max(a[,6])
maxlog1
a[maxlog1,]
for(i in 1:5){
cat(mean(a[draw,i]),'\n')
}
