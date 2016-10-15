init=Balpha(dt)
dt=diff(a1[,1])
taus=(a1[,2])
ma=updateBrowian(dt,taus,init[4,],init[1,],
              init[2,],init[3,],rep(0,length(taus)),
              1.5,0.25,0.3,3000,25000,1,5000,1)
plot(a1[,1],a1[,3],type="l",ylim=c(0,1))
lines(a1[,1],ma[4,],type="l")
which.max(abs(q-a[,3]))
q[1535]=0.8
a[1535,3]

likelihood_3(dt,taus,,rep(0,length(taus)),1.5,0.25,0.2,2000,35000)

lines(a[,1],q,lty=1,col="blue",lwd=1)
lines(a[,1],ma[4,],lty=1,col="red",lwd=1)
for(i in 1950:2000){
  lines(a[,1],ma$path[i,],lty=1,col="green",lwd=0.5)
}

ma$path[1,]-ma$path[100,]

pa=OUpath(dt,1,1)
up=updatext(dt,taus,a[,3],pa,100,0.01,
         1.5,0.25,0.2,1000,25000,0.1,1000)
plot(x=c(1,length(dt)),c(0,0),type="l",ylim=c(-1,1))
for(i in 950:1000)
{lines(up,col="red",lwd=0.1)}
lines(up[1,],col="blue",lwd=1)

updatext(datasets[[1]]$dt,datasets[[1]]$taus,datasets[[1]]$alpha,xts[[1]]$xt,0.,,1.5,0.25,0.3,3000,25000,0.5,100)
UpdateOUPara(rep(0,1526)+rnorm(1526,0,0.001),datasets[[1]]$dt,100,1,1,0.5,1000)
