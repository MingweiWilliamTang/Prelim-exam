ag_res1=agg_cont_diff1(datasets,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.5),3000,1,1)
ag_res1=agg_two_states(datasets,c(0.01,0.1,10,0.01,0.001,1),3000,1,1)
idx=seq(2000,4000,seq=3)

hist(ag_res1$r1[2500:3000,1],breaks=15,main="posterior sample of a",xlim=c(1.2,1.8)); abline(v=1.5,lwd=2,col="red")
hist(ag_res1$r1[2500:3000,2],breaks=12,main="posterior sample of b",xlim=c(0.22,0.27));abline(v=0.25,lwd=2,col="red")
hist(ag_res1$r1[2500:3000,3],main="posterior sample of pi_1");abline(v=0.3,lwd=2,col="red")
hist(ag_res1$r1[2500:3000,4],main="posterior sample k",freq=T,ylim=c(0,120));abline(v=3000,lwd=2,col="red")
hist(ag_res1$r1[500:1500,5],main="posterior sample of A0");abline(v=25000,lwd=2,col="red")


ag_likelihood(datasets,alphas,xts,1.4,0.25,0.3,3500,25000)

plot(c(0,cumsum(datasets[[1]]$dt)),ylim=c(0,1),datasets[[1]]$alpha,type="l")
lines(c(0,cumsum(datasets[[1]]$dt)),ag_res1$r3[[1]]$alpha,col="red")

plot(c(0,cumsum(datasets[[1]]$dt)),ag_res1$r2[[1]]$xt,type="l")
mean(ag_res1$r1[1500:3000,8])
ag_res1[,8]


agg_log_path(datasets,ag_res1$r2,250,5)
agg_path(datasets,xts,250,0.0001)


##############
par(mfrow=c(2,3))
ag_res_twostates=agg_two_states(datasets1,c(0.01,0.1,10,0.01,0.001,1),3000,1,1)
hist(ag_res_twostates$r1[2500:3000,1],breaks=5,main="posterior sample of a",xlab="a",xlim=c(1.4,1.55)); abline(v=1.5,lwd=2,col="red")
hist(ag_res_twostates$r1[2500:3000,2],breaks=5,main="posterior sample of b",xlab="b",xlim=c(0.248,0.252));abline(v=0.25,lwd=2,col="red")
hist(ag_res_twostates$r1[2500:3000,3],main="posterior sample of pi_1",xlab="pi1");abline(v=0.3,lwd=2,col="red")
hist(ag_res_twostates$r1[2500:3000,4],main="posterior sample k",freq=T,xlab="k",ylim=c(0,160));abline(v=3000,lwd=2,col="red")
hist(ag_res_twostates$r1[2500:3000,5],main="posterior sample of A0",xlab="A0");abline(v=25000,lwd=2,col="red")

hist(ag_res_twostates$r1[2500:3000,7],main="posterior sample of xi");


ag_res_diff=agg_cont_diff(datasets1,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.5),3000,1,1)
par(mfrow=c(2,3))
hist(ag_res_diff$r1[2500:3000,1],breaks=5,main="posterior sample of a",xlim=c(1.4,1.55),xlab="a"); abline(v=1.5,lwd=2,col="red")
hist(ag_res_diff$r1[2500:3000,2],breaks=5,main="posterior sample of b",xlim=c(0.244,0.255),xlab="b");abline(v=0.25,lwd=2,col="red")
hist(ag_res_diff$r1[2500:3000,3],main="posterior sample of pi_1",xlab="pi1");abline(v=0.3,lwd=2,col="red")
hist(ag_res_diff$r1[2500:3000,4],main="posterior sample k",freq=T,ylim=c(0,160),xlab="k");abline(v=3000,lwd=2,col="red")
hist(ag_res_diff$r1[2500:3000,5],main="posterior sample of A0",xlab="A0");abline(v=25000,lwd=2,col="red")
hist(ag_res_diff$r1[2500:3000,7],breaks=5,main="posterior sample of xi",xlim=c(0,0.1),xlab="xi")

hist(ag_res_twostates$r1[2500:3000,7],main="posterior sample of xi");

hist(ag_res_diff$r1[2500:3000,7])



ag_twostates_real=agg_cont_diff1(datasets2,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.1),3000,1,1)
hist(ag_twostates_real$r1[2500:3000,1],breaks=15,main="posterior sample of a")
hist(ag_twostates_real$r1[2500:3000,2],breaks=5,main="posterior sample of b")
hist(ag_twostates_real$r1[2500:3000,3],main="posterior sample of pi_1")
hist(ag_twostates_real$r1[2500:3000,4],main="posterior sample k")
hist(ag_twostates_real$r1[2500:3000,5],main="posterior sample of A0")
hist(ag_twostates_real$r1[2500:3000,7])

ag_twostates_real$r1[2500:3000,10]

plot(c(0,cumsum(datasets2[[6]]$dt)),ag_twostates_real$r2[[6]]$xt,type="l")


ag_twostates_real2=agg_cont_diff1(datasets2,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.5),3000,1,1)
hist(ag_twostates_real2$r1[2500:3000,1],breaks=15,main="posterior sample of a")
hist(ag_twostates_real2$r1[2500:3000,2],breaks=5,main="posterior sample of b")
hist(ag_twostates_real2$r1[2500:3000,3],main="posterior sample of pi_1")
hist(ag_twostates_real2$r1[2500:3000,4],main="posterior sample k")
hist(ag_twostates_real2$r1[2500:3000,5],main="posterior sample of A0")
hist(ag_twostates_real2$r1[2500:3000,7])
mean(ag_twostates_real2$r1[2000:3000,10])
plotrealhist(ag_twostates_real2)
#ag_twostates_real3=agg_two_states(datasets1,c(0.01,0.1,10,0.01,0.001,1),3000,1,1)

par(mfrow=c(2,3))
ag_twostates_real3=agg_cont_diff1(datasets2,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.5),3000,1,1)
hist(ag_twostates_real4$r1[2500:3000,1],xlab="a",breaks=15,main="posterior sample of a")
hist(ag_twostates_real4$r1[2500:3000,2],xlab="b",breaks=5,main="posterior sample of b")
hist(ag_twostates_real4$r1[2500:3000,3],xlab="pi1",main="posterior sample of pi_1")
hist(ag_twostates_real4$r1[2500:3000,4],xlab="k",main="posterior sample k")
hist(ag_twostates_real4$r1[2500:3000,5],xlab="A0",main="posterior sample of A0")
hist(ag_twostates_real3$r1[2500:3000,7],xlab="xi",main="posterior sample of xi",xlim=c(0,0.4))
mean(ag_twostates_real4$r1[2000:3000,8])
mean(ag_twostates_real4$r1[2000:3000,9])

plotrealhist(ag_twostates_real3)
ag_twostates_real4=agg_two_states(datasets2,c(0.01,0.1,10,0.01,0.001,1),3000,1,1)
ag_twostates_real3=agg_cont_diff1(datasets2,c(0.01,0.1,10,0.01,0.001,1,0.1,1,0.1),3000,1,1)

plotrealhist=function(ag_twostates_real4){
  par(mfrow=c(2,3))
  hist(ag_twostates_real4$r1[2500:3000,1],xlab="a",breaks=15,main="posterior sample of a")
  hist(ag_twostates_real4$r1[2500:3000,2],xlab="b",breaks=5,main="posterior sample of b")
  hist(ag_twostates_real4$r1[2500:3000,3],xlab="pi1",main="posterior sample of pi_1")
  hist(ag_twostates_real4$r1[2500:3000,4],xlab="k",main="posterior sample k")
  hist(ag_twostates_real4$r1[2500:3000,5],xlab="A0",main="posterior sample of A0")
}
