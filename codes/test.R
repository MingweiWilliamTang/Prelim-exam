# test
ag_likelihood0(datasets,alphas,1.5,0.25,0.2,3000,25000)
ag_likelihood(datasets,alphas,xts,1.5,0.25,0.3,3000,25000)

agg_log_path(datasets,xts,1,0.01)
res=updateTheta(datasets,alphas,xts,c(0.5,0.5,10,0.01,0.001),2,1,0.3,1000,28000,iter = 500)


ag_likelihood0(datasets,alphas,0.21,0.01,0.99,3000,50000)
ag_likelihood(datasets,alphas,xts,1.5,0.25,0.3,3000,28000)
ag_likelihood(datasets,alphas,ag_res1$r2,1.5,0.25,0.3,3000,28000)

ag_likelihood0(datasets,alphas,1,0.725,0.2,600,35000)
ag_likelihood(datasets,alphas,xts,1.5,0.25,0.2,2000,25000)


set.seed(518)
updateTheta(datasets[1],alphas[1],xts[1],c(0.5,0.5,10,0.01,0.001),1.5,0.25,0.2,2000,25000,iter = 500)
thetaupdate2(datasets[[1]]$dt,datasets[[1]]$taus,alphas[[1]]$alpha,
             xts[[1]]$xt,c(0.5,0.5,10,0.01,0.001),2,1,0.3,1000,30000,1000)

likelihood_3(datasets[[1]]$dt,datasets[[1]]$taus,alphas[[1]]$alpha,
             xts[[1]]$xt,2,0.96,0.54,2230,40000)

ag_likelihood(datasets[1],alphas[1],xts[1],1.5,0.25,0.4,2000,35000)


