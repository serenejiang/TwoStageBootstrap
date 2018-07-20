# test
source('dataSimulation.R')
source('powerMethods.R')

########### power plot 1: power ~ effect size (fixed sample size)
plot.power.boot = function(N.x, N.y, n.x, n.y, B, m, alpha.qs, mean.x, mean.y, sd.x, sd.y, dist.x, dist.y, stat.test, repeats){

	### descrption of parameters ###
	# N.x, N.y: sample size in each group of simulated data
	# n.x, n.y: subsample size in each group
	# B: number of bootstrap
	# m: number of repetition on simulating data
	# alpha.qs: type I error rate
	# mean.x, mean.y: mean of each group
	# sd.x, sd.y: sd of each group
	# dist.x, dist.y: distribution of each group in simulated data
	# stat.test: desired test to use for Kleiman, Colantuoni, Peddada's method
	# repeats: number of Monte Carlo simulation in Peddada's method


	eff.size = abs(mean.x - mean.y)/sd.x # assume sd.x = sd.y

	# calculate theoretical power for normal distribution
	pwr.normal = power.normal(eff=eff.size, alpha.qs=alpha.qs, size.x=n.x, size.y=n.y)
	pwr.t = power.t(mu.x=mean.x, mu.y=mean.y, sigma.x=sd.x, sigma.y=sd.y, size.x=n.x, size.y=n.y, alpha.qs=alpha.qs)
	pwr.t.R = power.t.R(mu.x=mean.x, mu.y=mean.y, sigma.x=sd.x, sigma.y=sd.y, size.x=n.x, size.y=n.y, alpha.qs=alpha.qs)


	# estimate power using different bootstrap methods
	pwr.TwoStageBootstrap.m = NULL
	pwr.Kleinman.m = NULL
	pwr.Colantuoni.m = NULL
	pwr.Peddada.m = NULL

	for (i in 1:m){
		dat = data.simulation(family.x=dist.x, family.y=dist.y, mu.x=mean.x, mu.y=mean.y, sigma.x=sd.x, sigma.y=sd.y, size.x=N.x, size.y=N.y)

		pwr.TwoStageBootstrap = NULL
		pwr.Kleinman = NULL
		pwr.Colantuoni = NULL
		pwr.Peddada = NULL

		for (b in B){
			pwr.TwoStageBootstrap = cbind(pwr.TwoStageBootstrap, power.TwoStageBootstrap(group.x=dat$X, group.y=dat$Y,
												  subsample.x=n.x, subsample.y=n.y, numboot1=b, numboot2=b, alpha.qs=alpha.qs))
			pwr.Kleinman = cbind(pwr.Kleinman, power.Kleinman(test=stat.test, group.x = dat$X, group.y = dat$Y, 
			                                      subsample.x=n.x, subsample.y=n.y, numboot=b, alpha.qs=alpha.qs))
			pwr.Colantuoni = cbind(pwr.Colantuoni, power.Colantuoni(test=stat.test, group.x = dat$X, group.y = dat$Y, 
			                                      subsample.x=n.x, subsample.y=n.y, numboot=b, alpha.qs=alpha.qs))
			pwr.Peddada = cbind(pwr.Peddada, power.Peddada(test=stat.test, group.x=dat$X, group.y=dat$Y, 
												  subsample.x=n.x, subsample.y=n.y, numboot=b, numrepeat=repeats, alpha.qs=alpha.qs))
		}

		pwr.TwoStageBootstrap.m = rbind(pwr.TwoStageBootstrap.m, pwr.TwoStageBootstrap)
		pwr.Kleinman.m = rbind(pwr.Kleinman.m, pwr.Kleinman)
		pwr.Colantuoni.m = rbind(pwr.Colantuoni.m, pwr.Colantuoni)
		pwr.Peddada.m = rbind(pwr.Peddada.m, pwr.Peddada)
	}

	# mean 
	mean.TwoStageBootstrap = apply(pwr.TwoStageBootstrap.m, MARGIN=2, FUN=mean) # column mean
	mean.Kleinman = apply(pwr.Kleinman.m, MARGIN=2, FUN=mean)
	mean.Colantuoni = apply(pwr.Colantuoni.m, MARGIN=2, FUN=mean)
	mean.Peddada = apply(pwr.Peddada.m, MARGIN=2, FUN=mean)

	# median
	median.TwoStageBootstrap = apply(pwr.TwoStageBootstrap.m, MARGIN=2, FUN=median) 
	median.Kleinman = apply(pwr.Kleinman.m, MARGIN=2, FUN=median)
	median.Colantuoni = apply(pwr.Colantuoni.m, MARGIN=2, FUN=median)
	median.Peddada = apply(pwr.Peddada.m, MARGIN=2, FUN=median)

	# 5% quntile
	Q5.TwoStageBootstrap = apply(pwr.TwoStageBootstrap.m, MARGIN=2, FUN=quantile, probs=0.05)
	Q5.Kleinman = apply(pwr.Kleinman, MARGIN=2, FUN=quantile, probs=0.05)
	Q5.Colantuoni = apply(pwr.Colantuoni, MARGIN=2, FUN=quantile, probs=0.05)
	Q5.Peddada = apply(pwr.Peddada, MARGIN=2, FUN=quantile, probs=0.05)

	# 95% quntile
	Q95.TwoStageBootstrap = apply(pwr.TwoStageBootstrap.m, MARGIN=2, FUN=quantile, probs=0.95)
	Q95.Kleinman = apply(pwr.Kleinman, MARGIN=2, FUN=quantile, probs=0.95)
	Q95.Colantuoni = apply(pwr.Colantuoni, MARGIN=2, FUN=quantile, probs=0.95)
	Q95.Peddada = apply(pwr.Peddada, MARGIN=2, FUN=quantile, probs=0.95)

	### export figures
	figure_name = paste(paste(dist.x, dist.y, N.x, N.y, sep='_'), 'pdf', sep='.')
	pdf(paste('figures/power_boot', figure_name, sep='_'))

	### plot figures
	plot(B, mean.TwoStageBootstrap, type='l', xlab='number of bootstrap', ylab='bootstrap estimated power', 
	     ylim=c(0,1), main=paste('power ~ #bootstrap', 'N.x', N.x, dist.x, 'eff.size', eff.size, 'm', m, sep='.'), col='orange', lwd=3) 
	lines(B, median.TwoStageBootstrap, lty=2, col='orange', lwd=2)
	lines(B, Q5.TwoStageBootstrap, lty=3, col='orange')
	lines(B, Q95.TwoStageBootstrap, lty=4, col='orange')

	lines(B, rep(pwr.normal, length(B)), type='l', col='red', lwd=3)
	lines(B, rep(pwr.t, length(B)), type='l', col='purple', lwd=3)
	lines(B, rep(pwr.t.R, length(B)), type='l', col='black', lwd=3)


	lines(B, mean.Kleinman, type='l', col='blue', lwd=3)
	lines(B, median.Kleinman, lty=2, col='blue', lwd=2)
	lines(B, Q5.Kleinman, lty=3, col='blue')
	lines(B, Q95.Kleinman, lty=4, col='blue')

	lines(B, mean.Colantuoni, type='l', col='green', lwd=3)
	lines(B, median.Colantuoni, lty=2, col='green', lwd=2)
	lines(B, Q5.Colantuoni, lty=3, col='green')
	lines(B, Q95.Colantuoni, lty=4, col='green')

	lines(B, mean.Peddada, type='l', col='pink', lwd=3)
	lines(B, median.Peddada, lty=2, col='pink', lwd=2)
	lines(B, Q5.Peddada, lty=3, col='pink')
	lines(B, Q95.Peddada, lty=4, col='pink')
	legend(10000, 0.7, c('pwr.normal', 'pwr.t', 'pwr.t.R', 'pwr.TwoStageBootstrap', 'pwr.Kleinman', 'pwr.Colantuoni', 'pwr.Peddada'), 
		col=c('red', 'purple', 'black', 'orange', 'blue', 'green', 'pink'), lty=c(1,1,1,1,1,1), lwd=c(3,3,3,3,3,3), bty='n', cex=0.8)
	dev.off()
	
}

# ### simple testing
plot.power.boot(N.x=100, N.y=100, n.x=100, n.y=100, B=seq(100, 2000, by=300), m=10, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=1, sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=10)



### power ~ # bootstrap with different original sample size
############ normal distribution ########
## N = 50
plot.power.boot(N.x=50, N.y=50, n.x=50, n.y=50, B=seq(500, 20000, by=500), m=1000, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=1, sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)


## N = 100
plot.power.boot(N.x=100, N.y=100, n.x=100, n.y=100, B=seq(500, 20000, by=500), m=1000, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=1, sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)

## N = 200
plot.power.boot(N.x=200, N.y=200, n.x=200, n.y=200, B=seq(500, 20000, by=500), m=1000, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=1, sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)

## N = 400
plot.power.boot(N.x=400, N.y=400, n.x=400, n.y=400, B=seq(500, 20000, by=500), m=1000, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=1, sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)


# ############ uniform distribution ########
# plot.power.boot(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=50, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='uniform', dist.y='uniform', stat.test='wilcoxon', repeats=10)

# plot.power.boot(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=500, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='laplace', dist.y='laplace', stat.test='wilcoxon', repeats=10)

# plot.power.boot(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=500, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='chi-squared', dist.y='chi-squared', stat.test='wilcoxon', repeats=10)






