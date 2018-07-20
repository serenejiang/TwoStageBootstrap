source('dataSimulation.R')
source('powerMethods.R')

########### power plot 1: power ~ effect size (fixed sample size)
plot.power.eff = function(N.x, N.y, n.x, n.y, B, m, alpha.qs, mean.x, mean.y, sd.x, sd.y, dist.x, dist.y, stat.test, repeats){

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


	eff.size = abs(mean.x - mean.y)/sd.x # assume sd.x = sd.y, also size.x = size.y

	# calculate theoretical power for normal and t distribution
	pwr.normal = power.normal(eff=eff.size, alpha.qs=alpha.qs, size.x=n.x, size.y=n.y)
	pwr.t = power.t(mu.x=mean.x, mu.y=mean.y, sigma.x=sd.x, sigma.y=sd.y, size.x=n.x, size.y=n.y, alpha.qs=alpha.qs)
	pwr.t.R = power.t.R(mu.x=mean.x, mu.y=mean.y, sigma.x=sd.x, sigma.y=sd.y, size.x=n.x, size.y=n.y, alpha.qs=alpha.qs)



	# estimate power using different bootstrap methods
	pwr.TwoStageBootstrap.m = NULL
	pwr.Kleinman.m = NULL
	pwr.Colantuoni.m = NULL
	pwr.Peddada.m = NULL

	for (i in 1:m){
		pwr.TwoStageBootstrap = NULL
		pwr.Kleinman = NULL
		pwr.Colantuoni = NULL
		pwr.Peddada = NULL

		for (i in mean.y){
			dat = data.simulation(family.x=dist.x, family.y=dist.y, mu.x=mean.x, mu.y=i, sigma.x=sd.x, sigma.y=sd.y, size.x=N.x, size.y=N.y)

			pwr.TwoStageBootstrap = cbind(pwr.TwoStageBootstrap, power.TwoStageBootstrap(group.x = dat$X, group.y = dat$Y, 
			                                      subsample.x=n.x, subsample.y=n.y, numboot1=B, numboot2=B, alpha.qs=alpha.qs))
			pwr.Kleinman = cbind(pwr.Kleinman, power.Kleinman(test=stat.test, group.x = dat$X, group.y = dat$Y, 
			                                      subsample.x=n.x, subsample.y=n.y, numboot=B, alpha.qs=alpha.qs))
			pwr.Colantuoni = cbind(pwr.Colantuoni, power.Colantuoni(test=stat.test, group.x = dat$X, group.y = dat$Y, 
			                                      subsample.x=n.x, subsample.y=n.y, numboot=B, alpha.qs=alpha.qs))
			pwr.Peddada = cbind(pwr.Peddada, power.Peddada(test=stat.test, group.x=dat$X, group.y=dat$Y, 
												  subsample.x=n.x, subsample.y=n.y, numboot=B, numrepeat=repeats, alpha.qs=alpha.qs))
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
	pdf(paste('figures/power_effSize', figure_name, sep='_'))

	### plot figures
	plot(eff.size, mean.TwoStageBootstrap, type='l', xlab='effect size', ylab='bootstrap estimated power', 
	     ylim=c(0,1), main=paste('power ~ effect size', 'N.x', N.x, dist.x, 'B', B, 'm', m, sep='.'), col='orange', lwd=3) 
	lines(eff.size, median.TwoStageBootstrap, lty=2, col='orange', lwd=2)
	lines(eff.size, Q5.TwoStageBootstrap, lty=3, col='orange')
	lines(eff.size, Q95.TwoStageBootstrap, lty=4, col='orange')

	lines(eff.size, pwr.normal, type='l', col='red', lwd=3)
	lines(eff.size, pwr.t, type='l', col='purple', lwd=3)
	lines(eff.size, pwr.t.R, type='l', col='black', lwd=3)


	lines(eff.size, mean.Kleinman, type='l', col='blue', lwd=3)
	lines(eff.size, median.Kleinman, lty=2, col='blue', lwd=2)
	lines(eff.size, Q5.Kleinman, lty=3, col='blue')
	lines(eff.size, Q95.Kleinman, lty=4, col='blue')

	lines(eff.size, mean.Colantuoni, type='l', col='green', lwd=3)
	lines(eff.size, median.Colantuoni, lty=2, col='green', lwd=2)
	lines(eff.size, Q5.Colantuoni, lty=3, col='green')
	lines(eff.size, Q95.Colantuoni, lty=4, col='green')

	lines(eff.size, mean.Peddada, type='l', col='pink', lwd=3)
	lines(eff.size, median.Peddada, lty=2, col='pink', lwd=2)
	lines(eff.size, Q5.Peddada, lty=3, col='pink')
	lines(eff.size, Q95.Peddada, lty=4, col='pink')
	legend(0.6, 0.3, c('pwr.normal', 'pwr.t', 'pwr.t.R', 'pwr.TwoStageBootstrap', 'pwr.Kleinman', 'pwr.Colantuoni', 'pwr.Peddada'), 
		col=c('red', 'purple', 'black', 'orange', 'blue', 'green', 'pink'), lty=c(1,1,1,1,1,1), lwd=c(3,3,3,3,3,3), bty='n', cex=0.8)
	dev.off()
	
}

### simple testing
plot.power.eff(N.x=10, N.y=10, n.x=10, n.y=10, B=100, m=10, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=10)

plot.power.eff(N.x=100, N.y=100, n.x=100, n.y=100, B=100, m=10, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=10)



### power ~ effect size with different original sample size
############ normal distribution ########
## N = 50
plot.power.eff(N.x=50, N.y=50, n.x=50, n.y=50, B=1000, m=500, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)


## N = 100
plot.power.eff(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=500, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)

## N = 200
plot.power.eff(N.x=200, N.y=200, n.x=200, n.y=200, B=1000, m=500, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)

## N = 400
plot.power.eff(N.x=400, N.y=400, n.x=400, n.y=400, B=1000, m=500, alpha.qs=0.05, 
			   				   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
			   				   dist.x='normal', dist.y='normal', stat.test='t.test', repeats=100)


# ############ uniform distribution ########
# plot.power.eff(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=50, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='uniform', dist.y='uniform', stat.test='wilcoxon', repeats=10)

# plot.power.eff(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=500, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='laplace', dist.y='laplace', stat.test='wilcoxon', repeats=10)

# plot.power.eff(N.x=100, N.y=100, n.x=100, n.y=100, B=1000, m=500, alpha.qs=0.05, 
# 			   mean.x=0.2, mean.y=seq(0.2, 2, by=0.05), sd.x=2, sd.y=2, 
# 			   dist.x='chi-squared', dist.y='chi-squared', stat.test='wilcoxon', repeats=10)






