######### different power calculation methods #######
## power.normal: theoreitcal power calculation for normally distribute data with equal variance
## power.Colantuoni: power calculation with single bootstrap 
## power.Kleinman: power calculation with single bootstrap under the alternative
## power.Peddada: power calculation with bootstrap and monte carlo simulation under the null
## power.DoubleBootstrap: our proposed method with double boostrap under the alternative


#### note: for Colantuoni, Keleinman, Peddada's method, seems that any test would work



#### theoretical power calculation for normally distributed data (assume same variance)
power.normal = function(eff, alpha.qs, size.x, size.y){
	lambda = sqrt((size.x * size.y) / (size.x + size.x)) * eff
	pwr = pnorm(lambda - qnorm(1-alpha.qs/2)) + pnorm(qnorm(alpha.qs/2) - lambda)
	return(pwr)
}

#### theoretical power calculation for two sample t test (assume equal variance)
power.t = function(mu.x, mu.y, sigma.x, sigma.y, size.x, size.y, alpha.qs){ 
	sigma.pooled = sqrt(((size.x - 1) * sigma.x^2 + (size.y - 1)*sigma.y^2) / (size.x + size.y - 2))
	df = size.x + size.y - 2
	lambda = (mu.x - mu.y) / (sigma.pooled * sqrt(1/size.x + 1/size.y))
	critical.high= qt(p=1-alpha.qs/2, df=df)
	critical.low= qt(p=alpha.qs/2, df=df)
	pwr = pt(q=critical.high, df=df, ncp=lambda, lower.tail=FALSE) + pt(q=critical.low, df=df, ncp=lambda, lower.tail=TRUE)
	return(pwr)
}

##### using existing R function pwr.t.test()
library(pwr)
power.t.R = function(mu.x, mu.y, sigma.x, sigma.y, size.x, size.y, alpha.qs){
	sigma.pooled = sqrt(((size.x - 1) * sigma.x^2 + (size.y - 1)*sigma.y^2) / (size.x + size.y - 2))
	cohen.d = (mu.x - mu.y) / sigma.pooled
	pwr = pwr.t.test(d=cohen.d, sig.level=alpha.qs, type='two.sample', alternative='two.sided', n=size.x)$power
	return(pwr)
}


##### power estimation from Colantuoni (https://atomicules.co.uk/2013/08/06/Calculating-Sample-Size-for-Non-Normal-Distributions-Using-Bootstrapping-in-R.html)
power.Colantuoni = function(test, group.x, group.y, subsample.x, subsample.y, numboot, alpha.qs){
	pvalue = NULL
	for (b in 1:numboot){
		X.boot = sample(group.x, size=subsample.x, replace=T) # seems that subsample.x = subsample.y
		Y.boot = sample(group.y, size=subsample.y, replace=T)
		if (test == 't.test'){
			p.val = t.test(X.boot, Y.boot)$p.value
		} else if (test == 'wilcoxon'){
			p.val = wilcox.test(X.boot, Y.boot)$p.value
		}
	pvalue = cbind(pvalue, p.val)	
	}

	power = sum(pvalue <= alpha.qs)/numboot
	return(power)
}

###### power estimation from Kleinman
power.Kleinman = function(test, group.x, group.y, subsample.x, subsample.y, numboot, alpha.qs){
	pvalue = NULL
	for (b in 1:numboot){
		X.boot = sample(group.x, size=subsample.x, replace=T) 
		Y.boot = sample(group.y, size=subsample.y, replace=T)
		Y.shift = Y.boot + 5
		if (test == 't.test'){
			p.val = t.test(X.boot, Y.shift)$p.value
		} else if (test == 'wilcoxon'){
			p.val = wilcox.test(X.boot, Y.shift)$p.value
		}
	pvalue = cbind(pvalue, p.val)	
	}

	power = sum(pvalue <= alpha.qs)/numboot
	return(power)
}



##### power estimation from Peddada
power.Peddada = function(test, group.x, group.y, subsample.x, subsample.y, numboot, numrepeat, alpha.qs){
	pvalue = NULL
	for (i in 1:numrepeat){
		X.sub = sample(group.x, size=subsample.x, replace=T) # Peddada states subsample.x = subsample.y
		Y.sub = sample(group.y, size=subsample.y, replace=T)
		if (test == 't.test'){
			T.sub = t.test(X.sub, Y.sub)$statistic
		} else if (test == 'wilcoxon'){
			T.sub = wilcox.test(X.sub, Y.sub)$statistic
		}

		# bootstrap using subsampled data by combining two groups
		T.boot = NULL
		for (b in 1:numboot){
			X.boot = sample(c(X.sub, Y.sub), size=subsample.x, replace=T)
			Y.boot = sample(c(X.sub, Y.sub), size=subsample.y, replace=T)
			if (test == 't.test'){
				T.boot.tmp = t.test(X.boot, Y.boot)$statistic
			} else if (test == 'wilcoxon'){
				T.boot.tmp = wilcox.test(X.boot, Y.boot)$statistic
			}
		T.boot = cbind(T.boot, T.boot.tmp)
		}

		# empirical p-value
		pvalue.emp = sum(T.sub > T.boot)/numboot # why not use observed T but need T.sub
		pvalue = cbind(pvalue, pvalue.emp)
	}	

	# estimate power
	pwr = sum(pvalue < alpha.qs) / numrepeat  # isn't such p-values calculated under the null but not under the alternative
	return(pwr)

}	


##### new method: two-stage bootstrap power estimation
power.TwoStageBootstrap = function(group.x, group.y, subsample.x, subsample.y, numboot1, numboot2, alpha.qs){

	#### under H0: mu.x = mu.y
	# center the mean of each group to be zero
	X.center = group.x - mean(group.x)
	Y.center = group.y - mean(group.y)

	# bootstrap from centered data to calculate test statistics
	T.stat.H0 = NULL
	for (b1 in 1:numboot1) {
		X.boot.H0 = sample(X.center, size=subsample.x, replace=T) # subsample.x is the new sample size we want from group X
		Y.boot.H0 = sample(Y.center, size=subsample.y, replace=T)
		T.stat.H0 = cbind(T.stat.H0, t.test(X.boot.H0, Y.boot.H0)$statistic)
	}

	cutoff.low.H0 = quantile(T.stat.H0, probs=alpha.qs/2)
	cutoff.high.H0 = quantile(T.stat.H0, probs=1-alpha.qs/2)

	#### under H1: mu.x != mu.y
	# bootstrap from observed data to calculate test statistics
	T.stat.H1 = NULL
	for (b2 in 1:numboot2) {
		X.boot.H1 = sample(group.x, size=subsample.x, replace=T) # subsample.x is the new sample size we want from group X
		Y.boot.H1 = sample(group.y, size=subsample.y, replace=T)
		T.stat.H1 = cbind(T.stat.H1, t.test(X.boot.H1, Y.boot.H1)$statistic)
	}

	# estimate power
	power = sum((T.stat.H1 > cutoff.high.H0) | (T.stat.H1 < cutoff.low.H0))/numboot2
	return(power)

}





