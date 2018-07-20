#### simulate data from location-scale family
library(rmutil) # for rplace()


## function to generate simulated data for group X and group Y
data.simulation = function(family.x, family.y, mu.x, mu.y, sigma.x, sigma.y, size.x, size.y){
  # family.x: distribution family for group X
  # family.y: distribution family for gorup Y
  # mu.x: mean of group X
  # mu.y: mean of group Y
  # sigma.x: variance of group X
  # sigma.y: variance of group Y

  ## general function for generate g
  data.g = function(family, size){
    # family: specifiy distribution family ('normal', 'uniform', 'laplace', 'chi-squared')
    # size: number of samples

      if (family == 'normal'){
      # standard normal with mean 0 and sd 1
      g = rnorm(n=size, mean=0, sd=1) 
    } else if (family == 'uniform'){
      # uniform with mean 0 and sd 1
      g = runif(n=size, min=-sqrt(3), max=sqrt(3)) 
    } else if (family == 'laplace'){
      # laplace with mean 0 and sd 1
      g = rlaplace(n=size, m=0, s=1/sqrt(2)) 
    } else if (family == 'chi-squared'){
      # Chi-squared with mean 0 and sd 1
      df_chi =3
      g = (rchisq(n=size, df=df_chi, ncp=0) - df_chi) /sqrt(2*df_chi)
    }
    return(g)
  }
  
  ## generate data for each group
  data.X = mu.x + sigma.x * data.g(family=family.x, size=size.x)
  data.Y = mu.y + sigma.y * data.g(family=family.y, size=size.y)

  ## output simulated data: store as list objects to allow for different length of X and Y
  data = list(X=data.X, Y=data.Y)
  return(data)
}

