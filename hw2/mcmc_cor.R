## mcmc function for estimating correlation in problem 4.4
library(coda)
mcmc_cor <- function(data, iterations,burn_in, rho = 0.1, thinning= 1, size = 0.1){
	set.seed(12345)

	ll <- function(rho, data){ ## compute loglikelihood of data given rho
		N = nrow(data)
		#lls <- apply(data, 1, ll_helper, rho)
		X = data[,1]
		Y = data[,2]
		lls = X^2 + Y^2 - 2*rho*X*Y
		ll <- -1/(2*(1-rho^2)) * sum(lls) - N * log(2*pi*sqrt(1-rho^2))
		return(ll)
	}

	prior <- function(rho){
		return(1/(pi*(1-rho^2)^(1/2)))
	}

	target_ratio <-function(rho_, rho, data) {  # compute ratio of target posterior(rho_)/posterior(rho)
	 	ll = ll(rho, data)
	 	ll_ = ll(rho_, data)
	  	ratio = exp(ll_ - ll) * prior(rho_)/prior(rho)
	 	return(ratio)
	}

	accept <- function(rho_, rho, data){ ## compute the acceptance probability
		return(min(1, target_ratio(rho_, rho, data)))
	}

	rw_sampler <- function(rho, size){ ## get new rho_ from old rho
		return(runif(1, rho-size, rho+size))
	}

	## computation

	n0=max(burn_in %/% thinning,1)
	n1=iterations %/% thinning
	mc_output=matrix(nrow=n1-n0+1,ncol=1)       #Make space for the output

	for (i in 0:n1) {
  		if (i>=n0) mc_output[i-n0+1,]=rho
  		for (j in 1:thinning) {
    		rho_ = rw_sampler(rho, size)
    		alpha= accept(rho_, rho, data)
    		if ( runif(1) < alpha ) rho=rho_}}
    x<-mcmc(mc_output,start=n0*thinning,thin=thinning)
	return(x)
}


running_var_mean <- function(mc_out){
  var_running = c()
  mean_running = c()
  for(i in 1:length(mc_out)){
    var_running <- c(var_running, var(mc_out[1:i]))
    mean_running <- c(var_running, mean(mc_out[1:i]))
  }
  return(list(var = var_running, mean = mean_running))
}












