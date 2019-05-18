## functions implememnting gibbs sampling for exercise problems
library(MASS)
library(rmutil)
library(e1071)
gibbs_graphical <- function(Theta, X = c(1,-1,-1,1,1),niter = 1000){
	set.seed(123)
	## helper function
	full_conditional <- function(Theta, X,idx){
	  fc = sum(Theta[,-idx] %*% X[-idx])
		fc = exp(2*X[idx]*fc)
		fc = fc/(1+fc)
		return(fc)
	}

	## begin gibbs sampling
	Xs = matrix(X, nrow = 1)
	Y = X ## prev
	d = length(X)
	for(t in 1:niter){
		for(i in 1:d){
			fc = full_conditional(Theta, Y, i)
			Y[i] = 2*rbinom(1,1,fc) - 1 ## get -1, or 1 based on full conditional
		}
		Xs = rbind(Xs, Y)
	}
	return(Xs)
}

flat2matrix <- function(theta, d = 5){
  Theta = matrix(replicate(d*d, 0), nrow = d)
  Theta[upper.tri(Theta)] = theta
  Theta[lower.tri(Theta)] = theta
  return(Theta)
}

likeli_log <- function(data, theta){
  ## Theta still needs to be a matrix !!
  Theta <- flat2matrix(theta)
  n = nrow(data)
  d = ncol(data)
  ## compute Z
  #browser()
  data_z = bincombinations(d) ## "data" to compute Z
  Zs = rowSums((data_z %*% Theta) * data_z) ## of size 2^d
  Z = sum(exp(Zs))
  #browser()
  ## compute ll for all data together
  XOX = data %*% Theta 
  XOX = XOX * data ##now rowsum is lik for each sample
  lls = rowSums(XOX)
  ll = sum(lls - log(Z))
  return(ll)
}

prior_log <- function(theta, lam){
  return(sum(dlaplace(theta,0,lam, log = T)))
}

posterir_log <- function(data, theta, lam){
  return(likeli_log(data, theta) + prior_log(theta, lam))
}

accept <- function(data, theta_, theta,lam){
  prob = exp(posterir_log(data, theta_, lam) - posterir_log(data, theta, lam))
  return(min(1, prob))
}

rw_sampler <- function(theta, s, d = 5){ ## get new rho_ from old rho
  ## need to replace hard coded d
  theta_new = theta + mvrnorm(1, mu = replicate(d, 0), s*diag(d))
  return(theta_new)
}

MH_graphical <- function(data, theta, s = 0.1, lam = 0.2, niter = 1000,thinning = 1){
	## helper functions
	## computation
	mc_output=matrix(nrow=niter,ncol=length(theta))       #Make space for the output

	for (i in 1:niter) {
  		mc_output[i,]=theta
  		for (j in 1:thinning) {
    		theta_ = rw_sampler(theta, s)
    		alpha= accept(data,theta_, theta, lam)
    		# print(alpha)
    		# if(is.nan(alpha)){
    		#   return(theta)
    		# }
    		if ( runif(1) < alpha ) theta=theta_}}
	return(mc_output)
}




## test
Theta = matrix(c(0, 0.5, 0, 0.5, 0, 
                 0.5, 0, 0.5, 0, 0.5,
                 0, 0.5, 0, 0.5, 0, 
                 0.5, 0, 0.5, 0, 0,
                 0, 0.5, 0, 0, 0), nrow = 5)
burn_in = 5000
N = 1000
start = proc.time()
out = gibbs_graphical(Theta, niter = N + burn_in)
cat(sprintf("runtime %s\n", (proc.time() - start)[[3]]))
samples = out[(burn_in + 1): nrow(out),]

d = 5
theta = replicate(d*(d-1)/2, 0)
#theta = test
theta_chain = MH_graphical(samples, theta, s = 0.1, lam = 0.2, niter = 10000,thinning = 1)

for(i in 1:10){
  plot(theta_chain[,i])
}
  











