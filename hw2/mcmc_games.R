## functions for mcmc samplers for problem 2, 3
library(coda)
mcmc_is <- function(iterations,burn_in,sd = 10, skill = c(0,0,0,0,0,0), thinning= 1){
	set.seed(12345)
	teams=as.matrix(read.table("teams.txt"))
	outcomes=read.table("outcomes.txt")
	g<-function(x) { 1/(1+exp(-x)) }

	pi<-function(skill) {  #Calculate the prior*likelihood, up to a constant.
	  x=rowSums(teams %*% skill)  #Calculate the parameter to the function g.
	  exp(-sum(skill**2)/8)*prod(ifelse(outcomes==1,g(x),1-g(x)))
	}

	n0=max(burn_in %/% thinning,1)
	n1=iterations %/% thinning
	mc_output=matrix(nrow=n1-n0+1,ncol=6)       #Make space for the output

	for (i in 0:n1) {
  		if (i>=n0) mc_output[i-n0+1,]=skill
  		for (j in 1:thinning) {
    		skill_=rnorm(6,0,sd)
    		alpha=min(1,pi(skill_)/pi(skill)*exp(sum(skill_**2-skill**2)/(2*sd**2)))
    		if ( runif(1) < alpha ) skill=skill_}}
    x<-mcmc(mc_output,start=n0*thinning,thin=thinning)
	return(x)
}

mcmc_rw <- function(iterations,burn_in,skill = c(0,0,0,0,0,0), thinning= 1){
	set.seed(12345)
	teams=as.matrix(read.table("teams.txt"))
	outcomes=read.table("outcomes.txt")
	g<-function(x) { 1/(1+exp(-x)) }

	pi<-function(skill) {  #Calculate the prior*likelihood, up to a constant.
	  x=rowSums(teams %*% skill)  #Calculate the parameter to the function g.
	  exp(-sum(skill**2)/8)*prod(ifelse(outcomes==1,g(x),1-g(x)))
	}

	n0=max(burn_in %/% thinning,1)
	n1=iterations %/% thinning
	mc_output=matrix(nrow=n1-n0+1,ncol=6)       #Make space for the output

	for (i in 0:n1) {
	if (i>=n0) mc_output[i-n0+1,]=skill
	for (j in 1:thinning) {
		idx = sample(1:6,1)
		skill_= skill
		skill_[idx] = skill_[idx] + rnorm(1,0,1)  

		#alpha=min(1,pi(skill_)/pi(skill))
		alpha=min(1,pi(skill_)/pi(skill))
		if ( runif(1) < alpha ) skill=skill_}}
	x<-mcmc(mc_output,start=n0*thinning,thin=thinning)
	return(x)
}

mc_vis <- function(x){
	plot(x)
	summary(x)
}


# mc_SamplePath <- function(mc_out){
# 	SamplePath = apply(mc_out, MARGIN = 2, cumsum) ## shape (niter,n_var)
# 	n = nrow(SamplePath)
# 	k = ncol(SamplePath)
# 	for(i in 1:k){
# 		p = plot(1:n, mc_out[,i]/(1:n), ylab = paste0("var ", i), xlab = "iteration", main = paste0("Sample Path of var ", i))
# 		print(p)
# 	}

# }













