#Compare perfomance of 6 tennis players in doubles matches.
teams=as.matrix(read.table("teams.txt"))
outcomes=read.table("outcomes.txt")
g<-function(x) { 1/(1+exp(-x)) }

pi<-function(skill) {  #Calculate the prior*likelihood, up to a constant.
  x=rowSums(teams %*% skill)  #Calculate the parameter to the function g.
  exp(-sum(skill**2)/8)*prod(ifelse(outcomes==1,g(x),1-g(x)))
}

skill=c(0,0,0,0,0,0)                        #Initial value of the Markov chain.
thinning=1                                  #The amount of thinning to do
iterations=1e+04                             #The number of iterations
burn_in=5000                                   #The burn-in length
n0=max(burn_in %/% thinning,1)
n1=iterations %/% thinning
mc_output=matrix(nrow=n1-n0+1,ncol=6)       #Make space for the output

for (i in 0:n1) {
  if (i>=n0) mc_output[i-n0+1,]=skill
  for (j in 1:thinning) {
    skill_=rnorm(6,0,10)
    alpha=min(1,pi(skill_)/pi(skill)*exp(sum(skill_**2-skill**2)/(2*10**2)))
    if ( runif(1) < alpha ) skill=skill_}}

library(coda)
png("mcmc_p2_%1d.png")
x<-mcmc(mc_output,start=n0*thinning,thin=thinning)
plot(x)
summary(x)

