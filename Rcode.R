n = 100
xbar = 0
gamma_rate = 0
mu = 5
for(i in 1:n) {
  xbar = xbar + Assign4Data[i]
  gamma_rate = gamma_rate + 0.5 * (Assign4Data[i] - mu)^2
}
xbar = as.numeric(xbar)/n
gamma_rate = as.numeric(gamma_rate)
tau = 0.25

#-------------------------------------------------------------------------------------------------------------------------------------
# The Gibbs sampling function
gibbs.f2 = function(mu0, tau0, m){
  mu.seq = tau.seq = rep(-1, m+1)
  mu.seq[1] = mu0
  tau.seq[1] = tau0
  for(j in 2:(m+1)) {
    mu.seq[j] = rnorm(1, xbar, sqrt(1/(n*tau.seq[j-1])))
    tau.seq[j] = rgamma(1, n/2, gamma_rate)
  }
  result = list(mu = mu.seq, tau = tau.seq)
  result
}

#m = number of samples of (mu, tau) values.
m = 1000
#Perform Gibbs sampling
gibbs1 = gibbs.f2(7, 1, m)
gibbs2 = gibbs.f2(0.5, 0.1, m)


plot(1:(m+1), gibbs1$mu, type='l', col='blue')
points(1:(m+1), gibbs2$mu, type='l', col='red')


plot(1:(m+1), gibbs1$tau, type='l', col='blue')
points(1:(m+1), gibbs2$tau, type='l', col='red')

#-------------------------------------------------------------------------------------------------------------------------------------

# Plots for marginal posterior distributions
par(mfrow=c(2,2))
plot(density(gibbs1$mu), ylab = "posterior dist", xlab = "mu1", main="")
plot(density(gibbs2$mu), ylab = "posterior dist", xlab = "mu2", main="")
plot(density(gibbs1$tau), ylab = "posterior dist", xlab = "tau1", main="")
plot(density(gibbs2$tau), ylab = "posterior dist", xlab = "tau2", main="")

# Marginal posterior means
mean(gibbs1$mu)
mean(gibbs2$mu)
mean(gibbs1$tau)
mean(gibbs2$tau)

#Fidning credible interval percentage for mu and tau
cdfmu = ecdf(gibbs1$mu)
cdfmu(c(4.77, 5.42))

cdftau = ecdf(gibbs1$tau)
cdftau(c(0.199,0.312))
#-------------------------------------------------------------------------------------------------------------------------------------
# Convert the data into a vector
data_list = c()
for(i in 1:100) {
  data_list = c(data_list, as.numeric(Assign4Data[i]))
}

# Determining the Likelihood Function
l.likelihood = function(parameters) {
  mu = parameters[1]
  tau = parameters[2]
  
  singlelikelihood = dnorm(data_list, mean=mu, sd = sqrt(1/tau), log=T)
  return(sum(singlelikelihood))
}

# Determining the Prior
prior = function(parameters) {
  return(log(1/parameters[2]))
}

# Determining the Posterior
posterior = function(parameters) {
  return(l.likelihood(parameters) + prior(parameters))
}

# Determining the Proposal Function
prop.f = function(parameters){
  tau.n = rgamma(1, 5*parameters[2], 5)
  mu.n = rnorm(1, parameters[1], sqrt(tau.n))
  return(c(mu.n, tau.n))
}

# The Metropolis-Hastings Algorithm
run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1, 2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = prop.f(chain[i,])
    probability = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probability){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

#m = Number of iterations
m = 10000
startvalue = c(7, 1)
startvalue2 = c(0.5, 0.1)
chain = run_metropolis_MCMC(startvalue, m)
chain2 = run_metropolis_MCMC(startvalue2, m)


plot(1:(m+1), chain[,1], ylab = "mu", type='l', col='blue')
points(1:(m+1), chain2[,1], ylab = "mu", type='l', col='red')
plot(1:(m+1), chain[,2], ylab = "tau", type='l', col='blue')
points(1:(m+1), chain2[,2], ylab = "tau", type='l', col='red')

#-------------------------------------------------------------------------------------------------------------------------------------
# Plots for marginal posterior distributions
par(mfrow=c(2,2))
plot(density(chain[,1]), ylab = "posterior dist", xlab = "mu1", main="")
plot(density(chain2[,1]), ylab = "posterior dist", xlab = "mu2", main="")
plot(density(chain[,2]), ylab = "posterior dist", xlab = "tau1", main="")
plot(density(chain2[,2]), ylab = "posterior dist", xlab = "tau2", main="")

# Marginal posterior means
mean(chain[,1])
mean(chain2[,1])
mean(chain[,2])
mean(chain2[,2])

#Find credible interval percentage for mu and tau
cdfmu = ecdf(chain[,1])
cdfmu(c(4.77, 5.41))

cdftau = ecdf(chain[,2])
cdftau(c(0.188,0.302))
