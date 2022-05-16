library(ggplot2)
library(coda)
library(rjags)
library(R2WinBUGS)

################################################################################
################################################################################
# Step O - Data prep 
################################################################################
################################################################################
rm(list=ls())
cannabis <- read.csv("cannabis.txt",sep="",header = TRUE)

################################################################################
################################################################################
# Question 3
################################################################################
################################################################################
## Functions for the Metropolis algorithm
log.f.2 <- function(pi,N,xplus) {
  (N-xplus)*log(1-pi)+xplus*log(1+35*pi)
}

metropolis <- function(pi0,M, sd.prop, lpost, N, xplus) {
  pi = numeric(M)
  pi[1]=pi0
  
  n.accept = 0
  for(i in 2:M) {
    pi.prop = pi[i-1]+rnorm(1,0,sd.prop)
    prop = min(1,exp(lpost(pi.prop,N,xplus) - lpost(pi[i-1],N,xplus)))
    accept = (runif(1) <= prop)
    if(accept) {
      n.accept = n.accept+1
      pi[i] = pi.prop
    } else {
      pi[i] = pi[i-1]
    }
  }
  
  accept.rate = round(n.accept/(M-1),2)
  cat("Acceptance rate: ",accept.rate,"\n")
  
  return(pi)
}

## Running of the metropolis algo on the full population
M <- 10000
pi = metropolis(
  pi0=0.5,M,sd=0.028,
  lpost=log.f.2,N=nrow(cannabis),xplus=sum(cannabis$y))

## Convergence
## acceptance rate is 40% : OK (univariate)
geweke.diag(mcmc(pi))
geweke.plot(mcmc(pi),nbins=100)

## summary
traceplot(mcmc(pi))
quantile(pi[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi[200:M]))
effectiveSize(mcmc(pi))

################################################################################
################################################################################
# Question 4
################################################################################
################################################################################
## Metropolis for the male subset
can_m =subset(cannabis,male==1)
pi_m = metropolis(
  pi0=0.5,M,sd.prop=0.030,
  lpost=log.f.2,N=nrow(can_m),xplus=sum(can_m$y))

## Convergence
geweke.diag(mcmc(pi_m))
geweke.plot(mcmc(pi_m),nbins=100)

## summary
traceplot(mcmc(pi_m))
quantile(pi_m[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi_m[200:M]))
effectiveSize(mcmc(pi_m))

## Metropolis for the female subset
can_f=subset(cannabis,male==0)
pi_f = metropolis(
  pi0=0.5,M,sd.prop=0.0250,
  lpost=log.f.2,N=nrow(can_f),xplus=sum(can_f$y))

## Convergence
geweke.diag(mcmc(pi_f))
geweke.plot(mcmc(pi_f),nbins=100)

## summary
traceplot(mcmc(pi_f))
quantile(pi_f[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi_f[200:M]))
effectiveSize(mcmc(pi_f))

## What about the delta
delta = pi_m[200:M] - pi_f[200:M]

## convergence
geweke.diag(mcmc(delta))
geweke.plot(mcmc(delta))

## summary
traceplot(mcmc(delta))
densplot(mcmc(delta))
HPDinterval(mcmc(delta))

################################################################################
################################################################################
# Question 5
################################################################################
################################################################################
model <- function() {
  #sample
  for(i in 1:n) {
    y[i] ~ dbin((1+35*pi[i])/36,1)
    logit(pi[i]) <- ((alpha0 + alpha1*(age[i]-40))*male[i]) + ((beta0 + beta1*(age[i]-40))*(1-male[i]))
    
    #logit(pi[i]) <- alpha0 + alpha1*(age[i]-40)
  }
  
  # prior
  alpha0 ~ dnorm(0,1E-6)
  alpha1 ~ dnorm(0,1E-6)
  beta0  ~ dnorm(0,1E-6)
  beta1  ~ dnorm(0,1E-6)
  
  # function to monitor
  delta <- alpha1 - beta1
  
  # predict the odds for males who are 25 year old
  etha_m_25 <- alpha0 + alpha1*(25-40)
}


write.model(model,"cannabis_2.bug")


cannabis.list <- with(cannabis,
                      list(y=y,
                           age=age,
                           male=male,
                           n=nrow(cannabis)))

inits.list <- list(list(
  alpha0=0,
  alpha1=0,
  beta0=0,
  beta1 =0),
  list(
    alpha0=-4,
    alpha1=-1,
    beta0=-5,
    beta1=-1
  )
)

params <- c("alpha0","alpha1","beta0","beta1","delta","etha_m_25")

cannabis.model <- jags.model(
  file="cannabis_2.bug",
  data=cannabis.list,
  inits = inits.list ,
  n.chains = length(inits.list))

update(cannabis.model,1000)

out = coda.samples(model=cannabis.model,
                   variable.names = params,
                   n.iter=10)

out.matrix = as.matrix(out)

effectiveSize(out)
#plot(out)
summary(out)
HPDinterval(out)

#hist(out.matrix[,"delta"],freq=F)
