library(rjags)
library(coda)
library(R2WinBUGS)

rm(list=ls())
set.seed(20220501)

cannabis = read.csv("cannabis.txt",sep="",header=T)


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
  beta1 =0)
  )

params <- c("alpha0","alpha1","beta0","beta1","delta","etha_m_25")
#params <- c("alpha0","alpha1")

cannabis.model <- jags.model(
  file="cannabis_2.bug",
  data=cannabis.list,
  inits = inits.list ,
  n.chains = length(inits.list))

update(cannabis.model,1000)

out = coda.samples(model=cannabis.model,
                   variable.names = params,
                   n.iter=10000)

out.matrix = as.matrix(out)

effectiveSize(out)
#plot(out)
summary(out)
HPDinterval(out)

#hist(out.matrix[,"delta"],freq=F)

geweke.diag(mcmc(out.matrix[,"beta1"]))
geweke.plot(mcmc(out.matrix[,"beta1"]),nbin=15)

traceplot(mcmc(out.matrix[,"beta1"]))
#gelman.diag(list(out[1],out[2]))


