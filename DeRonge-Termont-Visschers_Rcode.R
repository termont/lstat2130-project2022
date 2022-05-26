library(ggplot2)
library(coda)
library(rjags)
library(R2WinBUGS)

################################################################################
################################################################################
# Step O - Data preprocessing 
################################################################################
################################################################################
rm(list=ls())
cannabis <- read.csv("cannabis.txt",sep="",header = TRUE)
set.seed(20220501)

################################################################################
################################################################################
# Question 3
################################################################################
################################################################################

# Question 3.A

##The logarithm of the posterior distribution found in the previous question
log.f.2 <- function(pi,N,yplus) {
  (N-yplus)*log(1-pi)+yplus*log(1+35*pi)
}

## The Metropolis algorithm
metropolis <- function(pi0,M, sd.prop, lpost, N, yplus) {
  pi = numeric(M)
  pi[1]=pi0
  
  n.accept = 0
  for(i in 2:M) {
    pi.prop = pi[i-1]+rnorm(1,0,sd.prop)
    prop = min(1,exp(lpost(pi.prop,N,yplus) - lpost(pi[i-1],N,yplus)))
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
## acceptance rate in the range 0.3-0.5 is accepted

## first chain
M <- 10000
pi = metropolis(
  pi0=0.5,M,sd=0.02,
  lpost=log.f.2,N=nrow(cannabis),yplus=sum(cannabis$y))

## second chain
pi_2 = metropolis(
  pi0=0.75,M,sd=0.02,
  lpost=log.f.2,N=nrow(cannabis),yplus=sum(cannabis$y))

## The effectif size
effectiveSize(mcmc(pi[200:M]))

## Convergence
traceplot(mcmc(pi))

##1. the Gelman-Rubin diagnostic
gelman.diag(list(mcmc(pi),mcmc(pi_2)))

##2. the Geweke diagnostic
geweke.diag(mcmc(pi))
geweke.plot(mcmc(pi),nbins=100)

# Question 3.B
quantile(pi[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi[200:M]))

# Question 3.C
nrow(subset(data.frame(pi[200:M]),pi>=0.1))/(M-200)

################################################################################
################################################################################
# Question 4
################################################################################
################################################################################

# 1. Men

## Metropolis for the male subset
can_m =subset(cannabis,male==1)
pi_m = metropolis(
  pi0=0.5,M,sd.prop=0.033,
  lpost=log.f.2,N=nrow(can_m),yplus=sum(can_m$y))

## The effectif size
effectiveSize(mcmc(pi_m[200:M]))

## Convergence
traceplot(mcmc(pi_m))

## The Geweke diagnostic
geweke.diag(mcmc(pi_m))
geweke.plot(mcmc(pi_m),nbins=100)

## 95% credible interval for pi_m
quantile(pi_m[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi_m[200:M]))

## The posterior probability that among 20-59 years old men, 
##there is at least 10% of cannabis’ users 
mean(pi_m>0.1)

#2. Women

## Metropolis for the female subset
can_f=subset(cannabis,male==0)
pi_f = metropolis(
  pi0=0.5,M,sd.prop=0.024,
  lpost=log.f.2,N=nrow(can_f),yplus=sum(can_f$y))

## The effectif size
effectiveSize(mcmc(pi_f[200:M]))

## Convergence
traceplot(mcmc(pi_f))

## The Geweke diagnostic
geweke.diag(mcmc(pi_f))
geweke.plot(mcmc(pi_f),nbins=100)

## 95% credible interval for pi_f
quantile(pi_f[200:M],c(0.025,0.975))
HPDinterval(mcmc(pi_f[200:M]))

## The posterior probability that among 20-59 years old women, 
##there is at least 10% of cannabis’ users 
mean(pi_f>0.1)

# 3.Delta

##the value of pi for men
f.pi_m = data.frame(pi_m)
f.pi_m$seq = seq(1,length(pi_m))

##the value of pi for women
f.pi_f = data.frame(pi_f)
f.pi_f$seq = seq(1,length(pi_f))

##plots both the sample of the value of pi for men
## and the value of pi for women
ggplot() +
  theme_bw() +
  theme(panel.border = element_blank()) +
  geom_point(data=f.pi_m,aes(x=seq,y=pi_m),color="blue") +
  geom_point(data=f.pi_f,aes(x=seq,y=pi_f),color="red") +
  geom_vline(aes(xintercept=200),color="black",size=1) +
  labs(
    title="MCMC chain for pi_m and pi_f",
    subtitle="pi_m in blue, pi_f in red",
    y="pi_m \ pi_f")

## The plausible values of delta
delta = pi_m[200:M] - pi_f[200:M]

## Convergence
traceplot(mcmc(delta))

## The Geweke diagnostic
geweke.diag(mcmc(delta))
geweke.plot(mcmc(delta))

## The density of delta
densplot(mcmc(delta))

## Histogram of delta
plot(pi_m,ylim=c(0,0.5))
points(pi_f,col="red")
hist(delta)

## 95% credible interval for delta
HPDinterval(mcmc(delta))

## The probability that delta is bigger than 0
nrow(subset(data.frame(delta[200:M]),delta>0))/(M-200)


################################################################################
################################################################################
# Question 5
################################################################################
################################################################################

# Question 5.A

## Logistic Regression
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

# Create two chains with differences values
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
                   n.iter=20000)

out.matrix = as.matrix(out)

## Question 5.B

summary(out)
HPDinterval(out)

## Question 5.C

## mean etha_m_25
mean(out.matrix[,"etha_m_25"])

## sd etha_m_25
sd(out.matrix[,"etha_m_25"])

## The distribution of the posterior probability that 
##a 25 year old man is a recent cannabis users 
pi_m_25 = 1/(1 + exp(-out.matrix[,"etha_m_25"]))
hist(pi_m_25,freq=F)
