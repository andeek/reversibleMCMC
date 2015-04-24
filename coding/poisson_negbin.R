############################# Poisson Negative-Binomial Example ################################
library(ggplot2)

### Data ###
soccer <- read.csv("data/soccer.csv")

####### Set up functions for RJMCMC Sampler #######
### Posterior Distribution Functions ###
pk1 <- function(y, theta, alphal, betal){
  prod(dpois(y, theta))*dgamma(theta, shape=alphal, rate=betal)
}

pk2 <- function(y, theta, alphal, betal, alphak, betak){
  r <- 1/theta[2]
  prod(dnbinom(y, size=r, mean=theta[1]))*dgamma(theta[1], shape=alphal, rate=betal)*dgamma(theta[2], shape=alphak, rate=betak)
}

### Draw Lambda ###
qlam_k2<-function(lambda, y, kappa, alpha, beta){
  lambda^(sum(y) + alpha - 1)*(1+kappa*lambda)^(-sum(y)-1/kappa)*exp(-beta*lambda)
}

mh_lam_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for lambda
  lambda_star<-rgamma(1, a, b)
  r <- min(q_lam_k2(lambda_star, y, kappa, alpha, beta)/q_lam_k2(lambda, y, kappa, alpha, beta)*dgamma(lambda, a, b)/dgamma(lambda_star, a, b),1)
  return(ifelse(runif(1) <= r, z_star, z))
}

sim_lambda <- function(lambda, y, alpha, beta, kappa=NULL, k, a, b){
  stopifnot(k %in% c(1,2))
  n <- length(y)
  if(k=1){
    lnew <- rgamma(1, shape=alpha + sum(y), rate=beta + n) 
  }else{
    lnew <- mh_lam_k2(lambda, y, kappa, alpha, beta, a, b)
  }
  return(lnew)
}

### Draw Kappa ###
qkap_k2 <- function(lambda, y, kappa, alpha, beta){
  n <- length(y)
  gamma(1/kappa)^(-n)*prod(gamma(1/n + y))*kappa^(sum(y) + alpha - 1)*exp(-beta*kappa)*(1+kappa*lambda)^(-(sum(y) + 1/kappa))
}

mh_kap_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for kappa
  kappa_star<-rgamma(1, a, b)
  r <- min(q_kap_k2(lambda_star, y, kappa, alpha, beta)/q_kap_k2(lambda, y, kappa, alpha, beta)*dgamma(kappa, a, b)/dgamma(kappa_star, a, b),1)
  return(ifelse(runif(1) <= r, z_star, z))
}

sim_kappa <- function(lambda, y, alpha, beta, kappa, a, b){
  stopifnot(k == 2)
  n <- length(y)
  mh_kap_k2(lambda, y, kappa, alpha, beta, a, b)
}

################## RJMCMC SAMPLER ################
gibbs_sampler <- function(y, lambda0, kappa0, k0, p=0.5, mu, sigma, alphal, betal, alphak, betak, al, bl, ak, bk, mc.iter = 1000){
  # y is the data
  # lambda0, kappa0, and k0 are initial values
  # p is prior probability of model 1, set to 0.5 as default
  # mu and sigma are fixed parameters of Normal in the between-model step 
  # alphal, betal and alphak, betak are Gamma hyperparameters for priors on lambda and kappa, respectively
  # al, bl and ak, bk are Gamma parameters for the proposal distributions in their respective M-H algorithms
  
  
  
    # RJMCMC step to move between models
    if{k=1}{
      u <- rnorm(1, 0, sigma)
      theta <- lambda
      theta_new <- c(lambda, mu*exp(u))
      accept <- min((1-p)*pk2(y, theta_new, alphal, betal, alphak, betak)/(p*pk1(y, theta, alphal, betal))*1/dnorm(u, 0, sigma)*theta_new[2],1)
      if(runif(1) < accept){
        k <- 2
        kappa_new <- theta_new[2]
        lambda_new <- sim_lambda(lambda, y, kappa_new, alphal, betal, k, al, bl)
      }else{
        kappa_new <- kappa
        lambda_new <- sim_lambda(lambda, y, alphal, betal, k, al, bl){
      }
    }else if{k=2}{
      theta <- c(lambda, kappa)
      theta_new <- lambda
      accept <- min(p*pk1(y, theta_new, alphal, betal)/((1-p)*pk2(y, theta, alphal, betal, alphak, betak))*dnorm(log(theta[2]/mu), 0, sigma)*1/theta[2],1)
      if(runif(1) < accept){
        k <- 1
        kappa_new <- kappa
        lambda_new <- sim_lambda(lambda, y, alphal, betal, k, al, bl)
      }else{
        kappa_new <- sim_kappa(lambda, y, alphak, betak, kappa, ak, bk)
        lambda_new <- sim_lambda(lambda, y, kappa_new, alphal, betal, k, al, bl)
      }
    }
  
}

