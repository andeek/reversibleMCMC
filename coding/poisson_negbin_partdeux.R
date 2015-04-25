############################# Poisson Negative-Binomial Example ################################
library(ggplot2)

### Data ###
soccer <- read.csv("data/soccer.csv")

####### Set up functions for RJMCMC Sampler #######
### Posterior Distribution Functions ###
pk1 <- function(y, theta, alphal, betal){
  sum(log(dpois(y, theta))) + log(dgamma(theta, shape=alphal, rate=betal))
}

pk2 <- function(y, theta, alphal, betal, alphak, betak){
  r <- 1/theta[2]
  sum(log(dnbinom(y, size=r, mu=theta[1]))) + log(dgamma(theta[1], shape=alphal, rate=betal)) + log(dgamma(theta[2], shape=alphak, rate=betak))
}

### Draw Lambda ###
qlam_k2<-function(lambda, y, kappa, alpha, beta){
  log(lambda)*(sum(y) + alpha - 1) - (sum(y)+1/kappa)*log(1+kappa*lambda) - beta*lambda
}

mh_lam_k2 <- function(lambda, y, kappa, alpha, beta, sl){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for lambda
  a <- lambda^2/sl
  b <- lambda/sl
  lambda_star<-rgamma(1, a, b)
  as <- lambda_star^2/sl
  bs <- lambda_star/sl
  r <- qlam_k2(lambda_star, y, kappa, alpha, beta) - qlam_k2(lambda, y, kappa, alpha, beta) + log(dgamma(lambda, as, bs)) - log(dgamma(lambda_star, a, b))
  return(ifelse(log(runif(1)) <= r, lambda_star, lambda))
}

sim_lambda <- function(lambda, y, alpha, beta, kappa=NULL, k, sl){
  stopifnot(k %in% c(1,2))
  n <- length(y)
  if(k==1){
    lnew <- rgamma(1, shape=alpha + sum(y), rate=beta + n) 
  }else{
    lnew <- mh_lam_k2(lambda, y, kappa, alpha, beta, sl)
  }
  return(lnew)
}

### Draw Kappa ###
qkap_k2 <- function(lambda, y, kappa, alpha, beta){
  n <- length(y)
  -n*lgamma(1/kappa) + sum(lgamma(1/kappa + y)) + (sum(y) + alpha - 1)*log(kappa) - beta*kappa - (sum(y) + 1/kappa)*log(1+kappa*lambda)
}

mh_kap_k2 <- function(lambda, y, kappa, alpha, beta, sk){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for kappa
  a <- kappa^2/sk
  b <- kappa/sk
  kappa_star<-rgamma(1, a, b)
  as <- kappa_star^2/sk
  bs <- kappa_star/sk
  r <- qkap_k2(lambda, y, kappa_star, alpha, beta) - qkap_k2(lambda, y, kappa, alpha, beta) + log(dgamma(kappa, as, bs)) - log(dgamma(kappa_star, a, b))
  return(ifelse(log(runif(1)) <= r, kappa_star, kappa))
}

sim_kappa <- function(lambda, y, alpha, beta, kappa, sk){
  n <- length(y)
  mh_kap_k2(lambda, y, kappa, alpha, beta, sk)
}

################## RJMCMC SAMPLER ################
rjmcmc_sampler <- function(y, lambda0, kappa0, k0=1, p=0.5, mu, sigma, alphal, betal, 
                           alphak, betak, sl, sk, tune=TRUE, mc.iter = 1000){
  # y is the data
  # lambda0, kappa0, and k0 are initial values
  # p is prior probability of model 1, set to 0.5 as default
  # mu and sigma are fixed parameters of Normal in the between-model step 
  # alphal, betal and alphak, betak are Gamma hyperparameters for priors on lambda and kappa, respectively
  # al, bl and ak, bk are Gamma parameters for the proposal distributions in their respective M-H algorithms
  
  # initialize data frame to save chains
  theta_save <- as.data.frame(matrix(NA, ncol=3, nrow=mc.iter+1))
  names(theta_save)<-c("lambda", "kappa", "k")
  sigmas <- as.data.frame(matrix(NA, ncol=2, nrow=mc.iter+1))
  names(sigmas) <- c("lambda", "kappa")
  
  # store initial values
  theta_save[1,] <- c(lambda0, kappa0, k0)
  sigmas[1,] <- c(sl, sk)
  lambda <- lambda0
  kappa <- kappa0
  k <- k0
  
  for(i in 1:mc.iter){
    # Reversible JUMP MCMC to move between models
    if(k==1){
      u <- rnorm(1, 0, sigma)
      theta <- lambda
      theta_new <- c(lambda, mu*exp(u))
      accept <- log(1-p) + pk2(y, theta_new, alphal, betal, alphak, betak) - log(p) - pk1(y, theta, alphal, betal) - log(dnorm(u, 0, sigma)) + log(theta_new[2])
      if(log(runif(1)) < accept){
        k <- 2
        # within model moves
        kappa_new <- sim_kappa(lambda, y, alphak, betak, theta_new[2], sk)
        if(tune){
          if(kappa_new==theta_new[2]){
            sk <- sk/1.1
          }else{
            sk <- sk*1.1
          }
        }
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, sl)
        if(tune){
          if(lambda_new==lambda){
            sl <- sl/1.1
          }else{
            sl <- sl*1.1
          }
        }
      }else{
        # within model moves
        kappa_new <- NA
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, sl)
        if(tune){
          if(lambda_new==lambda){
            sl <- sl/1.1
          }else{
            sl <- sl*1.1
          }
        }
      }
    }else{
      theta <- c(lambda, kappa)
      theta_new <- lambda
      accept <- log(p) + pk1(y, theta_new, alphal, betal) - log(1-p) - pk2(y, theta, alphal, betal, alphak, betak) + log(dnorm(log(theta[2]/mu), 0, sigma)) - log(theta[2])
      if(log(runif(1)) < accept){
        k <- 1
        # within model moves
        kappa_new <- NA
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, sl)
        if(tune){
          if(lambda_new==lambda){
            sl <- sl/1.1
          }else{
            sl <- sl*1.1
          }
        }
      }else{
        # within model moves
        kappa_new <- sim_kappa(lambda, y, alphak, betak, kappa, sk) 
        if(tune){
          if(kappa_new==kappa){
            sk <- sk/1.1
          }else{
            sk <- sk*1.1
          }
        }
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, sl)
        if(tune){
          if(lambda_new==lambda){
            sl <- sl/1.1
          }else{
            sl <- sl*1.1
          }
        }
      }
    }
    theta_save[i+1,] <- c(lambda_new, kappa_new, k)
    sigmas[i+1,] <- c(sl, sk)
    lambda <- lambda_new
    kappa <- kappa_new
    cat("\r", i)
  }
  return(theta_save)
}

##### Try it out #####
y <- soccer$TotalGoals
alphal <- 25
betal <- 10
alphak <- 1
betak <- 10
mu = 0.015
sigma <- 1.5
sl <- 10
sk <- 1

test <- rjmcmc_sampler(soccer$TotalGoals, lambda0=2, kappa0=1, k0=1, mu=mu, sigma=sigma, 
                       alphal=alphal, betal=betal, alphak=alphak, betak=betak, sl=sl, sk=sk, 
                       tune=TRUE, mc.iter=10000)

testk1 <- subset(test, k==1)
testk2 <- subset(test, k==2)


qplot(data=testk1, x=lambda, geom="density")
qplot(data=testk2, x=lambda, geom="density")
qplot(data=testk2, x=kappa, geom="density")

qplot(x=1:nrow(testk2), testk2$kappa, geom="line")
