############################# Poisson Negative-Binomial Example ################################
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
  n <- length(y)
  log(lambda)*(sum(y) + alpha - 1) - (sum(y)+n/kappa)*log(1+kappa*lambda) - beta*lambda 
}

mh_lam_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for lambda
  lambda_star<-rgamma(1, a, b)
  r <- qlam_k2(lambda_star, y, kappa, alpha, beta) - qlam_k2(lambda, y, kappa, alpha, beta) + log(dgamma(lambda, a, b)) - log(dgamma(lambda_star, a, b))
  return(ifelse(log(runif(1)) <= r, lambda_star, lambda))
}

sim_lambda <- function(lambda, y, alpha, beta, kappa=NULL, k, a, b){
  stopifnot(k %in% c(1,2))
  n <- length(y)
  if(k==1){
    lnew <- rgamma(1, shape=alpha + sum(y), rate=beta + n) 
  }else{
    lnew <- mh_lam_k2(lambda, y, kappa, alpha, beta, a, b)
  }
  return(lnew)
}

### Draw Kappa ###
qkap_k2 <- function(lambda, y, kappa, alpha, beta){
  n <- length(y)
  -n*lgamma(1/kappa) + sum(lgamma(1/kappa + y)) + (sum(y) + alpha - 1)*log(kappa) - beta*kappa - (sum(y) + n/kappa)*log(1+kappa*lambda)
}

mh_kap_k2 <- function(lambda, y, kappa, alpha, beta, a, b){
  # alpha, beta are hyperparameters or priors
  # a, b are parameters of Gamma proposal distribution for kappa
  kappa_star<-rgamma(1, a, b)
  r <- qkap_k2(lambda, y, kappa_star, alpha, beta) - qkap_k2(lambda, y, kappa, alpha, beta) + log(dgamma(kappa, a, b)) - log(dgamma(kappa_star, a, b))
  return(ifelse(log(runif(1)) <= r, kappa_star, kappa))
}

sim_kappa <- function(lambda, y, alpha, beta, kappa, a, b){
  n <- length(y)
  mh_kap_k2(lambda, y, kappa, alpha, beta, a, b)
}

################## RJMCMC SAMPLER ##################
rjmcmc_sampler <- function(y, lambda0, kappa0, k0=1, p=0.5, mu, sigma, alphal, betal, alphak, betak, al, bl, ak, bk, mc.iter = 1000){
  # y is the data
  # lambda0, kappa0, and k0 are initial values
  # p is prior probability of model 1, set to 0.5 as default
  # mu and sigma are fixed parameters of Normal in the between-model step 
  # alphal, betal and alphak, betak are Gamma hyperparameters for priors on lambda and kappa, respectively
  # al, bl and ak, bk are Gamma parameters for the proposal distributions in their respective M-H algorithms
  
  # initialize data frame to save chains
  theta_save <- as.data.frame(matrix(NA, ncol=3, nrow=mc.iter+1))
  names(theta_save)<-c("lambda", "kappa", "k")
  
  # store initial values
  theta_save[1,] <- c(lambda0, kappa0, k0)
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
        kappa_new <- sim_kappa(lambda, y, alphak, betak, theta_new[2], ak, bk)
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, al, bl)
      }else{
        # within model moves
        kappa_new <- NA
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, al, bl)
      }
    }else{
      theta <- c(lambda, kappa)
      theta_new <- lambda
      accept <- log(p) + pk1(y, theta_new, alphal, betal) - log(1-p) - pk2(y, theta, alphal, betal, alphak, betak) + log(dnorm(log(theta[2]/mu), 0, sigma)) - log(theta[2])
      if(log(runif(1)) < accept){
        k <- 1
        # within model moves
        kappa_new <- kappa
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa=NULL, k, al, bl)
      }else{
        # within model moves
        kappa_new <- sim_kappa(lambda, y, alphak, betak, kappa, ak, bk) 
        lambda_new <- sim_lambda(lambda, y, alphal, betal, kappa_new, k, al, bl)
      }
    }
    theta_save[i+1,] <- c(lambda_new, kappa_new, k)
    lambda <- lambda_new
    kappa <- kappa_new
    cat("\r", i)
  }
  return(theta_save)
}

##### Try it out #####
y <- soccer$TotalGoals

# Parameters on priors
alphal <- 25
betal <- 10
alphak <- 1
betak <- 10

# Parameters for between-model step
mu = 0.015
sigma <- 1.5

# Proposal parameters for MH algorithms for lambda and kappa, k=2
al <- 30
bl <- 15
ak <- 2
bk <- 10

# Run MCMC
ptm <- proc.time()
test <- rjmcmc_sampler(soccer$TotalGoals, lambda0=2, kappa0=1, k0=1, mu=mu, sigma=sigma, 
                       alphal=alphal, betal=betal, alphak=alphak, betak=betak, al=al, bl=bl, 
                       ak=ak, bk=bk, mc.iter=50000)
time <- proc.time() - ptm



