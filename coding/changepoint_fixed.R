#### Libraries ####
library(ggplot2)
library(tidyr)
library(dplyr)

###### Try with fixed (true) steps ######
truefunction<-function(x){
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  temp <- 0
  for(i in 1:11) {
    temp <- temp + h[i]/2 * (1 + sign(x - t[i]))
  }
  return(temp)
}
## true and noisy data ##
n<-512
x<-(0:(n-1))/n
f<-truefunction(x)
set.seed(0401)
y=f+rnorm(f)/3
qplot(x=x, y=y, geom="point") + geom_line(aes(y=f), colour="red")

#### Number of steps = 11 ####
B <- 12 # B is number of steps as in HW 7 description
# i.e. need to estimate b_1,...,b_{B-1}, and f_1,...,f_B

### function to bin data ###
bsets <- function(b, x, y){
  m <- findInterval(x, b)
  y.sets <- split(y, m)
  n.sets <- as.numeric(sapply(y.sets, length))
  return(list(y.sets=y.sets, n.sets=n.sets))
}

### posterior distribution function ###
posterior <- function(f, sets, b, B=12, N=512, sigma, kappa0, a0, b0){
  # f is vector of B estimated heights
  # sets is output from bsets function
  # b is a vector of changepoints
  # B is total number of changepoints
  # N is number of observations
  # sigma^2 is variance of data normal
  # sigma^2/kappa0 is variance of normal priors on f
  # a0, b0 are parameters of inv-gamma prior on sigma
   
  ys <- sets$y.sets
  ns <- sets$n.sets
  sse <- sapply(1:B, function(i) sum((ys[[i]]-f[i])^2))
  tau <- sigma/sqrt(kappa0)
  
  # loglikelihood
  - sum(ns)*log(sqrt(2*pi)*sigma) - 1/(2*sigma^2)*sum(sse) - B*log(sqrt(2*pi*tau)) - 
    1/(2*tau^2)*sum(f^2) + lfactorial(N-B) - lfactorial(N) - log(dgamma(sigma, a0, b0))
}

## full conditional for f_j ##
sample_fj <- function(y, n, sigma, kappa0){
  # y is y's in j_th segment
  # n is num y's in j_th segment
  # sigma is data variance
  # kappa is variance scaling factor in f_j prior
  
  v <- sigma^2/(n + kappa0)
  m <- v*n/sigma^2*mean(y)
  rnorm(1, m, sqrt(v))
}

## full conditional for sigma ##
sample_sigma <- function(y.sets, ns, a0, b0, kappa0, N=512, B=12){
  # y.sets are observation sets
  # ns number of observations in each set
  # a0, b0 are hyperparameters of inv-gamma
  # kappa0 is variance scaling factor in f_j prior
  
  ssy <- 1/2*sum(sapply(y.sets, function(u) sum((u-mean(u))^2)))
  bm <- ns*kappa0/(2*(ns+kappa0))
  ssb <- sum(bm*sapply(y.sets, function(u) mean(u^2)))
  
  a1 <- a0 + N/2
  b1 <- b0 + ssy + ssb
  1/rgamma(1, a1, b1)
}

move <- function(b, b.sets, x, y, B=12, N=512, sigma, kappa0, a0, b0){
  pos <- which(x %in% b)
  b_star <- b
  i <- sample(1:length(b), 1)
  big_pos <- c(2, pos, N-1)
  int <- big_pos[c(i, i+2)]
  int2 <- c(int[1]+1, int[2]-1)
  i_star <- sample(int2[1]:int2[2], 1)
  b_new <- x[i_star]
  b_star[i] <- b_new
  b_star <- sort(b_star)
  new.sets <- bsets(b_star, x, y)
  accept <- posterior(f, new.sets, b_star, B=B, N=N, sigma, kappa0, a0,b0) - posterior(f, b.sets, b, B=B, N=N, sigma, kappa0, a0, b0)
  if(log(runif(1)) < accept){
    return(list(b=b_star, sets=new.sets))
  }else{
    return(list(b=b, sets=b.sets))
  }
}

move2 <- function(b, b.sets, x, y, B=12, N=512, sigma, kappa0, a0, b0){
  pos <- which(x %in% b)
  b_star <- b
  i <- sample(1:length(b), 1)
  N2 <- c(2:(N-1))[!2:(N-1) %in% pos]
  i_star <- sample(N2, 1)
  b_new <- x[i_star]
  b_star[i] <- b_new
  b_star <- sort(b_star)
  new.sets <- bsets(b_star, x, y)
  accept <- posterior(f, new.sets, b_star, B=B, N=N, sigma, kappa0, a0,b0) - posterior(f, b.sets, b, B=B, N=N, sigma, kappa0, a0, b0)
  if(log(runif(1)) < accept){
    return(list(b=b_star, sets=new.sets))
  }else{
    return(list(b=b, sets=b.sets))
  }
}



####### Sampler #######
fixed_sampler <- function(x, y, b, B, f, N, sigma, kappa0, a0, b0, mc.iter){
  b_save <- matrix(NA, ncol=length(b), nrow=mc.iter+1)
  f_save <- matrix(NA, ncol=B, nrow=mc.iter+1)
  sigma_save <- rep(NA, nrow=mc.iter)
  pos_save <- matrix(NA, ncol=length(b), nrow=mc.iter+1)

  b_save[1,] <- b
  f_save[1,] <- f
  sigma_save[1] <- sigma
  pos_save[1,] <- which(x %in% b)
  
  for(i in 1:mc.iter){
    b.sets <- bsets(b, x, y)
    
    ## check for move 
    check_newpos <- move(b, b.sets, x, y, B=B, N=N, sigma, kappa0, a0, b0)
    sets_new <- check_newpos$sets
    b_new <- check_newpos$b

    ## simulate fs
    f_new <- sapply(1:B, function(i) sample_fj(sets_new$y.sets[[i]],sets_new$n.sets[i], sigma, kappa0))
    
    ## simulate sigma
    sigma_new <- sample_sigma(sets_new$y.sets, sets_new$n.sets, a0, b0, kappa0, N=N, B=B)
      
    b_save[i+1,] <- b <- b_new
    f_save[i+1,] <- f <- f_new
    sigma_save[i+1] <- sigma <- sigma_new
    pos_save[i+1,] <- which(x %in% b) 
    cat("\r",i)
  }
  return(list(b=b_save, f=f_save, sigma=sigma_save, pos=pos_save))
}
fixed_sampler2 <- function(x, y, b, B, f, N, sigma, kappa0, a0, b0, mc.iter){
  b_save <- matrix(NA, ncol=length(b), nrow=mc.iter+1)
  f_save <- matrix(NA, ncol=B, nrow=mc.iter+1)
  sigma_save <- rep(NA, nrow=mc.iter)
  pos_save <- matrix(NA, ncol=length(b), nrow=mc.iter+1)
  
  b_save[1,] <- b
  f_save[1,] <- f
  sigma_save[1] <- sigma
  pos_save[1,] <- which(x %in% b)
  
  for(i in 1:mc.iter){
    b.sets <- bsets(b, x, y)
    
    ## check for move 
    check_newpos <- move2(b, b.sets, x, y, B=B, N=N, sigma, kappa0, a0, b0)
    sets_new <- check_newpos$sets
    b_new <- check_newpos$b
    
    ## simulate fs
    f_new <- sapply(1:B, function(i) sample_fj(sets_new$y.sets[[i]],sets_new$n.sets[i], sigma, kappa0))
    
    ## simulate sigma
    sigma_new <- sample_sigma(sets_new$y.sets, sets_new$n.sets, a0, b0, kappa0, N=N, B=B)
    
    b_save[i+1,] <- b <- b_new
    f_save[i+1,] <- f <- f_new
    sigma_save[i+1] <- sigma <- sigma_new
    pos_save[i+1,] <- which(x %in% b) 
    cat("\r",i)
  }
  return(list(b=b_save, f=f_save, sigma=sigma_save, pos=pos_save))
}

b.init <- sort(sample(x[-c(1,512)], 11))
f.init <- rnorm(12,0,2)
sigma.init <- 1
kappa0 <- 1
a0 <- 2
b0 <- 0.01

test<-fixed_sampler(x=x, y=y, b=b.init, B=12, f=f.init, N=512, 
              sigma=sigma.init, kappa0=kappa0, a0=a0, b0=b0, mc.iter=10000)

test2<-fixed_sampler2(x=x, y=y, b=b.init, B=12, f=f.init, N=512, 
                    sigma=sigma.init, kappa0=kappa0, a0=a0, b0=b0, mc.iter=50000)


qplot(x=1:50001, y=test2$sigma, geom="line")

bs <- test2$b[-c(1:5000),]
bs.m <- bs %>% data.frame %>% gather()
fs <- test2$f[-c(1:5000),]
fs.m <- fs %>% data.frame %>% gather()
qplot(data=bs.m, x=value) + facet_wrap(~key, ncol=3)
qplot(data=fs.m, x=value) + facet_wrap(~key, ncol=3)

fhats <- apply(fs, 2, median)
bhats <- apply(bs, 2, median)

qplot(x=x, y=y, geom="point") + geom_line(aes(y=f), colour="red") + geom_vline(aes(xintercept=bhats), colour="blue")

qplot(x=x, y=y, geom="point") + geom_line(aes(y=f), colour="red") + geom_hline(aes(yintercept=fhats), colour="blue")

setshat <- bsets(bhats, x,y)
yhats <- sapply(1:12, function(i) setshat$y.sets[[i]] + fhats[i] - setshat$y.sets[[i]])
qplot(x=x, y=y, geom="point") + geom_line(aes(y=f), colour="red") + geom_line(aes(y=unlist(yhats)), colour="blue")

# disaster number function m(.)
count.disasters <- function(s1, s2) {
  return(sum(coal >= s1 & coal < s2))
}

k0 <- 1
s0 <- c(1851,1963)
h0 <- c(1.7)

lambda = 3

# value of Gamma(a,b)
a=1
b=200/365.24


# accept probability function for Height move

alpha.hmove <- function(s,h,h_tilde,j) {
  m <- count.disasters(s[j],s[j+1])
  log_likelihood_ratio <- (h[j] - h_tilde) * (s[j+1] - s[j]) + m * (log(h_tilde) - log(h[j]))
  log_prior_ratio <- a * (log(h_tilde) - log(h[j])) - b * (h_tilde - h[j])
  log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
  return(exp(log_pi_ratio))
  
}


# accept probability function for Position move

alpha.smove <- function(s,h,s_tilde,j) {
  m1 <- count.disasters(s[j-1],s[j])
  m1_tilde <- count.disasters(s[j-1],s_tilde)
  m2 <- count.disasters(s[j],s[j+1])
  m2_tilde <- count.disasters(s_tilde,s[j+1])
  log_likelihood_ratio <- (h[j] - h[j-1]) * (s_tilde - s[j]) + (m1_tilde - m1) * log(h[j-1]) + (m2_tilde - m2) * log(h[j])
  log_prior_ratio <- log(s_tilde - s[j-1]) + log(s[j+1] - s_tilde) - log(s[j] - s[j-1]) - log(s[j+1] - s[j])
  log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
  return(exp(log_pi_ratio))
  
}


# accept probability function for birth of step

alpha.birth <- function(s,h,s_star,h1_prime,h2_prime,j,k) {
  m1_prime <- count.disasters(s[j],s_star)
  m2_prime <- count.disasters(s_star,s[j+1])
  mj <- m1_prime + m2_prime
  log_likelihood_ratio <- (
    - h1_prime * (s_star - s[j]) - h2_prime * (s[j+1] - s_star)
    + m1_prime * log(h1_prime) + m2_prime * log(h2_prime)
    + h[j] * (s[j+1] - s[j]) - mj * log(h[j])
  )
  log_prior_ratio <- (
    log(2) + log(lambda) + log(2*k+1)
    - log(s[k+1] - s[1]) + log(s_star - s[j]) + log (s[j+1] - s_star) - log(s[j+1] - s[j])
    + a * log(b) - log(gamma(a)) + (a-1) * (log(h1_prime) + log(h2_prime)) - a * log(h[j])
    - b * (h1_prime + h2_prime - h[j])
  )
  log_proposal_ratio <- log(death(k)) - log(birth(k-1)) - log(k)
  log_Jacobian <- 2 * log(h1_prime + h2_prime)
  log_pi_ratio <- log_likelihood_ratio + log_prior_ratio + log_proposal_ratio + log_Jacobian
  return(exp(log_pi_ratio))
}



# accept probability function for death of step

alpha.death <- function(s,h,hj_prime,j,k) {
  m1 <- count.disasters(s[j],s[j+1])
  m2 <- count.disasters(s[j+1],s[j+2])
  mj_prime <- m1 + m2
  log_likelihood_ratio <- (
    - hj_prime * (s[j+2] - s[j]) + mj_prime * log(hj_prime)
    + h[j] * (s[j+1] - s[j]) + h[j+1] * (s[j+2] - s[j+1])
    - m1 * log(h[j]) - m2 * log(h[j+1])
  )
  log_prior_ratio <- (
    - log(2) - log (lambda) - log(2*k-1)
    + log(s[k+1] - s[1]) - log(s[j+1] - s[j]) - log (s[j+2] - s[j+1]) + log(s[j+2] - s[j])
    - a * log(b) + log(gamma(a)) - (a-1) * (log(h[j]) + log(h[j+1])) + a * log(hj_prime)
    + b * (h[j] + h[j+1] - hj_prime)
  )
  log_proposal_ratio <- + log(birth(k-2)) + log(k-1) -log(death(k-1))
  log_Jacobian <- - 2 * log(h[j] + h[j+1])
  log_pi_ratio <- log_likelihood_ratio + log_prior_ratio + log_proposal_ratio + log_Jacobian
  return(exp(log_pi_ratio))
  
}


birth <- function(change_points){return(3.6/7 * min(1,lambda/(change_points+1)))}
death <- function(change_points){return(3.6/7 * min(1,change_points/lambda))}

Move <- function(n, h, s, k) {
  H <- list()
  S <- list()
  K <- vector(length=n)
  K[1] <- k
  H[[1]] <- h
  S[[1]] <- s
  
  for (i in 2:n){
    position_prob <- ifelse(k<=1, 0, 0.5 * (1 - birth(k-1) - death(k-1)))
    #if (k<=1) {position_prob <- 0} else {position_prob <- 0.5 * (1 - birth(k-1) - death(k-1))}
    
    height_prob <- 1 - birth(k-1) - death(k-1) - position_prob
    
    type <- runif(1)
    
    if (type>1-height_prob) {                       #from 1 to k is for Height h_1 to h_k
      j <- sample(1:k,size=1)
      u <- runif(1,-0.5,0.5)
      h_tilde <- h[j] * exp(u)
      U <- runif (1)
      if (U < alpha.hmove(s,h,h_tilde,j)) {
        h[j] <- h_tilde
      }
    }
    
    if (type<=1-height_prob && type>1-height_prob-position_prob) {                     #from 2 to k is for Position s_2 to s_k
      j <- sample(1:(k-1), size=1) + 1
      s_tilde <- runif(1,s[j-1],s[j+1])
      U <- runif (1)
      if (U < alpha.smove(s,h,s_tilde,j)) {
        s[j] <- s_tilde
      }
    }
    
    if (type>=birth(k-1) && type<=birth(k-1)+death(k-1)) {                    #from 1 to k-1 is for death of steps d_2 to d_k
      j <- sample(1:(k-1), size=1)
      r <- (s[j+2] - s[j+1]) / (s[j+2] - s[j])
      hj_prime <- h[j]^(1-r) * h[j+1]^r                    #exp((1-r) * log(h[j]) + r * log(h[j+1]))
      U <- runif(1)
      if (U < alpha.death(s,h,hj_prime,j,k)){
        k <- k - 1
        h[j] <- hj_prime
        h <- h[-(j+1)]
        s <- s[-(j+1)]
      }
    }
    
    if (type<=birth(k-1)){                                             #from 1 to k is for birth of steps b_1 to b_k
      j <- sample(1:k,size=1)
      s_star <- runif(1,s[j],s[j+1])
      u <- runif(1)
      r <- exp(log(s[j+1] - s_star) - log(s[j+1] - s[j]))
      h1_prime = h[j] * exp(r * (log(u) - log(1-u)))
      h2_prime = h[j] * exp((1-r) * (log(1-u) - log(u)))
      U <- runif(1)
      if (U < alpha.birth(s,h,s_star,h1_prime,h2_prime,j,k)){
        s <- c(s[1:j], s_star, s[(j+1):(k+1)])
        if (j > 1) {
          left <- h[1:(j-1)]
        } else {
          left <- c()
        }
        if (j < k) {
          right <- h[(j+1):k]
        } else {
          right <- c()
        }
        h <- c(left, h1_prime, h2_prime, right)
        k <- k + 1
      }
    }
    
    K[i] <- k
    H[[i]] <- h
    S[[i]] <- s
  } 
  return(list(h=H,k=K,s=S))
}