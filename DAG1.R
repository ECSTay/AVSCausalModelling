

library(truncnorm)
library(cmdstanr)
library(rethinking)
library(bayesplot)
options(mc.cores = parallel::detectCores())

N = 2000
Asim <- rtruncnorm(N,mean = 50, sd = 15)
Asim<-as.integer(Asim)
dummy <- function(Asim) {if (Asim < 50) {A<- 0} else {A<-1}}

hist(Asim)
print(summary(Asim))

A <- lapply(Asim, dummy)
A <- unlist(A)


response <- function(A) {if (A > 0) {RP <- rbern(1, 0.7)} else {RP <- rbern(1, 0.3)}}

RP <- lapply(A, response)
R <- unlist(RP)
dat <- data.frame(A,R)


mean(dat[dat$A == 0,]$R) #P(R = 1|A = 0)

mean(dat[dat$A == 1,]$R) #P(R = 1|A = 1)


DAG1 <- "
data {
  int<lower=0> N;
  int R[N];
  int A[N];
}


parameters {
  real a;
  real bAR;
  real<lower=0> sigma;
}

transformed parameters{
  vector[N] p;

  for(i in 1:N) p[i] = inv_logit(a + bAR*A[i]);
 
}

model {
  sigma ~ exponential( 1 );
    bAR ~ normal( 0, 2);
    a ~ normal( 0.36 , sigma );
   
    R ~ bernoulli( p );
}

"

dat = list(R=R,A=A,N=N)

n_chains <- 2
init <- function() {
  list(a = abs(rnorm(1, mean = 0, sd = 1.0)),
       bAR = abs(rnorm(1, mean = 0, sd = 2.0)),
      
  )
}

file <- file.path("C:/Users/ETay/Documents/DAG1/DAG1.stan")
mod <- cmdstan_model(file)
mod$print()

#write_stan_file(DAG1)
#mod <- cmdstan_model(DAG1)
#mod <- cmdstan_model("C:/Users/ETay/AppData/Local/Temp/Rtmp4KOlUx/model_59c9f36c0697b9d0d82b4762c26f9e44.stan")

fit <- mod$sample(
  data = dat, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 0 
)
fit$summary()


