
library(truncnorm)
library(cmdstanr)
library(rethinking)
library(bayesplot)
options(mc.cores = parallel::detectCores())

N = 2000
Asim <-abs(rtruncnorm(N,mean = 43.5, sd = 18.6))

Asim<-as.integer(Asim)
dummy <- function(Asim) {if (Asim < 50) {A<-0} else {A<-1}}

hist(Asim)
print(summary(Asim))

A <- lapply(Asim, dummy)
A <- unlist(A)


reaction <- function(A) {if (A > 0) {S <- rbern(1, 0.4)} else {S <- rbern(1, 0.6)}}

S <- lapply(A, reaction)
S <- unlist(S)

dat <- data.frame(A,S)

response <- function(dat) 
{A = dat[1]
S = dat[2]
if( A > 0 & S > 0) {R <- rbern(1, 0.85)}

else if( A > 0 & S < 1 )  {R <- rbern(1, 0.4)} 

else if( A < 1 & S > 0 ) {R <- rbern(1, 0.6)}

else  {R <- rbern(1, 0.3)} #(A < 1 & S < 1)

return(R)

}

R <- apply(dat, 1 ,response)
R <- unlist(R)

dat$R <- R


mean(dat[dat$A == 0,]$R) #P(R = 1|A = 0)

mean(dat[dat$A == 1,]$R) #P(R = 1|A = 1)

mean(dat[dat$A == 0,]$S) #P(S = 1|A = 0)

mean(dat[dat$A == 1,]$S) #P(S = 1|A = 1)

mean(dat[dat$A == 0 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 0)

mean(dat[dat$A == 0 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 0)

mean(dat[dat$A == 1 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 1)

mean(dat[dat$A == 1 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 1)




DAG2 <- "
data {

  int<lower=0> N;
  int R[N];
  int A[N];
  int S[N];
  
}


parameters{

    real aq;
    real bAR;
    real bSR;
    real<lower=0> sigmaq;
    real ap;
    real bAS;
    real<lower=0> sigmap;
    
}
transformed parameters{

  vector[N] p;
  vector[N] q;
  for(i in 1:N) p[i] = inv_logit(ap + bAS*A[i]);
  for(i in 1:N) q[i] = inv_logit(aq + bAR*A[i] + bSR*S[i]);
  
}

model{
    
    //priors
    sigmap ~ exponential( 1 );
    bAS ~ normal(0, 2);
    ap~ normal(0, sigmap);
    
    //likelihood
    S ~ bernoulli(p);
    
    //priors
    sigmaq ~ exponential(1);
    bSR ~ normal(0 , 2 );
    bAR ~ normal(0 , 2);
    aq ~ normal( 0 , sigmaq );
    
    //likelihood
    R ~ bernoulli( q );
}

generated quantities {
  //generate predictions
  real pA0 = inv_logit(ap);  //pA0 = P(S=1|A=0)
  real pA1 = inv_logit(ap + bAS); //pA1 = P(S=1|A=1)
  real qA0S0 = inv_logit(aq); //P(R=1|A=0,S=0)
  real qA0S1 = inv_logit(aq + bSR); //P(R=1|A=0,S=1)
  real qA1S0 = inv_logit(aq + bAR); //P(R=1|A=1,S=0)
  real qA1S1 = inv_logit(aq + bAR + bSR); //P(R=1|A=1,S=1)
}

"

dat = list(R=R,A=A,N=N,S=S)

n_chains <- 2
init <- function() {
  list(aq = abs(rnorm(1, mean = 0, sd = 1.0)),
       bAR = abs(rnorm(1, mean = 0, sd = 2.0)),
       bSR ~ normal(0 , 2 ),
       sigmap ~ exponential( 1 ),
       bAS ~ normal(0, 2),
       ap ~ normal(0, 1)
       
  )
}

file <- file.path("C:/Users/ETay/Documents/DAG1/DAG2.stan")
mod <- cmdstan_model(file)
mod$print()

#write_stan_file(DAG2)
#mod <- cmdstan_model(DAG2)



fit <- mod$sample(
  data = dat, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 0 
)
fit$summary()

