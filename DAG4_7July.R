
library(truncnorm)
library(cmdstanr)
library(rethinking)
library(bayesplot)
library(tidyverse)
library(posterior)
library(pROC)
options(mc.cores = parallel::detectCores())

N = 4000
Asim <-abs(rtruncnorm(N,mean = 43.5, sd = 18.6))

Asim<-as.integer(Asim)
dummy <- function(Asim) {if (Asim < 50) {A<-0} else {A<-1}}

hist(Asim)
print(summary(Asim))

A <- lapply(Asim, dummy)
A <- unlist(A)


reaction <- function(A) {if (A > 0) {S <- rbinom(1,1,0.4)} else {S <- rbinom(1,1, 0.6)}}

S <- lapply(A, reaction)
S <- unlist(S)
dat <- data.frame(A,S)

response <- function(dat) 
{A = dat[1]
S = dat[2]
if( A > 0 & S > 0) {R <- rbinom(1,1, 0.85)}

else if( A > 0 & S < 1 )  {R <- rbinom(1,1, 0.40)} 

else if( A < 1 & S > 0 ) {R <- rbinom(1, 1, 0.60)}

else  {R <- rbinom(1,1, 0.15)} 

return(R)

}

R <- apply(dat, 1 ,response)
R <- unlist(R)
dat$R <- R

###Simulating M

seek <- function(dat) 
{A = dat[1]
S = dat[2]
if( A > 0 & S > 0) {M <- rbinom(1,1, 0.15)} 

else if( A > 0 & S < 1 )  {M <- rbinom(1,1, 0.05)} 

else if( A < 1 & S > 0 ) {M <- rbinom(1,1, 0.05)}

else  {M <- rbinom(1,1, 0.005)} #(A < 1 & S < 1) 

return(M)

}

M <- apply(dat, 1 ,seek)
M <- unlist(M)
dat$M <- M

reportMA <- function(dat)
{A = dat[1]
S = dat[2]
R = dat[3]
M = dat[4]
if (R > 0 & M > 0 ) {D <- rbinom(1,1,0.999)}# R = 1, M = 1

else if( R > 0 & M < 1 )  {D <- rbinom(1,1,0.001)} #R = 1, M = 0

else  {D <- 0} #R = 0, M = 0

return(D)
}

D <- apply(dat, 1, reportMA)
D <- unlist(D)
dat$D <- D
#####################################################

mean(dat[dat$A == 0,]$S) #P(S = 1|A = 0)

mean(dat[dat$A == 1,]$S) #P(S = 1|A = 1)

mean(dat[dat$A == 0 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 1)

mean(dat[dat$A == 0 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 0)

mean(dat[dat$A == 1 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 1)

mean(dat[dat$A == 1 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 1)

mean(dat[dat$A == 0 & dat$S == 0,]$M) #P(M = 1|S = 0 A = 0)

mean(dat[dat$A == 0 & dat$S == 1,]$M) #P(M = 1|S = 1 A = 0)

mean(dat[dat$A == 1 & dat$S == 0,]$M) #P(M = 1|S = 0 A = 1)

mean(dat[dat$A == 1 & dat$S == 1,]$M) #P(M = 1|A = 1,S = 1)

mean(dat[dat$R == 1 & dat$M == 1,]$D) #P(D = 1|R = 1,M = 1)

dat_R <- dat %>% filter(R == 1)


##Totals of variables
#
mean(dat$A == 1)

mean(dat$S == 1)

mean(dat_R$S == 1)

mean(dat$R == 1)

mean(dat$M == 1)

mean(dat$D == 1)

mean(dat_R$D == 1)


#“int R[N]” use “array[N] int R”

DAG4 <- "
data{
  
  int<lower=0> N;
  array[N] int R;
  array[N] int A;
  array[N] int S;
  array[N] int M;
  array[N] int D;
  
}
parameters{
  real ap;
  real bAR;
  real bSR;
  real<lower=0> sigmap;
  real aq;
  real bAS;
  real<lower=0> sigmaq;
  real ag;
  real bAM;
  real bSM;
  real<lower=0> sigmag;
  real ah;
  real bMD;
 
  real<lower=0> sigmah;
}
transformed parameters{
  vector[N] p;
  vector[N] q;
  vector[N] g;
  vector[N] h;

  for(i in 1:N) p[i] = inv_logit(ap + bAR*A[i] + bSR*S[i]);
  for(i in 1:N) q[i] = inv_logit(aq + bAS*A[i]);
  for(i in 1:N) g[i] = inv_logit(ag + bAM*A[i] + bSM*S[i]);
  for(i in 1:N) h[i] = inv_logit(ah + bMD*M[i]);
}
model{

  //priors
  sigmap ~ exponential(1);
  bSR ~ normal( 0, 2);
  bAR ~ normal( 0, 2);
  ap ~ normal( 0 , sigmap );

  //likelihood
  R ~ bernoulli( p );

  //priors
  sigmaq ~ exponential( 1 ); 
  bAS ~ normal( 0, 2);
  aq ~ normal(0, sigmaq);

  //likelihood
  S ~ bernoulli(q);

  //priors
  sigmag ~ exponential(1);
  bSM ~ normal( 0, 2);
  bAM ~ normal( 0, 2);
  ag ~ normal( 0 , sigmag );

  //likelihood
  M ~ bernoulli( g );

  //priors
  sigmah ~ exponential(1);
  bMD ~ normal( 0, 2);
  ah ~ normal( 0, sigmah );

  //likelihood
  for(i in 1:N) if(R[i] == 1) D[i] ~ bernoulli(h[i]);

}

generated quantities {
  //generate predictions
  real qA0 = inv_logit(aq);  //pA0 = P(S=1|A=0)
  real qA1 = inv_logit(aq + bAS); //pA1 = P(S=1|A=1)
  real pA0S0 = inv_logit(ap); //P(R=1|A=0,S=0)
  real pA0S1 = inv_logit(ap + bSR); //P(R=1|A=0,S=1)
  real pA1S0 = inv_logit(ap + bAR); //P(R=1|A=1,S=0)
  real pA1S1 = inv_logit(ap + bAR + bSR); //P(R=1|A=1,S=1)
  real gA0S0 = inv_logit(ag); //P(R=1|A=0,S=0)
  real gA0S1 = inv_logit(ag + bSM); //P(R=1|A=0,S=1)
  real gA1S0 = inv_logit(ag + bAM); //P(R=1|A=1,S=0)
  real gA1S1 = inv_logit(ag + bAM + bSM); //P(R=1|A=1,S=1)
  real hR1M1 = inv_logit(ah + bMD); //P(D=1|R=1,M=1)

  //prior predictive check

  real<lower=0> sigsim = exponential_rng(1);
  real apsim = normal_rng(0,sigsim);
  real bARsim = normal_rng(0,2);
  real bSRsim = normal_rng(0,2);
  real agsim = normal_rng(0,sigsim);
  real bAMsim = normal_rng(0,2);
  real bSMsim = normal_rng(0,2);
  real ahsim = normal_rng(0,sigsim);
  real bMDsim = normal_rng(0,2);

  vector[N] psim;
  vector[N] Rsim;
  vector[N] gsim;
  vector[N] Msim;
  vector[N] hsim;
  vector[N] Dsim;

  for ( i in 1:N ) {
    psim[i] = apsim + bARsim * A[i] + bSRsim * S[i];
    psim[i] = inv_logit(psim[i]);
    Rsim[i] = bernoulli_rng(psim[i]);


  }
  for ( i in 1:N ) {

    gsim[i] = agsim + bAMsim * A[i] + bSMsim * S[i];
    gsim[i] = inv_logit(gsim[i]);
    Msim[i] = bernoulli_rng(gsim[i]);

  }
  for ( i in 1:N ) {

    hsim[i] = ahsim + bMDsim * S[i];
    hsim[i] = inv_logit(hsim[i]);
    Dsim[i] = bernoulli_rng(hsim[i]);

  }
}

"

dat = list(R=R,A=A,N=N,S=S, M=M, D=D)

n_chains <- 2
init <- function() 
  list(aq = abs(rnorm(1, mean = log(0.5), sd = 1.0)),
       bAR = abs(rnorm(1, mean = 0, sd = 2.0)),
       bSR ~ normal(0 , 2 ),
       sigmap ~ exponential( 1 ),
       bAS ~ normal(0, 2),
       ap ~ normal(log(0.5), 1),
       bSM ~ normal( 0, 2),
       bAM ~ normal( 0, 2),
       ag ~ normal( log(0.5) , 1 ),
       sigmag ~ exponential(1),
       sigmah ~ exponential(1),
       bMD ~ normal( 0, 2),
       ah ~ normal( log(0.5), 1 )
       
  )


file <- file.path("C:/Users/ETay/Documents/DAG1/DAG4_5July.stan")
mod <- cmdstan_model(file)
mod$print()

#write_stan_file(DAG2)
#mod <- cmdstan_model(DAG2)



fit <- mod$sample(
  data = dat, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 2,
  refresh = 0,
  init = list(list(aq = abs(rnorm(1, mean = log(0.5), sd = 1.0))),
         list(bAR = abs(rnorm(1, mean = 0, sd = 2.0))),
         list(bSR ~ normal(0 , 2 )),
         list(sigmap ~ exponential( 1 )),
         list(bAS ~ normal(0, 2)),
         list(ap ~ normal(log(0.5), 1)),
         list(bSM ~ normal( 0, 2)),
         list(bAM ~ normal( 0, 2)),
         list(ag ~ normal( log(0.5) , 1 )),
         list(sigmag ~ exponential(1)),
         list(sigmah ~ exponential(1)),
         list(bMD ~ normal( 0, 2)),
         list(ah ~ normal( log(0.5), 1 ))
         
    )
)



#get posterior draws

draws <- fit$draws("h")
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fit$draws(c("alpha", "beta")))

#plot posterior using bayesplot


mcmc_hist(fit$draws("h"))

fit$print(c("ap", "bAR", "bSR", "aq", "bAS", "ag","ah", "bAM", "bSM", 
              "bMD", 
              "qA0", "qA1","pA0S0","pA0S1","pA1S1","pA1S0",
              "gA0S0", "gA0S1", "gA1S0", "gA1S1", "hR1M1"),max_rows = 21)


h_summ <- summary(draws)
h_summ[,2]
review_h <- cbind(h_summ[,2], dat$D, dat$R, dat$M)
mean(review_h[which(dat$R == 1),1]) #this should be close to 0
mean(review_h[which(dat$M == 1),1]) #this number should be close to 1.0 and R ==1 in the which
mean(review_h[which(dat$D == 1),1])
mean(review_h[which(dat$R == 1 & dat$M ==1),1])
#? how to get the h probability - is it based on R or based on R and M together

mm_draws <- fit$draws("Msim")
Msim2 <- mm_draws[1:4000]
roc(dat$M ~ Msim2, plot = TRUE, print.auc = TRUE)

dd_draws <- fit$draws("Dsim")
Dsim2 <- dd_draws[1:4000]
roc(dat$D ~ Dsim2, plot = TRUE, print.auc = TRUE)
str(dd_draws)
