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
  //for(i in 1:N) h[i] = inv_logit(ah + bMD*M[i]);
  for(i in 1:N) h[i] = R[i]*(inv_logit(ah + bMD*M[i]));
 
}
model{

  //priors
  sigmap ~ exponential(1);
  bSR ~ normal( 0, 2);
  bAR ~ normal( 0, 2);
  ap ~ normal( log(0.5) , sigmap );

  //likelihood
  R ~ bernoulli( p );

  //priors
  sigmaq ~ exponential( 1 ); 
  bAS ~ normal( 0, 2);
  aq ~ normal(log(0.5), sigmaq);

  //likelihood
  S ~ bernoulli(q);

  //priors
  sigmag ~ exponential(1);
  bSM ~ normal( 0, 2);
  bAM ~ normal( 0, 2);
  ag ~ normal( log(0.5) , sigmag );

  //likelihood
  M ~ bernoulli( g );

  //priors
  sigmah ~ exponential(1);
  bMD ~ normal( 0, 2);
  ah ~ normal( log(0.5), sigmah );

  //likelihood
  
  //for(i in 1:N) D[i] ~ bernoulli(h[i]); try this one as well ( if gives lower value for M = 0 )
  
  
  for(i in 1:N) D[i] ~ bernoulli(h[i]);
  
 // for(i in 1:N) D[i] ~ bernoulli(h[i]);
  
 // R[i] = 1 ? D[i]~bernoulli(h[i]) : h[i] = 0; #? try this one as well

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
  real apsim = normal_rng(log(0.5),sigsim);
  real bARsim = normal_rng(0,2);
  real bSRsim = normal_rng(0,2);
  real agsim = normal_rng(log(0.5),sigsim);
  real bAMsim = normal_rng(0,2);
  real bSMsim = normal_rng(0,2);
  real ahsim = normal_rng(log(0.5),sigsim);
  real bMDsim = normal_rng(0,2);

  vector[N] psim;
  vector[N] Rsim;
  vector[N] gsim;
  vector[N] Msim;
  vector[N] hsim;
  vector[N] Dsim;

  for ( i in 1:N ) {
    psim[i] = inv_logit(apsim + bARsim * A[i] + bSRsim * S[i]);
    Rsim[i] = bernoulli_rng(psim[i]);


  }
  for ( i in 1:N ) {

    gsim[i] = inv_logit(agsim + bAMsim * A[i] + bSMsim * S[i]);
    Msim[i] = bernoulli_rng(gsim[i]);

  }
  for ( i in 1:N ) {

    hsim[i] = R[i]*(inv_logit(ahsim + bMDsim * S[i]));
    Dsim[i] = bernoulli_rng(hsim[i]);

  }
}
