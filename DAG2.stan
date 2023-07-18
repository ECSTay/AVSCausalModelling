
data {

  int<lower=0> N;
  array[N] int R;
  array[N] int A;
  array[N] int S;
  
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



