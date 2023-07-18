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
