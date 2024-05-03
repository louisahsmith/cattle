data {
  int<lower=0> J1; // indices with only one observation per study
  int<lower=0> J2; // indices with multiple observations per study
  array[J1] int<lower=1, upper=J1 + J2> ii1;
  array[J2] int<lower=1, upper=J1 + J2> ii2;
  array[J1] real<lower=0> y1;      // observed PFOS
  array[J2] real<lower=0> alpha2;  // observed mean PFOS
  array[J2] real<lower=0> sigma2;  // standard deviations of the PFOS measurements
}

transformed data {
  int<lower=0> J = J1 + J2;
}

parameters {
  real<lower=0> beta;               // overall mean
  real<lower=0> tau;       // deviation of study-specific means
  array[J2] real<lower=0> y2;       //
  array[J1] real<lower=0> alpha1;   // 
  array[J1] real<lower=0> sigma1;   //
}

transformed parameters {
  array[J] real<lower=0> alpha;
  array[J] real<lower=0> sigma;
  array[J] real<lower=0> y;
  alpha[ii1] = alpha1;
  alpha[ii2] = alpha2;
  sigma[ii1] = sigma1;
  sigma[ii2] = sigma2;
  y[ii1] = y1;
  y[ii2] = y2;
}

model {
  y ~ normal(alpha, sigma);
  alpha ~ normal(beta, tau);
  sigma ~ cauchy(0, 5);
  beta ~ normal(0, 10);
  tau ~ cauchy(0, 5);
}
