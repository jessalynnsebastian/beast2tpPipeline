// This model is adapted from Isaac Goldstein, see:
// https://github.com/igoldsteinh/kopanyo_tp_code/blob/main/R/misclass_logistic_regression.stan

//This model is for logistic regression with misclassified responses on some covariates
//here z is the observed response, possibly misclassified, y is the true response
//we assume there is no dependence between misclassification and the covariate
//we assume that sensitivity and specificity of misclassification are known
data {
  //number of observations
  int<lower=0> N;
  //number of covariates
  int<lower=0> K;
  //misclassified response
  int<lower=0, upper=1> z[N];
  //covariates (INCLUDING intercept)
  matrix[N,K+1] X;
  //prob(z=0 | y=0)
  real<lower=0> specificity;
  //prob(z = 1 | y=1)
  real<lower=0> sensitivity;
}

parameters {
  vector[K+1] beta;
}

model {
  vector[N]  prob_y1;
  vector[N]  prob_y0;
  vector[N]  p;
  
  // Priors on coefficients - std normal
  beta ~ std_normal();
  
  for (n in 1:N) {
    prob_y1[n] = inv_logit(dot_product(beta, X[n,]));
    prob_y0[n] = 1-prob_y1[n];
    p[n] = sensitivity*prob_y1[n] + (1-specificity)*prob_y0[n];

  }

  z ~ bernoulli(p);
}
