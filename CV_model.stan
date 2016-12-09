# The model specification:
model_string <- "
data {
  int<lower=0> N;
  real y[N];
  real meanY;
  real sdY;
}
transformed data {
  real unifLo;
  real unifHi;
  real normalSigma;
  unifLo = sdY / 100;
  unifHi = sdY * 100;
  normalSigma = sdY * 100;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> mu;
}
model {
  // If we had some prior information on the variance,
  // maybe it would make sence to make it follow a chisq distribution:
  // sigma ~ chi_square(1);

  // For now use a vague sigma prior:
  sigma ~ uniform(unifLo, unifHi);
  mu ~ normal(meanY, normalSigma);
  y ~ lognormal(mu, sigma);
}
generated quantities {
  real cv;
  cv = sigma / mu;
}"

# Make some test data:
library(perbbNovoMisc)
y <- draw_one_sampleCV_data(n_data1 = 30, mu1 = 1, std1 = 1, simdist = 'lognorm')$X1

# Running the model:
library(rstan)
modelfit <- stan(model_code=model_string, data=list(N=length(y), y=y, meanY=mean(y), sdY=sd(y)),
                 pars=c("mu", "sigma", "cv"), chains=3, iter=30000, warmup=10000)

# Exploring run convergence and estimates:
library("shinystan")
my_sso <- launch_shinystan(modelfit)