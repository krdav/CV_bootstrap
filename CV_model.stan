# The model specification:
model_string <- "
data {
  int<lower=0> N;
  real y[N];
  real<lower=0> unifLo;
  real<lower=0> unifHi;
}
parameters {
  real<lower=0> sigma_log;
  real mu_log;
}
model {
  // Use vague priors:
  sigma_log ~ uniform(unifLo, unifHi);
  mu_log ~ uniform(-unifHi, unifHi);
  y ~ lognormal(mu_log, sigma_log);
}
generated quantities {
  // Generate mean and variance
  // from the lognorm variables to calculate the CV:
  real<lower=0>  cv;
  real<lower=0>  mu;
  real<lower=0>  sigma;
  mu = exp(mu_log + sigma_log / 2);
  sigma = exp(sigma_log - 1) * exp(2 * mu_log + sigma_log);
  // cv = (1.0 + 1.0 / (4.0 * N)) * sigma / mu; // Same as below
  cv = (1.0 + 1.0 / (4.0 * N)) * sqrt(exp(sigma_log) - 1);
}"

# Make some test data:
library(perbbNovoMisc)
y <- draw_one_sampleCV_data(n_data1 = 100, mu1 = 1, std1 = 1, simdist = 'lognorm')$X1

# Running the model:
library(rstan)
meanY <- mean(y)
sdY <- sd(y)
meanlogY <- log(meanY) - 0.5 * log(sdY^2 / meanY^2 + 1)
sdlogY <- sqrt(log(sdY^2 / meanY^2 + 1))
unifLo <- sdlogY / 100
unifHi <- sdlogY * 100

data_list <- list(N=length(y), y=y, unifLo=unifLo, unifHi=unifHi)
modelfit <- stan(model_code=model_string, data=data_list,
                 pars=c("mu_log", "sigma_log", "mu", "sigma", "cv"),
                 chains=3, iter=30000, warmup=10000)

# Exploring run convergence and estimates:
library("shinystan")
my_sso <- launch_shinystan(modelfit)