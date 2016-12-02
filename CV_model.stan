# The model specification
model_string <- "
data {
    int<lower=0> N;
    real y[N];
}

parameters {
    real<lower=0> sigma;
    real<lower=0> mu;
}

model {
    mu ~ normal(1, 1);
    sigma ~ normal(0, 1);
    y ~ lognormal(mu,sigma);
}

generated quantities {
    real cv;
    cv = sigma / mu;
}"

# Get data from somewhere: y = vector of datavalues

# Running the model
library(rstan)
modelfit = stan(model_code=model_string, data=list(N=length(y), y=y), pars=c("mu", "sigma", "cv"), chains=3, iter=30000, warmup=10000)

# Exploring run convergence and estimates
library("shinystan")
my_sso <- launch_shinystan(modelfit)
