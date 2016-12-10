# The model specification:
model_string <- "
data {
int<lower=0> N1;
real y1[N1];
real<lower=0> sdY1;
real<lower=0> unifLo1;
real<lower=0> unifHi1;

int<lower=0> N2;
real y2[N2];
real<lower=0> sdY2;
real<lower=0> unifLo2;
real<lower=0> unifHi2;
}
parameters {
real<lower=0> sigma1;
real<lower=0> mu1;

real<lower=0> sigma2;
real<lower=0> mu2;
}
model {
sigma1 ~ uniform(unifLo1, unifHi1);
mu1 ~ uniform(0, unifHi1);
y1 ~ normal(mu1, sigma1);

sigma2 ~ uniform(unifLo2, unifHi2);
mu2 ~ uniform(0, unifHi2);
y2 ~ normal(mu2, sigma2);
}
generated quantities {
real<lower=0>  cv1;
real<lower=0>  cv2;
real cvdiff;

cv1 = (1.0 + 1.0 / (4.0 * N1)) * sigma1 / mu1;
cv2 = (1.0 + 1.0 / (4.0 * N2)) * sigma2 / mu2;

cvdiff = cv1 - cv2;
}"

# Make some test data:
library(perbbNovoMisc)
y1 <- draw_one_sampleCV_data(n_data1 = 16, mu1 = 1, std1 = 1, simdist = 'norm')$X1
y2 <- draw_one_sampleCV_data(n_data1 = 30, mu1 = 1, std1 = 0.2, simdist = 'norm')$X1


# Running the model:
library(rstan)
sdY1 <- sd(y1)
unifLo1 <- sdY1 / 100
unifHi1 <- sdY1 * 100

sdY2 <- sd(y2)
unifLo2 <- sdY2 / 100
unifHi2 <- sdY2 * 100



data_list <- list(N1=length(y1), y1=y1, sdY1=sdY1, unifLo1=unifLo1, unifHi1=unifHi1,
                  N2=length(y2), y2=y2, sdY2=sdY2, unifLo2=unifLo2, unifHi2=unifHi2)
modelfit <- stan(model_code=model_string, data=data_list,
                 pars=c("mu1", "sigma1", "cv1", "mu2", "sigma2", "cv2", "cvdiff"),
                 chains=3, iter=30000, warmup=10000)

# Exploring run convergence and estimates:
library("shinystan")
my_sso <- launch_shinystan(modelfit)