
# Fits a Components of Variance Model

source("http://www.stat.uchicago.edu/~pmcc/courses/regress.R")

# Basic : Assume each patient is iid and take advantage of the block structure

# Structure of the Data
	# Assume 2 data matrices
	# Table 1 = Survival Information
		# Column 1 = id
		# Column 2 = censored indicator
		# Column 3 = survival time
		# Columns 4-n = Baseline Covariates for 
	# Table 2 = Longitudinal Observations
		# Column 1 = id (match Column 1 from Table 1)
		# Column 2 = Observation time (on the same scale as the survival time)
		# Column 3 = Observed Value
		# Column 4-n = Time-Dependent Observed Covariates


# Basic Example: 100 patients with a linear revival model, and exponential survival model 

# No Censoring
beta_results = rep(0,0)
sigma_results = rep(0,0)

for(k in 1:25) {
n_patients = 100

# Table 1 ~ Survival Data is Exponential with rate parameter lambda = 5
	# Assume No Censoring

lambda = 1/5 # Patients survive an avg of 5 years

Table_1 = matrix(nrow = n_patients, ncol = 3) 

Table_1[,1] = seq(1,n_patients)
Table_1[,2] = rep(0, n_patients)
Table_1[,3] = rexp(n_patients, rate = lambda)

Table_1 = data.frame(Table_1)
names(Table_1) = c("id", "cens", "survival")

lambda_hat_uncens = n_patients/sum(Table_1$survival)

# Table 2 ~ Patients are observed once a year, assumed a linear revival model with exponential covariance between observations
Table_2 = rep(0,0)

alpha = 0
beta = 1
sigmasq_0 = 1
sigmasq_1 = 0
expcov <- function(x,y) {exp(-abs(x-y))}

for(i in 1:n_patients) {
	pat_surv = Table_1$survival[Table_1$id == i]
	obs_times = seq(0,floor(pat_surv))  # Observation Times for Patient i

	revival = pat_surv - obs_times
	pat_mean = alpha + beta*revival
	covariance = sigmasq_0 * diag(length(obs_times)) + sigmasq_1 * outer(revival, revival, FUN = 'expcov')

	C = t(chol(covariance))  # cholesky decomposition
	pat_obs = pat_mean + C%*%rnorm(floor(pat_surv)+1)

	Table_2 = rbind(Table_2, cbind(rep(i, floor(pat_surv)+1), obs_times, pat_obs))
}

Table_2 = data.frame(Table_2)

names(Table_2) = c('id', 'obs_times', 'obs')


# We compute the revival times for the uncensored individuals and create a table of only uncensored data

revival_model <- function(Table_1, Table_2) {
  censored <- function(x) {Table_1$cens[Table_1$id == x]}
  survival_time <- function(x) {Table_1$survival[Table_1$id == x]}

  cens = as.logical(lapply(Table_2$id, censored))

  uncens_Table_2 = Table_2[!cens,]

  surv = as.numeric(lapply(uncens_Table_2$id, survival_time))

  revival = surv - uncens_Table_2$obs_times

  revival_table = cbind(uncens_Table_2, revival)

  # Compute Proper Covariance Functions

  Patient <- outer(revival_table$id, revival_table$id, "==")  # Patient Indicator Matrix

  cov_lambda <- 1; Patient.ds <- exp(-abs(outer(revival,revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

  
  # Compute the Uncensored Results
  const = rep(1,length(revival))
  model <- regress(revival_table$obs ~ revival, ~ Patient.ds, kernel = const, pos = c(1))

  return(model)

}

model <- revival_model(Table_1, Table_2)

#summary(model)

#  Now assume a Censoring model that has a rate such that
#  the probabily of censoring is 30%.  So lambda_c = p * lambda / (1-p)

p = 0.5
cens_lambda = lambda*p / (1-p)

cens_time = rexp(n_patients, rate = cens_lambda)
surv_time = Table_1$survival
#cens_time = rep(quantile(surv_time, p), length(surv_time)) # This is if we want just right censoring at a particular value

Table_1$cens = as.numeric(cens_time < Table_1$survival)
Table_1$survival = apply(cbind(cens_time, surv_time),1, min)

model2 <- revival_model(Table_1, Table_2)

#summary(model2)

# Compute the Survival Model

library(survival)

#survmodel = coxph(Surv(Table_1$surv,Table_1$cens)~const)

#fit <- survfit(Surv(Table_1$surv,Table_1$cens)~1) 

#fit$surv

# For now, I'm fitting a simple mean estimate of the exponential rate function

lambda_hat = sum(!Table_1$cens)/ sum(Table_1$survival)

# Predict Survival Times Given The Results for the Censored Observations

exp_survtime <- function(i,Table_1, Table_2, model) {
  obs_times = Table_2$obs_times[Table_2$id == i]
  obs = Table_2$obs[Table_2$id == i]
  cens = Table_1$cens[Table_1$id == i]
  surv = Table_1$survival[Table_1$id == i]

  precens_times = obs_times[which(obs_times < surv)]
  y = obs[which(obs_times < surv)]

  cov = model$sigma[1] * diag(length(y)) + model$sigma[2] * exp(-abs(outer(precens_times, precens_times, "-")))

  likelihood <- function(t) {
    revival_t = t - precens_times
    mean_t = model$beta[1] + model$beta[2]*revival_t
    lik_t = dexp(t, rate = lambda_hat)*(2*pi)^(-length(precens_times)/2)*det(solve(cov))^(-1/2)*exp(-(y-mean_t)%*%solve(cov)%*%(y-mean_t)/2)
    return(lik_t)
  }
  
  integrand <- function(t) { likelihood(t)*t}
  
  if(cens == 1) {
    int <- integrate(Vectorize(integrand), surv, surv+30)
    normalizing_constant <- integrate(Vectorize(likelihood), surv, surv+30)
  }
  
  if(cens == 0) {
    int <- integrate(Vectorize(integrand), max(precens_times), max(precens_times)+30)
    normalizing_constant <- integrate(Vectorize(likelihood), max(precens_times), max(precens_times)+30)
  }
  
  expected_survival_time = int$value/normalizing_constant$value
  return(expected_survival_time)
}

expsurvival_time_hat = vector(length = length(Table_1$id))

for (i in Table_1$id) {
  expsurvival_time_hat[i] = exp_survtime(i,Table_1, Table_2, model2)
}

#rmse = sqrt(mean((expsurvival_time_hat[Table_1$cens == 0] - surv_time[Table_1$cens == 0])^2))

# Update Table_1 To Reflect These Estimates

Table_1_exp = Table_1
Table_2_exp = Table_2

tolerance = 0.1
iter = 1
orig_censored = Table_1_exp$cens
Table_1_exp$cens = 0
rmse = 10

while(iter < 100) {
Table_1_exp$survival[orig_censored == 1] <- expsurvival_time_hat[orig_censored == 1]
model3 <- revival_model(Table_1_exp, Table_2_exp)
#summary(model3)

# No iteration right now - seems to not converge! #
old = expsurvival_time_hat

for (i in Table_1$id) {
  expsurvival_time_hat[i] = exp_survtime(i,Table_1, Table_2, model3)
}

rmse_new = sqrt(mean(c((expsurvival_time_hat[orig_censored == 0] - Table_1$surv[orig_censored == 0])^2,(expsurvival_time_hat[orig_censored == 1] - old[orig_censored == 1])^2)))

if(abs(rmse - rmse_new) < tolerance) {break}

rmse = rmse_new

iter = iter + 1
}

# Model Results

beta_results = cbind(beta_results, c(model$beta, model2$beta, model3$beta))

sigma_results = cbind(sigma_results, c(model$sigma, model2$sigma, model3$sigma))

}

mean_beta = apply(beta_results,1,mean)
sd_beta = apply(beta_results,1,sd)

bias_beta =  rep(c(alpha,beta),3) - mean_beta

mean_sigma = apply(sigma_results,1,mean)
sd_sigma  = apply(sigma_results,1,sd)

bias_sigma = rep(c(sigmasq_0,sigmasq_1),3) - mean_sigma

# Embed inside a function
# Include Prediction and Iteration Capacity

# Do a Imputation Version  #

imp_survtime <- function(i,Table_1, Table_2, model) {
  obs_times = Table_2$obs_times[Table_2$id == i]
  obs = Table_2$obs[Table_2$id == i]
  cens = Table_1$cens[Table_1$id == i]
  surv = Table_1$survival[Table_1$id == i]
  
  precens_times = obs_times[which(obs_times < surv)]
  y = obs[which(obs_times < surv)]
  
  cov = model$sigma[1] * diag(length(y)) + model$sigma[2] * exp(-abs(outer(precens_times, precens_times, "-")))
  
  likelihood <- function(t) {
    revival_t = t - precens_times
    mean_t = model$beta[1] + model$beta[2]*revival_t
    lik_t = dexp(t, rate = lambda_hat)*(2*pi)^(-length(precens_times)/2)*det(solve(cov))^(-1/2)*exp(-(y-mean_t)%*%solve(cov)%*%(y-mean_t)/2)
    return(lik_t)
  }
  
  
  possible_t = seq(surv, surv+30, .05)
  normalizing_constant <- integrate(Vectorize(likelihood), possible_t[1], possible_t[length(possible_t)])
  
  int_likelihood <- function(final_time) {
    int <- integrate(Vectorize(likelihood), possible_t[1], final_time)
    cdf_at_finaltime = int$value/normalizing_constant$value
    return(cdf_at_finaltime)
  }
  # Probably can do a basic search through the area or find a function that does this on the line
  
  cdf = as.numeric(lapply(possible_t, int_likelihood))
  
  rand_uniform = runif(1)
  return(possible_t[abs(cdf-rand_uniform) == min(abs(cdf - rand_uniform))])  
  
}


imp_survival_time_hat = vector(length = length(Table_1$id))
Table_1_imp = Table_1
Table_2_imp = Table_2

for (i in Table_1$id) {
  imp_survival_time_hat[i] = imp_survtime(i,Table_1_imp, Table_2_imp, model2)
}

orig_censored = Table_1_imp$cens
Table_1_imp$cens = 0

Table_1_imp$survival[orig_censored == 1] <- imp_survival_time_hat[orig_censored == 1]

model4 <- revival_model(Table_1_imp, Table_2_imp)
summary(model4)
  