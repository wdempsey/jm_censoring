# Basic Example: 100 patients with a linear revival model, and exponential survival model 

# beta_results = rep(0,0)
# sigma_results = rep(0,0)

# for(k in 1:25) {
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
sigmasq_0 = 10
sigmasq_1 = 5
expcov <- function(x,y) {exp(-abs(x-y))}

for(i in 1:n_patients) {
	pat_surv = Table_1$survival[Table_1$id == i]
	observation_times = seq(0,floor(pat_surv))  # Observation Times for Patient i

	rev = pat_surv - observation_times
	pat_mean = alpha + beta*rev
	covariance = sigmasq_0 * diag(length(observation_times)) + sigmasq_1 * outer(rev, rev, FUN = 'expcov')

	C = t(chol(covariance))  # cholesky decomposition
	pat_obs = pat_mean + C%*%rnorm(floor(pat_surv)+1)

	Table_2 = rbind(Table_2, cbind(rep(i, floor(pat_surv)+1), observation_times, pat_obs))
}

Table_2 = data.frame(Table_2)

names(Table_2) = c('id', 'obs_times', 'obs')

### Build Revival Column

survival <- function(x) {Table_1$survival[Table_1$id == x]}

revival = as.numeric(lapply(Table_2$id, survival)) - Table_2$obs_times 

### Compute the Covariance Models

Patient <- outer(Table_2$id, Table_2$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1; Patient.ds <- exp(-abs(outer(revival,revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

### Formula

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

mean_formula <- obs ~ revival
cov_formula <- ~ Patient.ds

model1 <- regress(Table_2$obs~revival, ~Patient.ds)

#  Now assume a Censoring model that has a rate such that
#  the probabily of censoring is p.  So lambda_c = p * lambda / (1-p)

p = 0.5
cens_lambda = lambda*p / (1-p)

cens_time = rexp(n_patients, rate = cens_lambda)
surv_time = Table_1$survival
#cens_time = rep(quantile(surv_time, p), length(surv_time)) # This is if we want just right censoring at a particular value

Table_1_cens = Table_1
Table_1_cens$cens = as.numeric(cens_time < Table_1$survival)
Table_1_cens$survival = apply(cbind(cens_time, surv_time),1, min)
Table_2_cens = Table_2

### Build Revival Table

Table_2_cens$revival = as.numeric(lapply(Table_2_cens$id, survival)) - Table_2_cens$obs_times 

### Compute the Covariance Models
Patient_cens <- outer(Table_2_cens$id, Table_2_cens$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.0; Patient_cens.ds <- exp(-abs(outer(Table_2_cens $revival, Table_2_cens$revival, "-"))/cov_lambda) *Patient_cens # Patient Specific Exponential Covariance Matrix

model2 <- regress(Table_2_cens$obs ~ revival, ~Patient_cens.ds, kernel = 1)

# Create the Mean and Covariance Functions

pat_table = Table_2[Table_2$id == 94,]  
  
X_1 <- function(pat_table) {
  const = rep(1, dim(pat_table)[1])
}

X_2 <- function(t, pat_table) {
  # Returns the revival vector
  revival = t- pat_table$obs_times
  return(revival)
}

Sigma_calc <- function(cov_params, pat_table) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  lambda = cov_params[3]
  return( sigmasq_0 * diag(length(pat_table$obs_times)) + sigmasq_1 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/lambda))
}

# Initialization

mean_params <- model2$beta
cov_params <- c(model2$sigma,lambda)
theta <- (sum(Table_1_cens$cens))/sum(Table_1_cens$survival)

source('MLE_censoring.R')

rev <- revival_model(Table_1_cens, Table_2_cens, X_1, X_2, Sigma_calc, mean_params, cov_params, theta, fixed = FALSE)

mle = rev$mle
Hess = rev$hess

write.table(mle, 'sim_mle')
write.table(Hess, 'sim_hess')



