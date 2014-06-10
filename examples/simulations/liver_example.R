# Basic Example: 200 patients with the model based on the fit from the uncensored records for the Prothrombin, and exponential survival model 

library(mail)

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

compl_beta <- rep(0,0)
compl_sigma <- rep(0,0)

cens_beta <- rep(0,0)
cens_sigma <- rep(0,0)

rev_mle <- rep(0,0)

# for(k in 1:11) {
n_patients = 80

# Table 1 ~ Survival Data is Exponential with rate parameter lambda = 6
	# Assume No Censoring

lambda = 1/6 # Patients survive an avg of 6 years

Table_1 = matrix(nrow = n_patients, ncol = 4) 

Table_1[,1] = seq(1,n_patients)
Table_1[,2] = rep(0, n_patients)
Table_1[,3] = rexp(n_patients, rate = lambda)
Table_1[,4] = c(rep('control', n_patients/2), rep('prednisone', n_patients/2))

Table_1 = data.frame(Table_1)
names(Table_1) = c("id", "cens", "survival", "treatment")
Table_1$survival = as.numeric(Table_1$survival)

lambda_hat_uncens = n_patients/sum(Table_1$survival)

# Table 2 ~ Patients are observed once a year, assumed a linear revival model with exponential covariance between observations
Table_2 = rep(0,0)
delta = 1/365

intercept = 63.46
ctrl = 2.41
pred = 13.55
survival_coef = 1.75
revival_coef = -2.12
log_rev_coef = 4.66
pat_ar = 179.03
pat_spec = 208.86
meas_error = 210.89

expcov <- function(x,y) {exp(-abs(x-y)/1.67)}

fixed_apt_schedule = c(seq(0,1,1/3),seq(2,max(Table_1$survival),1))

for(i in 1:n_patients) {
	pat_surv = Table_1$survival[Table_1$id == i]

	observation_times = fixed_apt_schedule[fixed_apt_schedule< pat_surv]  # Observation Times for Patient i

	rev = pat_surv - observation_times
  
  if (Table_1$treatment[Table_1$id == i] == "control") {
    treat_effect = c(0,rep(ctrl, length(observation_times)-1))
  } else {
    treat_effect = c(0,rep(pred, length(observation_times)-1))    
  }
  
	pat_mean = intercept + treat_effect + survival_coef * pat_surv + revival_coef * rev + log_rev_coef * log(rev + delta)
	covariance = meas_error * diag(length(observation_times)) + pat_ar * outer(rev, rev, FUN = 'expcov') + pat_spec

	C = t(chol(covariance))  # cholesky decomposition
	pat_obs = pat_mean + C%*%rnorm(length(pat_mean))

	Table_2 = rbind(Table_2, cbind(rep(i, length(pat_obs)), observation_times, pat_obs))
}

Table_2 = data.frame(Table_2)

names(Table_2) = c('id', 'obs_times', 'obs')

### Build Revival and Treatment Columns
treatment <- function(x) {Table_1$treatment[Table_1$id == x]}
survival <- function(x) {Table_1$survival[Table_1$id == x]}

Table_2$treat = as.numeric(unlist(lapply(Table_2$id, treatment)))

Table_2$treat[Table_2$obs_times == 0] = 0

Table_2$treat = as.factor(Table_2$treat)

levels(Table_2$treat) = c('null', 'control', 'prednisone') 

Table_2$survival = unlist(lapply(Table_2$id, survival))
Table_2$revival =  Table_2$survival - Table_2$obs_times 

### Compute the Covariance Models

Patient <- outer(Table_2$id, Table_2$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.67; Patient.ds <- exp(-abs(outer(Table_2$revival,Table_2$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

### Formula

model1 <- regress(Table_2$obs ~ Table_2$treat + Table_2$survival + Table_2$revival  + log(Table_2$revival + delta), ~Patient + Patient.ds)

compl_beta <- cbind(compl_beta, model1$beta)
compl_sigma <- cbind(compl_sigma, model1$sigma)

#  Now assume a Censoring model that has a rate such that
#  the probabily of censoring is p.  So lambda_c = p * lambda / (1-p)

p = 1-292/488
cens_lambda = lambda*p / (1-p)

cens_time = rexp(n_patients, rate = cens_lambda)
surv_time = Table_1$survival
#cens_time = rep(quantile(surv_time, p), length(surv_time)) # This is if we want just right censoring at a particular value

Table_1_cens = Table_1
Table_1_cens$cens = as.numeric(cens_time < Table_1$survival)
Table_1_cens$survival = apply(cbind(cens_time, surv_time),1, min)

censored <- function(x) {Table_1_cens$cens[Table_1$id == x]}

Table_2_cens = Table_2[!as.logical(lapply(Table_2$id,censored)),]

### Compute the Covariance Models
Patient_cens <- outer(Table_2_cens$id, Table_2_cens$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.67; Patient_cens.ds <- exp(-abs(outer(Table_2_cens$revival, Table_2_cens$revival, "-"))/cov_lambda) *Patient_cens # Patient Specific Exponential Covariance Matrix

model2 <- regress(Table_2_cens$obs ~ Table_2_cens$treat + Table_2_cens$survival + Table_2_cens$revival  + log(Table_2_cens$revival + delta), ~Patient_cens + Patient_cens.ds, kernel = 0)

summary(model2)

cens_beta <- cbind(cens_beta, model2$beta)
cens_sigma <- cbind(cens_sigma, model2$sigma)

# Create the Mean and Covariance Functions
survival_cens_time <- function(x) {Table_1_cens$survival[Table_1_cens$id == x]}


Table_2_rev = Table_2[-which(Table_2$obs_times > as.numeric(lapply(Table_2$id,survival_cens_time))),]
  
Cov <- function(t, pat_table) {
  # Returns the revival vector
  revival = t- pat_table$obs_times
  return(model.matrix(~pat_table$treat + pat_table$survival + pat_table$revival  + log(pat_table$revival + delta)))
}

Sigma_calc <- function(cov_params, pat_table) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  sigmasq_2 = cov_params[3]
#   lambda = cov_params[3]
  lambda <- 1.67
  return( sigmasq_0 * diag(length(pat_table$obs_times)) + sigmasq_1 + sigmasq_2 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/lambda))
}

### Calculate Gradient of the Log-Likelihood ###
K_1 <- function(pat_table) {
	return(diag(dim(pat_table)[1]))
}

K_2 <- function(pat_table) {
  return(matrix(1, nrow = dim(pat_table)[1], ncol = dim(pat_table)[1]))
}

K_3 <- function(pat_table) {
	 return(exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/1))
}

K = list(K_1, K_2, K_3)

# Initialization

mean_params <- model2$beta
cov_params <- c(model2$sigma)

library(survival)

weibull_fit = survreg(Surv(Table_1_cens$survival,!Table_1_cens$cens)~1)

theta = list('k' = 1/weibull_fit$scale, 'lambda' = exp(weibull_fit$coef))

gamma = 0

source('../../weibull_code/fit.weibullph.R')

params = list('mean_params' = mean_params, 'cov_params' = cov_params, 'theta' = theta, 'gamma' = gamma)

rev <- fit.weibullph(Table_1_cens, Table_2_rev, Cov, Sigma_calc, K, params, control = list('fixed' = FALSE))

rev_fixed <- fit.weibullph(Table_1_cens, Table_2_rev, Cov, Sigma_calc, K, params, control = list('fixed' = TRUE))

write.table(rev$mle, './output/lin_mle', append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(rev$hess, '/output/lin_hess', append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(rev$conv, '/output/lin_conv', append = TRUE, row.names = FALSE, col.names = FALSE)

write.table(rev_fixed$mle, '/output/lin_fixed_mle', append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(rev_fixed$hess, '/output/lin_fixed_hess', append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(rev_fixed$conv, '/output/lin_fixed_conv', append = TRUE, row.names = FALSE, col.names = FALSE)

# ### Imputation Model ###

# source('Imputation_censoring.R')

# my_model <- function(table1, table2) {
	# ### Build the Table Mean Functions For the Model
	
	# survival <- function(x) {table1$survival[table1$id == x]}
	
	# table2$revival = as.numeric(lapply(table2$id, survival)) - table2$obs_times 
	
	# ### Compute the Covariance Models
	# Patient_cens <- outer(table2$id, table2$id, "==")  # Patient Indicator Matrix

	# cov_lambda <- 1.0; Patient_cens.ds <- exp(-abs(outer(table2$revival, table2$revival, "-"))/cov_lambda) *Patient_cens # Patient Specific Exponential Covariance Matrix

	# model <- regress(table2$obs ~ table2$revival, ~Patient_cens.ds, kernel = 1)
	
	# return(model)
	
# }

# M = 20

# rev_imp <- revival_model_imputation(Table_1_cens, Table_2_rev, X_1, X_2, Sigma_calc, my_model, M, mean_params, cov_params, theta, fixed = FALSE)


# results = cbind(rev$mle, c(rev_fixed$mle,theta), c(rev_imp$mle,theta))

# if (rev$conv == 0 & rev_fixed$conv == 0) {
	# write.table(results,"results", append = TRUE)

	# sendmail(recipient, subject="Notification from R", message="Simulation Successful!")
	
# }
# else {
	# sendmail(recipient, subject="Notification from R", message="Simulation Failure! Get 'em Next Time, Kid.")
# }


