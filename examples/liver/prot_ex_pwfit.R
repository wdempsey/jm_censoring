## Turn Command Line Arguments On

args <- commandArgs(TRUE)

library(mail)

# Create the Correct Tables for Prothrombin Data

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

Liver_Data = read.table('LiverData.csv', sep = ',', header = TRUE)

summary(Liver_Data)

Table_2 = Liver_Data[,c(1,2,3,4)]

Table_2$treatment[Table_2$time == 0] = -1

Table_2$treatment = as.factor(Table_2$treatment)

levels(Table_2$treatment) = c('null', 'control', 'treatment')

names(Table_2)[3] = 'obs_times'
names(Table_2)[2] = 'obs'

id = as.numeric(levels(as.factor(Liver_Data$id)))

Table_1 = data.frame(id)

id_check <- function(x) {which(Liver_Data$id == x)[1]}

Table_1$cens = 1-Liver_Data$cens[as.numeric(lapply(Table_1$id, id_check))]
Table_1$survival = Liver_Data$survival[as.numeric(lapply(Table_1$id, id_check))]

# Create the Mean and Covariance Functions

pat_table = Table_2[Table_2$id == 220,]  
  
Cov <- function(t, pat_table) {
  # Returns the revival vector
  # and the log(s+delta)
  # and treatment vectors
  delta = 1/365
  revival = t- pat_table$obs_times
  log_rev = log(revival + delta)
  treatment = pat_table$treatment
  T = rep(t, length(revival))
  form=~treatment+T+revival+log_rev
  return(model.matrix(form))
}


Cov_int <- function(t, pat_table) {
  # Returns the revival vector
  # and the log(s+delta)
  # and treatment vectors
  delta = 1/365
  revival = t- pat_table$obs_times
  log_rev = log(revival + delta)
  treatment = pat_table$treatment
  T = rep(t, length(revival))
  form=~treatment+T+revival+log_rev+log_rev*treatment
  return(model.matrix(form))
}


Sigma_calc <- function(cov_params, pat_table) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  sigmasq_2 = cov_params[3]
  lambda = 1.67
  
  return( sigmasq_0 * diag(length(pat_table$obs_times)) + sigmasq_1 * outer(rep(1,length(pat_table$obs_times)),rep(1,length(pat_table$obs_times)))
          + sigmasq_2 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/lambda))
}


## Initialization 
delta = 1/365
censored <- function(x) {Table_1$cens[Table_1$id == x]}
survival <- function(x) {Table_1$survival[Table_1$id == x]}

Table_2$cens = as.numeric(lapply(Table_2$id, censored))
Table_2$survival = as.numeric(lapply(Table_2$id, survival))
Table_2$revival =  Table_2$survival - Table_2$obs_times
Table_2$logrev = log(Table_2$revival+delta)


Table_2_uncens = Table_2[Table_2$cens==0,c(1:4,6:8)]
Table_1_uncens = Table_1[Table_1$cens==0,]

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

### Compute the Covariance Models

Patient <- outer(Table_2_uncens$id, Table_2_uncens$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.67
  
Patient.ds <- exp(-abs(outer(Table_2_uncens$revival,Table_2_uncens$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$treatment+Table_2_uncens$survival+Table_2_uncens$revival+Table_2_uncens$logrev, ~Patient + Patient.ds, kernel = 0)

int_baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$treatment+Table_2_uncens$survival+Table_2_uncens$revival+Table_2_uncens$logrev+Table_2_uncens$logrev*Table_2_uncens$treatment, ~Patient + Patient.ds, kernel = 0)

mean_params <- baseline_model$beta
cov_params <- c(baseline_model$sigma)

gamma = 0

source('../../piecewise_code/fit.piecewiseph.R')

### Running the Code on Only Censored Individuals ##

Table_2_cens = Table_2[Table_2$cens==1,c(1:4,6:8)]
Table_1_cens = Table_1[Table_1$cens==1,]

### Calculate Gradient of the Log-Likelihood ###
K_1 <- function(pat_table) {
	return(diag(dim(pat_table)[1]))
}

K_2 <- function(pat_table) {
	return(outer(rep(1,length(pat_table$obs_times)),rep(1,length(pat_table$obs_times))))
}

K_3 <- function(pat_table) {
	return(exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/1.67))
}

K = list(K_1, K_2, K_3)

recipient = "dempsey.walter@gmail.com"

if(args[1] == 'cens') {
	print('Computing Censored Only Model')
	
	params = list('mean_params' = mean_params, 'cov_params' = cov_params, 'gamma' = gamma)

	rev_mod_fixed <- fit.piecewiseph(Table_1_cens, Table_2_cens, Cov, Sigma_calc, K, params, control = list(fixed = TRUE))
	
   	mle_fixed = rev_mod_fixed$mle
	hess_fixed = rev_mod_fixed$hess
    conv_fixed = rev_mod_fixed$conv
    
   	write.table(mle_fixed, 'piecewise_prot_mle_cens_fixed')
	write.table(hess_fixed, 'piecewise_prot_hess_cens_fixed')
    write.table(conv_fixed, 'piecewise_prot_conv_cens_fixed')


    sendmail(recipient, subject="Notification from R", message="Censored Only Model Calculation finished!")

}

if(args[1] == 'int') {
	print('Computing Interaction Censored Only Model')
	
	mean_params_int <- int_baseline_model$beta
	cov_params_int <- int_baseline_model$sigma
	
	params = list('mean_params' = mean_params_int, 'cov_params' = cov_params_int, 'gamma' = gamma)

	rev_mod_fixed <- fit.piecewiseph(Table_1_cens, Table_2_cens, Cov_int, Sigma_calc, K, params, control = list(fixed = TRUE))
	
   	mle_fixed = rev_mod_fixed$mle
	hess_fixed = rev_mod_fixed$hess
    conv_fixed = rev_mod_fixed$conv
    
   	write.table(mle_fixed, 'piecewise_prot_mle_cens_int_fixed')
	write.table(hess_fixed, 'piecewise_prot_hess_cens_int_fixed')
    write.table(conv_fixed, 'piecewise_prot_conv_cens_int_fixed')

    sendmail(recipient, subject="Notification from R", message="Censored Only Model Interaction Calculation finished!")

}

if(args[1] == 'uncens') {
	print('Computing UnCensored Only Model')

	params = list('mean_params' = mean_params, 'cov_params' = cov_params, 'gamma' = gamma)
	
	rev_mod_fixed <- fit.piecewiseph(Table_1_uncens, Table_2_uncens, Cov, Sigma_calc, K, params, control = list(fixed = TRUE))
    
   	mle_fixed = rev_mod_fixed$mle
	hess_fixed = rev_mod_fixed$hess
    conv_fixed = rev_mod_fixed$conv
    
   	write.table(mle_fixed, 'piecewise_prot_mle_uncens_fixed')
	write.table(hess_fixed, 'piecewise_prot_hess_uncens_fixed')
    write.table(conv_fixed, 'piecewise_prot_conv_uncens_fixed')

    sendmail(recipient, subject="Notification from R", message="Uncensored Only Model Calculation finished!")
}


if(args[1] == 'compl') {
	print('Computing Complete Model')

	params = list('mean_params' = mean_params, 'cov_params' = cov_params, 'theta' = theta, 'gamma' = gamma)

	rev_mod <- fit.piecewiseph(Table_1, Table_2, Cov, Sigma_calc, K, params, control = list(fixed = FALSE))
	
	rev_mod_fixed <- fit.piecewiseph(Table_1, Table_2, Cov, Sigma_calc, K, params, control = list(fixed = TRUE))
	
	mle = rev_mod$mle
	hess = rev_mod$hess
    conv = rev_mod$conv
    
   	mle_fixed = rev_mod_fixed$mle
	hess_fixed = rev_mod_fixed$hess
    conv_fixed = rev_mod_fixed$conv
    
	write.table(mle, 'piecewise_prot_mle_uncens')
	write.table(hess, 'piecewise_prot_hess_uncens')
    write.table(conv, 'piecewise_prot_conv_uncens')
    
   	write.table(mle_fixed, 'piecewise_prot_mle_uncens_fixed')
	write.table(hess_fixed, 'piecewise_prot_hess_uncens_fixed')
    write.table(conv_fixed, 'piecewise_prot_conv_uncens_fixed')

    sendmail(recipient, subject="Notification from R", message="Complete Model Calculation finished!")
}