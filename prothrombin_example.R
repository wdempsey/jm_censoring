## Turn Command Line Arguments On

args <- commandArgs(TRUE)
library(sendmailR)

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

Table_1$cens = Liver_Data$cens[as.numeric(lapply(Table_1$id, id_check))]
Table_1$survival = Liver_Data$survival[as.numeric(lapply(Table_1$id, id_check))]

# Create the Mean and Covariance Functions

pat_table = Table_2[Table_2$id == 220,]  
  
X_1 <- function(pat_table) {
  const = rep(1, dim(pat_table)[1])
  return(const)
}

X_2 <- function(t, pat_table) {
  # Returns the revival vector
  # and the log(s+delta)
  # and treatment vectors
  delta = 1
  revival = t- pat_table$obs_times
  log_rev = log(revival + delta)
  treatment = pat_table$treatment
  T = rep(t, length(revival))
  form=~treatment+T+revival+log_rev
  return(model.matrix(form)[,2:6])
}

Sigma_calc <- function(cov_params, pat_table) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  sigmasq_2 = cov_params[3]
  lambda = cov_params[4]
  
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


Table_2_uncens = Table_2[Table_2$cens==1,c(1:4,6:8)]
Table_1_uncens = Table_1[Table_1$cens==1,]

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

### Compute the Covariance Models

Patient <- outer(Table_2_uncens$id, Table_2_uncens$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.67
  
Patient.ds <- exp(-abs(outer(Table_2_uncens$revival,Table_2_uncens$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$treatment+Table_2_uncens$survival+Table_2_uncens$revival+Table_2_uncens$logrev, ~Patient + Patient.ds)
  
summary(baseline_model)

mean_params <- baseline_model$beta
cov_params <- c(baseline_model$sigma,cov_lambda)
theta <- (sum(Table_1$cens))/sum(Table_1$survival)

### Running the Code on Only Censored Individuals ##

Table_2_cens = Table_2[Table_2$cens==0,c(1:4,6:8)]
Table_1_cens = Table_1[Table_1$cens==0,]

source('MLE_censoring.R')

if(args[1] == 'cens') {
	print('Computing Censored Only Model')

	rev_mod <- revival_model(Table_1_cens, Table_2_cens, X_1, 	X_2, Sigma_calc, mean_params, cov_params, theta)
	
	mle = rev_mod$mle
	hess = rev_mod$hess

	write.table(mle, 'prot_mle')
	write.table(hess, 'prot_hess')

	from = sprintf("<admirR@\\%s>", Sys.info()[4])
	to = "<dempsey.walter@gmail.com>"
	subject <- "Completed Fitting Censored Model"
	body <- list("Check The Output")
	sendmail(from, to, subject, body, control=list(smtpServer="ASPMX.L.GOOGLE.COM"))

}

if(args[1] == 'uncens') {
	print('Computing UnCensored Only Model')

	rev_mod <- revival_model(Table_1_uncens, Table_2_uncens, X_1, 	X_2, Sigma_calc, mean_params, cov_params, theta)
	
	mle = rev_mod$mle
	hess = rev_mod$hess

	write.table(mle, 'prot_mle_uncens')
	write.table(hess, 'prot_hess_uncens')

	from = sprintf("<admirR@\\%s>", Sys.info()[4])
	to = "<dempsey.walter@gmail.com>"
	subject <- "Completed Fitting UnCensored Model"
	body <- list("Check The Output")
	sendmail(from, to, subject, body, control=list(smtpServer="ASPMX.L.GOOGLE.COM"))

}


if(args[1] == 'compl') {
	print('Computing Complete Model')

	rev_mod_compl <- revival_model(Table_1, Table_2, X_1, X_2, Sigma_calc, mean_params, cov_params, theta)

	mle_compl = rev_mod_compl$mle
	hess_compl = rev_mod_compl$hess

	write.table(mle_compl, 'prot_mle_comp')
	write.table(hess, 'prot_hess_comp')
	
	from = sprintf("<admirR@\\%s>", Sys.info()[4])
	to = "<dempsey.walter@gmail.com>"
	subject <- "Completed Fitting Complete Model"
	body <- list("Check The Output")
	sendmail(from, to, subject, body, control=list(smtpServer="ASPMX.L.GOOGLE.COM"))

}

