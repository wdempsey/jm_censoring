# Create the Correct Tables for Prothrombin Data

setwd('/Users/walterdempsey/Documents/stat/research/joint_models/dual-alignment-project/code/prednisone')

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
  return(model.matrix(form)[2:5,])
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
cov_params <- c(baseline_model$sigma,1)
theta <- (sum(Table_1$cens))/sum(Table_1$survival)

### Running the Code on Only Censored Individuals ##

Table_2_cens = Table_2[Table_2$cens==0,c(1:4,6:8)]
Table_1_cens = Table_1[Table_1$cens==0,]

source('/Users/walterdempsey/Documents/stat/research/joint_models/censoring/jm_censoring/MLE_censoring_gen.R')

revival_model(Table_1_cens, Table_2_cens, X_1, X_2, Sigma_calc, mean_params, cov_params, theta)
