setwd("/Users/walterdempsey/Documents/stat/research/joint_models/censoring/jm_censoring/")

# Example Data #

Table_1 = read.table('Table_1_cens')
Table_2 = read.table('Table_2_cens')


# Check the Data for obs_time > survival#
#throw_away = rep(0,0)

#for(i in 1:(dim(Table_2)[1])) {
#	throw_away[i] = as.numeric(Table_2$obs_time[i] > Table_1$survival[Table_1$id == Table_2$id[i]])
	
#}

#Table_2_cens = Table_2[!throw_away,]

#write.table(Table_2_cens, 'Table_2_cens')

Table_1, Table_2, X_1,X_2, Cov, Sigma_calc, mean_params, cov_params, theta

# Generate the T-Dependent Covariates

X_1 <- function(pat_table) {
	return(rep(1, dim(pat_table)[1]))	
}

X_2 <- function(t, pat_table) {
	# Returns the revival vector
	return(t- pat_table$obs_times)
}

# Total Covariates

Cov <- function(t, pat_table) {
	return(cbind(X_1(pat_table), X_2(t, pat_table)))
}



Sigma_calc <- function(sigmasq_0, sigmasq_1, lambda, pat_table) {
	# Computes the Covariance Matrix for a Patient
	
	return( sigmasq_0 * diag(length(pat_table$obs_times)) 
	+ sigmasq_1 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/lambda))
	
}

g <- function( beta, sigmasq_0, sigmasq_1, lambda, pat_table) {
	# Provides a function of the survival time, t,
	# for the likelihood Y | T
	
	g2 <- function(t) {
		k = dim(pat_table)[1]
		X = Cov(t, pat_table)
		Sigma = Sigma_calc(sigmasq_0, sigmasq_1, lambda, pat_table)
		mu = X%*%beta
		y = pat_table$obs
		return( (2*pi)^(-k/2)*det(Sigma)^(-1/2)*exp(-t(y - mu)%*% solve(Sigma) %*% (y - mu)/2) )
	}
	
	return(g2)	
}

f <- function(theta) { 
	# Assume an Exponential Model with rate parameter theta
	
	f2 <- function(t) {return(theta * exp(-theta * t))}
	
	return(f2)
	
}

h <- function(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table) {
	# Returns a function that returns value of the joint density 
	# at T = t given the parameter values
	
	h2 <- function(t) {
		t_dens = f(theta)
		y_dens = g(beta, sigmasq_0, sigmasq_1, lambda, pat_table) 
		return( t_dens(t) * y_dens(t) )
	}
	
	return(h2)
	
}


logint <- function(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,c) {
	# Returns the log of the integral from censoring to infty 
	# of the joint density of (Y,t) 
	
	h_vec = Vectorize(h(beta,sigmasq_0,sigmasq_1,lambda,theta,pat_table))
	return(log(integrate(h_vec, c, Inf)$value))
}

### Need -log(1-F(c))

cdf_T <- function(theta, x) {
	dens = f(theta)
	return(integrate(Vectorize(dens), 0,x)$value)
}

log_cdf_T <- function( theta, c) {
	return(-log(1-cdf_T(theta,c)))
}

log_lik_cens <- function(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,c) {
	# Log-likelhood for Censored Patients
	return( log_cdf_T (theta,c) + logint(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,c))
}

log_lik_uncens <- function(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,T) {
	# Log-Likelihood For Uncensored Patients
	return(log( f(theta)(T) ) + log ( g(beta, sigmasq_0, sigmasq_1, lambda, pat_table)(T)))
}

### Complete Log-Likelihood

log_lik <- function(beta, sigmasq_0, sigmasq_1, lambda, theta, table1 = Table_1, table2 = Table_2) {
	llik = 0
	for (pat in table1$id) {
		pat_table = table2[table2$id == pat,]
		
		if (table1$cens[table1$id == pat] == 1) {
			c = table1$survival[table1$id == pat]
			llik = llik + log_lik_cens(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,c)
		}
		if (table1$cens[table1$id == pat] == 0) {
			T = table1$survival[table1$id == pat]
			llik = llik + log_lik_uncens(beta, sigmasq_0, sigmasq_1, lambda, theta, pat_table,T)
		}
	}	
	return(llik)
}

log_lik_vector <- function(params, table1 = Table_1, table2 = Table_2) {
	beta = params[1:2]
	sigmasq_0 = params[3]
	sigmasq_1 = params[4]
	lambda = params[5]
	theta = params[6]
	return(log_lik(beta,sigmasq_0, sigmasq_1, lambda,theta, table1, table2))
}

# Initialize the Model

censored <- function(x) {Table_1$cens[Table_1$id == x]}
survival <- function(x) {Table_1$survival[Table_1$id == x]}

Table_2$cens = as.numeric(lapply(Table_2$id, censored))
Table_2$revival = Table_2$obs_times - as.numeric(lapply(Table_2$id, survival))

Table_2_uncens = Table_2[!Table_2$cens,c(1:3,5)]
Table_1_uncens = Table_1[!Table_1$cens,]

source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

### Compute the Covariance Models

Patient <- outer(Table_2_uncens$id, Table_2_uncens$id, "==")  # Patient Indicator Matrix

baseline_llik <- function(lambda) {
	cov_lambda <- lambda
	
	Patient.ds <- exp(-abs(outer(Table_2_uncens$revival,Table_2_uncens$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

	baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$revival, ~Patient.ds, pos = c(1), kernel = 1)

	return(baseline_model$llik)

}

cov_lambda <- 1

Patient.ds <- exp(-abs(outer(Table_2_uncens$revival,Table_2_uncens$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$revival, ~Patient.ds, pos = c(1), kernel = 1)

summary(baseline_model)

# Fit from the Initialization Given By Model_Cens
library(numDeriv)

beta = as.vector(baseline_model$beta)
sigmasq_0 = baseline_model$sigma[1]
sigmasq_1 = baseline_model$sigma[2]
theta = (sum(Table_1$cens))/sum(Table_1$survival)

eps = 0.005
delta = 0.001
max_eps = 0.05

step_size = seq(eps,max_eps,delta)
precision = 0.01
x_old = c(beta, max(sigmasq_0,0.05), max(sigmasq_1,0.05), cov_lambda, theta)
x_old

for(iter in 1:100) {
	fprime = grad(log_lik_vector,x_old)
	
	# Line Search	
	x_new = x_old + eps*c(fprime[1:5],0)
#	x_new = x_old + eps*fprime
#	grad_step = eps
#	while (grad_step < max_eps) { 
#		if(log_lik_vector(x_new) > log_lik_vector(x_old) ) {
#			x_new = x_new + delta*c(fprime[1:5],0)
#			grad_step = grad_step + delta
#		}
#		else {break}
#	}
#	x_new
#	max(abs(x_new - x_old))
	if(max(abs(x_new - x_old)) < precision) {break}
	print(max(abs(x_new-x_old)))
	print(iter)
	x_old = x_new
}
x_new
mle_est = x_new

# Optim Function

inits = c(beta, max(sigmasq_0,0.05), max(sigmasq_1,0.05), cov_lambda, theta)

log_lik_vector_2 <- function(params, table1 = Table_1, table2 = Table_2) {
	return(-log_lik_vector(params, table1, table2))
}

op_llik <- optim(inits, log_lik_vector_2,lower = c(-Inf, -Inf, 0, 0, 0, 0), upper = c(Inf, Inf, Inf, Inf, 5))

par_values <- op_llik$par

par_values <- c(-0.4952341,  1.0780908,  0.0000000, 12.1815841,  0.4898579,  0.3466315)



# Compute the Hessian at the MLE estimate

obs_inf <- function(log_lik_vector, mle_est) {
	return(hessian(log_lik_vector, mle_est))	
}

mle_cov = solve(obs_inf(log_lik_vector, mle_est))

nonzero <- function(x) { return(max(x,0))}

std_error <- function(est_cov) {
	eig_vals = eigen(-est_cov)$value
	eig_vectors = eigen(-est_cov)$vectors

	D = diag(as.numeric(lapply(eig_vals,nonzero)))
	U = eig_vectors

	Est_Param_Cov <- U%*%D%*%t(U)

	est_std_error = sqrt(diag(Est_Param_Cov))
	
	return(est_std_error)

}

cbind(mle_est,std_error(mle_cov))

cbind(par_values, std_error(par_cov))

