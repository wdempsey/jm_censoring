revival_model_imputation <- function(Table_1, Table_2, X_1, X_2, Sigma_calc, model, M, mean_params, cov_params, theta, fixed = FALSE) {
	
	# Fits a Components of Variance Model

	# Basic : Assume each patient is iid and take advantage of the block structure

	# Structure of the Data
		# Assume 2 data matrices
		# Table 1 = Survival Information
			# Column 1 = id (called 'id')
			# Column 2 = censored indicator (called 'cens')
			# Column 3 = survival time (called 'survival')
		# Table 2 = Longitudinal Observations
			# Column 1 = id (match Column 1 from Table 1)
			# Column 2 = Observation time (on the same scale as the survival time) (called 'obs_time')
			# Column 3 = Observed Value (called 'obs')
			# Column 4-n = Observed Covariates to Be Used in Fitting Model
			
X_1 <- function(pat_table) {
  const = rep(1, dim(pat_table)[1])
  return(const)
}

X_2 <- function(t, pat_table) {
  # Returns the revival vector
  revival = t- pat_table$obs_times
  return(revival)
}

Sigma_calc <- function(cov_params, pat_table) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
#   lambda = cov_params[3]
  lambda <- 1
  return( sigmasq_0 * diag(length(pat_table$obs_times)) + sigmasq_1 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/lambda))
}

Cov <- function(t, pat_table) {
  if (dim(pat_table)[1] == 1) {
    return(c(X_1(pat_table), X_2(t, pat_table)))
  }
  else {
    return(cbind(X_1(pat_table), X_2(t, pat_table)))
  }
}

g <- function( mean_params, cov_params, pat_table) {
	# Provides a function of the survival time, t,
	# for the likelihood Y | T
	
	g2 <- function(t) {
		k = dim(pat_table)[1]
		X = Cov(t, pat_table)
		Sigma = Sigma_calc(cov_params, pat_table)
		mu = X%*%mean_params
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

h <- function(mean_params, cov_params, theta, pat_table) {
	# Returns a function that returns value of the joint density 
	# at T = t given the parameter values
	
	h2 <- function(t) {
		t_dens = f(theta)
		y_dens = g(mean_params, cov_params, pat_table) 
		return( t_dens(t) * y_dens(t) )
	}
	
	return(h2)
	
}

#### Imputation  ####

impute_T <- function(mean_params, cov_params, theta, pat, table1, table2) {
	pat_table = table2[table2$id == pat,]

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params, cov_params, theta, pat_table))

	eval_T = seq(c, c+10, by = 0.1)

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)

	diff = abs(runif(1) - cumsum(norm_probs))

	return(eval_T[which(diff == min(diff))])
	
}

impute_table <- function(mean_params, cov_params, theta, table1, table2) {
	num_pats = dim(table1)[1]
	num_cols = dim(table2)[2]
	table1_impute = table1
	table1_impute$cens = 0
	for (pat in 1:num_pats) {
		if (table1$cens[table1$id == pat] == 1) {
			table1_impute$survival[pat] = impute_T(mean_params, cov_params, theta, pat, table1, table2)
		}
	}
	return(table1_impute)
}

### Imputation Output ####
imp_mean_results = matrix(nrow = length(mean_params), ncol = M)
imp_cov_results = matrix(nrow = length(cov_params), ncol = M)

imp_mean_results.stderr = matrix(nrow = length(mean_params), ncol = M)
imp_cov_results.stderr = matrix(nrow = length(cov_params), ncol = M)


for (i in 1:M) {
	imp_table = impute_table(mean_params, cov_params, theta, Table_1, Table_2)
	imp_model <- model(imp_table,Table_2)

	imp_mean_params = imp_model$beta
	imp_cov_params = imp_model$sigma
	
	imp_mean_params.stderr = imp_model$beta.se
	imp_cov_params.stderr = sqrt(diag(imp_model$sigma.cov))
	
	imp_mean_results[,i] = imp_mean_params
	imp_cov_results[,i] = imp_cov_params

	imp_mean_results.stderr[,i] = imp_mean_params.stderr
	imp_cov_results.stderr[,i] = imp_cov_params.stderr
}

imp_mean <- apply(imp_mean_results,1,mean)
imp_cov <- apply(imp_cov_results,1,mean)

W_mean = apply(imp_mean_results.stderr^2,1,mean)
W_cov = apply(imp_cov_results.stderr^2,1,mean)

B_mean = (M+1)/(M*(M-1))*apply((imp_mean_results - imp_mean)^2,1,sum)
B_cov = (M+1)/(M*(M-1))*apply((imp_cov_results - imp_cov)^2,1,sum)

imp_variance_mean <- W_mean + B_mean
imp_variance_cov <- W_cov + B_cov

return(list("mle" = c(imp_mean,imp_cov), "stderr" = c(sqrt(imp_variance_mean), sqrt(imp_variance_cov))))

}