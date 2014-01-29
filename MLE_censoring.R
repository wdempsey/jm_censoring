revival_model <- function(Table_1, Table_2, X_1, X_2, Sigma_calc, mean_params, cov_params, theta, fixed = FALSE) {
	
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


logint <- function(mean_params, cov_params, theta, pat_table,c) {
	# Returns the log of the integral from censoring to infty 
	# of the joint density of (Y,t) 
	
	h_vec = Vectorize(h(mean_params, cov_params,theta,pat_table))
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

log_lik_cens <- function(mean_params, cov_params, theta, pat_table,c) {
	# Log-likelhood for Censored Patients
	return( log_cdf_T (theta,c) + logint(mean_params, cov_params, theta, pat_table,c))
}

log_lik_uncens <- function(mean_params, cov_params, theta, pat_table,T) {
	# Log-Likelihood For Uncensored Patients
	return(log( f(theta)(T) ) + log ( g(mean_params, cov_params, pat_table)(T)))
}

### Complete Log-Likelihood

log_lik <- function(mean_params, cov_params, theta, table1 = Table_1, table2 = Table_2) {
	llik = 0
	for (pat in table1$id) {
		pat_table = table2[table2$id == pat,]
		
		if (table1$cens[table1$id == pat] == 1) {
			c = table1$survival[table1$id == pat]
			llik = llik + log_lik_cens(mean_params, cov_params, theta, pat_table,c)
		}
		if (table1$cens[table1$id == pat] == 0) {
			T = table1$survival[table1$id == pat]
			llik = llik + log_lik_uncens(mean_params, cov_params, theta, pat_table,T)
		}
	}	
	return(llik)
}

log_lik_vector <- function(params, table1 = Table_1, table2 = Table_2) {
	mean_params = params[1:length(mean_params)]
	cov_params = params[(length(mean_params)+1):(length(params)-1)]
	theta = params[length(params)]
	return(-log_lik(mean_params, cov_params, theta, table1, table2))
}

log_lik_vector_fixed <- function(theta) {
  llik_theta <- function(params) {
  return(log_lik_vector(c(params,theta)))
  }
  return(llik_theta)  
}



### Calculate Gradient of the Log-Likelihood ###
K_1 <- function(pat_table) {
	return(diag(dim(pat_table)[1]))
}

K_2 <- function(pat_table) {
	 return(exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/1))
}

# do.call(paste("K_",2, sep = ""),list(sigmasq_1,pat_table))

expected_terms <- function(mean_params, cov_params, theta, pat, table1, table2) {	
	pat_table = table2[table2$id == pat,]

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params, cov_params, theta, pat_table))

	eval_T = seq(c, 30, by = 0.05)
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)
	K1 = K_1(pat_table)
	K2 = K_2(pat_table)

	quad_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%mean_params
		return(t(W_T)%*%solve(Sigma)%*%W_T)
	}	
	
	S_calc <- function(T) {
		if(dim(pat_table)[1] == 1) {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params
			return(X_T%*%solve(Sigma)%*%W_T)
		}
		else {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params
			return(t(X_T)%*%solve(Sigma)%*%W_T)
		}
	}
	
	K1_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%mean_params
		return(t(W_T)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%W_T)
	}		

	K2_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%mean_params
		return(t(W_T)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%W_T)
	}		

	quad = Vectorize(quad_calc)(eval_T)
	S = Vectorize(S_calc)(eval_T)
	K1_vec = Vectorize(K1_calc)(eval_T)
	K2_vec = Vectorize(K2_calc)(eval_T)

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)

	return(list("exp_T" = sum(eval_T*norm_probs),
	"exp_quad" = sum(quad*norm_probs), 
	"exp_S" = S%*%norm_probs, 
	"exp_K1" = sum(K1_vec*norm_probs),
	"exp_K2" = sum(K2_vec*norm_probs),
	"trace_K1" = sum(diag(Inv_Sigma%*%K1)),
	"trace_K2" = sum(diag(Inv_Sigma%*%K2)))
	)
	
}	


grad_calc <- function(mean_params, cov_params, theta, table1 = Table_1, table2 = Table_2) {
	grad_lambda = 0
	grad_beta = 0
	grad_sigma1 = 0
	grad_sigma2 = 0
	num_pats = dim(table1)[1]
	for (pat in 1:num_pats) {
		pat_table = table2[table2$id == pat, ]
		Sigma = Sigma_calc(cov_params, pat_table)
		Inv_Sigma = solve(Sigma)
			
		K1 = K_1(pat_table)
		K2 = K_2(pat_table)
			
		if(table1$cens[pat] == 0) {
			# Lambda
			T = table1$survival[pat]
			grad_lambda = grad_lambda + 1/theta - T
				
			# Beta
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params	
				
			if( dim(pat_table)[1] == 1) {
				grad_beta = grad_beta + (X_T)%*%Inv_Sigma%*%W_T							}
			else {
				grad_beta = grad_beta + t(X_T)%*%Inv_Sigma%*%W_T		
			}
								
			# Sigma			
			grad_sigma1 = grad_sigma1 +t(W_T)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%W_T/2 - sum(diag(Inv_Sigma%*%K1))/2
			grad_sigma2 = grad_sigma2 +t(W_T)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%W_T/2 - sum(diag(Inv_Sigma%*%K2))/2

		}
		if(table1$cens[pat] == 1) {
			c = table1$survival[pat]
			exp_terms = expected_terms(mean_params, cov_params, theta, pat, table1, table2) 

			# Lambda
			grad_lambda = grad_lambda + c + 1/theta - exp_terms$exp_T
				
			# Beta
			grad_beta = grad_beta + exp_terms$exp_S
				
			# Sigma			
			grad_sigma1 = grad_sigma1 - exp_terms$trace_K1/2 + exp_terms$exp_K1/2
			grad_sigma2 = grad_sigma2 - exp_terms$trace_K2/2 + exp_terms$exp_K2/2
		}
	}
	return(-c(grad_beta,grad_sigma1,grad_sigma2, grad_lambda))
}

grad_calc_vector <- function(params, table1 = Table_1, table2 = Table_2) {
	mean_params = params[1:length(mean_params)]
	cov_params = params[(length(mean_params)+1):(length(params)-1)]
	theta = params[length(params)]
	
	return(grad_calc(mean_params, cov_params, theta, table1, table2))
}

grad_calc_vector_fixed <- function(theta) {
	grad_theta <- function(params) {
		return(grad_calc_vector(c(params, theta)))
	}
	return(grad_calc)
}

log_lik_vector_fixed <- function(theta) {
  llik_theta <- function(params) {
  return(log_lik_vector(c(params,theta)))
  }
  return(llik_theta)  
}


# Fit from the Initialization Given By Model_Cens

print('Got to The Optimization Component')

#maxcens = max(Table_1$survival[Table_1$cens == 1])

#cdf_T(0.25, maxcens)

max_theta = 5*theta

if(fixed == TRUE) {
  print('Fixed Theta For Optimization')
  inits <- c(mean_params, cov_params)
  op_llik <- optim(inits, log_lik_vector_fixed(theta), lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), 0), upper = c(rep(Inf, length(mean_params) + length(cov_params)), 5))  
}
if(fixed == FALSE) {
  inits <- c(mean_params, cov_params, theta)
  op_llik <- optim(inits, log_lik_vector,grad_calc_vector, lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), 0), upper = c(rep(Inf, length(mean_params) + length(cov_params)), max_theta))  
}

print(op_llik$convergence)
print('Finished the Optimization Code')

mle_est <- op_llik$par

print(mle_est)

# Compute the Hessian at the MLE estimate

library(numDeriv)

print('Computing the Hessian')

obs_inf <- function(log_lik_vector, mle_est) {
	return(hessian(log_lik_vector, mle_est))	
}

if(fixed == TRUE) {
	Hess = obs_inf(log_lik_vector_fixed(theta), mle_est)	
}

if(fixed == FALSE) {
	Hess = obs_inf(log_lik_vector, mle_est)	
}


return(list("mle" = mle_est, "hess" = Hess, "conv" = op_llik$convergence))

}

