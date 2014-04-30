em_model <- function(Table_1, Table_2, X_1, X_2, Sigma_calc, mean_params, cov_params, theta, fixed = FALSE, max_iter = 1000) {
	
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


#### EM Algorithm  ####

# K_1 <- function(sigmasq_1, pat_table) {
	# return(sigmasq_1*diag(dim(pat_table)[1]))
# }

# K_2 <- function(sigmasq_2, pat_table) {
	 # return(sigmasq_2 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/1))
# }

# K = list(K_1,K_2)

num_varcomp = length(K)

expected_terms <- function(mean_params_old, cov_params_old, theta_old, cov_params, pat, table1, table2) {
		
	pat_table = table2[table2$id == pat,]
	cens = table1$cens[table1$id == pat]
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)
#	K1 = K_1(1, pat_table)
#	K2 = K_2(1, pat_table)
	
	quad_XSX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		if(dim(pat_table)[1] == 1) {
			return((X_T)%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(X_T)%*%Inv_Sigma%*%(X_T))
		}
	}	

	quad_YSX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		Y = pat_table$obs
		if(dim(pat_table)[1] == 1) {
			return((Y)%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(Y)%*%Inv_Sigma%*%(X_T))
		}
	}		
	
	if(cens == 1) {

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params_old, cov_params_old, theta_old, pat_table))

	eval_T = seq(c, c+10, by = 0.1)


	quad_XSX = Vectorize(quad_XSX_calc)(eval_T)
	quad_YSX = Vectorize(quad_YSX_calc)(eval_T)	

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)
	
	d = sqrt(dim(quad_XSX)[1])
	
	return(list("exp_T" = sum(eval_T*norm_probs),
	"exp_quad_XSX" = matrix(quad_XSX%*%norm_probs, nrow = d, ncol = d), 
	"exp_quad_YSX" = matrix(quad_YSX%*%norm_probs, nrow = 1, ncol = d)
	))

	}
	
	if (cens == 0) {
	T = table1$survival[table1$id == pat]
	
	quad_XSX = quad_XSX_calc(T)
	quad_YSX = quad_YSX_calc(T)
	
	return(list("exp_T" = T,
	"exp_quad_XSX" = quad_XSX, 
	"exp_quad_YSX" = quad_YSX
	))
		
	}
	
}	

info_exp <- function(mean_params_old, cov_params_old, theta_old, em_mean_params, cov_params, pat, table1, table2) {
		
	pat_table = table2[table2$id == pat,]
	cens = table1$cens[table1$id == pat]
	
#	K1 = K_1(1, pat_table)
#	K2 = K_2(1, pat_table)
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)	

	dW_dsigma <- function(j,l) {
		h <- function(T) {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%em_mean_params

			K_j = do.call(paste("K_",j, sep = ""), list(pat_table))
			K_l = do.call(paste("K_",l, sep = ""), list(pat_table))
			return( t(W_T)%*%Inv_Sigma%*%K_l%*%K_j%*%Inv_Sigma%*%W_T )
		}	
		return(h)
	}

	K_i_calc <- function(i) {
		quad_K_i <- function(T) {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params
			K_i = do.call(paste("K_",i, sep = ""), list(pat_table))
			return(t(W_T)%*%Inv_Sigma%*%K_i%*%Inv_Sigma%*%W_T)
		}			
	}
	
	if(cens == 1) {

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params_old, cov_params_old, theta_old, pat_table))

	eval_T = seq(c, c+20, by = 0.1)	
	
	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)
		
	exp_K = vector(length = num_varcomp)
	trace_K = vector(length = num_varcomp)
	dW_dsigma_ij = matrix(nrow = num_varcomp, ncol = num_varcomp)
	trace_K_ij = matrix(nrow = num_varcomp, ncol = num_varcomp)

	for(i in 1:num_varcomp) {
		quad_K_i = K_i_calc(i)	
		K_i_vec = Vectorize(quad_K_i)(eval_T)
		K_i = do.call(paste("K_",i, sep = ""), list(pat_table))

		exp_K[i] = sum(norm_probs*K_i_vec)
		trace_K[i] = sum(diag(Inv_Sigma%*%K_i))

		for(j in 1:num_varcomp) {
			quad_K_j = K_i_calc(j)	
			K_j_vec = Vectorize(quad_K_j)(eval_T)
			K_j = do.call(paste("K_",j, sep = ""), list(pat_table))
			
			dW_dsigma_ij[i,j] = sum(Vectorize(dW_dsigma(i,j))(eval_T)*norm_probs)
			trace_K_ij[i,j] = sum(diag(Inv_Sigma%*%K_i%*%Inv_Sigma%*%K_j))
		}
	}

	return(list(
	"trace_K" = trace_K,
	"dW_dsigma" = dW_dsigma_ij,
	"exp_K" = exp_K,
	"trace_K_ij" = trace_K_ij
	))

	}
	
	if (cens == 0) {
	
	T = table1$survival[table1$id == pat]
	
	exp_K = vector(length = num_varcomp)
	trace_K = vector(length = num_varcomp)
	dW_dsigma_ij = matrix(nrow = num_varcomp, ncol = num_varcomp)
	trace_K_ij = matrix(nrow = num_varcomp, ncol = num_varcomp)

	for(i in 1:num_varcomp) {
		quad_K_i = K_i_calc(i)	
		K_i_vec = Vectorize(quad_K_i)(T)
		K_i = do.call(paste("K_",i, sep = ""), list(pat_table))

		exp_K[i] = K_i_vec
		trace_K[i] = sum(diag(Inv_Sigma%*%K_i))

		for(j in 1:num_varcomp) {
			quad_K_j = K_i_calc(j)	
			K_j_vec = Vectorize(quad_K_j)(T)
			K_j = do.call(paste("K_",j, sep = ""), list(pat_table))
			
			dW_dsigma_ij[i,j] = dW_dsigma(i,j)(T)
			trace_K_ij[i,j] = sum(diag(Inv_Sigma%*%K_i%*%Inv_Sigma%*%K_j))
		}
	}

	return(list(
	"trace_K" = trace_K,
	"dW_dsigma" = dW_dsigma_ij,
	"exp_K" = exp_K,
	"trace_K_ij" = trace_K_ij
	))


	}
	
	

}

### EM Algorithm ###
em_mle <- function(mean_params,cov_params,theta, table1, table2, max_iter = 1000, fixed_theta = FALSE) {
	
	n = length(table1$id)
	
	### Initialize
	for(iter in 1:max_iter) { 
	
		mean_params_old = mean_params
		cov_params_old = cov_params
		theta_old = theta

		info = matrix(0,nrow = num_varcomp, ncol = num_varcomp)
		score = rep(0,num_varcomp)
		
		
		### First Calculate the New Covariance Parameters!
		
		for(pat in table1$id) {
			
			info_terms = info_exp(mean_params_old, cov_params_old, theta_old, mean_params, cov_params, pat, table1, table2)
		
			for (i in 1:num_varcomp){
				score[i] = score[i] + -info_terms$trace_K[i]/2 + info_terms$exp_K[i]/2	
				for(j in 1:num_varcomp){
					info[i,j] = info[i,j] + info_terms$trace_K_ij[i,j]/2 - info_terms$dW_dsigma[i,j]				
				}
			}
		}
		
		cov_params = cov_params_old-solve(info)%*%score
		cov_params[cov_params<0] = 0
		
		# Now calculate the new mean parameters and survival parameters
		# Note: Survival parameters are unaffected by the covariance calculation above
		
		quad_X = 0  # Quadaratic Term X^T S X
		quad_Y = 0  # Quadaratic Term Y^T S X 
		sum_exp_T = 0 # Sum of Expected Survival Terms
		
		for(pat in table1$id) {
		
			exp_terms = expected_terms(mean_params_old, cov_params_old, theta_old, cov_params, pat, table1, table2)
					
			quad_X = quad_X + exp_terms$exp_quad_XSX
			quad_Y = quad_Y + exp_terms$exp_quad_YSX
			
			sum_exp_T = sum_exp_T + exp_terms$exp_T
		}
		
		mean_params = solve(quad_X)%*%t(quad_Y)
		
		if(fixed_theta){
			theta = theta_old
		}
		if(!fixed_theta){
			theta = n*solve(sum_exp_T)
		}
	
		
		delta = c(mean_params_old, cov_params_old, theta_old)-c(mean_params, cov_params, theta)
		
		if(max(abs(delta)) < 0.01) {break}
	
		
		# results = rbind(results, c(mean_params, cov_params, theta))
		
	}
	
	if(iter == max_iter) {error = TRUE}
	else {error = FALSE}
	
	return(list(mean_params=mean_params,
		cov_params = cov_params,
		theta = theta,
		error = error,
		max_error = max(delta)
		))
	
}


output = em_mle(mean_params,cov_params,theta, Table_1, Table_2, max_iter, fixed)

return(list("estimates" = output))


}

