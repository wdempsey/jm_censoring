fit.weibullph <- function(Table_1, Table_2, Cov, Sigma_Calc, K, params, control) {
	
	# Fits a Components of Variance Model

	# Basic : Assume each patient is iid and take advantage of the block structure

	# Structure of the Data
		# Assume 2 data matrices (Save as Matrices Not data Frames)
		# Table 1 = Survival Information
			# Column 1 = id (called 'id')
			# Column 2 = censored indicator (called 'cens')
			# Column 3 = survival time (called 'survival')
		# Table 2 = Longitudinal Observations
			# Column 1 = id (match Column 1 from Table 1)
			# Column 2 = Observation time (on the same scale as the survival time) (called 'obs_time')
			# Column 3 = Observed Value (called 'obs')
			# Column 4-n = Observed Covariates to Be Used in Fitting Model
			
# Extract Parameters
mean_params = params$mean_params
cov_params = params$cov_params
theta = params$theta
gamma = params$gamma
num_varcomp = length(K)

# Fix W = 0 For Now
W = 0

list.params <- list(mean_params = mean_params, cov_params = cov_params, theta = theta, gamma = gamma)

params <- unlist(as.relistable(list.params))

# For Trial Period

# Cov <- function(t, pat_table) {
  # # Returns the revival vector
  # # and the log(s+delta)
  # # and treatment vectors
  # delta = 1/365
  # revival = t- pat_table$obs_times
  # log_rev = log(revival + delta)
  # treatment = pat_table$treatment
  # T = rep(t, length(revival))
  # form=~treatment+T+revival+log_rev
  # return(model.matrix(form))
# }

# Density Functions

log_g <- function(mean_params, cov_params, obs,pat_table) {
	# Provides a function of the survival time, t,
	# for the likelihood Y | T
	
	log_g_t <- function(t) {
		k = dim(pat_table)[1]
		X = Cov(t, pat_table)
		Sigma = Sigma_calc(cov_params, pat_table)
		mu = X%*%mean_params
		y = obs
		return( -k/2*log(2*pi) - as.numeric(determinant(Sigma)$modulus)/2 - t(y - mu)%*% solve(Sigma, (y - mu))/2 )
	}
	
	return(log_g_t)	
}

baseline_hazard <- function(theta) {
	# Returns Baseline Hazard Function
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	
	baseline_hazard_t <- function (t) {
		return(exp(log(theta$k) - theta$k * log(theta$lambda) + (theta$k-1)*log(t)))
	}
	
	return(baseline_hazard_t)
}

log_surv <- function(theta, gamma, W) {
	# Log Survival Distribution
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	# W is the baseline covariates
	
	log_surv_t <- function(t) {
		return(
			-integrate(baseline_hazard(theta),0,t)$value * exp(gamma%*%W)
		)	
	}	
	return(log_surv_t)
}

log_f <- function(theta, gamma,W) { 
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	# W is the baseline covariates
	
	log_f_t <- function(t) {
		
		return(log_surv(theta,gamma,W)(t) + gamma%*%W + log(baseline_hazard(theta)(t)) )
		
	}
	
	return(log_f_t)
	
}

h <- function(mean_params, cov_params, theta, gamma, W, obs, pat_table) {
	# Returns a function that returns value of the joint density 
	# at T = t given the parameter values
	
	h_t <- function(t) {
		t_dens = log_f(theta,gamma, W)
		y_dens = log_g(mean_params, cov_params, obs, pat_table) 
		return( exp(t_dens(t) + y_dens(t)))
	}
	
	return(h_t)
	
}

logint <- function(mean_params, cov_params, theta, gamma, W, obs, pat_table,c) {
	# Returns the log of the integral from censoring to infty 
	# of the joint density of (Y,t) 
	
	h_vec = Vectorize(h(mean_params, cov_params,theta,gamma, W, obs, pat_table))
	return(log(integrate(h_vec, c, c+30)$value))
}

log_lik_cens <- function(mean_params, cov_params, theta, gamma, W, obs, pat_table,c) {
	# Log-likelhood for Censored Patients
	return( logint(mean_params, cov_params, theta, gamma, W, obs, pat_table,c) ) #- log_surv(theta,gamma,W)(c) )
}

log_lik_uncens <- function(mean_params, cov_params, theta, gamma, W, obs, pat_table,T) {
	# Log-Likelihood For Uncensored Patients
	return(log_f(theta, gamma, W)(T) + log_g(mean_params, cov_params, obs, pat_table)(T) )
}

### Complete Log-Likelihood

log_lik <- function(params, table1 = Table_1, table2 = Table_2) {

    params <- relist(params, skeleton = list.params)
    
	mean_params = params$mean_params
	cov_params = params$cov_params
	theta = params$theta
	gamma = params$gamma

	llik = 0
	for (pat in table1$id) {
		pat_table = table2[table2$id == pat,]
		
		if (table1$cens[table1$id == pat] == 1) {
			c = table1$survival[table1$id == pat]
			llik = llik + log_lik_cens(mean_params, cov_params, theta, gamma, W, pat_table$obs, pat_table,c)
		}
		if (table1$cens[table1$id == pat] == 0) {
			T = table1$survival[table1$id == pat]
			llik = llik + log_lik_uncens(mean_params, cov_params, theta, gamma, W, pat_table$obs, pat_table,T)
		}
	}	
	return(-llik)
}

### Weibull Score Equations ##
shape_grad <- function(theta, gamma, W, c) {
		shape_grad_t <- function(t) {
			if (c != 0) {
				integral = exp(gamma%*%W) * ( (t/theta$lambda)^theta$k * log(t/theta$lambda)) # - (c/theta$lambda)^theta$k * log(c/theta$lambda)  )
			} else{
				integral = exp(gamma%*%W) * (t/theta$lambda)^theta$k * log(t/theta$lambda) 	
			}
			return(1/theta$k + log(t) - log(theta$lambda) - integral)
		}
		return(shape_grad_t)
	}
	
scale_grad <- function(theta, gamma, W, c) {
		scale_grad_t <- function(t) {
			integral = exp(gamma%*%W) * ((theta$k/theta$lambda)*(t/theta$lambda)^theta$k) # - (theta$k/theta$lambda)*(c/theta$lambda)^theta$k ) 
			return(-theta$k/theta$lambda + integral)
		}
		return(scale_grad_t)
	}

### Expected Terms
expected_terms <- function(params, pat_table, c) {    
	mean_params = params$mean_params
	cov_params = params$cov_params
	theta = params$theta
	gamma = params$gamma

	cond_dens = Vectorize(h(mean_params, cov_params, theta, gamma, W, pat_table$obs, pat_table))

	eval_T = seq(c,c+30, by = 0.05)
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)

	quad_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%mean_params
		return(t(W_T)%*%solve(Sigma)%*%W_T)
	}	
	
	S_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%mean_params
		return(t(X_T)%*%solve(Sigma)%*%W_T)
	}
	
	K_i_calc <- function(i) {
		quad_K_i <- function(T) {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params
			K = do.call(paste("K_",i, sep = ""), list(pat_table))
			return(t(W_T)%*%Inv_Sigma%*%K%*%Inv_Sigma%*%W_T)
		}			
	}
		
	shape = Vectorize(shape_grad(theta,gamma,W,c))(eval_T)
	scale = Vectorize(scale_grad(theta,gamma,W,c))(eval_T)
	
	quad = Vectorize(quad_calc)(eval_T)
	S = Vectorize(S_calc)(eval_T)

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)
	
	exp_K = list()
	trace_K = list()

	for(i in 1:num_varcomp) {
		quad_K_i = K_i_calc(i)	
		K_i_vec = Vectorize(quad_K_i)(eval_T)
		K = do.call(paste("K_",i, sep = ""), list(pat_table))

		exp_K[[i]] = sum(norm_probs*K_i_vec)
		trace_K[[i]] = sum(diag(Inv_Sigma%*%K))
	}

	return(list("exp_T" = sum(eval_T*norm_probs),
	"exp_quad" = sum(quad*norm_probs), 
	"exp_S" = S%*%norm_probs, 
	"exp_K" = exp_K,
	"trace_K" = trace_K,
	"shape" = sum(shape*norm_probs),
	"scale" = sum(scale*norm_probs)
	))
	
}	


grad_calc <- function(params, table1 = Table_1, table2 = Table_2) {

    params <- relist(params, skeleton = list.params)
    
	mean_params = params$mean_params
	cov_params = params$cov_params
	theta = params$theta
	gamma = params$gamma
	
	patients = table1$id

	grad_lambda = 0
	grad_k = 0
	grad_beta = 0
	grad_sigma = list()
	for(i in 1:num_varcomp) {	
		grad_sigma[[i]] = 0
	}
	num_pats = dim(table1)[1]
	for (pat in patients) {
		pat_table = table2[table2$id == pat, ]
		Sigma = Sigma_calc(cov_params, pat_table)
		Inv_Sigma = solve(Sigma)
			
		if(table1$cens[table1$id == pat] == 0) {
			# Lambda
			T = table1$survival[pat]
			
			grad_k = grad_k + shape_grad(theta,gamma, W, 0)(T)
			
			grad_lambda = grad_lambda + scale_grad(theta,gamma, W, 0)(T)
				
			# Beta
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params	
			
			grad_beta = grad_beta + t(X_T)%*%Inv_Sigma%*%W_T		
								
			# Sigma			
			for(i in 1:num_varcomp) {	
				K = do.call(paste("K_",i, sep = ""), list(pat_table))
				grad_sigma[[i]] = grad_sigma[[i]] + t(W_T)%*%Inv_Sigma%*%K%*%Inv_Sigma%*%W_T/2 - sum(diag(Inv_Sigma%*%K))/2
			}

		}
		if(table1$cens[table1$id == pat] == 1) {
			c = table1$survival[table1$id == pat]
			exp_terms = expected_terms(params,pat_table,c) 

			# Lambda
			grad_k = grad_k + exp_terms$shape
			
			grad_lambda = grad_lambda + exp_terms$scale
				
			# Beta
			grad_beta = grad_beta + exp_terms$exp_S
				
			# Sigma					
			for(i in 1:num_varcomp) {	
				grad_sigma[[i]] = grad_sigma[[i]] - exp_terms$trace_K[[i]]/2 + exp_terms$exp_K[[i]]/2
			}
		}
	}
	return(-c(grad_beta,unlist(grad_sigma), grad_k, grad_lambda,0))
}

### Fixed Theta Log-lik and Grad Calculations

list.params_fixed <- list(mean_params = mean_params, cov_params = cov_params, gamma = gamma)

params_fixed <- unlist(as.relistable(list.params_fixed))

log_lik_fixed <- function(theta, table1 = Table_1, table2 = Table_2) {
	
	log_lik_fixed_2 <- function(params_fixed) {
	    params_fixed <- relist(params_fixed, skeleton = list.params_fixed)
	    
	    list.params <- list(mean_params = params_fixed$mean_params, cov_params = params_fixed$cov_params, theta = theta, gamma = params_fixed$gamma)
	    
	    params <- unlist(as.relistable(list.params))
	    
	    return(log_lik(params, table1, table2))			
	}

	return( log_lik_fixed_2 )
}

grad_calc_fixed <- function(theta, table1 = Table_1, table2 = Table_2) {
		
	grad_calc_fixed_2 <- function(params_fixed) {
	    params_fixed <- relist(params_fixed, skeleton = list.params_fixed)
	    
	    list.params <- list(mean_params = params_fixed$mean_params, cov_params = params_fixed$cov_params, theta = theta, gamma = params_fixed$gamma)
	    
	    params <- unlist(as.relistable(list.params))
	    
	    num_longparams = length(params_fixed$mean_params)+length(params_fixed$cov_params)
	    
	    num_baseparams = length(params_fixed$gamma)
	    
	    num_totalparams = length(params)
	    
	    return(grad_calc(params, table1, table2)[c(1:num_longparams,(num_totalparams-num_baseparams+1):num_totalparams)]) 			
	}

	return( grad_calc_fixed_2 )

}

# Fit from the Initialization Given By Model_Cens

print('Got to The Optimization Component')

library(numDeriv)

max_k = 10*theta$k
max_lambda = 10*theta$lambda

min_k = theta$k/10
min_lambda = theta$lambda/10

if(control$fixed == FALSE) {
  print('Theta is Not Fixed For Optimization')
  inits <- params
  op_llik <- optim(inits, log_lik, grad_calc, lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), min_k, min_lambda, -Inf), upper = c(rep(Inf, length(mean_params) + length(cov_params)), max_k, max_lambda,Inf), hessian = TRUE, control = list(trace = 1, maxit = 500))  
}
if(control$fixed == TRUE) {
  print('Theta is Fixed For Optimization')
  inits <- params_fixed
  op_llik <- optim(inits, log_lik_fixed(theta), grad_calc_fixed(theta), lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), -Inf), upper = c(rep(Inf, length(mean_params) + length(cov_params)),Inf), hessian = TRUE, control = list(trace = 1, maxit = 500))  
}

print(op_llik$convergence)
print('Finished the Optimization Code')

# Return the MLE Estimates

return(list("mle" = op_llik$par, "hess" = op_llik$hessian, "conv" = op_llik$convergence))

}

