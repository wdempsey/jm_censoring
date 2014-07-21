fit.piecewiseph <- function(Table_1, Table_2, Cov, Sigma_Calc, K, params, control) {
	
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
gamma = params$gamma
num_varcomp = length(K)

# Number of Breaks
num_breaks = 7

censored = as.logical(Table_1$cens)
break_points = c(0,quantile(Table_1$survival[!censored],probs = seq(0,1, length.out = num_breaks+1))[-c(1,num_breaks+1)],Inf)

theta = vector(length = length(break_points)-1)

min_t <- function(t) {
  min_tx <- function(x) {
    return(min(x,t))
  }
  return(min_tx)
}

for(j in 2:(length(theta)+1)) {
  denom = sum(unlist(lapply(Table_1$survival,min_t(break_points[j]))) - unlist(lapply(Table_1$survival,min_t(break_points[j-1]))))
  numer = length(which(Table_1$survival < break_points[j] & Table_1$survival >= break_points[j-1] & Table_1$cens == 0))
  theta[j-1] = numer/denom
}

# Fix W = 0 For Now
W = 0

list.params <- list(mean_params = mean_params, cov_params = cov_params, theta = theta, gamma = gamma)

params <- unlist(as.relistable(list.params))

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

baseline_hazard <- function(theta, break_points) {
	# Returns Baseline Hazard Function
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	
	baseline_hazard_t <- function (t) {
	  return(theta[max(which(break_points <= t))])
	}
	
	return(baseline_hazard_t)
}

log_surv <- function(theta, gamma, W) {
	# Log Survival Distribution
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	# W is the baseline covariates
	
	log_surv_t <- function(t) {
		return(
			-integrate(Vectorize(baseline_hazard(theta, break_points)),0,t)$value * exp(gamma%*%W)
		)	
	}	
	return(log_surv_t)
}

log_f <- function(theta, gamma,W, break_points) { 
	# Assume a Weibull Distribution With Parameters k and lambda \in theta
	# W is the baseline covariates
	
	log_f_t <- function(t) {
		
		return(log_surv(theta,gamma,W)(t) + gamma%*%W + log(baseline_hazard(theta,break_points)(t)) )
		
	}
	
	return(log_f_t)
	
}

h <- function(mean_params, cov_params, theta, gamma, W, obs, pat_table) {
	# Returns a function that returns value of the joint density 
	# at T = t given the parameter values
	
	h_t <- function(t) {
		t_dens = log_f(theta,gamma, W, break_points)
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
	return(log_f(theta, gamma, W,break_points)(T) + log_g(mean_params, cov_params, obs, pat_table)(T) )
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
deriv_theta <- function(theta, gamma, W) {
	deriv_theta_t <- function(t) {
      
	    minimum_tx = unlist(lapply(break_points,min_t(t)))
	    indicator_t = vector(length = length(theta))
	    for (j in 2:length(break_points)) {
	        indicator_t[j-1] = ((break_points[j-1] <= t) & (break_points[j] > t))
	    }
      
		diff = minimum_tx[2:(length(minimum_tx))] - minimum_tx[1:(length(minimum_tx)-1)]
      
		integral = exp(gamma%*%W) * diff
			
		return(1/theta * indicator_t  - integral)
		}
	return(deriv_theta_t)
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
		
	theta_grad = Vectorize(deriv_theta(theta,gamma,W))(eval_T)
	
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
	"theta_grad" = theta_grad%*%norm_probs
	))
	
}	


grad_calc <- function(params, table1 = Table_1, table2 = Table_2) {

  params <- relist(params, skeleton = list.params)
    
	mean_params = params$mean_params
	cov_params = params$cov_params
	theta = params$theta
	gamma = params$gamma
	
	patients = table1$id

  grad_theta = 0
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
			
			grad_theta = grad_theta + deriv_theta(theta,gamma, W)(T)
							
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
			grad_theta = grad_theta + exp_terms$theta_grad
				
			# Beta
			grad_beta = grad_beta + exp_terms$exp_S
				
			# Sigma					
			for(i in 1:num_varcomp) {	
				grad_sigma[[i]] = grad_sigma[[i]] - exp_terms$trace_K[[i]]/2 + exp_terms$exp_K[[i]]/2
			}
		}
	}
	return(-c(grad_beta,unlist(grad_sigma), grad_theta,0))
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

max_theta = 10*theta
min_theta = theta/10

der = grad(log_lik, unlist(params))
der_calc = grad_calc(unlist(params))

if(control$fixed == FALSE) {
  print('Theta is Not Fixed For Optimization')
  inits <- params
  op_llik <- optim(inits, log_lik, grad_calc, lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), min_theta, -Inf), upper = c(rep(Inf, length(mean_params) + length(cov_params)), max_theta,Inf), hessian = TRUE, control = list(trace = 1, maxit = 500))  
}
if(control$fixed == TRUE) {
  print('Theta is Fixed For Optimization')
  inits <- params_fixed
  op_llik <- optim(inits, log_lik_fixed(theta), grad_calc_fixed(theta), lower = c(rep(-Inf,length(mean_params)), rep(0, length(cov_params)), -Inf), upper = c(rep(Inf, length(mean_params) + length(cov_params)),Inf), control = list(trace = 1, maxit = 500))  
}

print(op_llik$convergence)
print('Finished the Optimization Code')

# Return the MLE Estimates

return(list("mle" = op_llik$par, "hess" = op_llik$hessian, "conv" = op_llik$convergence))

}

