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

#### EM Algorithm  ####

K_1 <- function(sigmasq_1, pat_table) {
	return(sigmasq_1*diag(dim(pat_table)[1]))
}

K_2 <- function(sigmasq_2, pat_table) {
	 return(sigmasq_2 * exp(-abs(outer(pat_table$obs_times, pat_table$obs_times,"-"))/1))
}

K = list(K_1,K_2)

expected_terms <- function(mean_params_old, cov_params_old, theta_old, pat, table1, table2) {
		
	pat_table = table2[table2$id == pat,]
	cens = table1$cens[table1$id == pat]
	
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
	
	quad_XS_K1_SX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		if(dim(pat_table)[1] == 1) {
			return((X_T)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(X_T)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%(X_T))
		}
	}	

	quad_XS_K2_SX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		if(dim(pat_table)[1] == 1) {
			return((X_T)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(X_T)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%(X_T))
		}
	}	

	quad_YS_K1_SX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		Y = pat_table$obs
		if(dim(pat_table)[1] == 1) {
			return((Y)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(Y)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%(X_T))
		}
	}	
	
	quad_YS_K2_SX_calc <- function(T) {
		X_T = Cov(T, pat_table)
		Y = pat_table$obs
		if(dim(pat_table)[1] == 1) {
			return((Y)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%(X_T))
		}
		else{
			return(t(Y)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%(X_T))
		}
	}	
	

	
	
	if(cens == 1) {

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params_old, cov_params_old, theta_old, pat_table))

	eval_T = seq(c, c+10, by = 0.1)
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)
	K1 = K_1(cov_params[1], pat_table)
	K2 = K_2(cov_params[2], pat_table)


	quad_XSX = Vectorize(quad_XSX_calc)(eval_T)
	quad_YSX = Vectorize(quad_YSX_calc)(eval_T)	

	quad_XS_K1_SX = Vectorize(quad_XS_K1_SX_calc)(eval_T)	
	quad_XS_K2_SX = Vectorize(quad_XS_K2_SX_calc)(eval_T)	

	quad_YS_K1_SX = Vectorize(quad_YS_K1_SX_calc)(eval_T)	
	quad_YS_K2_SX = Vectorize(quad_YS_K2_SX_calc)(eval_T)	

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)
	
	d = sqrt(dim(quad_XSX)[1])
	
	return(list("exp_T" = sum(eval_T*norm_probs),
	"exp_quad_XSX" = matrix(quad_XSX%*%norm_probs, nrow = d, ncol = d), 
	"exp_quad_YSX" = matrix(quad_YSX%*%norm_probs, nrow = 1, ncol = d), 
	"exp_quad_XS_K1_SX" = matrix(quad_XS_K1_SX%*%norm_probs, nrow = d, ncol = d), 
	"exp_quad_XS_K2_SX" = matrix(quad_XS_K2_SX%*%norm_probs, nrow = d, ncol = d), 
	"exp_quad_YS_K1_SX" = matrix(quad_YS_K1_SX%*%norm_probs, nrow = 1, ncol = d), 
	"exp_quad_YS_K2_SX" = matrix(quad_YS_K2_SX%*%norm_probs, nrow = 1, ncol = d)
	))

	}
	
	if (cens == 0) {
	T = table1$survival[table1$id == pat]
	
	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)
	K1 = K_1(cov_params[1], pat_table)
	K2 = K_2(cov_params[2], pat_table)
	
	quad_XSX = quad_XSX_calc(T)
	quad_YSX = quad_YSX_calc(T)
	
	quad_XS_K1_SX = quad_XS_K1_SX_calc(T)
	quad_XS_K2_SX = quad_XS_K2_SX_calc(T)

	quad_YS_K1_SX = quad_YS_K1_SX_calc(T)
	quad_YS_K2_SX = quad_YS_K2_SX_calc(T)
	
	return(list("exp_T" = sum(eval_T*norm_probs),
	"exp_quad_XSX" = quad_XSX, 
	"exp_quad_YSX" = quad_YSX,
	"exp_quad_XS_K1_SX" = quad_XS_K1_SX, 
	"exp_quad_XS_K2_SX" = quad_XS_K2_SX, 
	"exp_quad_YS_K1_SX" = quad_YS_K1_SX, 
	"exp_quad_YS_K2_SX" = quad_YS_K2_SX
	))
		
	}
	
}	

info_exp <- function(mean_params_old, cov_params_old, theta_old, em_mean_params, pat, table1, table2) {
		
	pat_table = table2[table2$id == pat,]
	cens = table1$cens[table1$id == pat]

	dW_dsigma <- function(j,l) {
		h <- function(T) {
			X_T = Cov(T, pat_table)
			W_T = pat_table$obs - X_T%*%mean_params
			Y = pat_table$obs

			K_j = eval(parse(text=paste("K",j,sep="")))
			K_l = eval(parse(text=paste("K",l,sep="")))			
			dbeta_dsigma_l = eval(parse(text=paste("dbeta_dsigma",l,sep="")))
		
			term1 = (t(X_T%*%dbeta_dsigma_l)%*%Inv_Sigma%*%K_j + t(W_T)%*%Inv_Sigma%*%K_l%*%Inv_Sigma%*%K_j)%*%Inv_Sigma%*%W_T
	
			term2 = t(W_T)%*%Inv_Sigma%*%K_j%*%(Inv_Sigma%*%K_l%*%Inv_Sigma%*%W_T+Inv_Sigma%*%X_T%*%dbeta_dsigma_l)
			term1+term2	
			return(-(term1+term2))
		}	
		return(h)
	}

	K1_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%em_mean_params
		return(t(W_T)%*%Inv_Sigma%*%K1%*%Inv_Sigma%*%W_T)
	}		

	K2_calc <- function(T) {
		X_T = Cov(T, pat_table)
		W_T = pat_table$obs - X_T%*%em_mean_params
		return(t(W_T)%*%Inv_Sigma%*%K2%*%Inv_Sigma%*%W_T)
	}	


	Sigma = Sigma_calc(cov_params,pat_table)
	Inv_Sigma = solve(Sigma)
	K1 = K_1(cov_params[1], pat_table)
	K2 = K_2(cov_params[2], pat_table)
	

	dW_dsigma11 = dW_dsigma(1,1)
	dW_dsigma22 = dW_dsigma(2,2)
	dW_dsigma12 = dW_dsigma(1,2)	

	if(cens == 1) {

	c = table1$survival[table1$id == pat]

	cond_dens = Vectorize(h(mean_params_old, cov_params_old, theta_old, pat_table))

	eval_T = seq(c, c+10, by = 0.1)	
	
	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)
	
	K1_vec = Vectorize(K1_calc)(eval_T)
	K2_vec = Vectorize(K2_calc)(eval_T)

	
	return(list(
	"trace_K1_K1" = sum(diag(Inv_Sigma%*%K1%*%Inv_Sigma%*%K1)),
	"trace_K2_K2" = sum(diag(Inv_Sigma%*%K2%*%Inv_Sigma%*%K2)),	
	"trace_K1_K2" = sum(diag(Inv_Sigma%*%K1%*%Inv_Sigma%*%K2)),
	"dW_dsigma_1_1" = sum(Vectorize(dW_dsigma11)(eval_T)*norm_probs),
	"dW_dsigma_2_2" = sum(Vectorize(dW_dsigma22)(eval_T)*norm_probs),
	"dW_dsigma_1_2" = sum(Vectorize(dW_dsigma12)(eval_T)*norm_probs),
	"exp_K1" = sum(K1_vec*norm_probs),
	"exp_K2" = sum(K2_vec*norm_probs),
	"trace_K1" = sum(diag(Inv_Sigma%*%K1)),
	"trace_K2" = sum(diag(Inv_Sigma%*%K2))
	))

	}
	
	if (cens == 0) {
	T = table1$survival[table1$id == pat]
	
	K1_vec = K1_calc(T)
	K2_vec = K2_calc(T)

	
	return(list(
	"trace_K1_K1" = sum(diag(Inv_Sigma%*%K1%*%Inv_Sigma%*%K1)),
	"trace_K2_K2" = sum(diag(Inv_Sigma%*%K2%*%Inv_Sigma%*%K2)),	
	"trace_K1_K2" = sum(diag(Inv_Sigma%*%K1%*%Inv_Sigma%*%K2)),
	"dW_dsigma_1_1" = dW_dsigma11(T),
	"dW_dsigma_2_2" = dW_dsigma22(T),
	"dW_dsigma_1_2" = dW_dsigma12(T),
	"exp_K1" = sum(K1_vec),
	"exp_K2" = sum(K2_vec),
	"trace_K1" = sum(diag(Inv_Sigma%*%K1)),
	"trace_K2" = sum(diag(Inv_Sigma%*%K2))
	))
		
	}
	
	

}

### EM Algorithm ###

# Initial Estimates
theta_old = 1/mean( table1$survival[table1$cens == 0])
mean_params_old = model2$beta
cov_params_old = model2$sigma


# EM Algo #


YSX_sum = 0
XSX_sum = 0
num_pats = length(levels(as.factor(table1$id)))
Denominator = 0
Score1 = 0
Score2 = 0
total_quad_X_K1 = 0
total_quad_X_K2 = 0
total_quad_Y_K1 = 0
total_quad_Y_K2 = 0

for (pat in table1$id) {
	Denominator = Denominator + table1$survival[table1$id == pat] + (table1$cens[table1$id == pat] == 1)/theta_old
	
	exp_pat = expected_terms(mean_params_old, cov_params_old, theta_old, pat, table1, table2)
	
	YSX_sum = YSX_sum + exp_pat$exp_quad_YSX
	XSX_sum = XSX_sum + exp_pat$exp_quad_XSX	
	total_quad_X_K1 = total_quad_X_K1 + exp_pat$exp_quad_XS_K1_SX
	total_quad_X_K2 = total_quad_X_K2 + exp_pat$exp_quad_XS_K2_SX
	total_quad_Y_K1 = total_quad_Y_K1 + exp_pat$exp_quad_YS_K1_SX
	total_quad_Y_K2 = total_quad_Y_K2 + exp_pat$exp_quad_YS_K2_SX	
}

dbeta_dsigma1 = solve(XSX_sum,total_quad_K1)%*%solve(XSX_sum,t(YSX_sum)) - solve(XSX_sum,t(total_quad_Y_K1))

dbeta_dsigma2 = solve(XSX_sum,total_quad_K2)%*%solve(XSX_sum,t(YSX_sum)) - solve(XSX_sum,t(total_quad_Y_K2))

em_mean_params = solve(XSX_sum,t(YSX_sum))
em_theta = num_pats / Denominator

info = add_info = matrix(0,nrow = 2, ncol = 2)

for (pat in table1$id) {
	
	exp_pat = info_exp(mean_params_old, cov_params_old, theta_old, em_mean_params, pat, table1, table2)

	Score1 = Score1 + -exp_pat$trace_K1/2+exp_pat$exp_K1/2
	Score2 = Score2 + -exp_pat$trace_K2/2+exp_pat$exp_K2/2

	add_info[1,1] = exp_pat$trace_K1_K1/2-exp_pat$dW_dsigma_1_1/2
	add_info[1,2] = add_info[2,1] = exp_pat$trace_K1_K2/2-exp_pat$dW_dsigma_1_2/2
	add_info[2,2] = exp_pat$trace_K2_K2/2-exp_pat$dW_dsigma_2_2/2
	info = info + add_info
}

Score = c(Score1,Score2)

em_cov_params = cov_params_old + solve(info)%*%Score

em_cov_params = as.numeric(lapply(em_cov_params, function(x){max(x,0)}))

cbind(mean_params_old,em_mean_params)
cbind(cov_params_old,em_cov_params)
cbind(theta_old,em_theta)

max(c(abs(em_mean_params-mean_params_old),
abs(em_cov_params-cov_params_old),
abs(em_theta-theta_old)))

mean_params_old = em_mean_params
cov_params_old = em_cov_params
theta_old = em_theta
