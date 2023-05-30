##########################
# Functions	
##########################
simulate_from_orcl<-function(N,location,scale){
### What is it for? ###
# This is a proposal distribution to simulate the parameter values for 'linear_coef' from Cauchy(location,scale).

### In and Out ###
# Input: N - the Monte Carlo size; location - a p-dimensional vector; scale - a p*p matrix.  
# Output: a N*p matrix with rows containing linear_coef

	p<-length(location)
	simulation_of_cauchy<-matrix(rcauchy(N*p,location=0,scale=1),p,N)
	out<-t(scale%*%simulation_of_cauchy+location)
	return(out)
}

density_orcl<-function(x,location,scale){
### What is it for? ###
# Probability density of x simulated using 'simulate_from_orcl'

### In and Out ###
# Input: x - a N*p matrix with rows containing linear_coef; location - a p-dimensional vector; scale - a p*p matrix.  
# Output: a N*1 vector being the Cauchy density at rows of x

	simulation_of_cauchy<-solve(scale)%*%(t(x)-location) #p*N
	out<-c(apply(dcauchy(simulation_of_cauchy),2,prod)) #N*1
	return(out)
}


simulate_from_propsl<-function(N,subgroup_est,subgroup_ind,design_matr){
### What is it for? ###
# Simualte N values for the linear coefficients beta from the data-dependent r_n which is a uniformly weighted mixture of the Cauchy distribution centered at subgroup_est and scaled by sd matrix of subgroup_est, such as variance matrix of LSE.

### In and Out ###
# Input: N - the Monte Carlo size. subgroup_est - a K*p matrix containing the set of point estimates of observation subsets. subgroup_ind - a list containing the indices of the subgroup observations. design_matr - the n*p design matrix of the data
# Output: a N*p matrix with rows containing the parameter.
	
	K<-nrow(subgroup_est)
	param_p<-ncol(subgroup_est)
	mixture_means<-subgroup_est
	components_to_simulate<-sample(1:K,size=N,replace=T)
	out<-NULL
	for(j in 1:K){
		n_j<-sum(components_to_simulate==j)
		subgroup_desgn<-design_matr[subgroup_ind[[j]],]
		subgroup_scale<-sqrtmatrix(solve(t(subgroup_desgn)%*%subgroup_desgn))
		simulation_of_cauchy<-matrix(rcauchy(n_j*param_p,location=0,scale=1),param_p,n_j) 
		out<-rbind(out,t(subgroup_scale%*%simulation_of_cauchy+subgroup_est[j,]))
	}
	return(out)
}

density_propsl<-function(x,subgroup_est,subgroup_ind,design_matr){
### What is it for? ###
# Probability density of x simulated using 'simulate_from_propsl'

### In and Out ###
# Input: x - a N*p matrix with rows containing linear_coef; others are the same as those for 'simulate_from_propsl'
# Output: a N*1 vector being the Cauchy density at rows of x

	K<-nrow(subgroup_est)
	param_p<-ncol(subgroup_est)
	mixture_means<-subgroup_est
	out_matr<-matrix(0,N,K)
	for(j in 1:K){
		subgroup_desgn<-design_matr[subgroup_ind[[j]],]
		subgroup_scale_inv<-sqrtmatrix(t(subgroup_desgn)%*%subgroup_desgn)
		simulation_of_cauchy<-subgroup_scale_inv%*%(t(x)-subgroup_est[j,]) #p*N
		out_matr[,j]<-c(apply(dcauchy(simulation_of_cauchy),2,prod)) #N*1
	}
	out<-rowMeans(out_matr)
	return(out)
}

simulate_pseudo_datasets<-function(linear_coef,n,X){
### What is it for? ###
# Simulate N pseudo datasets from the Cauchy regression model using the design matrix X, each using a column of 'linear_coef' which is p-dimensional linear coefficients. Each dataset contains n independent observations i.i.d from X*beta+Cauchy error.

### Specification of Implementation ###
# 1. This function should not be used when N*n is large, say over 10^7. Because the set of data is stored in an N*n matrix.

### In and Out ###
# Input: linear_coef - a N*p matrix with columns being the linear coefficients. n - the number of independent observation in each dataset. X - the n*p design matrix
# Output: an N*n matrix with rows containing each dataset.

	N<-nrow(linear_coef)
	Cauchy_error<-matrix(rcauchy(N*n),n,N)
	pseudo_dataset<-X%*%t(linear_coef)+Cauchy_error # N*n
	pseudo_dataset<-t(pseudo_dataset)
	return(pseudo_dataset)
}

cal_summaries<-function(datasets,X){
### What is it for? ###
# Calculate the summary statistic for Cauchy regression model. The least square estimator is chosen as the summary statistic. All input datasets use the same design matrix.

### In and Out ###
# Input: datasets - an N*n matrix with rows containing each dataset; X - the n*p design matrix
# Output: an N*p matrix with cols containing the summary statistic
	
	N<-nrow(datasets) # the Monte Carlo size
	XTX_inv<-solve(t(X)%*%X)
	simulated_summaries<-t(XTX_inv%*%t(X)%*%t(datasets))
	return(simulated_summaries)
}

simulate_pseudo_summaries<-function(parameter_values,n,X){
### What is it for? ###
# Simulate pseudo summaries given a set of parameter values (linear_coef) Each dataset contains n independent observations i.i.d from X*beta+Cauchy error. This is used for 'simpleDR_in_paral' function.


### In and Out ###
# Input: parameter_values - a N*p matrix with rows being the linear coefficients. n - the number of independent observation in each dataset. X - the n*p design matrix. 
# Output: an N*p matrix with cols containing the summary statistic

	pseudo_dataset<-simulate_pseudo_datasets(parameter_values,n,X)
	pseudo_summaries<-cal_summaries(pseudo_dataset,X=X)
	return(pseudo_summaries)
}

ABC_ACC_results_summary<-function(results,type='ABC',w=NULL,true_params=NULL,param_est,CI_alpha=NULL,XTX){
### What is it for? ###
# Summaries the results of ABC/ACC algorithms by calculating the mean estimates, variance estimates, ESS of importance weights, two-sided credible/confidence bounds, one-sided credible/confidence bound, success/failure of covering true parameters by each interval, and the covering by two-dimension credible/confidence. 

### Specification of Implementation ###
# 1. For one-dimension covering, confidence bounds are calculated for each parameter.
# 2. For two-dimension covering, the pairs of parameters for which confidence regions are calculated need to be specified.

### In and Out ###
# Input: results - a list containing five elements: "parameters" N*p matrix, "summaries","acceptance_ind","tolerance" and "acceptance_rate", which usually is the output of 'rejABC'/'regABC' functions. type - either 'ABC' or 'ACC'. weighted/w - whether or not the ABC output are weighted, and the vector of importance weights. true_params - the vector of true parameter, used to calculate the success/failure of being coverred. XTX - used in Cauchy regression example for calculting the distance needed for confidence region
# Output: a list containing the above mentioned results.

	coverage<-0; wid_vol<-0
	N<-nrow(results$parameters)
	quant_vec<-rep(0,length(param_est))
	if(!is.null(w[1])){
		mean_vec<-wmean(results$parameters,w)
		cov_vec<-wvar(results$parameters,w)
		var_vec<-diag(cov_vec)

		IS_ESS<-ESS(w)/length(w)
		w_ks<-w/sum(w)*length(w)
		# for(param_i in 1:length(param_est)) quant_vec[param_i]<-wquant(results$parameters[,param_i],q=1-(1-CI_alpha)/2,weights=w) # For one-dimension
	}
	if(is.null(w[1])){
		mean_vec<-colMeans(results$parameters)
		cov_vec<-var(results$parameters)
		var_vec<-diag(cov_vec)
		IS_ESS<-1
		w_ks<-NULL
		# for(param_i in 1:length(param_est)) quant_vec[param_i]<-quantile(results$parameters[,param_i],probs=1-(1-CI_alpha)/2) # For one-dimension
	}
	if(type=='ABC'){
		# Change to multi-variate version
		param_twoD_coverred<-true_params
	}
	if(type=='ACC'){
		# Change to multi-variate version
		# param_twoD_coverred<-c(2*mean_vec[1]-true_params[1],true_params[2])
		param_twoD_coverred<-2*mean_vec-true_params
	}
	if(length(param_est)==2){
		# Change to multi-variate version
		coverage<-HDR_2d_coverred(param_twoD_coverred,x=results$parameters[,1],y=results$parameters[,2],alpha=CI_alpha,weights=w_ks)
		wid_vol<-pi*2*qf(1-CI_alpha,2,N-1)*det(sqrtmatrix(cov_vec)) # As if it's multivaraite normal
		coverage1<-0
	}
	if(length(param_est)>2){
		beta_dist<-function(beta_val,mean_vec,scale) return(t(beta_val-mean_vec)%*%scale%*%(beta_val-mean_vec)) # This is the depth function for the confidence region of regression linear coefficients
		XTX_inv<-solve(XTX)
		beta_disct_vec<-apply(results$parameters,1,beta_dist,mean_vec=mean_vec,scale=XTX_inv) 
		beta_dist_0<-beta_dist(true_params,mean_vec,XTX_inv)
		beta_dist_1<-beta_dist(param_twoD_coverred,mean_vec,XTX_inv)
		coverage<-(beta_dist_0<quantile(beta_disct_vec,CI_alpha))
		coverage1<-(beta_dist_1<quantile(beta_disct_vec,CI_alpha))
		wid_vol<-det(sqrtmatrix(cov_vec)) # As if it's multivaraite normal
	}
return(list(mean_vec=mean_vec,var_vec=var_vec,quant_vec=quant_vec,IS_ESS=IS_ESS,coverage=coverage,coverage1=coverage1,wid_vol=wid_vol))
} 


######################################
# Main function running experiment	
######################################
Example1_main<-function(test_size,parameter_true,tested_observations_all,design_matr,parameter_setting,minibatch=FALSE,N,platform,n_all,p_all,CI_alpha=0.05,divide=FALSE,divide_paral=FALSE,prior_choice='jeffery'){
### What is it for? ###
# Run the experiment in Cauchy regression example of Thorton, Li and Xie (2018) to compare ABC and ACC under different parameter/summary settings

### Specification of Implementation ###
# 1. Functions 'onesideCI_cal' only calculates confidence interval, so the coverage rate only works for univariate unknown parameter

### In and Out ###
# Input: test_size - number of different datasets tested in the experiment. parameter_true - the value of parameter with which the dataset is generated. tested_observations_all - the test_size*n matrix containing all datasets. design_matr - a list containing n*p design matrices for all test datasets. parameter_setting: 1 -- 5 covariates. minibatch - if TRUE, use the minibatch scheme to construct r_n; if FALSE, use the oracle proposal. N - Simulation size in ABC/ACC. platform - 'Mac' or 'Win'. n_all - data sizes to be tested. p_all - acceptance rates to be used in ABC/ACC. divide - if TRUE, dividing-recombining (DR) strategy will be applied to the simulation of summary statistics of pseudo datasets in order to limit the use of memory. divide_paral - if TRUE, each part in the DR implementation will run in parallel; be careful in the cost of overhead. 
# Output: A four-element list containing the posterior mean estimates, posterior variance estimates, coverage rate of the true parameter value in the test_size runs, and the tolerance values for all runs

	### Initialisation ###
	post_means_all<-list(0); post_variances_all<-list(0); post_coverage_all<-list(0); post_quant_all<-list(0)
	param_p<-length(parameter_true)
	tmp1<-matrix(0,test_size,param_p)

	for(method_i in 1:4){
		post_means_all[[method_i]]<-list(0)
		if(method_i==1) names(post_means_all)[method_i]<-'ISABC'
		if(method_i==2) names(post_means_all)[method_i]<-'ISregABC'
		if(method_i==3) names(post_means_all)[method_i]<-'rejACC'
		if(method_i==4) names(post_means_all)[method_i]<-'regACC'
		for(n_i in 1:length(n_all)){
			post_means_all[[method_i]][[n_i]]<-list(0)
			names(post_means_all[[method_i]])[n_i]<-paste0('n=',n_all[n_i])
			for(p_i in 1:length(p_all)){
				post_means_all[[method_i]][[n_i]][[p_i]]<-tmp1
				names(post_means_all[[method_i]][[n_i]])[p_i]<-paste0('acceptance rate=',p_all[p_i])
			}
		}
	}
	post_variances_all<-post_means_all
	post_quant_all<-post_means_all
	coverage_all<-list(list(0),list(0))
	for(method_i in 1:4) coverage_all[[method_i]]<-matrix(0,length(p_all),test_size)
	names(coverage_all)<-names(post_means_all)
	wid_vol_all<-coverage_all
	IS_ESSr_all<-coverage_all
	coverage1_all<-coverage_all

	##########################
	# Algorithms
	##########################
	### Things to be noted ###
	# 1. A N*n matrix 'simulated_observations' is stored and used to calculate 'simulated_summaries'. If N*n is larger than 10^7, this should be done by  dividing and recombining in parallel.
	# 2. avg_tolerance_all is defined by running the experiment once and finding out the average tolerance values when fixing the acceptance rate: avg_tolerance_all<-lapply(tolerance_all[[n_i]][[sumstat_i]],mean)

	if(minibatch){
		nobs<-ncol(tested_observations_all)
		group_no<-floor(nobs^(2/5))
		group_size<-floor(nobs/group_no)
		subgroup_indices<-divide(1:nobs,group_no)
		for(sub_i in 1:50) subgroup_indices[[group_no+sub_i]]<-sample(1:nobs,group_size,replace=T)
		group_no<-length(subgroup_indices)
	}

	tested_parameters<-matrix(rep(parameter_true,test_size),ncol=param_p,byrow=T) # a test_size*param_p matrix
	duration<-0
	for(n_i in 1:length(n_all)){
		n<-n_all[n_i] # the observation size
		tested_summaries_all<-cal_summaries(tested_observations_all,X=design_matr) #N*p

		for(data_i in 1:test_size){
			##########################
			#### Recrod start time ###
			start_time<-proc.time()
			##########################
			plot(0,parameter_setting,xlim=c(1,length(p_all)),ylim=c(1,5),xlab='p_i',ylab='parameter_setting',main=paste0('Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds, CI_alpha=',CI_alpha))
			data_i<-data_i
			tested_observations<-tested_observations_all[data_i,]
			tested_summary<-tested_summaries_all[data_i,]

			# Divide the dataset into subsets for the purpose of constructing initial distribution
			if(minibatch){
				subgroup_est<-matrix(0,group_no,param_p)
				for(group_i in 1:group_no){
					subgroup_obs<-tested_observations[subgroup_indices[[group_i]]]
					subgroup_desgn<-design_matr[subgroup_indices[[group_i]],]
					subgroup_XTX<-t(subgroup_desgn)%*%subgroup_desgn
					subgroup_est[group_i,]<-c(solve(subgroup_XTX)%*%t(subgroup_desgn)%*%subgroup_obs)
				} 		
			}

			# Simulate the parameter values and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
			param_est<-1:length(parameter_true)
			if(prior_choice=='jeffery') prior_densities<-rep(1,N)
			XTX<-t(design_matr)%*%design_matr

			### oracle proposal ###	
			if(!minibatch){
				group_size<-sqrt(n)
				proposal_location<-parameter_true
				proposal_scale<-solve(sqrtmatrix(XTX/nrow(design_matr)))/sqrt(group_size)
				simulated_parameters<-simulate_from_orcl(N,proposal_location,proposal_scale)
				proposal_densities<-density_orcl(simulated_parameters,proposal_location,proposal_scale)				
			}

			### mini-batch proposal ###
			if(minibatch){
				simulated_parameters<-simulate_from_propsl(N,subgroup_est,subgroup_indices,design_matr)
				proposal_densities<-density_propsl(simulated_parameters,subgroup_est,subgroup_indices,design_matr)
			}

			# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
			tested_value<-tested_parameters[data_i,]
			if(N*n>5*10^7) divide_counts<-N*n/(5*10^7)
			if(N*n<=5*10^7) divide<-FALSE
			if(divide) simulated_summaries<-simpleDR_in_paral(counts=divide_counts,x=simulated_parameters,FUN=simulate_pseudo_summaries,in_type='matrix',out_type='stacked_matrices',paral=divide_paral,platform=platform,packages=c('matrixStats'),n=n,X=design_matr)
			if(!divide){
				simulated_observations<-simulate_pseudo_datasets(simulated_parameters,n,X=design_matr) # a N*n matrix
				simulated_summaries<-cal_summaries(simulated_observations,X=design_matr) # a N*d matrix		
			}
			reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries)
			weights_of_sample<-prior_densities/proposal_densities
			for(p_i in 1:length(p_all)){
				points(p_i,parameter_setting,pch=20)
				tolerance_percentile<-p_all[p_i]

				###################################################
				# importance sampling regression ABC (ISregABC)
				###################################################
				method_i<-2
				results_ISregABC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample)
				results_ISregABC<-results_ISregABC[[1]] # This only deals with one observed dataset
				importance_weights<-results_ISregABC$weights
				results_summary<-ABC_ACC_results_summary(results_ISregABC,type='ABC',w=importance_weights,true_params=tested_value,param_est=param_est,CI_alpha=CI_alpha,XTX=XTX)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				coverage1_all[[method_i]][p_i,data_i]<-results_summary$coverage1
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol

				###################################################
				# regression ACC (regACC)
				###################################################
				method_i<-4
				results_regACC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
				results_regACC<-results_regACC[[1]]
				results_summary<-ABC_ACC_results_summary(results_regACC,type='ACC',w=NULL,true_params=tested_value,param_est=param_est,CI_alpha=CI_alpha,XTX=XTX)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				coverage1_all[[method_i]][p_i,data_i]<-results_summary$coverage1
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol
			}
		##########################
		#### Recrod running time ###
		duration<-signif(proc.time()-start_time,2)
		##########################	
		}
	}
	return(list(post_means_all=post_means_all,post_variances_all=post_variances_all,post_quant_all=post_quant_all,IS_ESSr_all=IS_ESSr_all,coverage_all=coverage_all,coverage1_all=coverage1_all,wid_vol_all=wid_vol_all))
}
