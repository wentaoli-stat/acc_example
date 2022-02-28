##########################
# Functions	
##########################
simulate_from_r1<-function(N,location,scale){
### What is it for? ###
# Simualte the parameter values for 'Cauchy_location' from the data-dependent distribution r1 that Cauchy_location~Cauchy(location,scale). 

### In and Out ###
# Input: N - the Monte Carlo size. 
# Output: a N*1 matrix with columns containing Cauchy_location.
	
	simulation_of_cauchy<-rcauchy(N,location=location,scale=scale)
	outert<-as.matrix(simulation_of_cauchy)
	return(out)
}

simulate_from_r2<-function(N,subgroup_est,nonnegative=FALSE){
### What is it for? ###
# Simualte N values for an individual parameter from the data-dependent distribution r2 which is a mixture normal of the kernel density estimate(KDE) formed by point estimates of observation subsets. If the parameter is nonnegative, then its logarithmic value is simulated from the mixture normal.

### In and Out ###
# Input: N - the Monte Carlo size. subgroup_est - a K*1 matrix containing the set of point estimates of observation subsets. nonnegative - indicate whether the logarithmic transformation is needed.
# Output: a N*1 matrix with column containing the targeted parameter.
	
	K<-length(subgroup_est)
	if(!nonnegative) mixture_normal_means<-subgroup_est
	if(nonnegative) mixture_normal_means<-log(subgroup_est)
	density_est<-density(mixture_normal_means,kernel='gaussian')
	bandwidth<-density_est$bw
	components_to_simulate<-sample(1:K,size=N,replace=T)
	normal_simulation<-rnorm(N)
	simulation_of_mixture<-normal_simulation*bandwidth+mixture_normal_means[components_to_simulate] # A length N vector		
	if(!nonnegative) simulation_of_parameter<-simulation_of_mixture
	if(nonnegative) simulation_of_parameter<-exp(simulation_of_mixture)
	out<-as.matrix(simulation_of_parameter)
	return(out)
}

density_r2<-function(x,subgroup_est,nonnegative=FALSE){
# Input: x - N*1 matrix of parameter values; subgroup_est - a K*1 matrix containing the set of point estimates of observation subsets. nonnegative - indicate whether the simulated sample is transformed from logrithm.
# Output: N*1 vector of proposal density

	K<-length(subgroup_est)
	if(!nonnegative) mixture_normal_means<-subgroup_est
	if(nonnegative) mixture_normal_means<-log(subgroup_est)
	density_est<-density(mixture_normal_means,kernel='gaussian')
	bandwidth<-density_est$bw
	if(!nonnegative){
		tmp<-sqdist_matr_matr(as.matrix(x/bandwidth),as.matrix(mixture_normal_means/bandwidth)) # N*K pairwise squared distance matrices between rows of x and mixture_normal_means
		tmp<-exp(-tmp/2) # N*K matrix of unnormalised normal density	
	}
	if(nonnegative){
		tmp<-sqdist_matr_matr(as.matrix(log(x)/bandwidth),as.matrix(mixture_normal_means/bandwidth)) 
		tmp<-exp(-tmp/2)/x
	}
	norm_constant<-1/(2*pi)^(1/2)/bandwidth
	tmp<-rowSums(tmp/K)*norm_constant # N*1 vector of mixture normal density		
	return(tmp)
}


simulate_pseudo_datasets<-function(parameter_values,n){
### What is it for? ###
# Simulate pseudo datasets given a set of parameter values (Cauchy_location,Cauchy_scale) Each dataset contains n observations i.i.d from Cauchy(Cauchy_location,Cauchy_scale).

### Specification of Implementation ###
# 1. This function should not be used when N*n is large, say over 10^7. Because the set of data is stored in an N*n matrix.
# 2. The set of data is created by multiplying the Cauchy scale parameter, tau, to the Cauchy noise and adding the Cauchy means.

### In and Out ###
# Input: parameter_values - a N*2 matrix with columns Cauchy_location and Cauchy_scale. n - the number of i.i.d observation in each dataset.
# Output: an N*n matrix with rows containing each dataset.
	
	N<-nrow(parameter_values) # the Monte Carlo size 
	Cauchy_location<-parameter_values[,1]
	Cauchy_scale<-parameter_values[,2]
	noise_in_observations<-matrix(rcauchy(N*n),nrow=N)
	pseudo_dataset<-noise_in_observations*Cauchy_scale+Cauchy_location
	return(pseudo_dataset)
}

cal_summaries<-function(datasets,setting,platform='Mac',paral=TRUE){
### What is it for? ###
# Calculate summary statistics for multiple given datasets. The summary statistics are selected according to Thorton, Li and Xie (2018), and depend on the parameter settings.

### Specification of Implementation ###
# 1. The function 'rowMedians' in <matrixStats> is used.

### In and Out ###
# Input: datasets - an N*n matrix with rows containing each dataset; setting -- choose from 1-4 the specification of which can be seen in the 'Example1_main' function for details.  function
# Output: an N*1 matrix with rows containing the summary statistic
	
	N<-nrow(datasets) # the Monte Carlo size
	
	if(setting==1){
		sample_medians<-rowMedians(datasets);
		simulated_summaries<-cbind(sample_medians=sample_medians)
	}
	if(setting==2){
		sample_means<-rowMeans(datasets);
		simulated_summaries<-cbind(sample_means=sample_means)
	}
	if(setting==3){
		sample_sds<-sqrt(rowVars(datasets))
		simulated_summaries<-cbind(sample_sds=sample_sds)
		# sample_mads<-rowWeightedMads(datasets);
		# simulated_summaries<-cbind(sample_mads=sample_mads)
	}
	if(setting==4){
		sample_means<-rowMeans(datasets);
		sample_mads<-rowWeightedMads(datasets);
		simulated_summaries<-cbind(sample_means=sample_means,sample_mads=sample_mads)	
	}
	return(simulated_summaries)
}

simulate_pseudo_summaries<-function(parameter_values,n,setting){
### What is it for? ###
# Simulate pseudo summaries given a set of parameter values (Cauchy_location,Cauchy_scale) Each dataset contains n observations i.i.d from Cauchy(Cauchy_location,Cauchy_scale). This is used for 'simpleDR_in_paral' function.


### In and Out ###
# Input: parameter_values - a N*2 matrix with columns Cauchy_location and Cauchy_scale. n - the number of i.i.d observation in each dataset. setting -- choose from 1-4; see 'Example1_main' function for details. 
# Output: an N*d matrix with rows containing each summary.

	pseudo_dataset<-simulate_pseudo_datasets(parameter_values,n)
	pseudo_summaries<-cal_summaries(pseudo_dataset,setting=setting,paral=FALSE)
	return(pseudo_summaries)
}

onesideCI_cal<-function(ACC_values,alpha,weights=NULL){
### What is it for? ###
# Calculate the 1-dimension one sided confidence interval using ACC sample. The interval (-\infty,x) covers the true parameter value with probbility alpha.

### Specification of Implementation ###
# 1. This calculates CI for one-dimension ACC_values

### In and Out ###
# Input: ACC_values - a N*p matrix. alpha - a scalar between 0 and 1. weights - a length N vector if the sample is weighted.
# Output: a scalar indicating the upper bound of the CI

	if(is.null(weights)) return(2*mean(ACC_values)-quantile(ACC_values,probs=1-alpha))
	if(!is.null(weights)) return(2*wmean(ACC_values,weights)-wquant(ACC_values,q=1-alpha,weights=weights))

}

test_percentile_cal<-function(ACC_parameter_values,theta0){
### What is it for? ###
# Calculate the proportion of 2\theta_hat-\theta_acc less than true parameter value, and check whether it is uniformly distributed

### In and Out ###
# Input: ACC_parameter_values - a length N vector. theta0 - true parameter value.
# Output: a scalar indicating the proportion

	mean(2*mean(ACC_parameter_values)-ACC_parameter_values<theta0)
}

ABC_ACC_results_summary<-function(results,type='ABC',w=NULL,true_params=NULL,param_est,CI_alpha=NULL,post_sample=NULL,method_i=0){
### What is it for? ###
# Summaries the results of ABC/ACC algorithms by calculating the mean estimates, variance estimates, ESS of importance weights, two-sided credible/confidence bounds, one-sided credible/confidence bound, success/failure of covering true parameters by each interval, and the covering by two-dimension credible/confidence. 

### Specification of Implementation ###
# 1. For one-dimension covering, confidence bounds are calculated for each parameter.
# 2. For two-dimension covering, the pairs of parameters for which confidence regions are calculated need to be specified.

### In and Out ###
# Input: results - a list containing five elements: "parameters", "summaries","acceptance_ind","tolerance" and "acceptance_rate", which usually is the output of 'rejABC'/'regABC' functions. type - either 'ABC' or 'ACC'. weighted/w - whether or not the ABC output are weighted, and the vector of importance weights. true_params - the vector of true parameter, used to calculate the success/failure of being coverred. oneD_covering/twoD_covering - whether or not the one-dimension/two-dimension credible/confidence regions should be calculated. twoD_pairs - a K*2 matrix specifying the confidence region of which pairs are calculated
# Output: a list containing the above mentioned results.

	coverage<-0
	N<-nrow(results$parameters)
	quant_vec<-rep(0,length(param_est))
	if(!is.null(w[1])){
		mean_vec<-wmean(results$parameters,w)
		cov_vec<-wvar(results$parameters,w)
		var_vec<-diag(cov_vec)
		for(param_i in 1:length(param_est)) quant_vec[param_i]<-wquant(results$parameters[,param_i],q=1-(1-CI_alpha)/2,weights=w)
		IS_ESS<-ESS(w)/length(w)
		w_ks<-w/sum(w)*length(w)
	}
	if(is.null(w[1])){
		mean_vec<-colMeans(results$parameters)
		cov_vec<-var(results$parameters)
		var_vec<-diag(cov_vec)
		for(param_i in 1:length(param_est)) quant_vec[param_i]<-quantile(results$parameters[,param_i],probs=1-(1-CI_alpha)/2)
		IS_ESS<-1
		w_ks<-NULL
	}
	if(type=='ABC'){
		oneD_bound_cal<-function(params,q,weights){
			if(!is.null(weights[1])) return(wquant(params,q=q,weights=weights))
			if(is.null(weights[1])) return(quantile(params,probs=q))
		}
		for(param_i in param_est){
			twoside_bounds<-c(oneD_bound_cal(results$parameters[,param_i],q=1-(1-CI_alpha)/2,weights=w),oneD_bound_cal(results$parameters[,param_i],q=(1-CI_alpha)/2,weights=w))
			twoside_upper<-max(twoside_bounds)
			twoside_lower<-min(twoside_bounds)
			post_ecdf<-ecdf(post_sample[,param_i])
			if(length(param_est)==1){
				coverage<-post_ecdf(twoside_upper)-post_ecdf(twoside_lower)
				wid_vol<-twoside_upper-twoside_lower
			}
			if(length(param_est)==2) coverage[param_i]<-post_ecdf(twoside_upper)-post_ecdf(twoside_lower)
		}
		if(length(param_est)==2){
			coverage<-mean(coverage)
			wid_vol<-pi*2/N*qf(1-CI_alpha,2,N-1)*det(sqrtmatrix(cov_vec)) # As if it's multivaraite normal
		}
	}
	if(type=='ACC'){
		if(length(param_est)==1){
			if(param_est==1){
			 	oneD_bound_cal<-function(params,q,weights){
					if(is.null(weights)) return(2*mean(params)-quantile(params,probs=q))
					if(!is.null(weights)) return(2*wmean(params,weights)-wquant(params,q=q,weights=weights))			 	
			 	}
			}
			if(param_est==2){		
				oneD_bound_cal<-function(params,q,weights){
					if(is.null(weights[1])) return(mean(params)^2/quantile(params,probs=q))
					if(!is.null(weights[1])) return(wmean(params,weights)^2/wquant(params,probs=q,weights=weights))
				}
			}	
			for(param_i in param_est){
				twoside_bounds<-c(oneD_bound_cal(results$parameters[,param_i],q=1-(1-CI_alpha)/2,weights=w),oneD_bound_cal(results$parameters[,param_i],q=(1-CI_alpha)/2,weights=w))
				twoside_upper<-max(twoside_bounds)
				twoside_lower<-min(twoside_bounds)
				coverage<-(true_params[param_i]<=twoside_upper)&(true_params[param_i]>=twoside_lower)
				wid_vol<-twoside_upper-twoside_lower
			}
		}
		if(length(param_est)==2){
			param_twoD_coverred<-c(2*mean_vec[1]-true_params[1],true_params[2])
			coverage<-HDR_2d_coverred(param_twoD_coverred,x=results$parameters[,1],y=results$parameters[,2],alpha=CI_alpha,weights=w_ks)
			wid_vol<-pi*2/N*qf(1-CI_alpha,2,N-1)*det(sqrtmatrix(cov_vec)) # As if it's multivaraite normal
		}
	}

# if(N==2500&(method_i==2||method_i==4)) browser()
# par(mfrow=c(6,3))
# c
# plot(x=results$parameters[,1],y=results$parameters[,2],pch=20,main=type)
# points(true_params[1],true_params[2],pch=21,col=2)

	# return(list(mean_vec=mean_vec,var_vec=var_vec,IS_ESS=IS_ESS,twoside_CI=twoside_CI,twoD_CR=twoD_CR))
	return(list(mean_vec=mean_vec,var_vec=var_vec,quant_vec=quant_vec,IS_ESS=IS_ESS,coverage=coverage,wid_vol=wid_vol))
} 


Process_all_output<-function(indices_all,output_all,param_est){
### What is it for? ###
# Summarise each part of output_all using some criterion function and list in a table. 

### Specification of Implementation ###
# 1. This is not a generic function.

### In and Out ###
# Input: indices_all - indices of runnings that are used to calculate the criterion. This is for the purpose of 'Monte_Carlo_stability_check'. output_all - a list of experiment variables and various outputs. 
# Output: a list

	n_all<-output_all$n; p_all<-output_all$p
	post_means_all<-output_all$post_means 
	post_variances_all<-output_all$post_variances 
	coverage<-output_all$coverage
	wid_vol<-output_all$wid_vol
	results<-output_all$results
	n_i<-1
	indices_groups<-divide(indices_all,group_No)
	for(group_i in 1:group_No){
		indices<-indices_groups[[group_i]]
		for(method_i in 1:number_methods){
			for(p_i in 1:length(p_all)){
				param_count<-1
				for(param_i in param_est){
					avg_coverage<-mean(coverage[[method_i]][p_i,indices])
					avg_wid_vol<-mean(wid_vol[[method_i]][p_i,indices])
					post_means<-post_means_all[[method_i]][[n_i]][[p_i]][indices,param_i]
					post_variances<-post_variances_all[[method_i]][[n_i]][[p_i]][indices,param_i]
					results[[group_i]][[param_count]][[p_i]][method_i,1]<-avg_coverage
					results[[group_i]][[param_count]][[p_i]][method_i,2]<-avg_wid_vol
					results[[group_i]][[param_count]][[p_i]][method_i,3]<-median(post_means)
					results[[group_i]][[param_count]][[p_i]][method_i,4]<-mad(post_means)
					results[[group_i]][[param_count]][[p_i]][method_i,5]<-median(post_variances)
					results[[group_i]][[param_count]][[p_i]][method_i,6]<-mad(post_variances)
					param_count<-param_count+1				
				}
			}
		}
	}
	return(results)
}


######################################
# Main function running experiment	
######################################
Example1_main<-function(test_size,parameter_true,tested_observations_all,parameter_setting,posterior_sample=NULL,N,platform,n_all,p_all,CI_alpha=0.05,divide=FALSE,divide_paral=FALSE,prior_choice='jeffery',prior_location=NULL,prior_scale=NULL,prior_df=NULL){
### What is it for? ###
# Run the experiment in example 1 of Thorton, Li and Xie (2018) to compare three algorithms (ABC, ACC with two initial distributions) under different parameter/summary settings

### Specification of Implementation ###
# 1. Functions 'onesideCI_cal' only calculates confidence interval, so the coverage rate only works for univariate unknown parameter

### In and Out ###
# Input: test_size - number of different datasets tested in the experiment. parameter_true - the value of parameter with which the dataset is generated. parameter_setting - three settings can be used: 1 -- location unknown and scale known with median as summary, 2 -- location unknown and scale known with mean as summary, 3 -- location known and scale unknown with median absolute deviation (MAD) as summary, 4 -- location known and scale unknown with Hodge-Lehman estimate as summary, 5 -- both unknown. N - Simulation size in ABC/ACC. platform - 'Mac' or 'Win'. n_all - data sizes to be tested. p_all - acceptance rates to be used in ABC/ACC. divide - if TRUE, dividing-recombining (DR) strategy will be applied to the simulation of summary statistics of pseudo datasets in order to limit the use of memory. divide_paral - if TRUE, each part in the DR implementation will run in parallel; be careful in the cost of overhead. 
# Output: A four-element list containing the posterior mean estimates, posterior variance estimates, coverage rate of the true parameter value in the test_size runs, and the tolerance values for all runs

	### Initialisation ###
	post_means_all<-list(0); post_variances_all<-list(0); post_coverage_all<-list(0); post_quant_all<-list(0)
	tmp1<-matrix(0,test_size,2)
	twoD_covering<-FALSE; oneD_covering<-TRUE

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

	##########################
	# Algorithms
	##########################
	### Things to be noted ###
	# 1. A N*n matrix 'simulated_observations' is stored and used to calculate 'simulated_summaries'. If N*n is larger than 10^7, this should be done by  dividing and recombining in parallel.
	# 2. avg_tolerance_all is defined by running the experiment once and finding out the average tolerance values when fixing the acceptance rate: avg_tolerance_all<-lapply(tolerance_all[[n_i]][[sumstat_i]],mean)

	tested_parameters<-matrix(rep(parameter_true,test_size),ncol=2,byrow=T) # a test_size*2 matrix
	duration<-0
	for(n_i in 1:length(n_all)){
		n<-n_all[n_i] # the observation size
		tested_summaries_all<-cal_summaries(tested_observations_all,setting=parameter_setting,platform=platform)

		for(data_i in 1:test_size){
			##########################
			#### Recrod start time ###
			start_time<-proc.time()
			##########################
			plot(0,parameter_setting,xlim=c(1,length(p_all)),ylim=c(1,5),xlab='p_i',ylab='parameter_setting',main=paste0('Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds, CI_alpha=',CI_alpha))
			data_i<-data_i
			tested_observations<-tested_observations_all[data_i,]
			tested_summary<-tested_summaries_all[data_i,]
			posterior_sample_i<-posterior_sample[[data_i]]

			# Divide the dataset into subsets for the purpose of constructing initial distribution
			group_no<-sqrt(n)
			subgroup<-divide(tested_observations,group_no)

			# Calculate point estimates for each sub-datasets
			# Simulate the parameter values and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
			if(parameter_setting==1){
				param_est<-1
				subgroup_medians<-unlist(lapply(subgroup,median)) # a N*1 subgroup_medians
				simulated_parameters<-cbind(simulate_from_r2(N,subgroup_medians),tested_parameters[data_i,2])
				proposal_densities<-density_r2(simulated_parameters[,1],subgroup_medians)
				if(prior_choice=='jeffery') prior_densities<-rep(1,N)
				if(prior_choice=='Cauchy') prior_densities<-dt((simulated_parameters[,1]-prior_location)/prior_scale,df=prior_df)
			}
			if(parameter_setting==2){
				param_est<-1
				subgroup_means<-unlist(lapply(subgroup,mean)) # a N*1 vector
				simulated_parameters<-cbind(simulate_from_r2(N,subgroup_means),tested_parameters[data_i,2])
				proposal_densities<-density_r2(simulated_parameters[,1],subgroup_means)
				if(prior_choice=='jeffery') prior_densities<-rep(1,N)
				if(prior_choice=='Cauchy') prior_densities<-dt((simulated_parameters[,1]-prior_location)/prior_scale,df=prior_df)
			}
			if(parameter_setting==3){
				param_est<-2
				subgroup_mads<-unlist(lapply(subgroup,mad)) # a N*1 subgroup standard deviation
				subgroup_mads<-subgroup_mads/1.4826 # The R function calculating MAD is the product of MAD and the constant 1.4826.
				simulated_parameters<-cbind(tested_parameters[data_i,1],simulate_from_r2(N,subgroup_mads,nonnegative=TRUE))
				proposal_densities<-density_r2(simulated_parameters[,2],subgroup_mads,nonnegative=TRUE)
				if(prior_choice=='jeffery') prior_densities<-1/simulated_parameters[,2]
				if(prior_choice=='Cauchy') prior_densities<-dt((log(simulated_parameters[,2])-log(prior_location))/prior_scale,df=prior_df)/simulated_parameters[,2]
			}
			if(parameter_setting==4){
				param_est<-1:2; twoD_covering<-TRUE
				subgroup_mean<-unlist(lapply(subgroup,mean)) # a N*1 subgroup_mean
				subgroup_mads<-unlist(lapply(subgroup,mad))/1.4826 # a N*1 subgroup standard deviation
				simulated_parameters<-cbind(simulate_from_r2(N,subgroup_mean),simulate_from_r2(N,subgroup_mads,nonnegative=TRUE))
				proposal_densities<-density_r2(simulated_parameters[,1],subgroup_mean)*density_r2(simulated_parameters[,2],subgroup_mads,nonnegative=TRUE)
				if(prior_choice!='jeffery') prior_densities<-dt((simulated_parameters[,1]-prior_location[1])/prior_scale[1],df=prior_df)*dt((log(simulated_parameters[,2])-log(prior_location[2]))/prior_scale[2],df=prior_df)/simulated_parameters[,2]
				if(prior_choice=='jeffery') prior_densities<-1/simulated_parameters[,2]
			}

			# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
			tested_value<-tested_parameters[data_i,]
			divide_counts<-N*n/(5*10^7)
			if(divide) simulated_summaries<-simpleDR_in_paral(counts=divide_counts,x=simulated_parameters,FUN=simulate_pseudo_summaries,in_type='matrix',out_type='stacked_matrices',paral=divide_paral,platform=platform,packages=c('matrixStats'),n=n,setting=parameter_setting)
			if(!divide){
				simulated_observations<-simulate_pseudo_datasets(simulated_parameters,n) # a N*n matrix
				simulated_summaries<-cal_summaries(simulated_observations,setting=parameter_setting,platform=platform) # a N*d matrix		
			}
			reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries)
			weights_of_sample<-prior_densities/proposal_densities
			for(p_i in 1:length(p_all)){
				points(p_i,parameter_setting,pch=20)
				tolerance_percentile<-p_all[p_i]

				#################################### 
				# importance sampling ABC (ISABC)
				####################################
				method_i<-1
				results_ISABC<-rejABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
				results_ISABC<-results_ISABC[[1]] # This only deals with one observed dataset
				importance_weights<-weights_of_sample[results_ISABC$acceptance_ind] 
				results_summary<-ABC_ACC_results_summary(results_ISABC,type='ABC',w=importance_weights,true_params=tested_value,param_est=param_est,CI_alpha=CI_alpha,post_sample=posterior_sample_i)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol

				###################################################
				# importance sampling regression ABC (ISregABC)
				###################################################
				method_i<-2
				results_ISregABC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample)
				results_ISregABC<-results_ISregABC[[1]] # This only deals with one observed dataset
				importance_weights<-results_ISregABC$weights
				results_summary<-ABC_ACC_results_summary(results_ISregABC,type='ABC',w=importance_weights,true_params=tested_value,param_est=param_est,CI_alpha=CI_alpha,post_sample=posterior_sample_i,method_i=method_i)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol

				###################################################
				# rejection ACC (rejACC)
				###################################################
				method_i<-3
				results_rejACC<-results_ISABC
				results_summary<-ABC_ACC_results_summary(results_rejACC,type='ACC',w=NULL,true_params=tested_value,,param_est=param_est,CI_alpha=CI_alpha)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol

				###################################################
				# regression ACC (regACC)
				###################################################
				method_i<-4
				results_regACC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
				results_regACC<-results_regACC[[1]]
				results_summary<-ABC_ACC_results_summary(results_regACC,type='ACC',w=NULL,true_params=tested_value,param_est=param_est,CI_alpha=CI_alpha,method_i=method_i)

				### This part is identical for all method_i ###
				IS_ESSr_all[[method_i]][p_i,data_i]<-results_summary$IS_ESS
				post_means_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$mean_vec
				post_variances_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$var_vec
				post_quant_all[[method_i]][[n_i]][[p_i]][data_i,]<-results_summary$quant_vec
				coverage_all[[method_i]][p_i,data_i]<-results_summary$coverage
				wid_vol_all[[method_i]][p_i,data_i]<-results_summary$wid_vol
			}
# browser()
# IS_ESSr_all[[2]][,data_i]
# c
		##########################
		#### Recrod running time ###
		duration<-signif(proc.time()-start_time,2)
		##########################	
		}
	}
	return(list(post_means_all=post_means_all,post_variances_all=post_variances_all,post_quant_all=post_quant_all,IS_ESSr_all=IS_ESSr_all,coverage_all=coverage_all,wid_vol_all=wid_vol_all))
}
