library(synlik)
library(matrixStats) # For function 'colVars'


##########################
# Functions	
##########################
point_estimate<-function(dataset,initPar,niter,nburn,priorFun,propCov,nsim = 500, multicore=FALSE,ncores_sim=detectCores() - 1,control=list()){
### What is it for? ### 
# This function calculates a rough point estimate of the model parameter for a given dataset. This is used to construct the KDE-type initial distribution for ABC/ACC. The parameter (logR,logSigma,logPhi) is estiamted. The prior distribution is transformed from the improper uniform distribution of (R,Sigma,Phi) over all positive values, hence the log density of the prior is logR+logSigma+logPhi. 

### In and Out ###
# Input: dataset - The dataset for which a point estimate is calculated. other options - see 'smcmc' function in 'synlik' package. multicore - this is for smcmc in 'synlik', don't use it.
# Output: an p*1 vector, where p is the parameter dimension

	ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),extraArgs = list("nObs" = 50, "nBurn" = 50)) # The extraArgs option can not use any expressions, so we change the value of 'nObs' after defining the synlik object
	ricker_sl@extraArgs$nObs=length(dataset)
	ricker_sl@data<-dataset
	ricker_sl@extraArgs$obsData<-dataset
	ricker_sl_postsample<-smcmc(ricker_sl, initPar =initPar,niter = niter,burn = nburn,priorFun =priorFun, propCov = propCov,targetRate=0.234, nsim = nsim, multicore=multicore,ncores=ncores_sim,control=control) 	

	# ### Debugging ### 
	# ricker_sl_postsample@accRate
	# par(mfrow=c(2,2))
	# plot(1:niter,ricker_sl_postsample@chains[,1],type='l',ylim=c(2,6))
	# abline(h=parameter_true[1],col=2,lty=2)
	# abline(h=mean(ricker_sl_postsample@chains[,1]),col=3,lty=2)
	# plot(1:niter,ricker_sl_postsample@chains[,2],type='l',ylim=c(-8,0))
	# abline(h=parameter_true[2],col=2,lty=2)
	# abline(h=mean(ricker_sl_postsample@chains[,2]),col=3,lty=2)
	# plot(1:niter,ricker_sl_postsample@chains[,3],type='l',ylim=c(1.5,3.5))
	# abline(h=parameter_true[3],col=2,lty=2)
	# abline(h=mean(ricker_sl_postsample@chains[,3]),col=3,lty=2)	

	return(colMeans(ricker_sl_postsample@chains))
}


ABC_PMC_Ricker<-function(group_i,subgroup,subgroup_indices,subgroup_param_est,ricker_sl,nobs,param_est_ind=1:3,N_PMC=10^3,bandwidth_percentile=0.1,iter_max=10,summary_normalising=TRUE){
### What is it for? ### 
# This function simulates from the posterior conditional on the summary statistics of a subset of the observed dataset, using ABC-PMC. The algorithm initiates from the KDE-type distribution formed by the rough point estimates from subsets of data. 

		subgroup_observations<-subgroup[[group_i]]
		subgroup_index<-subgroup_indices[[group_i]]
		ricker_sl@extraArgs$obsData<-subgroup_observations
		subgroup_summary<-c(ricker_sl@summaries(subgroup_observations,ricker_sl@extraArgs))			
		PMC_simulator<-function(parameter_values) simulate_pseudo_summaries(parameter_values,n=nobs,index=subgroup_index,d=length(subgroup_summary),obsData=subgroup_observations)
		prior_density<-function(parameter_values) rowProds(exp(parameter_values))
		prior_check<-function(parameter_values) rep(TRUE,nrow(parameter_values))
		proposal_func<-function(N,means,SDs,weights_normalised)	nonparamproposal(N,means,SDs,weights_normalised,check_prior=prior_check)
		PMC_ini_proposal<-function(N){
			simulated_parameters<-NULL
			proposal_densities<-1
			for(param_i in param_est_ind){
				simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,subgroup_param_est[,param_i]))
				proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],subgroup_param_est[,param_i])
			}				
			return(list(sample=simulated_parameters,densities=proposal_densities))
		}
		PMC_results<-ABC_PMC(observed_summary=subgroup_summary,simulator=PMC_simulator,initial_proposal_func=PMC_ini_proposal,prior_density=prior_density,proposal_func=proposal_func,N=N_PMC,bandwidth_percentile=bandwidth_percentile,iter_max=iter_max,summary_normalising=summary_normalising)
		component_sd<-sqrtmatrix(wvar(PMC_results$theta,PMC_results$weights))
		component_mean<-wmean(PMC_results$theta,PMC_results$weights)
		return(list(mean=component_mean,sd=component_sd))
}

simulate_pseudo_datasets<-function(parameter_values,n){
### What is it for? ### 
# Simulate summary statistics of pseudo datasets given a set of parameter values (logR,logSigma,logPhi). Each dataset contains length-n series of populations from Ricker model.

### Specification of Implementation ###

### In and Out ###
# Input: parameter_values - a N*3 matrix with columns (logR,logSigma,logPhi). n - data size in each dataset. d - dimension of summary statistic
# Output: an N*d matrix with rows containing each summary
	
	N<-nrow(parameter_values) # the Monte Carlo size 
	pseudo_dataset<-matrix(0,N,n)
	for(set_i in 1:N){
		ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = parameter_true,extraArgs = list("nObs" = 50, "nBurn" = 50))
		ricker_sl@param<-parameter_values[set_i,]
		ricker_sl@extraArgs$nObs<-n
		pseudo_dataset[set_i,]<-simulate(ricker_sl, nsim = 1)
	}
	return(pseudo_dataset)
}


simulate_pseudo_summaries<-function(parameter_values,n,index=NULL,d,obsData){
### What is it for? ### 
# Simulate summary statistics of pseudo datasets given a set of parameter values (logR,logSigma,logPhi). Each dataset contains length-n series of populations from Ricker model.

### Specification of Implementation ###

### In and Out ###
# Input: parameter_values - a N*3 matrix with columns (logR,logSigma,logPhi). n - data size in each dataset. d - dimension of summary statistic
# Output: an N*d matrix with rows containing each summary
	
	if(is.null(index)) index<-1:n
	N<-nrow(parameter_values) # the Monte Carlo size 
	pseudo_summaries<-matrix(0,N,d)
	for(set_i in 1:N){
		ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = parameter_true,extraArgs = list("nObs" = 50, "nBurn" = 50))
		ricker_sl@param<-parameter_values[set_i,]
		ricker_sl@extraArgs$nObs<-n
		pseudo_dataset<-simulate(ricker_sl, nsim = 1)
		pseudo_dataset<-pseudo_dataset[index]
		ricker_sl@extraArgs$obsData<-obsData
		if(sum(is.nan(pseudo_dataset))>0) pseudo_summaries[set_i,]<-rep(0,d) 
		if(sum(is.nan(pseudo_dataset))==0) pseudo_summaries[set_i,]<-ricker_sl@summaries(pseudo_dataset,ricker_sl@extraArgs)
	}
	return(pseudo_summaries)
}

simulate_from_r<-function(N,subgroup_est,nonnegative=FALSE){
### What is it for? ###
# Simualte N values for an individual parameter from the data-dependent distribution which is a mixture normal of the kernel density estimate(KDE) formed by point estimates of observation subsets. If the parameter is nonnegative, then its logarithmic value is simulated from the mixture normal.

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

density_r<-function(x,subgroup_est,nonnegative=FALSE){
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

density_prior<-function(parameter){
### What is it for? ###
# Calculate the prior density of the parameter (R,Sigma,Phi). R, Sigma and Phi follow improper uniform over all positive values all parameters are independent.

### In and Out ###
# Input: parameter - a N*3 matrix having parameter values 
# Output: a N*1 matrix having prior densities.
	R<-parameter[,1]; Sigma<-parameter[,2]; Phi<-parameter[,3]
	prior_dens<-(R>0)*(Sigma>0)*(Phi>0)
	return(prior_dens)
}

onesideCI_cal<-function(ACC_values,alpha,weights=NULL){
### What is it for? ###
# Calculate the 1-dimension one sided confidence interval using ACC sample. The interval (-\infty,x) covers the true parameter value with probbility alpha.

### Specification of Implementation ###
# 1. The current version only calculates CI for one-dimension ACC_values

### In and Out ###
# Input: ACC_values - a N*p matrix. alpha - a scalar between 0 and 1. weights - a length N vector if the sample is weighted.
# Output: a scalar indicating the upper bound of the CI

	if(is.null(weights[1])) return(2*mean(ACC_values)-quantile(ACC_values,probs=1-alpha))
	if(!is.null(weights[1])) return(2*wmean(ACC_values,weights)-wquant(ACC_values,q=1-alpha,weights=weights))
}

ABC_ACC_results_summary<-function(results,type='ABC',w=NULL,true_params=NULL,oneD_covering=FALSE,twoD_covering=FALSE,twoD_pairs=NULL,CI_alpha=NULL){
### What is it for? ###
# Summaries the results of ABC/ACC algorithms by calculating the mean estimates, variance estimates, ESS of importance weights, two-sided credible/confidence bounds, one-sided credible/confidence bound, success/failure of covering true parameters by each interval, and the covering by two-dimension credible/confidence. 

### Specification of Implementation ###
# 1. For one-dimension covering, confidence bounds are calculated for each parameter.
# 2. For two-dimension covering, the pairs of parameters for which confidence regions are calculated need to be specified.

### In and Out ###
# Input: results - a list containing five elements: "parameters", "summaries","acceptance_ind","tolerance" and "acceptance_rate", which usually is the output of 'rejABC'/'regABC' functions. type - either 'ABC' or 'ACC'. weighted/w - whether or not the ABC output are weighted, and the vector of importance weights. true_params - the vector of true parameter, used to calculate the success/failure of being coverred. oneD_covering/twoD_covering - whether or not the one-dimension/two-dimension credible/confidence regions should be calculated. twoD_pairs - a K*2 matrix specifying the confidence region of which pairs are calculated
# Output: a list containing the above mentioned results.

	IS_ESS<-0; twoside_CI<-0; oneside_upperCI<-0; oneside_lowerCI<-0; twoD_CR<-0; 
	wald_upper<-0; wald_lower<-0; wald_coverage<-0
	if(!is.null(w[1])){
		mean_vec<-wmean(results$parameters,w)
		var_vec<-wvar(results$parameters,w,diagonal=T)
		IS_ESS<-ESS(w)/length(w)
	}
	if(is.null(w[1])){
		mean_vec<-colMeans(results$parameters)
		var_vec<-colVars(results$parameters)
	}
	if(type=='ABC'){
		oneD_bound_cal<-function(params,q,weights){
			if(!is.null(weights[1])) return(wquant(params,q=q,weights=weights))
			if(is.null(weights[1])) return(quantile(params,q=q))
		}
		param_twoD_coverred<-true_params
	} 
	if(type=='ACC'){
		oneD_bound_cal<-function(params,q,weights) return(onesideCI_cal(params,alpha=q,weights))
		param_twoD_coverred<-2*mean_vec-true_params		
	}
	twoside_upper<-0; twoside_lower<-0; twoside_coverage<-0; oneside_upper<-0; oneside_uppercoverage<-0; oneside_lower<-0; oneside_lowercoverage<-0
	if(oneD_covering){
		for(param_i in 1:length(true_params)){
			twoside_bounds<-c(oneD_bound_cal(results$parameters[,param_i],q=1-(1-CI_alpha)/2,weights=w),oneD_bound_cal(results$parameters[,param_i],q=(1-CI_alpha)/2,weights=w))
			twoside_upper[param_i]<-max(twoside_bounds)
			twoside_lower[param_i]<-min(twoside_bounds)
			twoside_coverage[param_i]<-(true_params[param_i]<=twoside_upper[param_i])&(true_params[param_i]>=twoside_lower[param_i])
			oneside_upper[param_i]<-oneD_bound_cal(results$parameters[,param_i],q=CI_alpha,weights=w)
			oneside_uppercoverage[param_i]<-(true_params[param_i]<=oneside_upper[param_i])
			oneside_lower[param_i]<-oneD_bound_cal(results$parameters[,param_i],q=1-CI_alpha,weights=w)
			oneside_lowercoverage[param_i]<-(true_params[param_i]>=oneside_lower[param_i])
			wald_bounds<-c(mean_vec+qnorm(CI_alpha/2)*sqrt(var_vec),mean_vec-qnorm(CI_alpha/2)*sqrt(var_vec))
			wald_upper[param_i]<-max(wald_bounds)
			wald_lower[param_i]<-min(wald_bounds)
			wald_coverage[param_i]<-(true_params[param_i]<=wald_upper[param_i])&(true_params[param_i]>=wald_lower[param_i])
		}
		twoside_CI<-list(twoside_lower=twoside_lower,twoside_upper=twoside_upper,twoside_coverage=as.logical(twoside_coverage))
		oneside_upperCI<-list(oneside_upper=oneside_upper,oneside_coverage=as.logical(oneside_uppercoverage))
		oneside_lowerCI<-list(oneside_lower=oneside_lower,oneside_coverage=as.logical(oneside_lowercoverage))
		wald_CI<-list(wald_lower=wald_lower,wald_upper=wald_upper,wald_coverage=as.logical(wald_coverage))
	}
	if(twoD_covering){
		if(!is.null(w[1])) w_ks<-w/sum(w)*length(w)
		if(is.null(w[1])) w_ks<-NULL
		K<-nrow(twoD_pairs); twoD_CR<-list(0)
		for(twoD_i in 1:K){
			twoD_pair<-twoD_pairs[twoD_i,]
			twoD_CR[[twoD_i]]<-HDR_2d_coverred(param_twoD_coverred[twoD_pair],x=results$parameters[,twoD_pair[1]],y=results$parameters[,twoD_pair[2]],alpha=CI_alpha,weights=w_ks)
		}	
	}
	return(list(mean_vec=mean_vec,var_vec=var_vec,IS_ESS=IS_ESS,twoside_CI=twoside_CI,wald_CI=wald_CI,oneside_upperCI=oneside_upperCI,oneside_lowerCI=oneside_lowerCI,twoD_CR=twoD_CR))
}

Process_all_output<-function(indices,output_all,param_est,number_coverage_type){
### What is it for? ###
# Summarise each part of output_all using some criterion function and list in a table. 

### Specification of Implementation ###
# 1. This is not a generic function.

### In and Out ###
# Input: indices - indices of runnings that are used to calculate the criterion. This is for the purpose of 'Monte_Carlo_stability_check'. output_all - a list of experiment variables and various outputs. 
# Output: a list

	p_all<-output_all$p
	post_means_all<-output_all$post_means 
	post_variances_all<-output_all$post_variances 
	post_coverage_all<-output_all$post_coverage
	CIbounds_all<-output_all$CIbounds
	results_all<-output_all$results
	for(method_i in 1:number_methods){
		for(p_i in 1:length(p_all)){
			param_count<-1
			for(param_i in param_est){
				post_means<-post_means_all[[method_i]][[p_i]][indices,param_i]
				post_std<-sqrt(post_variances_all[[method_i]][[p_i]][indices,param_i])
				results_all[[param_count]][[1]][method_i,p_i]<-median(post_means)
				results_all[[param_count]][[2]][method_i,p_i]<-mad(post_means)
				results_all[[param_count]][[3]][method_i,p_i]<-median(post_std)
				results_all[[param_count]][[4]][method_i,p_i]<-mad(post_std)
				CIupper<-CIbounds_all$twosideupper[[param_i]][[method_i]][p_i,indices]
				CIlower<-CIbounds_all$twosidelower[[param_i]][[method_i]][p_i,indices]
				results_all[[length(param_est)+number_coverage_type+nrow(twoD_pairs)+param_count]][method_i,p_i]<-median(CIupper-CIlower)
				param_count<-param_count+1				
			}
			param_count<-param_count-1
			for(coverge_i in 1:number_coverage_type){
				post_coverage<-post_coverage_all[[coverge_i]][[method_i]][p_i,indices]
				results_all[[param_count+coverge_i]][method_i,p_i]<-mean(post_coverage)
			}
		}
	}
	return(results_all)
}


####################################################
# Simulation of reference tables and running 
# ABC/ACC algorithms in parallel
####################################################

Ricker_ABC_ACC_algorithms<-function(data_i,N,tested_observations_all,param_ests,p_all,CI_alpha,delete_summaries=-(2:3)){
### What is it for? ### 
# This function simulates sample of parameter values from the KDE-type proposal constructing by points in param_ests[[data_i]], and simulate summary statistic of pseudo dataset for each simulated parameter value using 'ricker_sl' object. Then the methods in comparison are run. 

### Specification of Implementation ###
# 1. 'ricker_sl' object is required as a global variable.

### In and Out ###
# Input: data_i - Used in param_ests[[data_i]]. N - MC sample size. tested_observations_all - A N*nobs matrix. param_ests - a list with test_size elements, and each element is a K*p matrix where p is the parameter dimension. p_all - the vector contains acceptance proportions to be tested. CI_alpha - the nominal level for confidence intervals. 
# Output: an p*1 vector.

	param_est_ind<-1:3
	post_means<-list(0); post_variances<-list(0);
	tmp1<-rep(0,length(param_est_ind))
	for(method_i in 1:6){
		post_means[[method_i]]<-list(0)
		if(method_i==1) names(post_means)[method_i]<-'ISABC'
		if(method_i==2) names(post_means)[method_i]<-'ISregABC'
		if(method_i==3) names(post_means)[method_i]<-'ISpi-ACC'
		if(method_i==4) names(post_means)[method_i]<-'ISpi-regACC'
		if(method_i==5) names(post_means)[method_i]<-'rejACC'
		if(method_i==6) names(post_means)[method_i]<-'regACC'
		for(p_i in 1:length(p_all)){
			post_means[[method_i]][[p_i]]<-tmp1
			names(post_means[[method_i]])[p_i]<-paste0('acceptance rate=',p_all[p_i])
		}
	}
	post_variances<-post_means
	post_coverage1D<-list(0)
	for(param_i in param_est_ind){
		post_coverage1D[[param_i]]<-post_means
		for(method_i in 1:6) post_coverage1D[[param_i]][[method_i]]<-rep(0,length(p_all))
	}
	names(post_coverage1D)<-names(parameter_true)
	wald_coverage<-post_coverage1D
	post_coverage2D<-list(0); twoD_pairs<-rbind(c(1,2),c(1,3))
	for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]]<-post_coverage1D[[1]]
	names(post_coverage2D)<-c('logR_and_logSigma','logR_and_logPhi')

	post_coverage_onesideupper<-post_coverage1D
	post_coverage_onesidelower<-post_coverage1D
	CIoneside_upper<-post_coverage1D
	CIoneside_lower<-post_coverage1D
	CItwoside_upper<-post_coverage1D
	CItwoside_lower<-post_coverage1D
	CIwald_upper<-post_coverage1D
	CIwald_lower<-post_coverage1D

	IS_ESSr<-post_coverage1D[[1]][[1]] # The effective sample size ratio, i.e. ESS/N0, of the importance sampling ABC
	param_coverage_type<-c('2D','1D') # Indicate what type of cofidence region should be calculated, '1D' for individual parameter, '2D' or '3D' for two or three parameters jointly.
	number_coverage_type<-6; type_name<-c(names(post_coverage2D),names(post_coverage1D))
	if('1D'%in%param_coverage_type) oneD_covering<-TRUE else oneD_covering<-FALSE
	if('2D'%in%param_coverage_type) twoD_covering<-TRUE else twoD_covering<-FALSE

	tested_observations<-tested_observations_all[data_i,]
	ricker_sl@extraArgs$obsData<-tested_observations
	tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))

	# Simulate the parameter values  and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
	simulated_parameters<-NULL
	proposal_densities<-1
	prior_densities<-1
	for(param_i in param_est_ind){
		simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,param_ests[[data_i]][,param_i]))
		proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],param_ests[[data_i]][,param_i])
		prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
	}
	# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
	simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
	tested_summary_less<-tested_summary[delete_summaries]
	reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries[,delete_summaries])
	weights_of_sample<-prior_densities/proposal_densities

	tested_value<-parameter_true
	for(p_i in 1:length(p_all)){
		tolerance_percentile<-p_all[p_i]
		#################################### 
		# importance sampling ABC (ISABC)
		####################################
		method_i<-1
		results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE)
		results_ISABC<-results_ISABC[[1]] # This only deals with one observed dataset
		importance_weights<-weights_of_sample[results_ISABC$acceptance_ind] 
		results_summary<-ABC_ACC_results_summary(results_ISABC,type='ABC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		IS_ESSr[p_i]<-results_summary$IS_ESS
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# importance sampling regression ABC (ISregABC)
		###################################################
		method_i<-2
		results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample)
		results_ISregABC<-results_ISregABC[[1]] # This only deals with one observed dataset
		importance_weights<-results_ISregABC$weights
		results_summary<-ABC_ACC_results_summary(results_ISregABC,type='ABC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]

# browser()

		###################################################
		# importance sampling prior-ACC (ISpi-ACC)
		###################################################
		method_i<-3
		results_ISpiACC<-results_ISABC
		importance_weights<-weights_of_sample[results_ISABC$acceptance_ind]
		results_summary<-ABC_ACC_results_summary(results_ISpiACC,type='ACC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# importance sampling prior-regression ACC (ISpi_regACC)
		###################################################
		method_i<-4
		results_ISpi_regACC<-results_ISregABC
		importance_weights<-results_ISregABC$weights
		results_summary<-ABC_ACC_results_summary(results_ISpi_regACC,type='ACC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# rejection ACC (rejACC)
		###################################################
		method_i<-5
		results_rejACC<-results_ISABC
		results_summary<-ABC_ACC_results_summary(results_rejACC,type='ACC',w=NULL,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# regression ACC (regACC)
		###################################################
		method_i<-6
		results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		results_regACC<-results_regACC[[1]]
		results_summary<-ABC_ACC_results_summary(results_regACC,type='ACC',w=NULL,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means[[method_i]][[p_i]]<-results_summary$mean_vec
		post_variances[[method_i]][[p_i]]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower[[param_i]][[method_i]][p_i]<-twoside_CI$twoside_lower[param_i]

			wald_CI<-results_summary$wald_CI
			wald_coverage[[param_i]][[method_i]][p_i]<-wald_CI$wald_coverage[param_i]
			CIwald_upper[[param_i]][[method_i]][p_i]<-wald_CI$wald_upper[param_i]
			CIwald_lower[[param_i]][[method_i]][p_i]<-wald_CI$wald_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper[[param_i]][[method_i]][p_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower[[param_i]][[method_i]][p_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D[[pair_i]][[method_i]][p_i]<-results_summary$twoD_CR[[pair_i]]
	}
	raw_results_all<-list(reference_tables=reference_tables,weights_of_sample=weights_of_sample,tested_summary_less=tested_summary_less,post_means=post_means,post_variances=post_variances,post_coverage=list(twoside=post_coverage1D,wald=wald_coverage,onesideupper=post_coverage_onesideupper,onesidelower=post_coverage_onesidelower,twoD=post_coverage2D),IS_ESSr=IS_ESSr,CIbounds=list(twosideupper=CItwoside_upper,twosidelower=CItwoside_lower,wald_upper=CIwald_upper,wald_lower=CIwald_lower,onesideupper=CIoneside_upper,onesidelower=CIoneside_lower))
	return(raw_results_all)
}

Example2_parametric_bootstrap<-function(tested_observations_all,parameter_true,N,nobs,CI_alpha=0.05,platform){
### What is it for? ###
# Run parametric bootstrap algorithm for example 1 of Thorton, Li and Xie (2018) under different parameter/summary settings. This version needs a quicker point estimate since the synthetic likelihood estimate uses MCMC and is too slow to calculate. Some alternatives are in https://github.com/umbertopicchini/pomp-ricker.

	p<-length(parameter_true)
	test_size<-nrow(tested_observations_all)
	post_coverage_all<-rep(0,test_size)
	CIupper_all<-post_coverage_all
	param_i<-1:3
	plot(0,1,xlim=c(1,test_size),ylim=c(0,2),xlab='data_i',main=paste0('Ricker logR and logSigma'))

	##########################
	# Algorithms
	##########################

	tested_parameters<-matrix(rep(parameter_true,test_size),ncol=p,byrow=T) # a test_size*3 matrix
	duration<-0
	for(data_i in 1:test_size){
		points(data_i,1,pch=20)
		##########################
		#### Recrod start time ###
		start_time<-proc.time()
		##########################
		tested_observations<-tested_observations_all[data_i,]
		system.time(obs_pt_estimate<-point_estimate(tested_observations,initPar = c(3.2, -1, 2.6),niter = niter,nburn = nburn,priorFun = function(input, ...) sum(input), propCov = diag(c(0.2, 0.35, 0.2))^2,control=list(verbose=TRUE)))
		tested_value<-tested_parameters[data_i,]
		pt_estimate_all<-matrix(0,N,p)
		for(N_i in 1:N){
			simulated_observations<-c(simulate_pseudo_datasets(t(obs_pt_estimate),n))
			pt_estimate_all[N_i,]<-point_estimate(simulated_observations,initPar = c(3.2, -1, 2.6),niter = niter,nburn = nburn,priorFun = function(input, ...) sum(input), propCov = diag(c(0.2, 0.35, 0.2))^2)
		}
		##########################
		#### Recrod running time ###
		duration<-signif(proc.time()-start_time,2)
		##########################			
		if(length(param_i)==1){
			upper<-quantile(pt_estimate_all[,param_i],1-CI_alpha)
			post_coverage_all[[n_i]][data_i]<-(tested_value<=upper)
			CIupper_all[[n_i]][data_i]<-upper					
		}
		if(length(param_i)>1){
			post_coverage_all[[n_i]][data_i]<-HDR_2d_coverred(tested_value,x=pt_estimate_all[,param_i[1]],y=pt_estimate_all[,param_i[2]],alpha=CI_alpha)
		}
	}
	return(list(post_coverage_all=post_coverage_all,CIupper_all=CIupper_all))
}

