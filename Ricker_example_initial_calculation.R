## Create for the Ricker example in Thorton, Li and Xie (2018) 
## This file is for simulating observed datasets and calculating point estimates for constructing initial ballpark distribution. 
## In this example, the point estimate is the synthetic likelihood estimate of Wood (2010). The calculation is very costly, since each estimate needs to generate an MCMC sample. Therefore the initial ballpark distribtuion is calculated beforehand and shared by all experiments.

library(matrixStats) # For function 'rowMads'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation

# setwd('/Users/wentao/bitbucket/')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\core-functions')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')

# setwd('/Users/wentao/bitbucket/ACC examples')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\acc-examples')
source('Ricker_example_functions.R')
# fix(smcmc) # The following lines in 'smcmc' need to be changed.
# Change 'tmp <- colMeans(mcmcSample[1:ii, ], na.rm = TRUE)' to 'tmp <- colMeans(mcmcSample[1:(ii-burn), ], na.rm = TRUE)''


### Setup for parallelisation ###
# platform<-'Mac'; ncores<-3
platform<-'Win'; ncores<-6

### Experiment variables ###
test_size<-150 # Number of runs in an experiment

# nobs<-500; group_no<-floor(sqrt(nobs)); rn_sub<-'half'; 
nobs<-500; group_no<-floor(nobs^(3/5)); rn_sub<-'threefifth'
# nobs<-50; group_size<-10; rn_sub<-'sequential_groupsize10'
if(nobs>100){ # When nobs is small, the number of non-overlapping subgroups is too small. Following codes are used to obtain overlapping subgroups.
	group_size<-floor(nobs/group_no); half_group_size<-floor(nobs/group_no/2)
	subgroup_indices<-divide(1:nobs,group_no); subgroup<-list(0)
}
if(nobs<=100){
	subgroup_starting_indices<-seq(1,(nobs-group_size+1),by=2)
	group_no<-length(subgroup_starting_indices)
	subgroup_indices<-list(0); subgroup<-list(0)
	for(group_i in 1:group_no) subgroup_indices[[group_i]]<-subgroup_starting_indices[group_i]:(subgroup_starting_indices[group_i]+group_size-1)
}

parameter_true<-c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)) # Follows the setting in Wood(2010)

### Random seed ###
set.seed(200)

### Observed datasets simulation ###
param_est_ind<-1:3
ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = parameter_true,extraArgs = list("nObs" = nobs, "nBurn" = 50))
tested_observations_all<-simulate(ricker_sl, nsim = test_size)
duration<-0
niter<-10000; nburn<-1000


######################################################
# Create r_n using rough point estimates
######################################################
subgroup_param_est<-list(0)
dev.new()
for(data_i in 90:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	tested_observations<-tested_observations_all[data_i,]
	ricker_sl@extraArgs$obsData<-tested_observations
	tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))

	# Divide the dataset into subsets for the purpose of constructing initial distribution
	for(group_i in 1:group_no) subgroup[[group_i]]<-tested_observations[subgroup_indices[[group_i]]]
	plot(0,data_i,xlim=c(1,group_no),ylim=c(data_i,data_i),xlab='group_i',ylab='data_i',main=paste0('r_n construction: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	# Calculate the point estimates for sub-datasets in order to construct KDE-type r_n 
	# With parallel_lapply in 6 cores and niter=1e4, running time for 22 subgroups is 520s
	# 
	tmp<-parallel_lapply(subgroup,point_estimate,ncores=ncores,platform=platform,packages=c('synlik','parallel'),initPar = c(3.2, -1, 2.6),niter = niter,nburn = nburn,priorFun = function(input, ...) sum(input), propCov = diag(c(0.2, 0.35, 0.2))^2, nsim = 500, multicore=FALSE,progress_bar=TRUE)
	subgroup_param_est[[data_i]]<-do_recombing(tmp,out_type='stacked_matrices')
	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	
}
rn_method<-'PT'

save.image(file=paste0('Ricker_SubgroupEst_',test_size,'rep_n',nobs,'.RData'))


####################################################################
# sampling from rough point estimates r_n (IN PARALLEL)
####################################################################
test_size<-150; N<-10^5
set.seed(200)
tmp<-parallel_lapply(1:test_size,Ricker_reference_table_sampling,ncores=ncores,platform=platform,packages=c('synlik','parallel'),N=N,tested_observations_all=tested_observations_all,param_ests=subgroup_param_est,progress_bar=TRUE)
tmp<-do_recombing(tmp,out_type='fixed_list_combined_list')
reference_tables_all<-tmp$reference_tables
weights_of_sample_all<-tmp$weights_of_sample
tested_summary_less_all<-tmp$tested_summary_less
rm(tmp); gc()
save.image('Ricker_initial_cal_150rep_N1e5_n50_rn_PT.RData')


######################################################
# r_n refined by ABC-PMC
######################################################
# rm(reference_tables_all); rm(weights_of_sample_all); rm(tested_summary_less_all); gc()
PMC_refine<-FALSE; N_PMC<-10^4; bandwidth_percentile<-0.05
subgroup_refine_est<-list(0)
tmp<-list(0)
options(error=recover)
set.seed(200)
dev.new()
for(data_i in 1:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################
	plot(0,data_i,xlim=c(1,group_no),ylim=c(data_i,data_i),xlab='group_i',ylab='data_i',main=paste0('r_n construction: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	tested_observations<-tested_observations_all[data_i,]
	for(group_i in 1:group_no) subgroup[[group_i]]<-tested_observations[subgroup_indices[[group_i]]]
	# for(group_i in 1:group_no) tmp<-ABC_PMC_Ricker(group_i,subgroup=subgroup,subgroup_indices=subgroup_indices,subgroup_param_est=subgroup_param_est[[data_i]],ricker_sl=ricker_sl,nobs=nobs,N_PMC=10^3,bandwidth_percentile=0.1,iter_max=10)
	tmp<-parallel_lapply(as.list(1:length(subgroup)),ABC_PMC_Ricker,ncores=ncores,platform=platform,packages=c('synlik','parallel','matrixStats'),subgroup=subgroup,subgroup_indices=subgroup_indices,subgroup_param_est=subgroup_param_est[[data_i]],ricker_sl=ricker_sl,nobs=nobs,N_PMC=10^3,bandwidth_percentile=0.1,iter_max=10,progress_bar=TRUE)
	subgroup_refine_est[[data_i]]<-do_recombing(tmp,out_type='fixed_list_combined_list')
	subgroup_refine_est[[data_i]]$mean<-do.call(rbind,subgroup_refine_est[[data_i]]$mean)
	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	
}
rn_method<-'refinedPT'

subgroup_means<-list(0)
for(data_i in 1:test_size){
	subgroup_ind<-sample(1:subgroup_length,floor(subgroup_length/2))
	subgroup_means[[data_i]]<-rbind(subgroup_param_est[[data_i]][subgroup_ind,],(subgroup_refine_est[[data_i]]$mean)[-subgroup_ind,])
	subgroup_means[[data_i]][,2]<-(subgroup_refine_est[[data_i]]$mean)[,2]
}

save.image(file=paste0('Ricker_SubgroupRefinedEst_',test_size,'rep_n',nobs,'.RData'))


######################################################
# Arguments on why we should use the refined estimates of logSigma and rough estimates of the other two
######################################################
load('Ricker_SubgroupRefinedEst_150rep_n50.RData')
# Compare the rough and refined point estimates in marginal densities. The conclusion is the refined one is better in logSigma, but worse in logR and logPhi; possible reason is the PMC only uses 10 iterations and hasn't converged. 
par_names<-names(parameter_true)
subgroup_param_est_all<-list(NULL,NULL,NULL)
subgroup_refine_est_all<-list(NULL,NULL,NULL)
subgroup_comb_est_all<-list(NULL,NULL,NULL)
for(param_i in 1:3){
	for(data_i in 1:1){	
		subgroup_param_est_all[[param_i]]<-c(subgroup_param_est_all[[param_i]],subgroup_param_est[[data_i]][,param_i])
		subgroup_refine_est_all[[param_i]]<-c(subgroup_refine_est_all[[param_i]],(subgroup_refine_est[[data_i]]$mean)[,param_i])
		subgroup_comb_est_all[[param_i]]<-c(subgroup_comb_est_all[[param_i]],(subgroup_refine_est[[data_i]]$mean)[,param_i],subgroup_param_est[[data_i]][,param_i])
	}
}
dev.new(); par(mfrow=c(2,2))
for(param_i in 1:3){
	draw_densities(list(rough=density(subgroup_param_est_all[[param_i]]),refined=density(subgroup_refine_est_all[[param_i]]),comb=density(subgroup_comb_est_all[[param_i]])),main=par_names[param_i],with_legend=T)
	abline(v=parameter_true[param_i],lty=2)
}

# Compare the rough and refined point estimates in bivariate plots; black points are these point estimates; red points are simulations from r using them; the conclusion is although the r_n proposal just uses marginal distributions and doesn't consider correlation, it looks OK. 
plot_par_i<-rbind(c(1,2),c(2,3),c(1,3))
par_names<-names(parameter_true)
for(data_i in 1:1){
	dev.new(); par(mfrow=c(3,2))
	simulated_parameters<-NULL
	simulated_parameters_refined<-NULL
	plot_min_max<-matrix(0,3,2)
	for(param_i in 1:3){
		simulated_parameters<-cbind(simulated_parameters,simulate_from_r(1000,subgroup_param_est[[data_i]][,param_i]))
		simulated_parameters_refined<-cbind(simulated_parameters_refined,simulate_from_r(1000,(subgroup_refine_est[[data_i]]$mean)[,param_i]))
		tmp<-c(simulated_parameters[,param_i],simulated_parameters_refined[,param_i])
		plot_min_max[param_i,]<-c(min(tmp),max(tmp))
	}
	for(plot_i in 1:3){
		plot(simulated_parameters[,plot_par_i[plot_i,]],col=2,pch=20,xlim=plot_min_max[plot_par_i[plot_i,1],],ylim=plot_min_max[plot_par_i[plot_i,2],],main=paste0('data_i=',data_i),xlab=par_names[plot_par_i[plot_i,1]],ylab=par_names[plot_par_i[plot_i,2]])
		points(subgroup_param_est[[data_i]][,plot_par_i[plot_i,]],pch=19)
		abline(v=parameter_true[plot_par_i[plot_i,1]],lty=2)
		abline(h=parameter_true[plot_par_i[plot_i,2]],lty=2)


		plot(simulated_parameters_refined[,plot_par_i[plot_i,]],col=2,pch=20,xlim=plot_min_max[plot_par_i[plot_i,1],],ylim=plot_min_max[plot_par_i[plot_i,2],],main='refined',xlab=par_names[plot_par_i[plot_i,1]],ylab=par_names[plot_par_i[plot_i,2]])
		points((subgroup_refine_est[[data_i]]$mean)[,plot_par_i[plot_i,]],pch=19)
		abline(v=parameter_true[plot_par_i[plot_i,1]],lty=2)
		abline(h=parameter_true[plot_par_i[plot_i,2]],lty=2)
	}
}


# ######################################################
# # sampling from refined point estimates r_n (IN PARALLEL)
# ######################################################
# test_size<-150; N<-10^5
# subgroup_refined_mean<-list(0)
# for(data_i in 1:test_size) subgroup_refined_mean[[data_i]]<-subgroup_refine_est[[data_i]]$mean
# set.seed(200)
# no_parts<-5; run_parts<-divide(1:test_size,no_parts); tmp<-NULL
# for(i in 1:no_parts){
# 	tmp<-c(tmp,parallel_lapply(run_parts[[i]],Ricker_reference_table_sampling,ncores=ncores,platform=platform,packages=c('synlik','parallel'),N=N,tested_observations_all=tested_observations_all,param_ests=subgroup_refined_mean,progress_bar=TRUE))
# }
# tmp<-do_recombing(tmp,out_type='fixed_list_combined_list')
# rn_method<-'refinedPT'
# reference_tables_all<-tmp$reference_tables
# weights_of_sample_all<-tmp$weights_of_sample
# tested_summary_less_all<-tmp$tested_summary_less
# rm(tmp); gc()
# save.image('Ricker_initial_cal_150rep_N1e5_n50_rn_refinedPT.RData')


# ###############################################################
# # sampling from rough point estimates r_n (SEQUENTIAL)
# ###############################################################
# test_size<-150; N<-10^5
# reference_tables_all<-list(0)
# weights_of_sample_all<-list(0)
# tested_summary_less_all<-list(0)
# set.seed(200)
# dev.new()
# for(data_i in 1:test_size){
# 	##########################
# 	#### Recrod start time ###
# 	start_time<-proc.time()
# 	##########################

# 	plot(data_i,data_i,xlim=c(1,test_size),ylim=c(1,test_size),xlab='data_i',ylab='data_i',main=paste0('Simulation: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
# 	# Simulate the parameter values  and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
# 	tested_observations<-tested_observations_all[data_i,]
# 	ricker_sl@extraArgs$obsData<-tested_observations
# 	tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))

# 	simulated_parameters<-NULL
# 	proposal_densities<-1
# 	prior_densities<-1
# 	for(param_i in param_est_ind){
# 		simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,subgroup_param_est[[data_i]][,param_i]))
# 		proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],subgroup_param_est[[data_i]][,param_i])
# 		prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
# 	}
# 	# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
# 	simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
# 	tested_summary_less<-tested_summary[-(2:3)]
# 	reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries[,-(2:3)])
# 	weights_of_sample<-prior_densities/proposal_densities

# 	reference_tables_all[[data_i]]<-reference_tables
# 	weights_of_sample_all[[data_i]]<-weights_of_sample
# 	tested_summary_less_all[[data_i]]<-tested_summary_less

# 	##########################
# 	#### Recrod running time ###
# 	duration<-signif(proc.time()-start_time,2)
# 	##########################	
# }

# save.image('Ricker_initial_cal_300rep_N1e5_n500_rn_PT.RData')



# ######################################################
# # sampling from refined point estimates r_n (SEQUENTIAL)
# ######################################################
# set.seed(300)
# reference_tables_all<-list(0)
# weights_of_sample_all<-list(0)
# tested_summary_less_all<-list(0)
# dev.new()
# for(data_i in 1:test_size){
# 	##########################
# 	#### Recrod start time ###
# 	start_time<-proc.time()
# 	##########################

# 	plot(data_i,data_i,xlim=c(1,test_size),ylim=c(1,test_size),xlab='data_i',ylab='data_i',main=paste0('Simulation: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
# 	tested_observations<-tested_observations_all[data_i,]
# 	ricker_sl@extraArgs$obsData<-tested_observations
# 	tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))
# 	# Simulate the parameter values  and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
# 	simulated_parameters<-NULL
# 	proposal_densities<-1
# 	prior_densities<-1
# 	for(param_i in param_est_ind){
# 		refined_mean<-subgroup_refine_est[[data_i]]$mean
# 		simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,refined_mean[,param_i]))
# 		proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],refined_mean[,param_i])
# 		prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
# 	}
# 	# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
# 	simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
# 	tested_summary_less<-tested_summary[-(2:3)]
# 	reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries[,-(2:3)])
# 	weights_of_sample<-prior_densities/proposal_densities

# 	reference_tables_all[[data_i]]<-reference_tables
# 	weights_of_sample_all[[data_i]]<-weights_of_sample
# 	tested_summary_less_all[[data_i]]<-tested_summary_less

# 	##########################
# 	#### Recrod running time ###
# 	duration<-signif(proc.time()-start_time,2)
# 	##########################	
# }
# rn_method<-'refinedPT'
# ######################################################
# # the sample of r_n using refined point estimates
# ######################################################

# save.image('Ricker_initial_cal_167rep_N1e5_n500_half_rn_refinedPT.RData')



########################################################################
# sampling from r_n using refined mixture distribution (SEQUENTIAL)
########################################################################
set.seed(400)
reference_tables_all<-list(0)
weights_of_sample_all<-list(0)
tested_summary_less_all<-list(0)
dev.new()
for(data_i in 1:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	plot(data_i,data_i,xlim=c(1,test_size),ylim=c(1,test_size),xlab='data_i',ylab='data_i',main=paste0('Simulation: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	tested_observations<-tested_observations_all[data_i,]
	ricker_sl@extraArgs$obsData<-tested_observations
	tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))
	# Simulate the parameter values  and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
	simulated_parameters<-NULL
	proposal_densities<-1; proposal_densities_matrix<-matrix(1,group_no,N)
	prior_densities<-1
	mixsample_component_ind <- sample(1:group_no,size=N,replace=TRUE)
	for(param_i in param_est_ind){
		refined_means<-(subgroup_refine_est[[data_i]]$mean)[,param_i]
		refined_sds<-unlist(lapply(subgroup_refine_est[[data_i]]$sd,function(x)x[param_i,param_i]))
		simulated_parameters<-cbind(simulated_parameters,rnorm(N,mean=refined_means[mixsample_component_ind],sd=refined_sds[mixsample_component_ind]))
		for(group_i in 1:group_no) proposal_densities_matrix[group_i,]<-proposal_densities_matrix[group_i,]*dnorm(simulated_parameters[,param_i],mean=refined_means[group_i],sd=refined_sds[group_i])
		prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
	}
	# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
	simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
	proposal_densities<-colSums(proposal_densities_matrix)/group_no
	tested_summary_less<-tested_summary[-(2:3)]
	reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries[,-(2:3)])
	weights_of_sample<-prior_densities/proposal_densities

	reference_tables_all[[data_i]]<-reference_tables
	weights_of_sample_all[[data_i]]<-weights_of_sample
	tested_summary_less_all[[data_i]]<-tested_summary_less

	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	
}
rn_method<-'refinedmixture'

save.image('Ricker_initial_cal_167rep_N1e5_n500_half_rn_refinedmixture.RData')


##########################
# Debugging
##########################

par(mfcol=c(3,5))
for(data_i in seq(1,test_size,length=5)){
	for(param_i in 1:3){
		point_est_density<-density(subgroup_param_est[[data_i]][,param_i])
		draw_densities(list(point_est=point_est_density),main=names(parameter_true)[param_i])
		abline(v=parameter_true[param_i],col=4,lty=2)
	}	
}

par(mfcol=c(3,5))
for(data_i in c(1,25,134,250,278)){
	refine_mixture_density<-list(0)
	for(param_i in 1:3){
		mixsample_component_ind <- sample(1:group_no,size=1e4,replace=TRUE)
		mus <- (subgroup_refine_est[[data_i]]$mean)[,param_i]
		sds <- unlist(lapply(subgroup_refine_est[[data_i]]$sd,function(x)x[param_i,param_i]))
		mixsample <- rnorm(n=1e4,mean=mus[mixsample_component_ind],sd=sds[mixsample_component_ind])
		refine_mixture_density[[param_i]]<-density(mixsample)
	}
	for(param_i in 1:3){
		refine_mean_density<-density((subgroup_refine_est[[data_i]]$mean)[,param_i])
		point_est_density<-density(subgroup_param_est[[data_i]][,param_i])
		draw_densities(list(point_est=point_est_density,refine_mean=refine_mean_density,refine_mixture=refine_mixture_density[[param_i]]),main=names(parameter_true)[param_i])
		abline(v=parameter_true[param_i],col=4,lty=2)
	}	
}

# The followings are run after lines 44
data_i<-20
tested_observations<-tested_observations_all[data_i,]
group_no<-floor(nobs^(3/5))
subgroup1<-divide(tested_observations,group_no)
system.time(tmp1<-parallel_lapply(subgroup1,point_estimate,ncores=ncores,platform=platform,packages=c('synlik','parallel'),initPar = c(3.2, -1, 2.6),niter = niter,nburn = nburn,priorFun = function(input, ...) sum(input), propCov = diag(c(0.2, 0.35, 0.2))^2, nsim = 500, multicore=FALSE))
subgroup_param_est_tmp1<-do_recombing(tmp1,out_type='stacked_matrices')

group_no<-floor(nobs^(1/2))
subgroup2<-divide(tested_observations,group_no)
system.time(tmp2<-parallel_lapply(subgroup2,point_estimate,ncores=ncores,platform=platform,packages=c('synlik','parallel'),initPar = c(3.2, -1, 2.6),niter = niter,nburn = nburn,priorFun = function(input, ...) sum(input), propCov = diag(c(0.2, 0.35, 0.2))^2, nsim = 500, multicore=FALSE))
subgroup_param_est_tmp2<-do_recombing(tmp2,out_type='stacked_matrices')
par(mfrow=c(2,2))
for(param_i in 1:3){
	tmp1_dens<-density(subgroup_param_est_tmp1[,param_i])
	tmp2_dens<-density(subgroup_param_est_tmp2[,param_i])
	draw_densities(list(threefifth_group_r=tmp1_dens,sqrt_group_r=tmp2_dens),main=names(parameter_true)[param_i])
}

