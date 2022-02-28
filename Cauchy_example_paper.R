## Create for examples in Thorton, Li and Xie (2018) 
library(matrixStats) # For function 'rowMads'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation
library(pbapply) # For function 'pblapply'
library(mcmc)

# setwd('/Users/wentao/bitbucket/core-functions')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\core-functions')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')


# setwd('/Users/wentao/bitbucket/acc-examples')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\acc-examples')
source('Cauchy_example_functions.R')

##########################
# Experiment Setting
##########################
### Setup for parallelisation ###
# platform<-'Mac'
platform<-'Win'

### Experiment variables ###
# test_size<-500; N<-5*10^5; 
test_size<-1000; N<-5e4; n_all<-400
p_all<-signif2(c(0.005,0.05,0.1,0.2,0.4))
parameter_true<-c(location=10,scale=0.55); # Data model is Cauchy(location_cubicrt,scale_sqr)

set.seed(100)
tested_observations_all<-simulate_pseudo_datasets(matrix(rep(parameter_true,test_size),ncol=2,byrow=T),n_all[1]) # a test_size*n matrix
tested_parameters<-matrix(rep(parameter_true,test_size),ncol=2,byrow=T) # a test_size*2 matrix



logposterior_Cauchy<-function(params,obs){ 
### What is it for? ###
# Output the logarithm of the unnormalized posterior density; use Jeffrey's prior

### In and Out ###
# Input: params - location and scale parameters. obs - a length-n vector containing the dataset
# Output: the unnormalised posterior density

	loglik_val<-sum(dcauchy(obs,location=params[1],scale=params[2],log=TRUE))
	if(is.na(loglik_val)) loglik_val<--Inf
	prior_val<-1/params[2]
	if(prior_val>0) logprior_val<-log(prior_val)
	if(prior_val<=0) logprior_val<--Inf
	return(loglik_val+logprior_val)
}

N_MCMC<-1e5; sample_posterior<-list(0); duration<-0
MCMC_accept<-0
plot(0,0,xlim=c(1,test_size),ylim=c(1,test_size))
for(data_i in 1:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	plot(data_i,data_i,xlim=c(1,test_size),ylim=c(1,test_size),pch=20,xlab='data_i',ylab='data_i',main=paste0('Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds'))
	MCMC_results<-metrop(logposterior_Cauchy,initial=parameter_true,scale=diag(c(0.03,0.02)),nbatch=N_MCMC,obs=tested_observations_all[data_i,])
	MCMC_accept[data_i]<-MCMC_results$accept
	sample_posterior[[data_i]]<-MCMC_results$batch[-(1:1e4),]

	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	
}

save(sample_posterior,file=paste0('Cauchy_true_posterior_sample_',test_size,'rep_n400.RData'))

load(paste0('Cauchy_true_posterior_sample_',test_size,'rep_n400.RData'))



##########################
# Experiments for 4.2 
##########################
# Parameter settings include: 1 -- location unknown and scale known with median as summary, 2 -- location unknown and scale known with mean as summary, 3 -- location known and scale unknown with median absolute deviation (MAD) as summary, 4 -- both unknown with mean and MAD as summary
# Parameter proposals include: 1 -- KDE of certain point estimates, e.g. median and MAD, 2 -- Cauchy distribution with given location and scale.
# Choices of prior include: Cauchy distribution with scale parameter=10; Jeffrey prior for location and scale family, i.e. uniform over all location and uniform over all log-scale.

###########################################################################
# Cauchy example: N1e5, test_size 1000 
# parameter_setting 2:4
# For prior_set=2: prior is Jeffery
# For prior_set=1: prior is Cauchy
###########################################################################


CI_alpha<-0.95; prior_set<-1
raw_results_all<-list(list(0),list(0),list(0),list(0))
dev.new()
for(i in 2:4){
	# if(i!=2) prior_set<-2
	# if(i==2) prior_set<-1
	if(i==1) {parameter_setting<-1; prior_location<-parameter_true[1]}
	if(i==2) {parameter_setting<-2; prior_location<-parameter_true[1]}
	if(i==3) {parameter_setting<-3; prior_location<-parameter_true[2]}
	if(i==4) {parameter_setting<-4; prior_location<-parameter_true}
	names(raw_results_all)[parameter_setting]<-paste0('setting_',parameter_setting)
	prior_scale<-rep(10,length(prior_location))

	##########################
	# Set Random Seed 
	##########################
	set.seed(100)
	
	if(prior_set==1){
		tmp<-Example1_main(test_size,parameter_true,tested_observations_all=tested_observations_all,parameter_setting,posterior_sample=sample_posterior,N=N,platform=platform,n_all=n_all,p_all=p_all,CI_alpha=CI_alpha,divide=TRUE,prior_choice='Cauchy',prior_location=prior_location,prior_scale=prior_scale,prior_df=1)
		raw_results_all[[parameter_setting]]<-tmp
	}
	if(prior_set==2){
		tmp<-Example1_main(test_size,parameter_true,tested_observations_all=tested_observations_all,parameter_setting,posterior_sample=sample_posterior,N=N,platform=platform,n_all=n_all,p_all=p_all,CI_alpha=CI_alpha,divide=TRUE)
		raw_results_all[[parameter_setting]]<-tmp
	}
}

# For acceptance rate = c(0.005,0.05,0.1,0.2,0.4)
signif2(rowMeans((raw_results_all[[i]]$coverage_all)[[2]])) # coverage of IS-ABC

signif2(rowMeans((raw_results_all[[i]]$coverage_all)[[4]])) # coverage of ACC

signif2(rowMeans((raw_results_all[[i]]$wid_vol_all)[[2]])) # volume of IS-ABC

signif2(rowMeans((raw_results_all[[i]]$wid_vol_all)[[4]])) # volume of ACC


###########################################################################
# CI_alpha<-0.95
# Cauchy example: N1e5, test_size 1000 
# parameter_setting 2:4
# For setting 1,2,4: prior is uniform
# For setting 3: prior is t4
###########################################################################


##########################
# Results
##########################
# Under each setting, the raw result is a list with four elements: the posterior mean estimates, posterior variance estimates, coverage rate of the true parameter value in the test_size runs, and the tolerance values for all runs. Each element is a list with two levels. Level one is for the data sizes n_all, level two is for the acceptance proportions p_all, and the basic element is a test_size*2 matrix, for moments estimate, or a length test_size vector.

number_methods<-4; n<-n_all; group_No<-10
coverage_rejABC<-list(0); coverage_regABC<-list(0); coverage_rejACC<-list(0); coverage_regACC<-list(0); 
wid_vol_rejABC<-list(0); wid_vol_regABC<-list(0); wid_vol_rejACC<-list(0); wid_vol_regACC<-list(0); 
results_all<-list(0,0,0,0) # Corresponds to four experiment settings
tmp_matrix<-matrix(0,number_methods,6)
colnames(tmp_matrix)<-c('coverage','wid_vol','mean_avg','mean_var','var_avg','var_var')
rownames(tmp_matrix)<-c('ISABC','ISregABC','rejACC','regACC')
for(i in 2:2){
	if(i==1) {parameter_setting<-1; param_est<-1}
	if(i==2) {parameter_setting<-2; param_est<-1}
	if(i==3) {parameter_setting<-3; param_est<-2}
	if(i==4) {parameter_setting<-4; param_est<-1:2}
	if(i!=2) prior_set<-2
	if(i==2) prior_set<-1
	names(results_all)[parameter_setting]<-paste0('setting_',parameter_setting)
	coverage_rejABC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	coverage_regABC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	coverage_rejACC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	coverage_regACC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	wid_vol_rejABC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	wid_vol_regABC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	wid_vol_rejACC[[parameter_setting]]<-matrix(0,group_No,length(p_all))
	wid_vol_regACC[[parameter_setting]]<-matrix(0,group_No,length(p_all))

	tmp<-list(0)
	for(group_i in 1:group_No){
		if(length(param_est)==1){
			tmp[[group_i]]<-list(list(0))
			names(tmp[[group_i]])<-names(parameter_true)[param_est]
		}
		if(length(param_est)==2){
			tmp[[group_i]]<-list(list(0),list(0))
			names(tmp[[group_i]])<-names(parameter_true)
		}
		names(tmp)[group_i]<-paste0('group',group_i)
		for(p_i in 1:length(p_all)){ 
			tmp[[group_i]][[1]][[p_i]]<-tmp_matrix
			names(tmp[[group_i]][[1]])[p_i]<-paste0('p= ',p_all[p_i])
		}
		if(length(param_est)==2) tmp[[group_i]][[2]]<-tmp[[group_i]][[1]]
	}
	tmp_raw<-raw_results_all[[parameter_setting]]
	output_all<-list(n=n_all,p=p_all,post_means=tmp_raw$post_means_all,quant_vec=tmp_raw$post_quant_all,post_variances=tmp_raw$post_variances_all,coverage=tmp_raw$coverage_all,wid_vol=tmp_raw$wid_vol_all,results=tmp)

# methods=c('ISABC','ISregABC','rejACC','regACC')
# p_all<-signif2(c(0.005,0.05,0.1,0.2,0.4))
	p_i<-4
par(mfrow=c(2,2))
plot(tmp_raw$post_quant_all[[2]][[1]][[p_i]],tmp_raw$post_quant_all[[4]][[1]][[p_i]],pch=20,xlab='ISregABC',ylab='regACC',main='CI_right_endpoint')
abline(coef=c(0,1),lty=2)
plot(tmp_raw$wid_vol_all[[2]][p_i,],tmp_raw$wid_vol_all[[4]][p_i,],pch=20,xlab='ISregABC',ylab='regACC',main='width')
abline(coef=c(0,1),lty=2)
plot(density(tmp_raw$post_means_all[[4]][[1]][[p_i]][,1]),xlab='regACC',main='mean')
lines(density(tmp_raw$post_means_all[[2]][[1]][[p_i]][,1]),xlab='ISregABC',col=2)

	results_all[[parameter_setting]]<-Process_all_output(1:test_size,output_all,param_est=param_est)
	for(group_i in 1:group_No){
		for(p_i in 1:length(p_all)){
			coverage_rejABC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][1,1]
			coverage_regABC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][2,1]
			coverage_rejACC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][3,1]
			coverage_regACC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][4,1]
			wid_vol_rejABC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][1,2]		
			wid_vol_regABC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][2,2]		
			wid_vol_rejACC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][3,2]		
			wid_vol_regACC[[parameter_setting]][group_i,p_i]<-results_all[[parameter_setting]][[group_i]][[1]][[p_i]][4,2]		
		}
	}
}

par(mfrow=c(2,1))
p_sel<-c(1,3,4)
p_plot<-rep(p_sel,c(group_No,group_No,group_No))
for(i in 2:3){
	plot(p_plot,c(coverage_regABC[[i]][,p_sel]),pch=20,ylim=c(0.7,1))
	points(p_plot,c(coverage_regACC[[i]][,p_sel]),pch=20,col=2)
	abline(h=CI_alpha,lty=2)
	legend('bottomright',c('ABC','ACC'),pch=20,col=c(1,2))
}



##########################
# Experiments for 4.1
##########################
n<-400
set.seed(100)
tested_observations_all<-simulate_pseudo_datasets(matrix(rep(parameter_true,test_size),ncol=2,byrow=T),n) # a test_size*n matrix
tested_parameters<-matrix(rep(parameter_true,test_size),ncol=2,byrow=T) # a test_size*2 matrix

data_i<-3; N<-1e4

logposterior_Cauchy_loc<-function(location,obs){
	logposterior_Cauchy(c(location,parameter_true[2]),obs)
}
MCMC_results<-metrop(logposterior_Cauchy_loc,initial=parameter_true[1],scale=0.01,nbatch=N_MCMC,obs=tested_observations_all[data_i,])
sample_posterior_loc<-MCMC_results$batch[-(1:1e5),]

# dev.new(); par(mfrow=c(4,4))

# for(data_i in 17:32){
# parameter_setting<-2

for(data_i in 1:test_size){
for(parameter_setting in 1){
dev.new()

tested_summaries_all<-cal_summaries(tested_observations_all,setting=parameter_setting,platform=platform)
tested_observations<-tested_observations_all[data_i,]
tested_summary<-tested_summaries_all[data_i,]
group_no<-sqrt(n)
subgroup<-divide(tested_observations,group_no)

set.seed(200)
if(parameter_setting==1){
	param_i<-1
	subgroup_medians<-unlist(lapply(subgroup,median)) # a N*1 subgroup_medians
	pt_estimates<-subgroup_medians^(1/3)
	simulated_parameters<-cbind(simulate_from_r2(N,pt_estimates),tested_parameters[data_i,2])
	proposal_densities<-density_r2(simulated_parameters[,1],pt_estimates)
	prior_densities<-rep(1,N)
	xlim<-c(9.7,10.3)
}
if(parameter_setting==2){
	param_i<-1
	subgroup_means<-unlist(lapply(subgroup,mean)) # a N*1 vector
	pt_estimates<-sign(subgroup_means)*abs(subgroup_means)^(1/3)
	simulated_parameters<-cbind(simulate_from_r2(N,pt_estimates),tested_parameters[data_i,2])
	proposal_densities<-density_r2(simulated_parameters[,1],pt_estimates)
	prior_densities<-rep(1,N)
	xlim<-c(6,14)
}
simulated_observations<-simulate_pseudo_datasets(simulated_parameters,n) # a N*n matrix
simulated_summaries<-cal_summaries(simulated_observations,setting=parameter_setting,platform=platform) # a N*d matrix
reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries)
weights_of_sample<-prior_densities/proposal_densities

# simulated_parameters_rejABC<-cbind(simulate_from_r1(N,location=tested_parameters[data_i,1],scale=parameter_true[2]),tested_parameters[data_i,2])
# simulated_observations_rejABC<-simulate_pseudo_datasets(simulated_parameters_rejABC,n) # a N*n matrix
# simulated_summaries_rejABC<-cal_summaries(simulated_observations_rejABC,setting=parameter_setting,platform=platform) # a N*d matrix		
# reference_tables_rejABC<-list(parameters=simulated_parameters_rejABC,summaries=simulated_summaries_rejABC)
# reference_tables<-reference_tables_rejABC
# prior_densities<-rep(1,N); 
# proposal_densities<-dcauchy(simulated_parameters_rejABC[,1],location=tested_parameters[data_i,1],scale=parameter_true[2])
# weights_of_sample<-prior_densities/proposal_densities
# plot(density(simulated_parameters_rejABC[,1]))


p_i<-9; tolerance_percentile<-p_all[p_i]
results_ISABC<-rejABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
results_ISABC<-results_ISABC[[1]] # This only deals with one observed dataset
IS_w<-weights_of_sample[results_ISABC$acceptance_ind]
IS_w<-IS_w/sum(IS_w)
results_ISregABC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample)
results_ISregABC<-results_ISregABC[[1]] # This only deals with one observed dataset
IS_w_reg<-results_ISregABC$weights/sum(results_ISregABC$weights)
results_rejACC<-results_ISABC
results_regACC<-regABC(reference_tables,t(tested_summary),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
results_regACC<-results_regACC[[1]]

# density_all<-list(posterior=density(sample_posterior_loc),ISABC=density(results_ISABC[[1]][,1]^3,weights=IS_w),ISregABC=density(results_ISregABC[[1]][,1]^3,weights=IS_w_reg),rejACC=density(results_rejACC[[1]]^3),regACC=density(results_regACC[[1]][,1]^3))

if(parameter_setting==2){
	density_all<-list(ISregABC=density(results_ISregABC[[1]][,1]^3,weights=IS_w_reg),regACC=density(results_regACC[[1]][,1]^3))
	draw_densities(density_all,color_vec=c(1,1),lty_vec=c(1,2),xlim=xlim,main='',xlab=expression(theta))
	lines(density(sample_posterior_loc^3),col=8)
}
if(parameter_setting==1){ 
	density_all<-list(posterior=density(sample_posterior_loc^3),ISregABC=density(results_ISregABC[[1]][,1]^3,weights=IS_w_reg),regACC=density(results_regACC[[1]][,1]^3))
	draw_densities(density_all,color_vec=c(8,1,1),lty_vec=c(1,1,2),xlim=xlim,main='',xlab=expression(theta))
}
}
}


