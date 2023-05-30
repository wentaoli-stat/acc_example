## Create for examples in Thorton, Li and Xie (2018) 
library(matrixStats) # For function 'rowMads'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation
library(pbapply) # For function 'pblapply'
library(mcmc)
library(MASS) # For kde2d

# setwd('/Users/wentao/bitbucket/core-functions')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\core-functions')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')


# setwd('/Users/wentao/bitbucket/acc-examples')
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\acc-examples')
source('Cauchy_regression_example_functions.R')

options(error=recover)

##########################
# Experiment Setting
##########################
### Setup for parallelisation ###
# platform<-'Mac'
platform<-'Win'

### Experiment variables ###
# test_size<-500; N<-5e4; n_all<-50 
test_size<-200; N<-1e6; n_all<-100
test_size<-200; N<-1e6; n_all<-200
p_all<-signif2(c(0.005,0.05,0.1))
param_p<-5
parameter_true<-1:param_p; # Data model is beta*rnorm(n*p)+Cauchy(0,1)
design_matr<-matrix(rnorm(n_all[1]*length(parameter_true)),n_all[1])


set.seed(100)
tested_observations_all<-simulate_pseudo_datasets(matrix(rep(parameter_true,test_size),ncol=length(parameter_true),byrow=T),n_all[1],design_matr) # a test_size*n matrix
tested_parameters<-matrix(rep(parameter_true,test_size),nrow=length(parameter_true)) # a test_size*p matrix


##########################
# Experiments 
##########################
# Parameter settings include: 1 -- param_p covariates
# Parameter proposal is 
# The prior is the Jeffary prior


CI_alpha<-0.95; prior_set<-1
raw_results_all<-list(list(0))
dev.new()
for(i in 1){
	parameter_setting<-i
	names(raw_results_all)[parameter_setting]<-paste0('setting_',parameter_setting)

	##########################
	# Set Random Seed 
	##########################
	set.seed(200)	
	tmp<-Example1_main(test_size,parameter_true,tested_observations_all=tested_observations_all,design_matr=design_matr,parameter_setting=i,minibatch=TRUE,N=N,platform=platform,n_all=n_all,p_all=p_all,CI_alpha=CI_alpha,divide=TRUE,prior_choice='jeffery')
	raw_results_all[[parameter_setting]]<-tmp
}


dev.new()
par(mfrow=c(2,3))
p_i_plot<-1; n_i<-1
for(param_i in 1:length(parameter_true)){
	param_name<-paste0('beta',param_i)
	density_ABC<-density(raw_results_all[[1]]$post_means_all[[2]][[n_i]][[p_i_plot]][,param_i])
	density_ACC<-density(raw_results_all[[1]]$post_means_all[[4]][[n_i]][[p_i_plot]][,param_i])
	density_all<-list(ABC=density_ABC,ACDC=density_ACC)
	draw_densities(density_all,main='Est_mean',with_legend=T,xlab=param_name,xlim=c(-10,10))
	abline(v=parameter_true[param_i],lty=2,col=3)
}

dev.new()
par(mfrow=c(2,3))
p_i_plot<-1; n_i<-1
for(param_i in 1:length(parameter_true)){
	param_name<-paste0('beta',param_i)
	density_ABC<-density(sqrt(raw_results_all[[1]]$post_variances_all[[2]][[n_i]][[p_i_plot]][,param_i]))
	density_ACC<-density(sqrt(raw_results_all[[1]]$post_variances_all[[4]][[n_i]][[p_i_plot]][,param_i]))
	density_all<-list(ABC=density_ABC,ACDC=density_ACC)
	draw_densities(density_all,main='Est_sd',with_legend=T,xlab=param_name,xlim=c(0,10))
}

signif2(rowMeans((raw_results_all[[i]]$coverage_all)[[2]])) # coverage of IS-ABC
signif2(rowMeans((raw_results_all[[i]]$coverage_all)[[4]])) # coverage of ACC
signif2(rowMeans((raw_results_all[[i]]$coverage1_all)[[4]])) # coverage of ACC

vol_ratio<-(raw_results_all[[i]]$wid_vol_all)[[4]]/(raw_results_all[[i]]$wid_vol_all)[[2]]
signif2(rowMedians(vol_ratio)) # Ratios of volume
signif2(rowMads((raw_results_all[[i]]$wid_vol_all)[[4]]/(raw_results_all[[i]]$wid_vol_all)[[2]])) # Ratios of volume


vol_acc<-(raw_results_all[[i]]$wid_vol_all)[[4]]
vol_abc<-(raw_results_all[[i]]$wid_vol_all)[[2]]
# cbind(vol_acc[2,],vol_abc[2,])
signif2(rowMedians(vol_abc)) # volume of IS-ABC
signif2(rowMedians(vol_acc)) # volume of ACC
signif2(rowMedians(vol_acc)/rowMedians(vol_abc))


