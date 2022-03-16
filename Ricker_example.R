## Create for examples in Thorton, Li and Xie (2018) 

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


##########################
# Experiment Setting
############ ##############
setwd('C:\\Users\\liwen\\OneDrive\\Documents\\acc-examples')
# setwd('/Users/wentao/Dropbox')
# load('Ricker_initial_cal_167rep_N1e5_n500_half_rn_PT.RData') # Load the observed datasets and the initial ballpark distribution for each dataset; Load the variables 'nobs', 

load('Ricker_SubgroupRefinedEst_150rep_n50.RData')
# load('Ricker_SubgroupRefinedEst_150rep_n500.RData')
param_ests<-list(0)

# for(data_i in 1:test_size)	param_ests[[data_i]]<-subgroup_refine_est[[data_i]]$mean
# rn_method<-'refinedPT'; rn_sub<-'sequential_groupsize10'

# param_ests<-subgroup_means
# rn_method<-'semi_refinedPT'; rn_sub<-'sequential_groupsize10'

for(data_i in 1:test_size){
 	param_ests[[data_i]]<-rbind(subgroup_param_est[[data_i]],subgroup_refine_est[[data_i]]$mean)
 	param_ests[[data_i]][,2]<-(subgroup_refine_est[[data_i]]$mean)[,2]
}
rn_method<-'comb_refinedPT'; rn_sub<-'sequential_groupsize10'

# load('Ricker_SubgroupEst_150rep_n50.RData')
# param_ests<-subgroup_param_est
# rn_method<-'PT'; rn_sub<-'sequential_groupsize10'

### Setup for parallelisation ###
# platform<-'Mac'; ncores<-3
platform<-'Win'; ncores<-6

### Experiment variables ###
test_size<-150; N<-10^6; 
# test_size<-150; N<-10^5; 
# test_size<-300; N<-10^5; 
# test_size<-500; N<-5*10^5
p_all<-c(0.005,0.01,0.05,0.1)
parameter_true<-c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)) 
ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = parameter_true,extraArgs = list("nObs" = nobs, "nBurn" = 50))
CI_alpha<-0.95
# CI_alpha<-0.90

# tmp<-Ricker_ABC_ACC_algorithms(1,N=N,tested_observations_all=tested_observations_all,param_ests=param_ests,p_all=p_all,CI_alpha=CI_alpha)

set.seed(200)
no_parts<-10; run_parts<-divide(1:test_size,no_parts)
for(i in 1:no_parts){
	tmp<-parallel_lapply(run_parts[[i]],Ricker_ABC_ACC_algorithms,ncores=ncores,platform=platform,packages=c('synlik','parallel','matrixStats','ks','spatstat'),N=N,tested_observations_all=tested_observations_all,param_ests=param_ests,p_all=p_all,CI_alpha=CI_alpha,progress_bar=TRUE)
	save(tmp,file=paste0('Ricker_results_',test_size,'rep_N',N,'_n',nobs,'_rn_',rn_method,'_',rn_sub,'_tmp',i,'.RData'))
	rm(tmp); gc()
}

tmp_all<-NULL
for(i in 1:no_parts){
	load(file=paste0('Ricker_results_',test_size,'rep_N',N,'_n',nobs,'_rn_',rn_method,'_',rn_sub,'_tmp',i,'.RData'))
	for(i in 1:length(tmp)){
		tmp[[i]]$reference_tables<-NULL; gc()
	}	
	tmp_all<-c(tmp_all,tmp)
}

# Reformulate the structure of tmp_all so it fits the analysis code 
post_means_all<-list(0); post_variances_all<-list(0);
tmp1<-matrix(0,test_size,length(param_est_ind))
for(method_i in 1:6){
	post_means_all[[method_i]]<-list(0)
	if(method_i==1) names(post_means_all)[method_i]<-'ISABC'
	if(method_i==2) names(post_means_all)[method_i]<-'ISregABC'
	if(method_i==3) names(post_means_all)[method_i]<-'ISpi-ACC'
	if(method_i==4) names(post_means_all)[method_i]<-'ISpi-regACC'
	if(method_i==5) names(post_means_all)[method_i]<-'rejACC'
	if(method_i==6) names(post_means_all)[method_i]<-'regACC'
	for(p_i in 1:length(p_all)){
		post_means_all[[method_i]][[p_i]]<-tmp1
		names(post_means_all[[method_i]])[p_i]<-paste0('acceptance rate=',p_all[p_i])
	}
}
post_variances_all<-post_means_all
post_coverage1D_all<-list(0)
for(param_i in param_est_ind){
	post_coverage1D_all[[param_i]]<-post_means_all
	for(method_i in 1:6) post_coverage1D_all[[param_i]][[method_i]]<-matrix(0,length(p_all),test_size)
}
names(post_coverage1D_all)<-names(parameter_true)
post_coverage2D_all<-list(0); twoD_pairs<-rbind(c(1,2),c(1,3))
for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]]<-post_coverage1D_all[[1]]
names(post_coverage2D_all)<-c('logR_and_logSigma','logR_and_logPhi')
param_coverage_type<-c('2D','1D') 
number_coverage_type<-5; type_name<-c(names(post_coverage2D_all),names(post_coverage1D_all))
if('1D'%in%param_coverage_type) oneD_covering<-TRUE else oneD_covering<-FALSE
if('2D'%in%param_coverage_type) twoD_covering<-TRUE else twoD_covering<-FALSE
post_coverage_onesideupper_all<-post_coverage1D_all
post_coverage_onesidelower_all<-post_coverage1D_all
wald_coverage<-post_coverage1D_all

CIoneside_upper_all<-post_coverage1D_all
CIoneside_lower_all<-post_coverage1D_all
CItwoside_upper_all<-post_coverage1D_all
CItwoside_lower_all<-post_coverage1D_all
CIwald_upper_all<-post_coverage1D_all
CIwald_lower_all<-post_coverage1D_all

IS_ESSr_all<-post_coverage1D_all[[1]][[1]]
reference_tables_all<-list(0); weights_of_sample_all<-list(0); tested_summary_less_all<-list(0)
for(data_i in 1:test_size){
	# reference_tables_all[[data_i]]<-tmp_all[[data_i]]$reference_tables
	weights_of_sample_all[[data_i]]<-tmp_all[[data_i]]$weights_of_sample
	tested_summary_less_all[[data_i]]<-tmp_all[[data_i]]$tested_summary_less
	
	post_coverage<-tmp_all[[data_i]]$post_coverage
	CIbounds<-tmp_all[[data_i]]$CIbounds
	for(p_i in 1:length(p_all)){
		IS_ESSr_all[p_i,data_i]<-tmp_all[[data_i]]$IS_ESSr[[p_i]]
		for(method_i in 1:6){
			post_means_all[[method_i]][[p_i]][data_i,]<-tmp_all[[data_i]]$post_means[[method_i]][[p_i]]
			post_variances_all[[method_i]][[p_i]][data_i,]<-tmp_all[[data_i]]$post_variances[[method_i]][[p_i]]
			for(param_i in 1:3){
				post_coverage1D_all[[param_i]][[method_i]][p_i,data_i]<-post_coverage$twoside[[param_i]][[method_i]][p_i]
				CItwoside_upper_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$twosideupper[[param_i]][[method_i]][p_i]
				CItwoside_lower_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$twosidelower[[param_i]][[method_i]][p_i]

				wald_coverage[[param_i]][[method_i]][p_i,data_i]<-post_coverage$wald[[param_i]][[method_i]][p_i]
				CIwald_upper_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$wald_upper[[param_i]][[method_i]][p_i]
				CIwald_lower_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$wald_lower[[param_i]][[method_i]][p_i]

				post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,data_i]<-post_coverage$onesideupper[[param_i]][[method_i]][p_i]
				CIoneside_upper_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$onesideupper[[param_i]][[method_i]][p_i]			
				post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,data_i]<-post_coverage$onesidelower[[param_i]][[method_i]][p_i]
				CIoneside_lower_all[[param_i]][[method_i]][p_i,data_i]<-CIbounds$onesidelower[[param_i]][[method_i]][p_i]	
			}
			for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,data_i]<-post_coverage$twoD[[pair_i]][[method_i]][p_i]
		}
	}	
}

rm(tmp); rm(tmp_all); gc()
save.image(file=paste0('Ricker_results_',test_size,'rep_N',N,'_n',nobs,'_rn_',rn_method,'_',rn_sub,'_alpha',CI_alpha*100,'.RData'))

p_i<-1
par(mfrow=c(2,2))
for(par_i in 1:3){
	density1<-density(post_variances_all[[6]][[p_i]][1:30,par_i])
	density2<-density(post_variances_all[[2]][[p_i]][1:30,par_i])
	draw_densities(list(density1,density2),main=names(parameter_true)[par_i])
	abline(v=var(post_means_all[[6]][[p_i]][1:30,par_i]),lty=2)
	abline(v=var(post_means_all[[2]][[p_i]][1:30,par_i]),lty=2,col=2)
	abline(v=median(post_variances_all[[6]][[p_i]][1:30,par_i]))
	abline(v=median(post_variances_all[[2]][[p_i]][1:30,par_i]),col=2)
	legend('topright',c('regACC','regABC'),lty=c(1,1),col=1:2)
}

post_coverage1D_all[[2]][[2]][,1:30]
post_coverage1D_all[[2]][[6]][,1:30]

var(post_means_all[[2]][[2]][1:30,2])
mean(post_variances_all[[2]][[2]][1:30,2])

var(post_means_all[[6]][[2]][1:30,2])
mean(post_variances_all[[6]][[2]][1:30,2])


##########################
# Results
##########################
# Under each setting, the raw result is a list with four elements: the posterior mean estimates, posterior variance estimates, coverage rate of the true parameter value in the test_size runs, and the tolerance values for all runs. Each element is a list with the acceptance proportions p_all, and the basic element is a test_size*2 matrix, for moments estimate, or a length test_size vector.

number_methods<-6 # 'ISABC', 'ISregABC', 'rejACC', 'regACC'
results_all<-list(list(0,0),list(0),list(0),list(0),list(0)) # Corresponds to five experiment settings
tmp_matrix<-matrix(0,number_methods,length(p_all))
colnames(tmp_matrix)<-paste0('p=',p_all)
rownames(tmp_matrix)<-c('ISABC','ISregABC','ISpi_ACC','ISpi_regACC','rejACC','regACC')
tmp<-list(0); param_est<-1:3; param_count<-1
for(param_i in param_est){
	tmp[[param_count]]<-list(tmp_matrix,tmp_matrix,tmp_matrix,tmp_matrix)
	if(param_i==1) names(tmp)[param_count]<-c('logR')
	if(param_i==2) names(tmp)[param_count]<-c('logSigma')
	if(param_i==3) names(tmp)[param_count]<-c('logPhi')
	names(tmp[[param_count]])<-c('mean_med','mean_mad','var_med','var_mad')
	param_count<-param_count+1
}
param_count<-param_count-1
for(coverge_i in 1:number_coverage_type){
	tmp[[param_count+coverge_i]]<-tmp_matrix
	names(tmp)[param_count+coverge_i]<-paste0(type_name[coverge_i],'_coverage')
}
for(coverge_i in 1:number_coverage_type){
	tmp[[param_count+number_coverage_type+coverge_i]]<-tmp_matrix
	names(tmp)[param_count+number_coverage_type+coverge_i]<-paste0(type_name[coverge_i],'_CIwidth')
}
output_all<-list(p=p_all,post_means=post_means_all,post_variances=post_variances_all,post_coverage=c(post_coverage2D_all,post_coverage1D_all),CIbounds=list(twosideupper=CItwoside_upper_all,twosidelower=CItwoside_lower_all,onesideupper=CIoneside_upper_all,onesidelower=CIoneside_lower_all),results_all=tmp)
results_all<-Process_all_output(1:test_size,output_all,param_est=param_est,number_coverage_type=number_coverage_type)


p_i_table<-c(1,2,3,4); method_adj<-c(2,4,6)
tmp_matrix<-matrix(0,nrow=2*length(p_i_table),ncol=8)
rownames(tmp_matrix)<-paste0('p=',rep(p_all[p_i_table],rep(2,length(p_i_table))),' ',rep(c('confidence','credible'),length(p_i_table)))
colnames(tmp_matrix)<-c('rej-ABC','width','IS-ABC','width','piACC','width','rACC','width')
tmp_matrix[2*p_i_table,5:8]<-NaN; 
results_table<-list(0)
for(coverge_i in 1:number_coverage_type){
	results_table[[coverge_i]]<-list(0)
	names(results_table)[coverge_i]<-paste0(type_name[coverge_i],'_coverage')
	results_table[[coverge_i]]<-tmp_matrix
	results_table[[coverge_i]][2*p_i_table-1,c(3,5,7)]<-signif(t(results_all[[length(param_est)+coverge_i]][method_adj,p_i_table]),3)
	results_table[[coverge_i]][2*p_i_table-1,c(4,6,8)]<-signif(t(results_all[[length(param_est)+number_coverage_type+coverge_i]][method_adj,p_i_table]),3)
}		
results_table


line_col<-c('black','black','red','red','green','green','blue','blue'); line_type<-c(1,1,1,1,1,1,3,3); line_width<-c(1,2,1,2,1,2,1,2); lines_names<-c('IS ABC/No Adj','IS ABC/Adj','IS piACC/No Adj','IS piACC/Adj','ACC/No Adj','ACC/Adj','True Value')
method_adj<-c(2,4,6); 
adj_line_col<-line_col[c(method_adj,(1:length(line_col))[-(1:6)])]; adj_line_type<-line_type[c(method_adj,(1:length(line_col))[-(1:6)])]; adj_line_width<-line_width[c(method_adj,(1:length(line_col))[-(1:6)])]; adj_lines_names<-lines_names[c(method_adj,(1:length(line_col))[-(1:6)])]
min_max_quantile<-1
# setEPS()
# postscript('Ricker_example.eps')#,width=7,height=0.8*7)
dev.new()
par(mfcol=c(3,4))
param_est<-1:3; param_num<-length(param_est)
for(coverage_i in 1:number_coverage_type){
	if(('2D'%in%param_coverage_type)&coverage_i==3) plot(1,1,type="n", axes=F, xlab="", ylab="")
	if(coverage_i==1) setting_name<-paste0(data_i,' rep / N ',N,' / ',type_name[coverage_i],' / ',rn_sub,' ',rn_method)
	if(coverage_i>1) setting_name<-type_name[coverage_i]
	yname<-'Coverage'
	coverage_matrix<-cbind(t(results_all[[param_num+coverage_i]]),CI_alpha)
	if(coverage_i==1) with_legend<-TRUE else with_legend<-FALSE
	draw_lines(p_all,coverage_matrix,line_col=line_col,line_type=line_type,line_width=line_width,main=setting_name,ylim=c(0.8,1),xname='Acceptance Proportion',yname=yname,with_legend=with_legend,line_names=lines_names)
}
with_legend<-FALSE
# ylim_matrix<-list(rbind(c(3.6,4),c(-2,-1),c(2.25,2.35)),0,rbind(c(0,0.14),c(0,0.6),c(0,0.04)))
ylim_matrix<-list(rbind(c(3,4.5),c(-2,0),c(2,3)),0,rbind(c(0,0.6),c(0,0.8),c(0,0.3)))
for(mean_or_std in c(1,3)){
	for(param_i in 1:param_num){
		results_plot<-results_all[[param_i]]
		param_name<-names(parameter_true)[param_i]
		est_med_matrix<-t(results_plot[[mean_or_std]][method_adj,])
		est_mad_matrix<-t(results_plot[[mean_or_std+1]][method_adj,])
		if(mean_or_std==1){
			yname<-'Estimated Mean'
			truth<-parameter_true[param_est[param_i]]
			est_med_matrix<-cbind(est_med_matrix,truth=truth)
			est_mad_matrix<-cbind(est_mad_matrix,0)
		}
		if(mean_or_std==3){
			yname<-'Estimated STD'
			est_med_matrix<-cbind(est_med_matrix,med_regACC_mad=results_plot[[mean_or_std-1]][6,])
			est_mad_matrix<-cbind(est_mad_matrix,0)
		}
		draw_lines(p_all,est_med_matrix,line_col=adj_line_col,line_type=adj_line_type,line_width=adj_line_width,main=param_name,xname='Acceptance Proportion',yname=yname,ylim=ylim_matrix[[mean_or_std]][param_i,],with_legend=with_legend,line_names=adj_lines_names,with_CI=TRUE,lines_sd_matrix=est_mad_matrix,min_max_quantile=min_max_quantile)
		# draw_lines(p_all,est_med_matrix,line_col=adj_line_col,line_type=adj_line_type,line_width=adj_line_width,main=param_name,xname='Acceptance Proportion',yname=yname,with_legend=with_legend,line_names=adj_lines_names,with_CI=TRUE,lines_sd_matrix=est_mad_matrix,min_max_quantile=min_max_quantile)
	}
}
# dev.off()

### Random seed ###
set.seed(200)

####################################################
# Experiment (SEQUENTIAL); need reference_tables_all,
# weights_of_sample_all and tested_summary_less_all 
# to be simulated beforehand
####################################################
### Initialisation ###
param_est_ind<-1:3
post_means_all<-list(0); post_variances_all<-list(0);
tmp1<-matrix(0,test_size,length(param_est_ind))
for(method_i in 1:6){
	post_means_all[[method_i]]<-list(0)
	if(method_i==1) names(post_means_all)[method_i]<-'ISABC'
	if(method_i==2) names(post_means_all)[method_i]<-'ISregABC'
	if(method_i==3) names(post_means_all)[method_i]<-'ISpi-ACC'
	if(method_i==4) names(post_means_all)[method_i]<-'ISpi-regACC'
	if(method_i==5) names(post_means_all)[method_i]<-'rejACC'
	if(method_i==6) names(post_means_all)[method_i]<-'regACC'
	for(p_i in 1:length(p_all)){
		post_means_all[[method_i]][[p_i]]<-tmp1
		names(post_means_all[[method_i]])[p_i]<-paste0('acceptance rate=',p_all[p_i])
	}
}
post_variances_all<-post_means_all
post_coverage1D_all<-list(0)
for(param_i in param_est_ind){
	post_coverage1D_all[[param_i]]<-post_means_all
	for(method_i in 1:6) post_coverage1D_all[[param_i]][[method_i]]<-matrix(0,length(p_all),test_size)
}
names(post_coverage1D_all)<-names(parameter_true)
post_coverage2D_all<-list(0); twoD_pairs<-rbind(c(1,2),c(1,3))
for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]]<-post_coverage1D_all[[1]]
names(post_coverage2D_all)<-c('logR_and_logSigma','logR_and_logPhi')

post_coverage_onesideupper_all<-post_coverage1D_all
post_coverage_onesidelower_all<-post_coverage1D_all
CIoneside_upper_all<-post_coverage1D_all
CIoneside_lower_all<-post_coverage1D_all
CItwoside_upper_all<-post_coverage1D_all
CItwoside_lower_all<-post_coverage1D_all

IS_ESSr_all<-post_coverage1D_all[[1]][[1]] # The effective sample size ratio, i.e. ESS/N0, of the importance sampling ABC
param_coverage_type<-c('2D','1D') # Indicate what type of cofidence region should be calculated, '1D' for individual parameter, '2D' or '3D' for two or three parameters jointly.
number_coverage_type<-5; type_name<-c(names(post_coverage2D_all),names(post_coverage1D_all))
if('1D'%in%param_coverage_type) oneD_covering<-TRUE else oneD_covering<-FALSE
if('2D'%in%param_coverage_type) twoD_covering<-TRUE else twoD_covering<-FALSE


### Algorithms ###
dev.new()
for(data_i in 1:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	plot(0,data_i,xlim=c(1,length(p_all)),ylim=c(data_i,data_i),xlab='p_i',ylab='data_i',main=paste0('Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	indice_i<-data_i

	tested_value<-parameter_true[param_est_ind]
	tested_observations<-tested_observations_all[indice_i,]
	reference_tables<-reference_tables_all[[indice_i]]
	weights_of_sample<-weights_of_sample_all[[indice_i]]
	tested_summary_less<-tested_summary_less_all[[indice_i]]

	for(p_i in 1:length(p_all)){ # for 50 p_i, running time is 130s
		points(p_i,indice_i,pch=20)
		tolerance_percentile<-p_all[p_i]
		#################################### 
		# importance sampling ABC (ISABC)
		####################################
		method_i<-1
		results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		results_ISABC<-results_ISABC[[1]] # This only deals with one observed dataset
		importance_weights<-weights_of_sample[results_ISABC$acceptance_ind] 
		results_summary<-ABC_ACC_results_summary(results_ISABC,type='ABC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		IS_ESSr_all[p_i,indice_i]<-results_summary$IS_ESS
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# importance sampling regression ABC (ISregABC)
		###################################################
		method_i<-2
		results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample)
		results_ISregABC<-results_ISregABC[[1]] # This only deals with one observed dataset
		importance_weights<-results_ISregABC$weights
		results_summary<-ABC_ACC_results_summary(results_ISregABC,type='ABC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# importance sampling prior-ACC (ISpi-ACC)
		###################################################
		method_i<-3
		results_ISpiACC<-results_ISABC
		importance_weights<-weights_of_sample[results_ISABC$acceptance_ind]
		results_summary<-ABC_ACC_results_summary(results_ISpiACC,type='ACC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)


		### This part is identical for all method_i ###
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# importance sampling prior-regression ACC (ISpi_regACC)
		###################################################
		method_i<-4
		results_ISpi_regACC<-results_ISregABC
		importance_weights<-results_ISregABC$weights
		results_summary<-ABC_ACC_results_summary(results_ISpi_regACC,type='ACC',w=importance_weights,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# rejection ACC (rejACC)
		###################################################
		method_i<-5
		results_rejACC<-results_ISABC
		results_summary<-ABC_ACC_results_summary(results_rejACC,type='ACC',w=NULL,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)


		### This part is identical for all method_i ###
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]

		###################################################
		# regression ACC (regACC)
		###################################################
		method_i<-6
		results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		results_regACC<-results_regACC[[1]]
		results_summary<-ABC_ACC_results_summary(results_regACC,type='ACC',w=NULL,true_params=tested_value,oneD_covering=oneD_covering,twoD_covering=twoD_covering,twoD_pairs=twoD_pairs,CI_alpha=CI_alpha)

		### This part is identical for all method_i ###
		post_means_all[[method_i]][[p_i]][indice_i,]<-results_summary$mean_vec
		post_variances_all[[method_i]][[p_i]][indice_i,]<-results_summary$var_vec
		for(param_i in param_est_ind){
			twoside_CI<-results_summary$twoside_CI
			post_coverage1D_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_coverage[param_i]
			CItwoside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_upper[param_i]
			CItwoside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-twoside_CI$twoside_lower[param_i]

			oneside_upperCI<-results_summary$oneside_upperCI
			post_coverage_onesideupper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_coverage[param_i]
			CIoneside_upper_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_upperCI$oneside_upper[param_i]			
			oneside_lowerCI<-results_summary$oneside_lowerCI
			post_coverage_onesidelower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_coverage[param_i]
			CIoneside_lower_all[[param_i]][[method_i]][p_i,indice_i]<-oneside_lowerCI$oneside_lower[param_i]	
		}
		for(pair_i in 1:nrow(twoD_pairs)) post_coverage2D_all[[pair_i]][[method_i]][p_i,indice_i]<-results_summary$twoD_CR[[pair_i]]
	}

	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	


######################################################
# Ricker example: N1e5, test_size 300, PT, 0.95 alpha
######################################################

}
# raw_results_all<-list(post_means_all=post_means_all,post_variances_all=post_variances_all,post_coverage_all=c(post_coverage2D_all,post_coverage1D_all),IS_ESSr_all=IS_ESSr_all,CIupper_all=CIupper_all)
raw_results_all<-list(post_means_all=post_means_all,post_variances_all=post_variances_all,post_coverage_all=list(twoside=post_coverage1D_all,onesideupper=post_coverage_onesideupper_all,onesidelower=post_coverage_onesidelower_all,twoD=post_coverage2D_all),IS_ESSr_all=IS_ESSr_all,CIbounds_all=list(twosideupper=CItwoside_upper_all,twosidelower=CItwoside_lower_all,onesideupper=CIoneside_upper_all,onesidelower=CIoneside_lower_all))


rm(reference_tables_all); rm(weights_of_sample_all); rm(tested_summary_less_all); gc()
save.image('Ricker_results_300rep_N1e5_n500_half_rn_PT_allcoverage_80.RData')



##########################
##########################
# Debugging
#
##########################
##########################


###################################################################################
# Run IS-ABC on 50 datasets, using r_n constructed from \sqrt{n} and n^{3/5}
# subgroups to check why the variance of \hat{\theta} given by the latter is 
# too small
###################################################################################

# load('Ricker_initialcal_300rep_N1e5_n500_half_rn_PT.RData')
# load('Ricker_initialcal_300rep_N1e5_n500_threefifth_rn_PT.RData')

test_size<-50; N<-10^5; 

reference_tables_all<-list(0)
weights_of_sample_all<-list(0)
tested_summary_less_all<-list(0)
set.seed(200)
dev.new()
for(data_i in 1:test_size){
	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	plot(data_i,data_i,xlim=c(1,test_size),ylim=c(1,test_size),xlab='data_i',ylab='data_i',main=paste0('Simulation: Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	# Simulate the parameter values  and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
	simulated_parameters<-NULL
	proposal_densities<-1
	prior_densities<-1
	for(param_i in param_est_ind){
		simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,subgroup_param_est[[data_i]][,param_i]))
		proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],subgroup_param_est[[data_i]][,param_i])
		prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
	}
	# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
	simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
	tested_summary_less<-tested_summary[-(1:3)]
	reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries[,-(1:3)])
	weights_of_sample<-prior_densities/proposal_densities

	reference_tables_all[[data_i]]<-reference_tables
	weights_of_sample_all[[data_i]]<-weights_of_sample
	tested_summary_less_all[[data_i]]<-tested_summary_less

	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	
}


load('tested_summary_less_all_n500.RData')


test_size<-50; N<-10^5; 
# setwd('/Users/wentao/bitbucket/')
setwd('D:\\bitbucket\\core-functions')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')
# setwd('/Users/wentao/bitbucket/ACC examples')
setwd('D:\\bitbucket\\acc-examples')
source('Ricker_example_functions.R')


platform<-'Win'; ncores<-6
p_all<-signif2(seq(0.001,0.5,length=20))
parameter_true<-c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)) 
CI_alpha<-0.90

set.seed(200)
param_est_ind<-1:3
post_means_all<-list(0); post_variances_all<-list(0);
tmp1<-matrix(0,test_size,length(param_est_ind))
for(method_i in 1:4){
	post_means_all[[method_i]]<-list(0)
	if(method_i==1) names(post_means_all)[method_i]<-'ISABC'
	if(method_i==2) names(post_means_all)[method_i]<-'ISregABC'
	if(method_i==3) names(post_means_all)[method_i]<-'rejACC'
	if(method_i==4) names(post_means_all)[method_i]<-'regACC'
	for(p_i in 1:length(p_all)){
		post_means_all[[method_i]][[p_i]]<-tmp1
		names(post_means_all[[method_i]])[p_i]<-paste0('acceptance rate=',p_all[p_i])
	}
}
post_variances_all<-post_means_all
post_coverage1D_all<-list(0)
for(param_i in param_est_ind){
	post_coverage1D_all[[param_i]]<-post_means_all
	for(method_i in 1:4) post_coverage1D_all[[param_i]][[method_i]]<-matrix(0,length(p_all),test_size)
}
names(post_coverage1D_all)<-names(parameter_true)

CIupper_all<-post_coverage1D_all
IS_ESSr_all<-post_coverage1D_all[[1]][[1]] # The effective sample size ratio, i.e. ESS/N0, of the importance sampling ABC
param_coverage_type<-c('1D') # Indicate what type of cofidence region should be calculated, '1D' for individual parameter, '2D' or '3D' for two or three parameters jointly.
number_coverage_type<-3; type_name<-c(names(post_coverage1D_all))

### Algorithms ###
dev.new()
for(data_i in 1:test_size){

	##########################
	#### Recrod start time ###
	start_time<-proc.time()
	##########################

	plot(0,data_i,xlim=c(1,length(p_all)),ylim=c(data_i,data_i),xlab='p_i',ylab='data_i',main=paste0('Run ',data_i,'/',test_size,'; last run costs ',duration[3],' seconds')) # Plot the progress of experiment, showing the index of running dataset, p_i and the running time of last dataset
	indice_i<-data_i

	tested_value<-parameter_true[param_est_ind]
	tested_observations<-tested_observations_all[indice_i,]
	reference_tables<-list(parameters=reference_tables_all[[indice_i]][[1]][1:N,],summaries=reference_tables_all[[indice_i]][[2]][1:N,])
	weights_of_sample<-weights_of_sample_all[[indice_i]][1:N]
	tested_summary_less<-tested_summary_less_all[[indice_i]]

	for(p_i in 1:length(p_all)){ # for 50 p_i, running time is 130s
		points(p_i,indice_i,pch=20)
		tolerance_percentile<-p_all[p_i]
		### importance sampling ABC (ISABC)###
		results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform) 
		weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind] # This only deals with one observed dataset
		IS_ESSr_all[p_i,indice_i]<-ESS(weights_ISABC)/length(weights_ISABC)
		post_means_all[[1]][[p_i]][indice_i,]<-Process_outputlist(results_ISABC,function(x) wmean(x$parameters,weights_ISABC),rbind)
		post_variances_all[[1]][[p_i]][indice_i,]<-Process_outputlist(results_ISABC,function(x) wvar(x$parameters,weights_ISABC,diagonal=T),rbind)
		if('1D'%in%param_coverage_type){
			for(param_i in param_est_ind){
				upper<-Process_outputlist(results_ISABC,function(x)wquant(x$parameters[,param_i],q=CI_alpha,weights=weights_ISABC),c)
				post_coverage1D_all[[param_i]][[1]][p_i,indice_i]<-(tested_value[param_i]<=upper)
				CIupper_all[[param_i]][[1]][p_i,indice_i]<-upper
			} 					
		}
		### importance sampling regression ABC (ISregABC)###
		results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
		post_means_all[[2]][[p_i]][indice_i,]<-Process_outputlist(results_ISregABC,function(x) wmean(x$parameters,x$weights),rbind)
		post_variances_all[[2]][[p_i]][indice_i,]<-Process_outputlist(results_ISregABC,function(x) wvar(x$parameters,x$weights,diagonal=T),rbind)
		if('1D'%in%param_coverage_type){
			for(param_i in param_est_ind){
				upper<-Process_outputlist(results_ISregABC,function(x)wquant(x$parameters[,param_i],q=CI_alpha,weights=x$weights),c)
				post_coverage1D_all[[param_i]][[2]][p_i,indice_i]<-(tested_value[param_i]<=upper)
				CIupper_all[[param_i]][[2]][p_i,indice_i]<-upper
			}
		}
		### rejection ACC (rejACC)###
		results_rejACC<-results_ISABC
		post_means_all[[3]][[p_i]][indice_i,]<-Process_outputlist(results_rejACC,function(x) colMeans(x$parameters),rbind)
		post_variances_all[[3]][[p_i]][indice_i,]<-Process_outputlist(results_rejACC,function(x) colVars(x$parameters),rbind)
		ACC_mean<-post_means_all[[3]][[p_i]][indice_i,]
		if('1D'%in%param_coverage_type){
			for(param_i in param_est_ind){
				upper<-Process_outputlist(results_rejACC,function(x)onesideCI_cal(x$parameters[,param_i],CI_alpha),c)
				post_coverage1D_all[[param_i]][[3]][p_i,indice_i]<-(tested_value[param_i]<=upper)
				CIupper_all[[param_i]][[3]][p_i,indice_i]<-upper
			}
		}

		### regression ACC (regACC)###
		results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		# results_regACC<-results_ISregABC
		post_means_all[[4]][[p_i]][indice_i,]<-Process_outputlist(results_regACC,function(x) colMeans(x$parameters),rbind)
		post_variances_all[[4]][[p_i]][indice_i,]<-Process_outputlist(results_regACC,function(x) colVars(x$parameters),rbind)
		ACC_mean<-post_means_all[[4]][[p_i]][indice_i,]
		if('1D'%in%param_coverage_type){
			for(param_i in param_est_ind){
				upper<-Process_outputlist(results_regACC,function(x)onesideCI_cal(x$parameters[,param_i],CI_alpha),c)
				post_coverage1D_all[[param_i]][[4]][p_i,indice_i]<-(tested_value[param_i]<=upper)
				CIupper_all[[param_i]][[4]][p_i,indice_i]<-upper
			}
		}
	}
	##########################
	#### Recrod running time ###
	duration<-signif(proc.time()-start_time,2)
	##########################	


######################################################
# Ricker example: N1e5, test_size 50, PT, half
######################################################

}
raw_results_all<-list(post_means_all=post_means_all,post_variances_all=post_variances_all,post_coverage_all=c(post_coverage1D_all),IS_ESSr_all=IS_ESSr_all,CIupper_all=CIupper_all)


number_methods<-4
tmp_matrix<-matrix(0,number_methods,length(p_all))
colnames(tmp_matrix)<-paste0('p=',p_all)
rownames(tmp_matrix)<-c('ISABC','ISregABC','rejACC','regACC')
tmp<-list(0); param_est<-1:3; param_count<-1
for(param_i in param_est){
	tmp[[param_count]]<-list(tmp_matrix,tmp_matrix,tmp_matrix,tmp_matrix)
	if(param_i==1) names(tmp)[param_count]<-c('logR')
	if(param_i==2) names(tmp)[param_count]<-c('logSigma')
	if(param_i==3) names(tmp)[param_count]<-c('logPhi')
	names(tmp[[param_count]])<-c('mean_med','mean_mad','var_med','var_mad')
	param_count<-param_count+1
}
param_count<-param_count-1
for(coverge_i in 1:number_coverage_type){
	tmp[[param_count+coverge_i]]<-tmp_matrix
	names(tmp)[param_count+coverge_i]<-paste0(type_name[coverge_i],'_coverage')
}
output_all<-list(p=p_all,post_means=raw_results_all$post_means_all,post_variances=raw_results_all$post_variances_all,post_coverage=raw_results_all$post_coverage_all,results=tmp)
results_all<-Process_all_output(1:50,output_all,param_est=param_est,number_coverage_type=number_coverage_type)



line_col<-c('black','black','red','red','blue','blue'); line_type<-c(1,1,1,1,3,3); line_width<-c(1,2,1,2,1,2); lines_names<-c('IS ABC/No Adj','IS ABC/Adj','ACC/No Adj','ACC/Adj','True Value')
min_max_quantile<-1
method_adj<-c(2,4); 
adj_line_col<-line_col[c(method_adj,(1:length(line_col))[-(1:4)])]; adj_line_type<-line_type[c(method_adj,(1:length(line_col))[-(1:4)])]; adj_line_width<-line_width[c(method_adj,(1:length(line_col))[-(1:4)])]; adj_lines_names<-lines_names[c(method_adj,(1:length(line_col))[-(1:4)])]
dev.new()
par(mfcol=c(3,3))
param_est<-1:3; param_num<-length(param_est)
for(coverage_i in 1:number_coverage_type){
	if(coverage_i==1) setting_name<-paste0(test_size,' rep / N ',N,' / ',rn_sub,' ',rn_method)
	if(coverage_i>1) setting_name<-type_name[coverage_i]
	yname<-'Coverage'
	coverage_matrix<-cbind(t(results_all[[param_num+coverage_i]]),CI_alpha)
	if(coverage_i==1) with_legend<-TRUE else with_legend<-FALSE
	draw_lines(p_all,coverage_matrix,line_col=line_col,line_type=line_type,line_width=line_width,main=setting_name,ylim=c(0.5,1),xname='Acceptance Proportion',yname=yname,with_legend=with_legend,line_names=lines_names)
}
ylim_matrix<-list(rbind(c(3.6,4),c(-2,-1),c(2.25,2.35)),0,rbind(c(0,0.14),c(0,0.6),c(0,0.04)))
with_legend<-FALSE
for(mean_or_std in c(1,3)){
	for(param_i in 1:param_num){
		results_plot<-results_all[[param_i]]
		param_name<-names(parameter_true)[param_i]
		est_med_matrix<-t(results_plot[[mean_or_std]][method_adj,])
		est_mad_matrix<-t(results_plot[[mean_or_std+1]][method_adj,])
		if(mean_or_std==1){
			yname<-'Estimated Mean'
			truth<-parameter_true[param_est[param_i]]
			est_med_matrix<-cbind(est_med_matrix,truth=truth)
			est_mad_matrix<-cbind(est_mad_matrix,0)
		}
		if(mean_or_std==3){
			yname<-'Estimated STD'
			est_med_matrix<-cbind(est_med_matrix,med_regACC_mad=results_plot[[mean_or_std-1]][2,])
			est_mad_matrix<-cbind(est_mad_matrix,0)
		}
		draw_lines(p_all,est_med_matrix,line_col=adj_line_col,line_type=adj_line_type,line_width=adj_line_width,main=param_name,xname='Acceptance Proportion',ylim=ylim_matrix[[mean_or_std]][param_i,],yname=yname,with_legend=with_legend,line_names=adj_lines_names,with_CI=TRUE,lines_sd_matrix=est_mad_matrix,min_max_quantile=min_max_quantile)
	}
}

rn_sub<-'half'
rn_method<-'PT'

setwd('D:\\bitbucket\\acc-examples')
rn_sub<-'half'
rn_method<-'refinedPT'
save.image('Ricker_initial_cal_167rep_N1e5_n500_half_rn_refinedPT.RData')

rn_sub<-'half'
rn_method<-'refinedmixture'
save.image('Ricker_initial_cal_167rep_N1e5_n500_half_rn_refinedmixture.RData')


setwd('D:\\bitbucket\\acc-examples')
rn_sub<-'threefifth'
rn_method<-'PT'
save.image('Ricker_initial_cal_300rep_N1e5_n500_threefifth_rn_PT.RData')

setwd('D:\\bitbucket\\acc-examples')
rn_sub<-'threefifth'
rn_method<-'refinedPT'
save.image('Ricker_initial_cal_300rep_N1e5_n500_threefifth_rn_refinedPT.RData')

setwd('D:\\bitbucket\\acc-examples')
rn_sub<-'threefifth'
rn_method<-'refinedmixture'
save.image('Ricker_initial_cal_300rep_N1e5_n500_threefifth_rn_refinedmixture.RData')


#############################################################
# Compare the density plots of r_n, refined r_n, pi_ACC output and rACC output
#############################################################


# save(parameter_true,tested_observations_all,subgroup_param_est,file='Ricker_PMC_tmp.RData')

platform<-'Win'; ncores<-6
### Experiment variables ###
test_size<-300 # Number of runs in an experiment
nobs<-500; N<-10^5
# nobs<-50 # Follows the data size in Wood(2010) where the 95% credible region only covers 85%
parameter_true<-c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)) # Follows the setting in Wood(2010)
PMC_refine<-FALSE; N_PMC<-10^4; bandwidth_percentile<-0.05
### Random seed ###
set.seed(200)
### Observed datasets simulation ###
param_est_ind<-1:3
ricker_sl<-synlik(simulator = rickerSimul,summaries = rickerStats,param = parameter_true,extraArgs = list("nObs" = nobs, "nBurn" = 50))
tested_observations_all<-simulate(ricker_sl, nsim = test_size)
duration<-0
niter<-10000; nburn<-1000
group_no<-floor(nobs^(3/5))
group_size<-floor(nobs/group_no); half_group_size<-floor(nobs/group_no/2)


# dev.new()
# par(mfrow=c(2,2))
data_test<-c(1,5,10,15); p_test<-45; param_test<-1; N_PMC<-10^4; bandwidth_percentile<-0.01
data_i<-1; group_i<-1

# for(data_i in data_test){

	tested_value<-parameter_true[param_est_ind]
	tested_observations<-tested_observations_all[data_i,]

	subgroup<-divide(tested_observations,group_no)
	subgroup_indices<-divide(1:nobs,group_no)
	tested_subgroup_param_est<-subgroup_param_est[[data_i]]
	component_sd_mean<-list(mean=matrix(0,length(subgroup),3),sd=list(0))
	dev.new()
	for(group_i in 1:length(subgroup)){
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
				simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,tested_subgroup_param_est[,param_i]))
				proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],tested_subgroup_param_est[,param_i])
			}				
			return(list(sample=simulated_parameters,densities=proposal_densities))
		}
		PMC_results<-ABC_PMC(observed_summary=subgroup_summary,simulator=PMC_simulator,initial_proposal_func=PMC_ini_proposal,prior_density=prior_density,proposal_func=proposal_func,N=N_PMC,bandwidth_percentile=bandwidth_percentile,iter_max=10,summary_normalising=TRUE)

		(component_sd_mean$sd)[[group_i]]<-sqrtmatrix(wvar(PMC_results$theta,PMC_results$weights))
		(component_sd_mean$mean)[group_i,]<-wmean(PMC_results$theta,PMC_results$weights)
	}

component_sd_mean<-subgroup_refine_est[[data_i]]

refine_mixture_density<-list(0)
mixsample_component_ind <- sample(1:group_no,size=1e4,replace=TRUE)
for(param_i in 1:3){
	refined_means<-(component_sd_mean$mean)[,param_i]
	refined_sds<-unlist(lapply(component_sd_mean$sd,function(x)x[param_i,param_i]))
	mixsample <- rnorm(n=1e4,mean=refined_means[mixsample_component_ind],sd=refined_sds[mixsample_component_ind])
	refine_mixture_density[[param_i]]<-density(mixsample)
}
dev.new()
par(mfrow=c(2,2))
for(param_i in 1:3){
	refine_mean_density<-density((component_sd_mean$mean)[,param_i])
	point_est_density<-density(subgroup_param_est[[data_i]][,param_i])
	draw_densities(list(point_est=point_est_density,refine_mean=refine_mean_density,refine_mixture=refine_mixture_density[[param_i]]),main=names(parameter_true)[param_i])
	abline(v=parameter_true[param_i],col=4,lty=2)
}

ricker_sl@extraArgs$obsData<-tested_observations
tested_summary<-c(ricker_sl@summaries(tested_observations,ricker_sl@extraArgs))
tested_summary_less<-tested_summary[-(2:3)]

simulated_parameters<-NULL
proposal_densities<-1
prior_densities<-1
for(param_i in param_est_ind){
	simulated_parameters<-cbind(simulated_parameters,simulate_from_r(N,subgroup_param_est[[data_i]][,param_i]))
	proposal_densities<-proposal_densities*density_r(simulated_parameters[,param_i],subgroup_param_est[[data_i]][,param_i])
	prior_densities<-prior_densities*exp(simulated_parameters[,param_i])
}
# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
simulated_summaries<-simulate_pseudo_summaries(simulated_parameters,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
reference_tables_ptest<-list(parameters=simulated_parameters,summaries=simulated_summaries[,-(2:3)])
weights_of_sample_ptest<-prior_densities/proposal_densities

refinemean_sample<-NULL
proposal_densities<-1
prior_densities<-1
for(param_i in param_est_ind){
	refinemean_sample<-cbind(refinemean_sample,simulate_from_r(N,(component_sd_mean$mean)[,param_i]))
	proposal_densities<-proposal_densities*density_r(refinemean_sample[,param_i],(component_sd_mean$mean)[,param_i])
	prior_densities<-prior_densities*exp(refinemean_sample[,param_i])
}
# Simulate pseudo datasets and summaries, and build up the reference table to run ABC/ACC algortihms
simulated_summaries<-simulate_pseudo_summaries(refinemean_sample,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
reference_tables_refineest<-list(parameters=refinemean_sample,summaries=simulated_summaries[,-(2:3)])
weights_of_sample_refineest<-prior_densities/proposal_densities


mixsample<-NULL
prior_densities<-1; proposal_densities_matrix<-matrix(1,group_no,N)
mixsample_component_ind <- sample(1:group_no,size=N,replace=TRUE)
for(param_i in 1:3){
	refined_means<-(component_sd_mean$mean)[,param_i]
	refined_sds<-unlist(lapply(component_sd_mean$sd,function(x)x[param_i,param_i]))
	mixsample <- cbind(mixsample,rnorm(n=N,mean=refined_means[mixsample_component_ind],sd=refined_sds[mixsample_component_ind]))
	for(group_i in 1:group_no) proposal_densities_matrix[group_i,]<-proposal_densities_matrix[group_i,]*dnorm(mixsample[,param_i],mean=mus[group_i],sd=sds[group_i])
	prior_densities<-prior_densities*exp(mixsample[,param_i])
}
simulated_summaries<-simulate_pseudo_summaries(mixsample,n=nobs,d=length(tested_summary),obsData=tested_observations) # a N*d matrix; for N=1e5, running time is 60s
proposal_densities<-colSums(proposal_densities_matrix)/group_no
reference_tables_mixest<-list(parameters=mixsample,summaries=simulated_summaries[,-(2:3)])
weights_of_sample_mixest<-prior_densities/proposal_densities

p_all<-signif2(seq(0.005,0.9,length=50)); CI_alpha<-0.9
p_test<-40; tolerance_percentile<-p_all[p_test]
xlim_all<-rbind(c(2.5,5.5),c(-6,2),c(1.8,2.9))
dev.new(); 
par(mfrow=c(3,3))	
for(test_i in 1:3){
	if(test_i==1) {reference_tables<-reference_tables_ptest; weights_of_sample<-weights_of_sample_ptest; r_name<-'point_est'}
	if(test_i==2) {reference_tables<-reference_tables_refineest; weights_of_sample<-weights_of_sample_refineest; r_name<-'refine_mean_est'}
	if(test_i==3) {reference_tables<-reference_tables_mixest; weights_of_sample<-weights_of_sample_mixest; r_name<-'mixture_est'}
	results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
	weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind]
	results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
	results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
	ABCpost_mean<-Process_outputlist(results_ISregABC,function(x) wmean(x$parameters,x$weights),rbind)
	regACC_mean<-Process_outputlist(results_regACC,function(x) colMeans(x$parameters),rbind)
	for(param_test in 1:3){
		initial_mean<-mean(reference_tables$parameters[,param_test])
		ABCpost_upper<-Process_outputlist(results_ISregABC,function(x)onesideCI_cal(x$parameters[,param_test],CI_alpha,weights=x$weights),c)
		regACC_upper<-Process_outputlist(results_regACC,function(x)onesideCI_cal(x$parameters[,param_test],CI_alpha),c)
		initial_density<-density(reference_tables$parameters[,param_test])
		ABCpost_density<-density(results_ISregABC[[1]]$parameters[,param_test],weights=weights_ISABC/sum(weights_ISABC))
		rACC_density<-density(results_regACC[[1]]$parameters[,param_test])
		densities_plot<-list(initial_density=initial_density,ABCpost_density=ABCpost_density,rACC_density=rACC_density)
		draw_densities(densities_plot,main=paste(names(parameter_true)[param_test],r_name),xlim=xlim_all[param_test,])
		abline(v=initial_mean,col=1,lty=2)
		abline(v=ABCpost_mean[param_test],col=2,lty=2)
		abline(v=regACC_mean[param_test],col=3,lty=2)
		abline(v=parameter_true[param_test],col=4,lty=2)
		abline(v=ABCpost_upper,col=2,lty=3,lwd=2)
		abline(v=regACC_upper,col=3,lty=3,lwd=2)
	}
}

}

rbind(pi_ACC_mean=signif(post_means_all[[4]][[p_test]][data_test,param_test],4),ACC_mean=signif(post_means_all[[6]][[p_test]][data_test,param_test],4))
rbind(pi_ACC=post_coverage1D_all[[1]][[4]][p_test,data_test],ACC=post_coverage1D_all[[1]][[6]][p_test,data_test])
rbind(pi_ACC=signif(CIupper_all[[1]][[4]][p_test,data_test],4),ACC=signif(CIupper_all[[1]][[6]][p_test,data_test],4))



#############################################################
# Check coverage rates that look weird
#############################################################

p_test<-40; param_test<-1
rbind(pi_ACC=signif(CIupper_all[[param_test]][[4]][p_test,],4),ACC=signif(CIupper_all[[param_test]][[6]][p_test,],4),true=signif(parameter_true[param_test]))
rbind(pi_ACC=post_coverage1D_all[[param_test]][[4]][p_test,],ACC=post_coverage1D_all[[param_test]][[6]][p_test,])
rbind(pi_ACC_mean=signif(post_means_all[[4]][[p_test]][,param_test],4),ACC_mean=signif(post_means_all[[6]][[p_test]][,param_test],4),pi_ACC_variance=signif(post_variances_all[[4]][[p_test]][,param_test],4),ACC_variance=signif(post_variances_all[[6]][[p_test]][,param_test],4))


(raw_results_all$post_coverage_all)[[3]][[6]][40,]
results_all[[6]][6,40]




data_i<-1; param_i<-1
reference_tables<-reference_tables_all[[data_i]]
weights_of_sample<-weights_of_sample_all[[data_i]]
tested_summary_less<-tested_summary_less_all[[data_i]]
results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind]
results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
par(mfrow=c(2,2))
for(param_i in 1:3){
	ABCpost_density<-density(results_ISregABC[[1]]$parameters[,param_i],weights=weights_ISABC/sum(weights_ISABC))
	draw_densities(list(output_dist=ABCpost_density),main=names(parameter_true)[param_i])
}

method_names<-c('ISregABC','regACC')
line_col<-rep(c('red','green'),4); line_type<-rep(c(2,2),4)
for(param_i in 1:3){
	dev.new()
	par(mfrow=c(2,2))
	for(p_i in c(1,5,10,15)){
		tolerance_percentile<-p_all[p_i]
		density_all<-NULL
		for(data_i in c(1,5,10,15)){
			reference_tables<-reference_tables_all[[data_i]]
			weights_of_sample<-weights_of_sample_all[[data_i]]
			tested_summary_less<-tested_summary_less_all[[data_i]]
			results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		 	weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind]
			results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
			results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
			ISpi_regACC_mean<-post_means_all[[4]][[p_i]][data_i,]
			regACC_mean<-post_means_all[[6]][[p_i]][data_i,]
			ABCpost_density<-density(results_ISregABC[[1]]$parameters[,param_i]-ISpi_regACC_mean[param_i],weights=weights_ISABC/sum(weights_ISABC))
			rACC_density<-density(results_regACC[[1]]$parameters[,param_i]-regACC_mean[param_i])
			density_all<-c(density_all,list(ABCpost_density,rACC_density))
		}
		draw_densities(density_all,color_vec=line_col,lty_vec=line_type,main=paste0(names(parameter_true)[param_i],' p=',p_all[p_i]))
		abline(v=0,col=3)
	}
}


tmp_raw$post_means_all[[method_i]][[p_i]][indices,param_i]

list(logR=(results_all$logR)[[1]],logSigma=(results_all$logSigma)[[1]],logPhi=(results_all$logPhi)[[1]])
list(logR=(results_all$logR)[[2]],logSigma=(results_all$logSigma)[[2]],logPhi=(results_all$logPhi)[[2]])
parameter_true : c(logR = 3.8, logSigma = -1.204, logPhi = 2.302) 


dev.new()
# xlim_all<-rbind(c(2.5,5.5),c(-6,2),c(1.8,2.9))
xlim=c(2.5,5)
data_test<-floor(seq(from=1, to=50,length=9)); param_test<-1
par(mfrow=c(3,3))
p_test<-3; tolerance_percentile<-p_all[p_test]
for(data_i in data_test){
	reference_tables<-reference_tables_all[[data_i]]
	weights_of_sample<-weights_of_sample_all[[data_i]]
	tested_summary_less<-tested_summary_less_all[[data_i]]
	results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
	weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind]
	results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
	results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
	ABCpost_mean<-Process_outputlist(results_ISregABC,function(x) wmean(x$parameters,x$weights),rbind)
	regACC_mean<-Process_outputlist(results_regACC,function(x) colMeans(x$parameters),rbind)
	ABCpost_rejmean<-Process_outputlist(results_ISABC,function(x) wmean(x$parameters,weights_ISABC),rbind)
	rejACC_mean<-Process_outputlist(results_ISABC,function(x) colMeans(x$parameters),rbind)
	for(param_i in param_test){
		initial_mean<-mean(reference_tables$parameters[,param_i])
		r_n_density<-density(reference_tables_all[[data_i]]$parameters[,param_i])
		ABCpost_density<-density(results_ISregABC[[1]]$parameters[,param_i],weights=weights_ISABC/sum(weights_ISABC))
		# regACC_density<-density(results_regACC[[1]]$parameters[,param_i])
		ABCpost_rejdensity<-density(results_ISABC[[1]]$parameters[,param_i],weights=weights_ISABC/sum(weights_ISABC))
		# rejACC_density<-density(results_ISABC[[1]]$parameters[,param_i])
		ABCpost_upper<-Process_outputlist(results_ISregABC,function(x)onesideCI_cal(x$parameters[,param_i],CI_alpha,weights=x$weights),c)
		regACC_upper<-Process_outputlist(results_regACC,function(x)onesideCI_cal(x$parameters[,param_i],CI_alpha),c)
		densities_plot<-list(r_n_density=r_n_density,ABCpost_density=ABCpost_density,ABCpost_rejdensity=ABCpost_rejdensity)
		# draw_densities(densities_plot,main=names(parameter_true)[param_i],xlim=xlim_all[param_i,])
		if(data_i==min(data_test)) with_legend<-TRUE else with_legend<-FALSE
		draw_densities(densities_plot,main=names(parameter_true)[param_i],xlim=xlim,with_legend=with_legend)
		abline(v=initial_mean,col=1,lty=2)
		abline(v=ABCpost_mean[param_i],col=2,lty=2)
		abline(v=ABCpost_rejmean[param_i],col=3,lty=2)
		# abline(v=regACC_mean[param_i],col=3,lty=2)
		abline(v=parameter_true[param_i],col=4,lty=2)			
		# abline(v=ABCpost_rejmean[param_i],col=2,lty=2,lwd=2)
		# abline(v=rejACC_mean[param_i],col=3,lty=2,lwd=2)
		# abline(v=ABCpost_upper,col=2,lty=3,lwd=2)
		# abline(v=regACC_upper,col=3,lty=3,lwd=2)
	}
}


reference_tables_all[[1]]$parameters[1:10,]
summary(weights_of_sample_all[[1]])
tested_summary_less_all[[1]]

#############################################################
# Check whether the covariance between each 
# pair of parameters estimated by ABC and rACC are similar. 
# 
# Result: They are very similar. 
#############################################################
library(ggtern)
library(MASS)
method_names<-c('ISregABC','regACC'); pairs<-list(c(1,2),c(1,3))
line_col<-c('black','green')
for(pair_i in 1:2){
	# quartz()
	param_names<-names(parameter_true)[pairs[[pair_i]]]
	dev.new()
	par(mfcol=c(3,3))
	for(p_i in c(5,15,30)){
		tolerance_percentile<-p_all[p_i]
		density_all<-NULL
		for(data_i in c(10,15,30)){
			setting_name<-paste0('p=',p_all[p_i],', data_i=',data_i)
			reference_tables<-reference_tables_all[[data_i]]
			weights_of_sample<-weights_of_sample_all[[data_i]]
			tested_summary_less<-tested_summary_less_all[[data_i]]
			results_ISABC<-rejABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
		 	weights_ISABC<-weights_of_sample[results_ISABC[[1]]$acceptance_ind]
			results_ISregABC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,weights_of_sample=weights_of_sample) 
			results_regACC<-regABC(reference_tables,t(tested_summary_less),tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)

			ISregABC_params<-results_ISregABC[[1]]$parameters[,pairs[[pair_i]]]
			regACC_params<-results_regACC[[1]]$parameters[,pairs[[pair_i]]]
			ISregABC_mean<-post_means_all[[2]][[p_i]][data_i,pairs[[pair_i]]]
			regACC_mean<-post_means_all[[6]][[p_i]][data_i,pairs[[pair_i]]]
			ABCpost_density<-kde2d.weighted(ISregABC_params[,1],ISregABC_params[,2],w=weights_ISABC/sum(weights_ISABC))
			rACC_density<-kde2d(regACC_params[,1],regACC_params[,2])
			contour(ABCpost_density,xlab=param_names[1],ylab=param_names[2],main=setting_name,col=line_col[1])
			abline(v=ISregABC_mean[1],col=line_col[1],lty=2)
			abline(h=ISregABC_mean[2],col=line_col[1],lty=2)
			contour(rACC_density,xlab=param_names[1],ylab=param_names[2],col=line_col[2],add=TRUE)
			abline(v=regACC_mean[1],col=line_col[2],lty=2)
			abline(h=regACC_mean[2],col=line_col[2],lty=2)
		}
	}
}
