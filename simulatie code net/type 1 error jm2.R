############# simulate random effects
library(mvtnorm)
library(nlme)
library(survival)

simulation_function_histo_slope<-function(n,ncontrol,fixed_intercept,fixed_slope,residual_variance, cov_matrix, alpha, baseline, scale, followups, followupscontrol, nsim, correction=0){
  bias_estimates_lmm<-data.frame(bias_b0=rep(0,nsim), bias_b1=rep(0,nsim), bias_varb0=rep(0, nsim), bias_varb1=rep(0, nsim), bias_cov=rep(0, nsim),
                                 bias_res = rep(0, nsim))
  
  bias_estimates_lmm_trt<-data.frame(bias_b0=rep(0,nsim), bias_b1=rep(0,nsim), bias_varb0=rep(0, nsim), bias_varb1=rep(0, nsim), bias_cov=rep(0, nsim),
                                     bias_res = rep(0, nsim), trt_effect = rep(0,nsim), wald_p_value = rep(0, nsim), std_error = rep(0, nsim))
  counter<-0
  total_iterations <- 0
  followups[1]<-correction
  p_values<-c()
  
  while(counter != nsim){
    #sample random effects
    random_intercepts_slopes<-rmvnorm(n, c(0,0), cov_matrix)
    
    ########simulate survival times and random dropout times
    #survival times
    tij<-c()
    #random dropout times
    rij<-c()
    #control group
    for(i in 1:ncontrol){
      #specification of hazard
      hazard<-function(t,base = baseline, rho = scale, intercept = (fixed_intercept + random_intercepts_slopes[i,1]), slope = (fixed_slope + random_intercepts_slopes[i,2])){
        rho*t^(rho-1)*exp(base + alpha*slope)}
      s<-runif(1)
      #S(t) minus random number on [0,1] as a function
      #this is inverse transform sampling
      prob_minus_rand<-Vectorize(function(x){exp(-integrate(hazard,correction,x)$value)-s})
      #find zero of the function to solve for t
      tij[i]<-tryCatch(uniroot(prob_minus_rand, c(min(followupscontrol), max(followupscontrol)))$root, error = function(e){max(followupscontrol)})
      
      lambda<- -log(0.9)/12 #lambda for survival probability (wrt random dropout) of 90% at 12 months
      s_dropout <- runif(1)#sample survival probability for random dropout
      #solve for t
      rij[i]<-log(s_dropout)/-lambda
    }
    #treatment group
    for(i in (ncontrol+1):n){
      
      hazard<-function(t,base = baseline, rho = scale, intercept = (fixed_intercept + random_intercepts_slopes[i,1]), slope = (fixed_slope + random_intercepts_slopes[i,2])){
        rho*t^(rho-1)*exp(base + alpha*slope)}
      s<-runif(1)
      
      prob_minus_rand<-Vectorize(function(x){exp(-integrate(hazard,correction,x)$value)-s})
      
      tij[i]<-tryCatch(uniroot(prob_minus_rand, c(min(followups), max(followups)))$root, error = function(e){max(followups)})
      
      lambda<- -log(0.9)/12
      s_dropout <- runif(1)
      rij[i]<-log(s_dropout)/-lambda
    }
    
    #Control group
    ######create the dataframe with fixed effects, random effects, survival times, random dropout times.
    ######patients that decease or dropout randomly not dropped yet
    longitudinal_data_long_control<-data.frame(id = rep(1:ncontrol, each = length(followupscontrol)), time = rep(followupscontrol, ncontrol), 
                                               random_int = rep(random_intercepts_slopes[,1][1:ncontrol], each = length(followupscontrol)),
                                               rand_sl = rep(random_intercepts_slopes[,2][1:ncontrol], each = length(followupscontrol)),
                                               fixed_int = rep(fixed_intercept, ncontrol*length(followupscontrol)), 
                                               fixed_sl = rep(fixed_slope, ncontrol*length(followupscontrol)), 
                                               meas_error = rnorm(ncontrol*length(followupscontrol), 0, sqrt(residual_variance)), 
                                               random_dropout = rep(rij[1:ncontrol], each = length(followupscontrol)), survival_time = rep(tij[1:ncontrol], each = length(followupscontrol)),
                                               treatment = rep(0, ncontrol*length(followupscontrol))) #treatment value for each followup and then for each person in contol
    
    longitudinal_data_long_control$observed_scores<-longitudinal_data_long_control$fixed_int+longitudinal_data_long_control$random_int+
      longitudinal_data_long_control$time*(longitudinal_data_long_control$fixed_sl+longitudinal_data_long_control$rand_sl)+longitudinal_data_long_control$meas_error
    #event if someone's survival time is smaller than random dropout time or time study
    longitudinal_data_long_control$event<-ifelse(longitudinal_data_long_control$survival_time < max(followupscontrol) & longitudinal_data_long_control$survival_time < 
                                                   longitudinal_data_long_control$random_dropout, 1, 0)
    #drop deceased
    #drop measurements whose survival time is lower than measurement time of patients that die; not randomly dropout (a month before death)
    longitudinal_deceased_gone_control<-longitudinal_data_long_control[!(longitudinal_data_long_control$time > longitudinal_data_long_control$survival_time-1 & 
                                                                           longitudinal_data_long_control$event==1),]
    #drop measurements whose random dropout is lower than measurement time of paients that dropout randomly
    longitudinal_deceased_random_gone_control <- longitudinal_deceased_gone_control[!(longitudinal_deceased_gone_control$time > longitudinal_deceased_gone_control$random_dropout & 
                                                                                        longitudinal_deceased_gone_control$event==0),]
    
    longitudinal_data_long_treatment<-data.frame(id = rep((ncontrol+1):n, each = length(followups)), time = rep(followups, (n-ncontrol)), 
                                                 random_int = rep(random_intercepts_slopes[,1][(ncontrol+1):n], each = length(followups)),
                                                 rand_sl = rep(random_intercepts_slopes[,2][(ncontrol+1):n], each = length(followups)),
                                                 fixed_int = rep(fixed_intercept, (n-ncontrol)*length(followups)), 
                                                 fixed_sl = rep(fixed_slope, (n-ncontrol)*length(followups)), 
                                                 meas_error = rnorm((n-ncontrol)*length(followups), 0, sqrt(residual_variance)), 
                                                 random_dropout = rep(rij[(ncontrol+1):n], each = length(followups)), survival_time = rep(tij[(ncontrol+1):n], each = length(followups)),
                                                 treatment = rep(1, (n-ncontrol)*length(followups))) #treatment value for each followup and then for each person in contol
    
    longitudinal_data_long_treatment$observed_scores<-longitudinal_data_long_treatment$fixed_int+longitudinal_data_long_treatment$random_int+
      longitudinal_data_long_treatment$time*(longitudinal_data_long_treatment$fixed_sl+longitudinal_data_long_treatment$rand_sl)+longitudinal_data_long_treatment$meas_error
    #event if someone's survival time is smaller than random dropout time or time study
    longitudinal_data_long_treatment$event<-ifelse(longitudinal_data_long_treatment$survival_time < max(followups) & longitudinal_data_long_treatment$survival_time < 
                                                     longitudinal_data_long_treatment$random_dropout, 1, 0)
    #drop deceased
    #drop measurements whose survival time is lower than measurement time of patients that die; not randomly dropout (a month before death)
    longitudinal_deceased_gone_treatment<-longitudinal_data_long_treatment[!(longitudinal_data_long_treatment$time > longitudinal_data_long_treatment$survival_time-1 & 
                                                                               longitudinal_data_long_treatment$event==1),]
    #drop measurements whose random dropout is lower than measurement time of paients that dropout randomly
    longitudinal_deceased_random_gone_treatment <- longitudinal_deceased_gone_treatment[!(longitudinal_deceased_gone_treatment$time > longitudinal_deceased_gone_treatment$random_dropout & 
                                                                                            longitudinal_deceased_gone_treatment$event==0),]
    
    
    longitudinal_deceased_random_gone<-rbind(longitudinal_deceased_random_gone_control, longitudinal_deceased_random_gone_treatment)
    
    #fit LMMs
    lmm<-try(lme(observed_scores ~ time, random = ~ time|id, data = longitudinal_deceased_random_gone, method = 'ML', 
                 control = lmeControl(opt = "optim")))
    lmm_trt<-try(lme(observed_scores ~ time + time:treatment, random = ~ time|id, data = longitudinal_deceased_random_gone, method = 'ML', 
                     control = lmeControl(opt = "optim")))
    
    
    if(is.list(lmm) & is.list(lmm_trt)){
      p_values<-c(p_values, anova(lmm,lmm_trt)[[9]][2])
      
      
      bias_estimates_lmm[(counter+1), 3]<-getVarCov(lmm)[1,1]-cov_matrix[1,1]
      bias_estimates_lmm[(counter+1), 4]<-getVarCov(lmm)[2,2]-cov_matrix[2,2]
      bias_estimates_lmm[(counter+1), 5]<-getVarCov(lmm)[1,2]-cov_matrix[2,1]
      bias_estimates_lmm[(counter+1), 6]<-as.numeric(VarCorr(lmm)[3,1]) - residual_variance
      bias_estimates_lmm[(counter+1), 1]<-lmm$coefficients$fixed[1]-fixed_intercept
      bias_estimates_lmm[(counter+1), 2]<-lmm$coefficients$fixed[2]-fixed_slope
      
      bias_estimates_lmm_trt[(counter+1), 3]<-getVarCov(lmm_trt)[1,1]-cov_matrix[1,1]
      bias_estimates_lmm_trt[(counter+1), 4]<-getVarCov(lmm_trt)[2,2]-cov_matrix[2,2]
      bias_estimates_lmm_trt[(counter+1), 5]<-getVarCov(lmm_trt)[1,2]-cov_matrix[2,1]
      bias_estimates_lmm_trt[(counter+1), 6]<-as.numeric(VarCorr(lmm_trt)[3,1]) - residual_variance
      bias_estimates_lmm_trt[(counter+1), 1]<-lmm_trt$coefficients$fixed[1]-fixed_intercept
      bias_estimates_lmm_trt[(counter+1), 2]<-lmm_trt$coefficients$fixed[2]-fixed_slope
      
      bias_estimates_lmm_trt[(counter+1), 7]<-summary(lmm_trt)$tTable[3,1]
      bias_estimates_lmm_trt[(counter+1), 8]<-summary(lmm_trt)$tTable[3,5]
      bias_estimates_lmm_trt[(counter+1), 9]<-summary(lmm_trt)$tTable[3,2]
      
      
      counter <- counter + 1
    }
    
    total_iterations <- total_iterations + 1
  }
  est_type1_error<-sum(p_values<=0.05)/nsim
  n_failed_iterations <- total_iterations - nsim
  return(list(n_failed_iterations, est_type1_error, p_values, bias_estimates_lmm_trt, bias_estimates_lmm, longitudinal_deceased_random_gone))
}

#data generating values

#all progressors
b0_all<-38
b1_all<- -1
alpha_all = -1.6
intercept_all = -12
scal_all = 3.1
cov_b0_b1_all<-0.73
var_b0_all<-29
var_b1_all<-0.62
var_e_all<-4
covar_mat_all<-matrix(c(var_b0_all,cov_b0_b1_all,cov_b0_b1_all,var_b1_all), 2, 2)

#follow-up schemes
monthly_halfyear<-c(0,1,2,3,4,5,6)
threemonthly_halfyear<-c(0, 3, 6)

monthly_year<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
threemonthly_year<-c(0, 3, 6, 9, 12)

monthly_oneandhalfyear<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
threemonthly_oneandhalfyear<-c(0, 3, 6, 9, 12, 15, 18)

#simulations
set.seed(199)
sim_diff_size_slope<-simulation_function_histo_slope(200,150, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_year,monthly_year, 100000)
set.seed(199)
sim_diff_duration_slope<-simulation_function_histo_slope(200,100, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_year,monthly_oneandhalfyear, 100000)
set.seed(199)
sim_diff_interval_slope<-simulation_function_histo_slope(200,100, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year,monthly_year, 100000)
set.seed(199)
sim_diff_size_duration_slope<-simulation_function_histo_slope(200,150, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_year,monthly_oneandhalfyear, 100000)
set.seed(199)
sim_diff_size_interval_slope<-simulation_function_histo_slope(200,150, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year,monthly_year, 100000)
set.seed(199)
sim_diff_duration_interval_slope<-simulation_function_histo_slope(200,100, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year,monthly_oneandhalfyear, 100000)
set.seed(199)
sim_diff_duration_interval_size_slope<-simulation_function_histo_slope(200,150, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year,monthly_oneandhalfyear, 100000)
set.seed(199)
sim_same_duration_interval_size_slope<-simulation_function_histo_slope(200,100, b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year,threemonthly_year, 100000)
