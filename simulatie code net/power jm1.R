library(mvtnorm)
library(nlme)

simulation_function_power<-function(n=0,fixed_intercept,fixed_slope,residual_variance, cov_matrix, alpha, baseline, scale, followups, nsim, correction=0){
  bias_estimates_lmm<-data.frame(bias_b0=rep(0,nsim), bias_b1=rep(0,nsim), bias_varb0=rep(0, nsim), bias_varb1=rep(0, nsim), bias_cov=rep(0, nsim),
                                 bias_res = rep(0, nsim))
  counter<-0
  total_iterations <- 0
  followups[1]<-correction
  p_values<-c()
  #calculate samole size needed for 0.8 power if alpha is 0.05
  if(n==0){
  followupterm<-sum((followups-mean(followups))^2)
  n_arm<-2*(1.96+0.842)^2*(cov_matrix[2,2]+residual_variance/followupterm)/(fixed_slope*0.3)^2  #sample size for 0.8 power and 0.05 alpha
  n<-ceiling(n_arm)*2} #round n/arm upwards
  trt_effect<-c()
  while(counter != nsim){
    ##simulated random effects
    random_intercepts_slopes<-rmvnorm(n, c(0,0), cov_matrix)
    
    ########simulate survival times and random dropout times
    #survival times
    tij<-c()
    #random dropout times
    rij<-c()
    #control group
    for(i in 1:(n/2)){
      #specification hazard function
      hazard<-function(t,base = baseline, rho = scale, intercept = (fixed_intercept + random_intercepts_slopes[i,1]), slope = (fixed_slope + random_intercepts_slopes[i,2])){
        rho*t^(rho-1)*exp(base + alpha*(intercept+t*slope))}
      s<-runif(1)# sample survival probability
      #S(t) minus random number on [0,1] as a function
      #this is inverse transform sampling
      
      prob_minus_rand<-Vectorize(function(x){exp(-integrate(hazard,correction,x)$value)-s})
      #find zero of the function to solve for t
      tij[i]<-tryCatch(uniroot(prob_minus_rand, c(min(followups), max(followups)))$root, error = function(e){max(followups)})
      
      lambda<- -log(0.9)/12 #lambda for survival probability (wrt random dropout) of 90% at 12 months
      s_dropout <- runif(1)#sample survival probability for random dropout
      #random dropout time
      rij[i]<-log(s_dropout)/-lambda #solve for t
    }
    #treatment group
    for(j in (n/2+1):n){
      
      hazard<-function(t,base = baseline, rho = scale, intercept = (fixed_intercept + random_intercepts_slopes[j,1]), slope = (fixed_slope*0.7 + random_intercepts_slopes[j,2])){
        rho*t^(rho-1)*exp(base + alpha*(intercept+t*slope))}
      s<-runif(1)
      
      prob_minus_rand<-Vectorize(function(x){exp(-integrate(hazard,correction,x)$value)-s})
      tij[j]<-tryCatch(uniroot(prob_minus_rand, c(min(followups), max(followups)))$root, error = function(e){max(followups)})
      
      lambda<- -log(0.9)/12
      s_dropout <- runif(1)
      rij[j]<-log(s_dropout)/-lambda
    }
    
    ######create the dataframe with fixed effects, random effects, survival times, random dropout times.
    ######patients that decease or dropout randomly not dropped yet
    longitudinal_data_long<-data.frame(id = rep(1:n, each = length(followups)), time = rep(followups, n), 
                                       random_int = rep(random_intercepts_slopes[,1], each = length(followups)),
                                       rand_sl = rep(random_intercepts_slopes[,2], each = length(followups)),
                                       fixed_int = rep(fixed_intercept, n*length(followups)), 
                                       fixed_sl = c(rep(fixed_slope, n*length(followups)/2), rep(fixed_slope*0.7, n*length(followups)/2)), 
                                       meas_error = rnorm(n*length(followups), 0, sqrt(residual_variance)), 
                                       random_dropout = rep(rij, each = length(followups)), survival_time = rep(tij, each = length(followups)),
                                       treatment = c(rep(rep(0,length(followups)), n/2), rep(rep(1,length(followups)), n/2)))
    longitudinal_data_long$observed_scores<-longitudinal_data_long$fixed_int+longitudinal_data_long$random_int+
      longitudinal_data_long$time*(longitudinal_data_long$fixed_sl+longitudinal_data_long$rand_sl)+longitudinal_data_long$meas_error
    #event if someone's survival time is smaller than random dropout time or time study
    longitudinal_data_long$event<-ifelse(longitudinal_data_long$survival_time < max(followups) & longitudinal_data_long$survival_time < 
                                           longitudinal_data_long$random_dropout, 1, 0)
    #drop deceased
    #drop measurements whose survival time is lower than measurement time of patients that die; not randomly dropout (a month before death)
    longitudinal_deceased_gone<-longitudinal_data_long[!(longitudinal_data_long$time > longitudinal_data_long$survival_time-1 & 
                                                           longitudinal_data_long$event==1),]
    #drop measurements whose random dropout is lower than measurement time of paients that dropout randomly
    longitudinal_deceased_random_gone <- longitudinal_deceased_gone[!(longitudinal_deceased_gone$time > longitudinal_deceased_gone$random_dropout & 
                                                                        longitudinal_deceased_gone$event==0),]
    #fit LMMs
    lmm<-try(lme(observed_scores ~ time, random = ~ time|id, data = longitudinal_deceased_random_gone, method = 'ML', 
                 control = lmeControl(opt = "optim")))
    lmm_trt<-try(lme(observed_scores ~ time + time:treatment, random = ~ time|id, data = longitudinal_deceased_random_gone, method = 'ML', 
                     control = lmeControl(opt = "optim")))
    
    
    if(is.list(lmm) & is.list(lmm_trt)){
      p_values<-c(p_values, anova(lmm,lmm_trt)[[9]][2])
      trt_effect<- c(trt_effect, summary(lmm_trt)$tTable[3,1])
      counter <- counter + 1
    }
    
    
    total_iterations <- total_iterations + 1
  }
  est_power<-sum(p_values<=0.05)/nsim#power estimate
  n_failed_iterations <- total_iterations - nsim
  return(list(n_failed_iterations, est_power, p_values, n, trt_effect))
}

#data generating values

#slow progressors
b0_slow <- 42
b1_slow<- -0.74
alpha_slow <- -0.14 
intercept_slow <- -1.9
scal_slow <- 1.4
cov_b0_b1_slow<- -0.01
var_b0_slow<- 14
var_b1_slow<- 0.4
var_e_slow<- 3.1
covar_mat_slow<-matrix(c(var_b0_slow,cov_b0_b1_slow,cov_b0_b1_slow,var_b1_slow), 2, 2)
#intermediate progressors
b0_interm<-37
b1_interm<- -1.1
alpha_interm <- -0.11 
intercept_interm <- -3.3
scal_interm <- 1.7
cov_b0_b1_interm<- -0.13
var_b0_interm<- 19
var_b1_interm<- 0.52
var_e_interm<- 4.2
covar_mat_interm<-matrix(c(var_b0_interm,cov_b0_b1_interm,cov_b0_b1_interm,var_b1_interm), 2, 2)
#fast progressors
b0_fast<- 33
b1_fast <- -1.5
alpha_fast = -0.14 
intercept_fast = -0.6
scal_fast <- 0.92
cov_b0_b1_fast<- -1.1
var_b0_fast<- 32
var_b1_fast<-0.92
var_e_fast<-6.2
covar_mat_fast<-matrix(c(var_b0_fast,cov_b0_b1_fast,cov_b0_b1_fast,var_b1_fast), 2, 2)
#all progressors
b0_all<-38
b1_all<- -1
alpha_all = -0.12 
intercept_all = -2.2
scal_all = 1.4
cov_b0_b1_all<-0.66
var_b0_all<-29
var_b1_all<-0.61
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
sim_one_six_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, monthly_halfyear, 10000)
set.seed(199)
sim_three_six_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, threemonthly_halfyear, 10000)
set.seed(199)
sim_one_twelve_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, monthly_year, 10000)
set.seed(199)
sim_three_twelve_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, threemonthly_year, 10000)
set.seed(199)
sim_one_eighteen_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, monthly_oneandhalfyear, 10000)
set.seed(199)
sim_three_eighteen_slow_power<-simulation_function_power(0,b0_slow, b1_slow, var_e_slow, covar_mat_slow, alpha_slow, intercept_slow, scal_slow, threemonthly_oneandhalfyear, 10000)
#intermediate
set.seed(199)
sim_one_six_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, monthly_halfyear, 10000)
set.seed(199)
sim_three_six_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, threemonthly_halfyear, 10000)
set.seed(199)
sim_one_twelve_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, monthly_year, 10000)
set.seed(199)
sim_three_twelve_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, threemonthly_year, 10000)
set.seed(199)
sim_one_eighteen_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, monthly_oneandhalfyear, 10000)
set.seed(199)
sim_three_eighteen_interm_power<-simulation_function_power(0,b0_interm, b1_interm, var_e_interm, covar_mat_interm, alpha_interm, intercept_interm, scal_interm, threemonthly_oneandhalfyear, 10000)
#fast
set.seed(199)
sim_one_six_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, monthly_halfyear, 10000, 0.01)#0.01, because if t=0, division by 0 in hazard
set.seed(199)
sim_three_six_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, threemonthly_halfyear, 10000, 0.01)
set.seed(199)
sim_one_twelve_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, monthly_year, 10000, 0.01)
set.seed(199)
sim_three_twelve_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, threemonthly_year, 10000, 0.01)
set.seed(199)
sim_one_eighteen_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, monthly_oneandhalfyear, 10000, 0.01)
set.seed(199)
sim_three_eighteen_fast_power<-simulation_function_power(0,b0_fast, b1_fast, var_e_fast, covar_mat_fast, alpha_fast, intercept_fast, scal_fast, threemonthly_oneandhalfyear, 10000, 0.01)
#all
set.seed(199)
sim_one_six_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_halfyear, 10000)
set.seed(199)
sim_three_six_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_halfyear, 10000)
set.seed(199)
sim_one_twelve_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_year, 10000)
set.seed(199)
sim_three_twelve_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_year, 10000)
set.seed(199)
sim_one_eighteen_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, monthly_oneandhalfyear,10000)
set.seed(199)
sim_three_eighteen_all_power<-simulation_function_power(0,b0_all, b1_all, var_e_all, covar_mat_all, alpha_all, intercept_all, scal_all, threemonthly_oneandhalfyear, 10000)
