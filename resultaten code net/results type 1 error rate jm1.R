lijst<-list(sim_same_duration_interval_size,sim_diff_size, sim_diff_duration,sim_diff_interval, sim_diff_duration_interval,  
              sim_diff_size_duration, sim_diff_size_interval, sim_diff_duration_interval_size)
treat_effects<-c()
errors<-c()
wald_p<-c()
type1<-c()
for(i in 1:8){
  treat_effects[i] <- mean(lijst[i][[1]][[4]]$trt_effect)
  errors[i] <- mean(lijst[i][[1]][[4]]$std_error)
  wald_p[i] <- mean(lijst[i][[1]][[4]]$wald_p_value)
  type1[i] <- lijst[i][[1]][[2]]
}

