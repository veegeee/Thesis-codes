##first load results into work space
lijst<-list(sim_one_six_slow, sim_three_six_slow, sim_one_twelve_slow, sim_three_twelve_slow, 
            sim_one_eighteen_slow, sim_three_eighteen_slow, sim_one_six_interm, sim_three_six_interm,
            sim_one_twelve_interm, sim_three_twelve_interm, sim_one_eighteen_interm, 
            sim_three_eighteen_interm, sim_one_six_fast, sim_three_six_fast, sim_one_twelve_fast,
            sim_three_twelve_fast, sim_one_eighteen_fast, sim_three_eighteen_fast, sim_one_six_all,
            sim_three_six_all, sim_one_twelve_all, sim_three_twelve_all, sim_one_eighteen_all,
            sim_three_eighteen_all, 
            sim_one_six_slow_small, sim_three_six_slow_small, sim_one_twelve_slow_small, 
            sim_three_twelve_slow_small, sim_one_eighteen_slow_small, sim_three_eighteen_slow_small, sim_one_six_interm_small, sim_three_six_interm_small,
            sim_one_twelve_interm_small, sim_three_twelve_interm_small, sim_one_eighteen_interm_small, 
            sim_three_eighteen_interm_small, sim_one_six_fast_small, sim_three_six_fast_small, sim_one_twelve_fast_small,
            sim_three_twelve_fast_small, sim_one_eighteen_fast_small, sim_three_eighteen_fast_small, sim_one_six_all_small,
            sim_three_six_all_small, sim_one_twelve_all_small, sim_three_twelve_all_small, sim_one_eighteen_all_small,
            sim_three_eighteen_all_small)

interval<-as.character(rep(rep(c(1,1,1,1,1,1,3,3,3,3,3,3), 12), 2)) #rep(vector, 2) want voor 2 sample sizen
duration<- as.factor(rep(rep(c(rep(6,12), rep(12,12), rep(18,12)),4), 2))
type<-rep(c(rep('slow', 36), rep('moderate',36), rep('fast',36), rep('all', 36)),2)
var_comp<-rep(rep(c('varb0', 'varb1', 'cov', 'b1','b2', 'resvar'), 24), 2) ##180
true<-rep(c(rep(c(14, 0.4, -0.01, -0.74,(0.3*0.74), 3.1), 6), rep(c(19, 0.52, -0.13, -1.1,(0.3*1.1), 4.2), 6), rep(c(32, 0.92, -1.1, -1.5,(0.3*1.5), 6.2), 6), rep(c(29, 0.61, 0.66, -1,(0.3*1), 4), 6)), 2)
sample_size<-as.character(c(rep(240, 144), rep(80, 144)))  
group<-as.character(c(rep(c(240,240,240,240,240,240,241,241,241,241,241,241), 12), rep(c(80,80,80,80,80,80,81,81,81,81,81,81), 12)))
results<-data.frame(Bias = rep(0, 288), lower = rep(0,288), upper = rep(0,288), Interval = interval, Duration = duration, 
                    type = type, var_comp = var_comp, true = true, Sample_size = sample_size, group=group, sim_se=rep(0, 288))
#the row in which result is stored
roww<-1
#ga door lijst met resultaatlijsten
for(i in lijst){
  #ga door kolommen van dataframe in lijst
  for(j in c(3,4,5,2,7,6)){
    #i[1] is accesses the dataframe
    results[roww,1] <- mean(i[[2]][,j])
    results[roww,2] <-as.numeric(quantile(i[[2]][,j], c(0.025, 0.975))[1])
    results[roww,3] <-as.numeric(quantile(i[[2]][,j], c(0.025, 0.975))[2])
    results[roww, 11] <- sd(i[[2]][,j])/100
    roww<-roww+1
    
  }
}
#for tables
results_table<-results[results$var_comp=='varb1' | results$var_comp=='b1' | results$var_comp=='b2' | results$var_comp=='resvar',]
results_table$Bias<-round(results_table$Bias,4)
results_table$sim_se<-round(results_table$sim_se,4)


results$observed<-results$true + results$Bias
results$percentage<-results$Bias / results$true * 100

results_slow<-results[results$type=='slow',]
results_inter<-results[results$type=='moderate',]
results_fast<-results[results$type=='fast',]
results_all<-results[results$type=='all',]


windowsFonts(tnr = windowsFont('Helvetica'))

library(ggplot2)
#slow progressors
results_slow_cov<-results_slow[results_slow$var_comp=='cov',]
p1<-ggplot(results_slow_cov, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point() + labs(caption = 'Cov = -0.01', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top',  )+ scale_linetype_discrete(name='Sample size')
results_slow_varb0<-results_slow[results_slow$var_comp=='varb0',]
p2<-ggplot(results_slow_varb0, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Var b0 = 14', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_slow_varb1<-results_slow[results_slow$var_comp=='varb1',]
p3<-ggplot(results_slow_varb1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = expression(sigma[b1]^2==0.40), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.05,0.001)+ scale_linetype_discrete(name='Sample size')

#intermediate progressors
results_inter_cov<-results_inter[results_inter$var_comp=='cov',]
p4<-ggplot(results_inter_cov, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Cov = -0.13', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_inter_varb0<-results_inter[results_inter$var_comp=='varb0',]
p5<-ggplot(results_inter_varb0, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Var b0 = 19', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_inter_varb1<-results_inter[results_inter$var_comp=='varb1',]
p6<-ggplot(results_inter_varb1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = expression(sigma[b1]^2==0.52), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top', text = element_text(family = 'tnr') )+ylim(-0.05,0.001)+ scale_linetype_discrete(name='Sample size')

#fast progressors
results_fast_cov<-results_fast[results_fast$var_comp=='cov',]
p7<-ggplot(results_fast_cov, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Cov = -1.1', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_fast_varb0<-results_fast[results_fast$var_comp=='varb0',]
p8<-ggplot(results_fast_varb0, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Var b0 = 32', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_fast_varb1<-results_fast[results_fast$var_comp=='varb1',]
p9<-ggplot(results_fast_varb1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[b1]^2==0.92), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.05,0.001)+ scale_linetype_discrete(name='Sample size')

#all progressors
results_all_cov<-results_all[results_all$var_comp=='cov',]
p10<-ggplot(results_all_cov, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Cov = 0.66', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_all_varb0<-results_all[results_all$var_comp=='varb0',]
p11<-ggplot(results_all_varb0, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption = 'Var b0 = 29', x='Duration')+theme(
  plot.caption = element_text(hjust = 0), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' )
results_all_varb1<-results_all[results_all$var_comp=='varb1',]
p12<-ggplot(results_all_varb1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[b1]^2==0.61), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.05,0.001)+ scale_linetype_discrete(name='Sample size')

#slow
results_slow_b1<-results_slow[results_slow$var_comp=='b1',]
p13<-ggplot(results_slow_b1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[1]==-0.74), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.001,0.07)+ scale_linetype_discrete(name='Sample size')
results_slow_res<-results_slow[results_slow$var_comp=='resvar',]
p14<-ggplot(results_slow_res, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[e]^2==3.1), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.012, 0.045)+ scale_linetype_discrete(name='Sample size')

results_slow_b2<-results_slow[results_slow$var_comp=='b2',]
p15<-ggplot(results_slow_b2, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[2]==0.22), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.018, 0.003)+ scale_linetype_discrete(name='Sample size')

#mod
results_inter_b1<-results_inter[results_inter$var_comp=='b1',]
p16<-ggplot(results_inter_b1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[1]==-1.1), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.001,0.07)+ scale_linetype_discrete(name='Sample size')
results_inter_res<-results_inter[results_inter$var_comp=='resvar',]
p17<-ggplot(results_inter_res, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[e]^2==4.2), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.012, 0.045)+ scale_linetype_discrete(name='Sample size')

results_inter_b2<-results_inter[results_inter$var_comp=='b2',]
p18<-ggplot(results_inter_b2, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[2]==0.33), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.018, 0.003)+ scale_linetype_discrete(name='Sample size')


#fast
results_fast_b1<-results_fast[results_fast$var_comp=='b1',]
p19<-ggplot(results_fast_b1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[1]==-1.5), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.001,0.07)+ scale_linetype_discrete(name='Sample size')
results_fast_res<-results_fast[results_fast$var_comp=='resvar',]
p20<-ggplot(results_fast_res, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[e]^2==6.2), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.012, 0.045)+ scale_linetype_discrete(name='Sample size')

results_fast_b2<-results_fast[results_fast$var_comp=='b2',]
p21<-ggplot(results_fast_b2, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[2]==0.45), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top', text = element_text(family = 'tnr') )+ylim(-0.018, 0.003)+ scale_linetype_discrete(name='Sample size')


#all
results_all_b1<-results_all[results_all$var_comp=='b1',]
p22<-ggplot(results_all_b1, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[1]==-1.0), x='Duration')+theme(
  plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.001,0.07)+ scale_linetype_discrete(name='Sample size')
results_all_res<-results_all[results_all$var_comp=='resvar',]
p23<-ggplot(results_all_res, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(sigma[e]^2==4.0), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.012, 0.045)+ scale_linetype_discrete(name='Sample size')

results_all_b2<-results_all[results_all$var_comp=='b2',]
p24<-ggplot(results_all_b2, aes(Duration, Bias, color=Interval, linetype=Sample_size, group = group)) + geom_line() + geom_point()+ labs(caption =expression(beta[2]==0.3), x='Duration')+theme(
  plot.caption = element_text(hjust = 0,size=9), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9), legend.title =element_text(size=9), legend.text = element_text(size=8),legend.position = 'top' , text = element_text(family = 'tnr'))+ylim(-0.018, 0.003)+ scale_linetype_discrete(name='Sample size')


p_lab0<-ggplot()+theme_void()
p_lab3<-ggplot()+theme_void()+annotate('text',x=1, y=1, label = expression(sigma[b1]^2),size=4)
p_lab4<-ggplot()+theme_void()+annotate('text',x=1, y=1, label = expression(beta[1]), size=4)
p_lab5<-ggplot()+theme_void()+annotate('text',x=1, y=1, label = expression(beta[2]),size = 4)

p_lab6<-ggplot()+theme_void()+annotate('text',x=1, y=1, label = expression(sigma[e]^2),size =4)

p_lab7<-ggplot()+theme_void()+geom_text(aes(0,0,label='Slow'), size=4)
p_lab8<-ggplot()+theme_void()+geom_text(aes(0,0,label='Moderate'), size=4)
p_lab9<-ggplot()+theme_void()+geom_text(aes(0,0,label='Fast'), size=4)
p_lab10<-ggplot()+theme_void()+geom_text(aes(0,0,label='All'), size=4)

#library(gridExtra)
library(ggpubr) #ggarrange

ho<-ggarrange(p_lab0,p_lab7, p_lab8, p_lab9, p_lab10, p_lab4, p13, p16, p19, p22, p_lab5, p15,p18,p21,p24,
              p_lab3, p3,p6,p9,p12,p_lab6,p14,p17,p20,p23, ncol=5, nrow=5, common.legend = TRUE, legend = 'bottom')
