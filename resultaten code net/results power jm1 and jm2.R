lijst<-list(sim_one_six_slow_power, sim_three_six_slow_power, sim_one_twelve_slow_power, sim_three_twelve_slow_power, 
            sim_one_eighteen_slow_power, sim_three_eighteen_slow_power, sim_one_six_interm_power, sim_three_six_interm_power,
            sim_one_twelve_interm_power, sim_three_twelve_interm_power, sim_one_eighteen_interm_power, 
            sim_three_eighteen_interm_power, sim_one_six_fast_power, sim_three_six_fast_power, sim_one_twelve_fast_power,
            sim_three_twelve_fast_power, sim_one_eighteen_fast_power, sim_three_eighteen_fast_power, sim_one_six_all_power,
            sim_three_six_all_power, sim_one_twelve_all_power, sim_three_twelve_all_power, sim_one_eighteen_all_power,
            sim_three_eighteen_all_power)
interval<-as.character(rep(c(1,3), 12)) #rep(vector, 2) want voor 2 sample sizen
duration<- as.factor(rep(c(rep(6,2), rep(12,2), rep(18,2)),4))
type<-c(rep('slow', 6), rep('moderate',6), rep('fast',6), rep('all', 6))
powers<-c()
counter<-1
for(i in lijst){
  powers[counter]<-i[[2]]
  counter<-counter + 1
}

resi_frame <- data.frame(Interval = interval, Duration = duration, type = type, Power = powers)

data_slow <- resi_frame[resi_frame$type=='slow',]
data_mod <- resi_frame[resi_frame$type=='moderate',]
data_fast <- resi_frame[resi_frame$type=='fast',]
data_all <- resi_frame[resi_frame$type=='all',]

library(ggplot2)
windowsFonts(tnr = windowsFont('Helvetica'))
p1<-ggplot(data_slow, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10), text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p2<-ggplot(data_mod, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p3<-ggplot(data_fast, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p4<-ggplot(data_all, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)


lijst2<-list(sim_one_six_slow_power_slope, sim_three_six_slow_power_slope, sim_one_twelve_slow_power_slope, sim_three_twelve_slow_power_slope, 
            sim_one_eighteen_slow_power_slope, sim_three_eighteen_slow_power_slope, sim_one_six_interm_power_slope, sim_three_six_interm_power_slope,
            sim_one_twelve_interm_power_slope, sim_three_twelve_interm_power_slope, sim_one_eighteen_interm_power_slope, 
            sim_three_eighteen_interm_power_slope, sim_one_six_fast_power_slope, sim_three_six_fast_power_slope, sim_one_twelve_fast_power_slope,
            sim_three_twelve_fast_power_slope, sim_one_eighteen_fast_power_slope, sim_three_eighteen_fast_power_slope, sim_one_six_all_power_slope,
            sim_three_six_all_power_slope, sim_one_twelve_all_power_slope, sim_three_twelve_all_power_slope, sim_one_eighteen_all_power_slope,
            sim_three_eighteen_all_power_slope)

interval2<-as.character(rep(c(1,3), 12)) #rep(vector, 2) want voor 2 sample sizen
duration2<- as.factor(rep(c(rep(6,2), rep(12,2), rep(18,2)),4))
type2<-c(rep('slow', 6), rep('moderate',6), rep('fast',6), rep('all', 6))
powers2<-c()
counter2<-1
for(i in lijst2){
  powers2[counter2]<-i[[2]]
  counter2<-counter2 + 1
}

resi_frame2 <- data.frame(Interval = interval2, Duration = duration2, type = type2, Power = powers2)

data_slow2 <- resi_frame2[resi_frame2$type=='slow',]
data_mod2 <- resi_frame2[resi_frame2$type=='moderate',]
data_fast2 <- resi_frame2[resi_frame2$type=='fast',]
data_all2 <- resi_frame2[resi_frame2$type=='all',]

library(ggplot2)
windowsFonts(tnr = windowsFont('Helvetica'))
p5<-ggplot(data_slow2, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10), text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p6<-ggplot(data_mod2, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p7<-ggplot(data_fast2, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)
p8<-ggplot(data_all2, aes(Duration, Power, color=Interval, group = Interval)) + geom_line() + geom_point()+
  theme(title = element_text(size = 10),text = element_text(family = 'tnr'))+ylim(0.72,0.81)

p_lab1<-ggplot()+theme_void()
p_lab2<-ggplot()+theme_void()+geom_text(aes(0,0,label='Slow'), family='tnr')
p_lab3<-ggplot()+theme_void()+geom_text(aes(0,0,label='Moderate'), family = 'tnr')
p_lab4<-ggplot()+theme_void()+geom_text(aes(0,0,label='Fast'), family = 'tnr')
p_lab5<-ggplot()+theme_void()+geom_text(aes(0,0,label='All'), family = 'tnr')
p_lab6<-ggplot()+theme_void()+geom_text(aes(0,0,label='JM1'), family = 'tnr')
p_lab7<-ggplot()+theme_void()+geom_text(aes(0,0,label='JM2'), family= 'tnr')


library(ggpubr)
hi<-ggarrange(p_lab1, p_lab6, p_lab7,p_lab2, p1, p5,p_lab3, p2, p6,p_lab4, p3,p7,p_lab5, p4,p8, ncol=3, nrow=5, common.legend = TRUE, legend = 'bottom')
